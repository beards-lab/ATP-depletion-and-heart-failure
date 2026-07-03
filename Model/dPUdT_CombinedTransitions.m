function [f, outputs, rates] = dPUdT_CombinedTransitions(t,PU,params)
% DPUDT_COMBINEDTRANSITIONS  ODE right-hand side for the cardiac cross-bridge model.
%
%   [F, OUTPUTS, RATES] = DPUDT_COMBINEDTRANSITIONS(T, PU, PARAMS) computes
%   the time derivatives of the state vector PU at time T given model
%   parameters PARAMS. Intended for use with ode15s via evaluateModel.
%
%   State vector layout:
%     PU = [p1(1..ss), p2(1..ss), [p3(1..ss),] P_SR, NP, SL, LSE, PD, [P_SRD,] [x_dash]]
%     where ss = params.ss is the number of strain bins.
%     p1, p2, [p3] — strain-discretized cross-bridge state probability distributions
%     P_SR          — super-relaxed (SRX) fraction
%     NP            — nucleotide phosphate (ADP·Pi) fraction
%     SL            — sarcomere length (um)
%     LSE           — serial elastic element length (um)
%     PD            — super-relaxed ADP state fraction (if UseSuperRelaxedADP)
%     P_SRD         — optional SRX-ADP state
%     x_dash        — Maxwell dashpot position (if UseMaxwellDashpot)
%
%   Inputs:
%     T       - current time (s)
%     PU      - state vector (see layout above)
%     PARAMS  - parameter struct from getParams; must include Vums (velocity, um/s)
%
%   Outputs:
%     F       - dPU/dt, time derivatives (same length as PU)
%     OUTPUTS - scalar diagnostics: [Force, F_active, F_passive, N_overlap,
%                 p1_0, p2_0, [p3_0,] p1_1, p2_1, [p3_1,] PT, F_Maxwell,
%                 f_lattice, f_saturation]
%     RATES   - transition rate scalars for post-processing
%
%   Behavior is controlled by boolean flags in PARAMS (UseOverlap,
%   UseSuperRelaxed, UseSerialStiffness, UseA1AttachmentKernel, etc.).
%
%   See also: evaluateModel, getParams

vel = params.Vums;

% TODO: profile computational bottlenecks before optimizing.
%   Run: profile on; RunBakersExp; profile viewer
%   Candidate hotspot: overlap/lattice factor recomputed on every ODE RHS call.
%   (SL changes continuously so true precomputation is not straightforward.)

% Named constants
MAX_RATE = 1e4;  % clamp for strain-dependent transition rates — prevents numerical blow-up during stiff phases
POP_RATE = 1e6;  % large rate used to forcibly drain cross-bridge state (UseA2Popping)

%% Unpack state vector
% Decompose State Variables from PU vector
ss = params.ss; % space size (length of the s for each of p1-p3)
dS = params.dS; % step size

p1 = PU(1:ss); p1(p1<0) = 0;
p2 = PU(ss+1:2*ss); p2(p2<0) = 0;
if params.NumberOfStates == 3
    p3 = PU(2*ss+1:3*ss); p3(p3<0) = 0;  % clamp negatives like p1/p2 (prevents 3-state solver stalls)
    Ns = 3; % number of states
else
    p3 = 0;
    Ns = 2; % number of states
end

P_SR = max(0, PU(Ns*ss+1));  % always read; transitions gated separately below

% if ~params.UseCa
    NP = 0;

% if ~params.UseSLInput
    SL = PU(Ns*ss + 3);
    dSL = vel;
% else
%     % too slow
%     if t >= params.datatable(end-1, 1)
%         % if the sim time is over the datatable length, hold the SL
%         SL = params.datatable(end, 2);
%         dSL = 0;
%     elseif t <= params.datatable(1, 1)
%         SL = params.datatable(1, 2);
%         dSL = 0;
%     else
%         % TODO make it faster in sorted list
%         % https://stackoverflow.com/questions/20166847/faster-version-of-find-for-sorted-vectors-matlab
%         i = find(params.datatable(:, 1) >= t,1,'First');    
%     %     i = min(length(params.datatable(:, 1))-1, i);
%         SL = params.datatable(i, 2);
%         if i == 1
%             dSL = 0;
%         else
%             dSL = (params.datatable(i, 2) - params.datatable(i-1, 2))/((params.datatable(i, 1) - params.datatable(i-1, 1)));
%         end
%     end
%     vel = dSL;
% end    

    
% P_SRs = 1 - P_SR;
LSE = PU(Ns*ss + 4); % length of the serial stiffness
PD = PU(Ns*ss + 5);

P_SRD = max(0, PU(Ns*ss+6));  % always read; transitions gated separately below

if params.UseMaxwellDashpot
    % Titin viscoelastic model: F_titin = F_elastic(SL_c) + X_visc
    % X_visc is a viscous stress deviation; decays to zero so F_titin decays toward F_elastic.
    % dX/dt = kSE_M * dSLc - X / eta_M   (computed after dLSEdt is known, below)
    X_visc = PU(Ns*ss+7);
    F_elastic_titin = params.k_pas * max((SL - LSE) - 1.51, 0)^params.gamma;
    F_Maxwell = F_elastic_titin*0 + X_visc;
    dx_dash_dt = 0; % placeholder; overwritten after dLSEdt
else
    X_visc = 0; F_Maxwell = 0; dx_dash_dt = 0;
end

% Registration availability (target-zone registration; the FV "shoulder").
% A_reg in [A0,1] is the single scalar state appended at the tail of PU
% (index Ns*ss+8) that gates attachment availability. Relaxation ODE in the
% assembly section; gate applied to ka_eff below. When the feature is off the
% state is absent and A_reg=1 (identity -> byte-compatible with the legacy model).
if params.UseRegistrationAvailability
    A_reg = min(1, max(params.A0, PU(Ns*ss + 8)));
else
    A_reg = 1;
end





%% Overlap, lattice spacing, and filament geometry corrections
% Sarcomere geometry
if params.UseOverlap
    L_thick = params.L_thick;% = 1.67; % Length of thick filament, um
    L_hbare = params.L_hbare;% = 0.10; % Length of bare region of thick filament, um
    L_thin = params.L_thin;  %= 1.20; % Length of thin filament, um

    % deltaR  = 0.010; % um    
    L_T_HS1 = min(L_thick*0.5, SL*0.5);
    L_T_HS2 = max((SL-LSE)*0.5 - ((SL-LSE)-L_thin),L_hbare*0.5);
    L_ov = L_T_HS1 - L_T_HS2; % Length of single overlap region
    N_overlap = (L_ov*2/(L_thick - L_hbare))^1;
    % N_overlap = min(1, 0.6 + (SL-1.8)); % this is just a function that Dan made up
else
    N_overlap = 1;
end

% Lattice spacing correction on attachment rate
% Wider lattice (shorter SL) reduces attachment; narrower lattice (longer SL) saturates at 1
if params.UseLatticeSpacing
    d10 = params.d10_ref * sqrt(params.SL_ref_lattice / SL);
    d_surface = d10 - params.R_thick - params.R_thin;
    excess = max(0, d_surface - params.d_optimal); % only penalize when lattice is wider than optimal
    f_lattice = exp(-(excess^2) / (2 * params.sigma_lattice^2));
else
    f_lattice = 1;
end

% strain - above the myosin heads is zero. 
% Negative - shorter, Positive - longer
% need to cut the change in two because half-sarcomere means half the speed
% and half the space change
if params.LegacyStrainFlipping
    s = params.s' + (-(SL - LSE) + params.LXBpivot)/2;
    s = flipud(-s);
else
    s = params.s - (-(SL - LSE) + params.LXBpivot)/2;
end

%% Force calculation
% sum of all probabilities
p1_0 = dS*sum(p1);
p2_0 = dS*sum(p2); % p2_1 = dS*sum((sign(s+params.dr).*abs(s+params.dr).^params.estiff).*p2);
p3_0 = dS*sum(p3);
p3_1 = dS*sum((s+params.drp3).*p3);

% Route B — p1 force offset: s1 = s + dr1 makes the PRE-stroke state p1 force-bearing
% (dr1>0 = pre-stroke force; dr1<0 = back-force). dr1=0 (default) recovers the original
% Σ(s·p1)≈0 exactly. Decouples force-duty (p1+p2) from the p2 detachment that welds ktr<->FV.
s1 = s + params.dr1;
if params.UseNegativeKstiff
    p1_1_pos = dS*sum(s1.*p1.*(s1>=0));
    p1_1_neg = dS*sum(s1.*p1.*(s1<=0));
    p2_1_pos = dS*sum((sign(s+params.dr).*abs(s+params.dr).^params.estiff).*p2.*(s >= -params.dr));
    p2_1_neg = dS*sum((sign(s+params.dr).*abs(s+params.dr).^params.estiff).*p2.*(s < -params.dr));
    F_active = params.kstiff3*(p3_1) + params.kstiff2*(p2_1_pos) + params.kstiff2_n*(p2_1_neg) + params.kstiff1*(p1_1_pos) + params.kstiff1_n*(p1_1_neg);
    p1_1 = p1_1_neg + p1_1_pos;
    p2_1 = p2_1_pos + p2_1_neg;
else
    p1_1 = dS*sum(s1.*p1);
    p2_1 = dS*sum((sign(s+params.dr).*abs(s+params.dr).^params.estiff).*p2);
    F_active = params.kstiff3*(p3_1) + params.kstiff2*(p2_1) + params.kstiff1*(p1_1);
    % F_active = max(params.kstiff2*(p2_1), 0) + max(params.kstiff1*(p1_1), 0);
end

% non-hydrolized ATP in non-super relaxed state
PT = max(0, 1*(1.0 - NP) - (p1_0 + p2_0 + p3_0 + PD + P_SR + P_SRD)); % unattached permissive fraction - 
    

if params.UseOverlapFactor
    N_overlap = (N_overlap - (p1_0 + p2_0))/(P_SR + P_SRD + PT + PD); % overlap factor
end    

F_passive = 0;
if params.UsePassive && ~isempty(params.LoadPassiveMat)
    F_passive = sigLookup(params.driveSig, t);
else
    if params.UsePassive && ~params.UseTitinIdentifiedPassive
        % Lsc0    = 1.51;
        % gamma = 7.5;
        F_passive = F_passive + params.k_pas*max(SL-LSE -  params.Lsc0, 0)^params.gamma; 
    elseif params.UsePassive && params.UseTitinIdentifiedPassive
        % identified from FitPassiveRampUp.m
        y = @(k_pas, x0, gamma, x) k_pas.*(x-x0).^gamma - 4*0 - x0*0 + 0.5e9.*(x-0.95).^13;    
        F_passive = y(0.4, -0.4, 7.9, (SL-LSE)/2);
        % F_passive = y(0.4, -0.4, 7.9, (SL-LSE+LSE)/2);
    end
    
    if params.UseTitinInterpolation
        F_passive = F_passive + max(0, interp1(params.TitinTable.Time, params.TitinTable.ForceV, ...
            min(max(t, params.TitinTable.Time(1)), params.TitinTable.Time(end)), "linear") - params.TitinTable.ForceV(end));
    end
    
    % H2: Dynamic passive — titin viscoelasticity during rapid restretch.
    % Titin is stiffer at high stretch rates; after slack it re-engages with a
    % transient force proportional to restretch velocity.
    if params.UseDynamicPassive && vel > 0
        F_passive = F_passive + params.c_titin_visc * vel;
    end
end

if params.UseMaxwellDashpot
    F_passive = F_passive + F_Maxwell;
end

% H3: Driving signal passive — time-varying force from preprocessed lookup.
% params.LoadPassiveMat must be set and sig loaded into params.driveSig by getParams.


% default
F_total = F_active + F_passive;

% F_total = max(0, F_total);


if params.FudgeVmax && t > -1
    Force = max(params.MaxSlackNegativeForce, params.kSE*LSE);
    if params.kSE*LSE >= 0
        velHS = (Force - F_total)/params.mu;
    else
        Force = params.MaxSlackNegativeForce;
        velHS = -(params.FudgeA*SL^2 +params.FudgeB*SL + params.FudgeC);
    end
    dLSEdt = vel - velHS;
else
    Force = sign(LSE)*abs(max(params.MaxSlackNegativeForce, params.kSE*LSE)).^params.ekSE;
    % Force = params.kSE*LSE*(LSE >= 0) + params.kSEn*LSE*(LSE < 0);
    if params.UseViscoelasticSE && vel > 0
        % Kelvin-Voigt SE: F_SE = kSE*LSE + c_SE*dLSEdt
        % Implicit solve: (mu + c_SE)*dLSEdt = mu*vel + F_total - Force
        % → dLSEdt = (mu*vel + F_total - Force) / (mu + c_SE)
        % Effect: SE resists rapid extension → more restretch strain reaches cross-bridges.
        % Inactive at vel≤0 so steady-state and slack phases are unaffected.
        if (Force - F_total) > 0
            dLSEdt = (params.mu     * vel + F_total - Force) / (params.mu     + params.c_SE_visc);
        else
            dLSEdt = (params.mu_neg * vel + F_total - Force) / (params.mu_neg + params.c_SE_visc);
        end
    else
        if (Force - F_total) > 0
            velHS = (Force - F_total)/params.mu;
        else
            velHS = (Force - F_total)/params.mu_neg;
        end
        dLSEdt = vel - velHS;
    end
end

Force = Force + params.mu2*vel;

if params.UseMaxwellDashpot
    dSLc = vel - dLSEdt;
    dx_dash_dt = params.kSE_M * max(0, dSLc) - X_visc / params.eta_M;
end

%% TRANSITIONS
% plotStateTransitionsFlag = true;
if params.justPlotStateTransitionsFlag
    % s = s - (s(end) - s(1))/2; 
    s = (-0.02:0.001:0.02)';
    F_total = linspace(-50, 80, 8);
    F_passive = linspace(-5, 10, 8);
    p1 = ones(size(s));
    p2 = ones(size(s));
    
    PD = 1;PT = 1;P_SR = 1;P_SRD = 1;
end

%% Transition rates
% quasi-equilibrium binding factor functions
% TODO move to evalModel for optim
% MgATP = params.MgATP;
% Pi = params.Pi;
% MgADP = params.MgADP;

g1 = 1; g2 = 1; f1 = 0; f2 = 1;

RTD = g2*params.kah*PT;
RDT = g2*params.kamh*PD;

% H4: Stretch activation — boost ka during rapid restretch.
if params.UseStretchActivation && vel > 0
    ka_eff = params.ka * (1 + params.k_SA * vel);
else
    ka_eff = params.ka;
end
% H5: 1-Gaussian velocity-dependent attachment modifier.
% At isometric (v=0): target zone congestion reduces available sites.
% At shortening velocities: old heads clear sites -> recovery toward full ka.
% Formula: f = 1 - amplitude * exp(-(|v| - center)^2 / (2*sigma^2))
if params.UseVelGaussAttachment
    v_abs = abs(velHS);
    f_vel_gauss = 1 - params.v_att_amplitude * exp(-(v_abs - params.v_att_center).^2 / (2 * params.v_att_sigma^2));
    ka_eff = ka_eff * max(0, f_vel_gauss);
end
% Registration-availability gate: multiply attachment by the slow availability
% state A_reg. Driven by the IMPOSED vel (not velHS) so SE velocity ringing cannot
% corrupt it -> ktr stays single-order. A_reg==1 when the feature is off (no change).
if params.UseRegistrationAvailability
    ka_eff = ka_eff * A_reg;
end
RD1 = ka_eff*PD*N_overlap*f_lattice; % to loosely attachment state

% Global (mean-field) binding-site occupancy saturation of attachment.
% The strain axis is not a site axis, so the only faithful occupancy proxy in
% this strain-distribution model is the total bound fraction P_bound. Attenuate
% the whole attachment flux as occupancy rises toward the physical ceiling
% P_bound_max (rat-cardiac max-Ca isometric ~0.12-0.15). p1_0/p2_0/p3_0 are the
% integrated bound fractions computed above.
if params.UseGlobalOccupancySaturation
    P_bound = p1_0 + p2_0 + p3_0;
    if strcmpi(params.OccupancyForm, 'langmuir')
        f_sat_global = 1 / (1 + P_bound / params.P_bound_max);
    else
        f_sat_global = max(0, 1 - P_bound / params.P_bound_max);
    end
    RD1 = RD1 * f_sat_global;
end

if params.UsePassiveForSR
    F_SR = F_passive;
else
    if params.useHalfActiveForSR
        F_SR  = F_active/2 + F_passive;
    else
    % default
        F_SR  = F_total;
    end
end

if params.UsePieceWiseStrainDep
    
    if params.UseFastPPval && isfield(params.cached, 'pwsd_brk')
        % Fast path: inline piecewise cubic evaluation — no ppval/unmkpp call overhead.
        % Breaks (row) and coefs (n_segs×4) were pre-unpacked in getParams Step 6b.
        % Segment search: count breaks strictly less than s (vectorised broadcast).
        % Broadcasting: s (ss×1) vs b(1:end-1) (1×k) → (ss×k) logical, then sum→(ss×1).

        % b stored as (k+1)×1 column → b(1:end-1)' is (1×k) row for broadcast
        % s (ss×1) >= b_row (1×k) → (ss×k); sum along dim 2 → (ss×1)
        % b(si) with b column and si (ss×1) always returns (ss×1)  ← key invariant

        b = params.cached.pwsdR1D_brk; c = params.cached.pwsdR1D_cof;
        si = max(1, min(size(c,1), sum(s >= b(1:end-1)', 2)));
        dx = s - b(si); R1D = params.kd  * p1 .* (((c(si,1).*dx + c(si,2)).*dx + c(si,3)).*dx + c(si,4));

        b = params.cached.pwsd_brk;    c = params.cached.pwsd_cof;
        si = max(1, min(size(c,1), sum(s >= b(1:end-1)', 2)));
        dx = s - b(si); R12 = params.k1  * p1 .* (((c(si,1).*dx + c(si,2)).*dx + c(si,3)).*dx + c(si,4));

        b = params.cached.pwsdR21_brk; c = params.cached.pwsdR21_cof;
        si = max(1, min(size(c,1), sum(s >= b(1:end-1)', 2)));
        dx = s - b(si); R21 = params.k_1 * p2 .* (((c(si,1).*dx + c(si,2)).*dx + c(si,3)).*dx + c(si,4));

        s2 = s + params.dr2 - params.dr;
        b = params.cached.pwsd2_brk;   c = params.cached.pwsd2_cof;
        si = max(1, min(size(c,1), sum(s2 >= b(1:end-1)', 2)));
        dx = s2 - b(si); R2 = params.k2  * p2 .* (((c(si,1).*dx + c(si,2)).*dx + c(si,3)).*dx + c(si,4));
    else
        % Original path: piecewise cubic (pchip) evaluation via ppval.
        R1D = params.kd*p1.*ppval(params.PieceWiseStrainDepR1D, s);
        R12 = params.k1*p1.*ppval(params.PieceWiseStrainDep, s);
        R21 = params.k_1.*p2.*ppval(params.PieceWiseStrainDepR21, s);
        R2  = params.k2*p2.*ppval(params.PieceWiseStrainDep2, s + params.dr2 - params.dr);
    end

elseif params.UseUniformTransitionFunc
    % the cycle goes: PT (ATP bound) <-> PD(ready) <-> P1 <-> P2 -> P3 -> PT
    % dPUdT_TransitionRates;
    sd = @(kx, alphaL, alphaR, dr,eL, eR) min(MAX_RATE, kx*(exp((alphaL*(s-dr)).^eL).*(s<dr) + exp((alphaR*(s-dr)).^eR).*(s>=dr)));
    R1D = p1.*sd(params.kd, params.alpha0_L, params.alpha0_R, params.dr0, 2, 2);
    
    R12 = p1.*sd(params.k1, params.alpha1, 0, params.dr1, 2, 2); % P1 to P2
    R21 = f1*p2.*sd(params.k_1, params.alpha_1, params.alpha_1, params.dr_1, 2, 2); % p2 to p1
    
    R2 = p2.*sd(params.k2, params.alpha2_L, params.alpha2_R, params.dr2, params.e2L, params.e2R);
else 
    % s_b = s;
    % s = s*( max(-1, -Force))
    strainDep = @(alpha, dr) min(MAX_RATE, exp((alpha*(s+dr)).^params.StrainExp));																		

    R12 = params.k1*p1.*exp(-params.alpha1*s); % P1 to P2
    R21 = f1*params.k_1*p2.*strainDep(params.alpha_1, params.dr_1); % p2 to p1

    if params.UseWDetachment
        R2 = params.k2_L*exp(-params.alpha2_L*(s - (params.dr2_L-1))) + params.k2*exp(-params.alpha2*(s - (params.dr2-1)).^2)  + params.k2_R*exp(params.alpha2_R*(s - (params.dr2_R-1)));
        R2 = min(MAX_RATE, R2.*p2);
    else
        if params.UseNegativeForceRip
            f = @(x,A,s) (x <= 1).* (1 + (A - 1).*(1 - x).^s) + (x > 1);
            R1D = params.kd*p1.*(strainDep(params.alpha0_L, params.dr0).*(s<= 0) ...
        + strainDep(params.alpha0_R, params.dr0).*(s> 0))*f(F_SR, 2, 3); %(exp(-params.alpha1*s)) + params.TK*(s>params.TK0).*s.*p1; % p1 to PU - detachment rate + tearing constant           
            kL = min(MAX_RATE, params.k2_L*((s+params.dr2_L)<=0).*(1 - exp(-(s+params.dr2_L)*params.alpha2_L)).^2)*f(F_SR, 2, 3);
         else
           R1D = params.kd*p1.*(strainDep(params.alpha0_L, params.dr0).*(s<= 0) ...
        + strainDep(params.alpha0_R, params.dr0).*(s> 0)); %(exp(-params.alpha1*s)) + params.TK*(s>params.TK0).*s.*p1; % p1 to PU - detachment rate + tearing constant            
            kL = min(MAX_RATE, params.k2_L*((s+params.dr2_L)<=0).*(1 - exp(-(s+params.dr2_L)*params.alpha2_L)).^2);
        end
        kR = max(0, params.k2_R*(s-params.dr2_R)).^params.alpha2_R; %.*(s>0.002);
        R2 = p2.*(params.k2 + kL + kR);
    end
end

R3 = params.k3*p3;
R3m = params.k3m*p3;

% to PT state directly
XB_Ripped = params.k2rip*p2.*min(1e9, max(0, exp(params.alphaRip*(s+params.dr3))));

if params.UseStrictDetachmentAt > 0
    strictArea = s > params.UseStrictDetachmentAt | s < -params.UseStrictDetachmentAt;
    R2(strictArea) = p2(strictArea)*(10000);
    R1D(strictArea) = p1(strictArea)*(10000);
end

% H3: Catch bonds — during rapid restretch, reduce detachment rates.
% Physiological basis: actomyosin bonds strengthen under rapid stretch
% (catch bond behaviour), producing the initial force spike. As velocity
% drops the catch effect disappears, bridges detach → valley, then
% reattachment → second peak.
if params.UseCatchBond && vel > 0
    catch_factor = max(0, 1 - params.k_catch_bond * vel);
    R1D = R1D * catch_factor;
    if ~params.UseCatchBondR1DOnly
        % Strain-limited catch: apply R2 suppression only below CatchBondStrainMax.
        % At high strain (slip region) bridges detach normally — fixes first-peak dSL timing.
        catch_mask = s <= params.CatchBondStrainMax;
        R2(catch_mask)  = R2(catch_mask)  * catch_factor;
    end
end

%% Super-relaxed (SRX/DRX) dynamics
% F_SR = (SL-LSE);
if params.UseSuperRelaxedADP
    RSRD2PD = params.kmsrd*exp(F_SR/params.sigma_srd1)*max(0, P_SRD);
    if params.UseLimitedSRDTransition
        RPD2SRD = params.Srd_max./(1+exp(params.Srd_n*(F_SR)))*PD;
    else
        RPD2SRD = params.ksrd*exp(-F_SR/params.sigma_srd2)*PD;
    end    
else
    RSRD2PD = 0;
    RPD2SRD = 0;
end

if params.UseSuperRelaxed     
    RSR2PT = params.kmsr*exp(F_SR/params.sigma1)*max(0, P_SR);
	% RSR2PT = params.kmsr*(1 + log(max(1, F_SR))*params.sigma1)*P_SR;
    RPT2SR = params.ksr0*exp(-F_SR/params.sigma2)*max(0, PT);
else 
    RSR2PT = 0;
    RPT2SR = 0;    
end


if params.UseSuperRelaxed && params.UseSuperRelaxedADP
    RSRD2SR = params.ksrd2sr*max(0, P_SRD);    
    RSR2SRD = params.ksr2srd*max(0, P_SR);
else
    RSRD2SR = 0;    
    RSR2SRD = 0;
end


dU_SRD = RSR2SRD - RSRD2SR - RSRD2PD  + RPD2SRD ;

if params.UseDirectSRXTransition
    dU_SR = -RSR2PT  + RPT2SR + sum(R2)*dS + RSRD2SR - RSR2SRD;
else    
    dU_SR = -RSR2PT  + RPT2SR + RSRD2SR - RSR2SRD;
end

if params.UseA2AttachmentShift

    % A2 hopping: strained p2 heads hop to neighboring actin site (d_actin apart)
    slope_over_dr = params.slope / params.dr;

    % 1. Hopping rate: ramps up beyond threshold strain
    RA_K = max(0, slope_over_dr * (s - params.s_threshold_R)) ...
         + max(0, -slope_over_dr * (s + params.s_threshold_L));
    dp2_RAm = p2 .* RA_K;

    % 2. Target positions (shifted toward center by one actin monomer)
    s_target = (s > params.s_threshold_R) .* (s - params.d_actin) ...
             + (s < -params.s_threshold_L) .* (s + params.d_actin);

    % 3. Bin indices via arithmetic (uniform grid, replaces discretize)
    L = max(1, min(ss-1, 1 + floor((s_target - s(1)) / dS)));
    R = L + 1;

    % 4. Linear interpolation weights
    w_R = (s_target - s(L)) ./ (s(R) - s(L));
    w_L = 1 - w_R;

    % 5. Accumulate hopping mass into target bins
    dp2_RAL = accumarray(L, dp2_RAm .* w_L, [ss, 1]);
    dp2_RAR = accumarray(R, dp2_RAm .* w_R, [ss, 1]);
else
    dp2_RAm = 0; dp2_RAL = 0; dp2_RAR = 0;
end


if params.UseA2Reattaching
    [s_i0_a2, s_i1_a2, s_i0k_a2] = attachmentPoint(s(1) + params.dr, params.dS, params.ss);
    RT2 = params.k_2*PT;
    dp2(s_i0_a2) = dp2(s_i0_a2) + s_i0k*(RT2/dS); % attachment
    dp2(s_i1_a2) = dp2(s_i1_a2) + (1-s_i0k_a2)*(RT2/dS); % attachment    
    % PT state is calculated as a complement of all states to 1, so we do
    % not have to specify the outflow
    % PT = ...
else
    RT2 = 0;
end

if params.UseA2MechanicalRecocking
    R2D = params.k2d*p2.*(s > params.drmr);
else
    R2D = 0;
end

if params.justPlotStateTransitionsFlag    
    plotStateTransitions;
    % plotStateTransitions_subpanels;
    error('Quitting after plotting states');
end
%% governing flows
% PT - calculated as complement of sum of all probabilities
% state 0: unattached, ATP-cocked

dPD = RSRD2PD - RPD2SRD + RTD - RDT - RD1 + sum(R1D)*dS + sum(R2D)*dS ;
% dPD = RSRD2PD - RPD2SRD + RTD - RDT - RD1 + sum(R1D);
dp1 = - R1D -  R12 + R21; % state 1: loosely attached, just sitting&waiting
dp2 = + R12 - R21  - R2 - XB_Ripped + R3m - R2D + dp2_RAR + dp2_RAL - dp2_RAm; % strongly attached, post-ratcheted: hydrolyzed ATP to ADP, producing Pi - ready to ratchet
dp3 = - R3 + R3m + R2;


% if ~params.UseMutualPairingAttachment && params.UseSpaceInterpolation
% Attachment
try
    % if params.A1AttachmentWidth <= dS
    Ax = params.A1AttachmentWidth;
    % Vernier / target zone modulation of attachment (per-bin or velocity-based)
    f_sat = ones(size(p1)); % default: no modulation
    if params.UseTargetZoneSaturation
        % [DEPRECATED as site occupancy — see report] Per-strain-bin attachment
        % modifier. dS-invariant: compare attached-head DENSITY (p1+p2+p3) [1/um]
        % to a density cap rho_attach_max, so the grid spacing dS cancels. (The
        % old form (p1+p2+p3)*dS / max_attached_per_bin was per-bin MASS and thus
        % grid-dependent.) NOTE: this penalizes heads sharing a strain value, not
        % heads sharing an actin site; use UseGlobalOccupancySaturation for the
        % physically faithful occupancy term.
        f_sat = max(0, 1 - (p1 + p2 + p3) / params.rho_attach_max);
    elseif params.UseVernierVelocity
        v_hs = abs(velHS);
        % f_sat_scalar = 1 + params.alpha_vernier * v_hs / (v_hs + params.v_ref_vernier);
        % f_sat = f_sat_scalar * ones(size(p1));
        f_sat = params.alpha_vernier + (1 - params.alpha_vernier)./(1+params.v_ref_vernier./v_hs);
    end

    if params.UseA1AttachmentKernel
        % Find nearest grid index of zero
        s_i0 = 1 + floor(-s(1)/dS);
        if s_i0 > params.cached.halfSpan && s_i0 <= length(dp1)-params.cached.halfSpan
            s_pos = s_i0 + params.cached.jj;
            dp1(s_pos) = dp1(s_pos) + (RD1) * params.cached.K0 .* f_sat(s_pos);
        else
            s_pos = s_i0 + params.cached.jj;
            vp = s_pos >= 1 & s_pos <= length(dp1);
            dp1(s_pos(vp)) = dp1(s_pos(vp)) + (RD1) * params.cached.K0(vp) .* f_sat(s_pos(vp));
        end

    elseif Ax == 0
            [s_i0, s_i1, s_i0k] = attachmentPoint(s(1), dS, ss);
            dp1(s_i0) = dp1(s_i0) + s_i0k*(RD1/dS)*f_sat(s_i0);
            dp1(s_i1) = dp1(s_i1) + (1-s_i0k)*(RD1/dS)*f_sat(s_i1);
    else
        s_0 = 1 + round(-s(1)/dS, 6);
        halfspan = ceil(Ax/dS);
        att = zeros(size(dp1));
        for jj = floor(max(1,s_0-halfspan)):ceil(min(length(dp1), s_0+halfspan))
            xj = s(jj);
            dist = abs(xj - 0.0);
            w = max(0,1 - dist/Ax);
            att(jj) = w;
        end
        att_norm = att/max(1, sum(att))/dS;
        dp1 = dp1 + att_norm .* f_sat .* RD1;
    end
catch e
    
    disp([num2str(t) ' : ' e.message]);
end

% if ~params.UseCa
    dNP = 0;

if params.UseA2Popping
    s_pop = s <= params.pop_s;
    % ignoring at this level
    % all probabilities are auto shifted to UT
    if any(s_pop==1) && t > 0
        dp2(s_pop) = -p2(s_pop)*POP_RATE;
        % dp1(s_pop) = -p1(s_pop)*1e6;
    end
end


% Scalar f_saturation for output: value at attachment point (s≈0)
if params.UseTargetZoneSaturation
    s_i0_out = max(1, min(ss, 1 + floor(-s(1)/dS)));
    f_saturation = f_sat(s_i0_out);
elseif params.UseVernierVelocity
    f_saturation = f_sat;
else
    f_saturation = 1;
end

%% Assemble ODE derivatives and outputs
if Ns == 2
    f = [dp1; dp2; dU_SR; dNP; dSL;dLSEdt;dPD;dU_SRD;dx_dash_dt];
    outputs = [Force, F_active, F_passive, N_overlap, p1_0, p2_0, p1_1, p2_1, PT, F_Maxwell, f_lattice, f_saturation];
elseif Ns == 3
    f = [dp1; dp2; dp3; dU_SR; dNP; dSL;dLSEdt;dPD;dU_SRD; dx_dash_dt];
    outputs = [Force, F_active, F_passive, N_overlap, p1_0, p2_0, p3_0, p1_1, p2_1, p3_1, PT, F_Maxwell, f_lattice, f_saturation];
end

% Registration-availability relaxation (FV shoulder). First-order approach to the
% velocity-set target A_inf with time constant tau_reg (~1/ktr -> adds no new slow
% mode; the isometric depression appears as a lower plateau, not a separate phase).
% Driven by the IMPOSED vel (params.Vums). Appended as the final state (Ns*ss+8);
% outputs/rates are deliberately left unchanged to avoid shifting downstream indices.
if params.UseRegistrationAvailability
    % Speed that drives availability. Default: |vel| (symmetric). Shortening-only
    % (max(0,-vel)) means restretch/lengthening does NOT boost attachment, so the
    % post-restretch undershoot is preserved while the FV (shortening) shoulder is kept.
    if isfield(params, 'RegAvailShorteningOnly') && params.RegAvailShorteningOnly
        vreg = max(0, -vel);
    else
        vreg = abs(vel);
    end
    A_inf  = params.A0 + (1 - params.A0) * vreg / (vreg + params.v_ref_reg);
    dA_reg = (A_inf - A_reg) / params.tau_reg;
    f = [f; dA_reg];
end

rates = [RTD, RDT, RD1, sum([R1D, R12,R21,R2, XB_Ripped], 1)*dS, RSR2PT, RPT2SR, RSRD2PD, RPD2SRD, RSR2SRD, RSRD2SR, RT2, sum(R2D)*dS];

if params.DryRun
    f = zeros(size(PU));
    % outputs = [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1];
    rates = zeros(size(rates));
    return;
end

%% breakpints
if any(isnan(f)) || t > 76 % && any(PD > 0)
    % P_SR < 0 % && dU_SR < 0 
    % || any(~isreal(f)) || t > 0.012 % || t > 0 && (p1_0 + p2_0 + PD + P_SR) > 1
    numberofthebeast = 66788;
    
end

% disp('oj')
end

function [s_i0, s_i1, s_i0k] = attachmentPoint(s0, dS, ss)
    % estimate the position of the actual index of zero strain
    % IMPORTANT: s MUST be around 0 somewhere!
    s_p0 = 1 + round(-s0/dS, 6); % strain position in space at 0
    
    if isnan(s_p0) || floor(s_p0) < 0
        % this can happen during numerical rottfinding iteration step, should be avoided in the result
        s_i0 = 1; % just hold on..
        s_i1 = 2;
        s_i0k = 1;
        % warning(sprintf('Out of bounds! At %0.6fs and SL %0.2fum, the s(1) was %0.2f', t, SL, s(1)));
    elseif floor(s_p0) == 0
        s_i0 = 1;
        s_i1 = 2;
        s_i0k = 1;
    elseif ceil(s_p0) == ss
        s_i0 = ss - 1;
        s_i1 = ss;
        s_i0k = 0;
    elseif ceil(s_p0) > ss
        s_i0 = ss - 1;
        s_i1 = ss;
        s_i0k = 0;
        % error(sprintf('Out of bounds! At %0.6fs and SL %0.2fum, the s(end) was %0.2f', t, SL, s(end)));
    else
        % if params.UseSpaceInterpolation
        s_i0 = floor(s_p0);
        s_i1 = ceil(s_p0);
        s_i0k = s_i1 - s_p0;
    end
end

