function [f, outputs, rates] = dPUdT_CombinedTransitions(t,PU,params)
% ODE function for the d/dt operator for the cross-bridge model of half-sarcomere.
%  first 2N-1 entries of PU represent p1(s,t)
%  second 2N-1 entries represent p2(s,t)
%  third-last entry is U_NR, the fraction of myosin heads in not-relaxed state
% then second last NP and last SL 

vel = params.Vums;

% Decompose State Variables from PU vector
ss = params.ss; % space size (length of the s for each of p1-p3)
dS = params.dS; % step size

p1 = PU(1:ss); p1(p1<0) = 0;
p2 = PU(ss+1:2*ss); p2(p2<0) = 0;
if length(PU) > 3*ss
    p3 = PU(2*ss+1:3*ss);
    Ns = 3; % number of states
else
    p3 = 0;
    Ns = 2; % number of states
end

if params.UseSuperRelaxed
    P_SR = PU(Ns*ss+1);
else
    P_SR = 0;
end

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

if params.UseSuperRelaxedADP
    P_SRD = PU(Ns*ss+6);
else
    P_SRD = 0;
end

% x_dash_n = 2.0;
if params.UseMaxwellDashpot
    x_dash = PU(Ns*ss+7); % current middle of the parallel visco-elasticity
    % 3. Maxwell Force (State-Dependent)
    F_Maxwell  = params.kSE_M * (SL - x_dash);
    
    % Update the internal dashpot state (Simple Euler or use your ODE solver)
    % dx_dash/dt = (Spring_Force) / Viscosity
    dx_dash_dt = F_Maxwell  / params.eta_M;

else
    x_dash = 0; F_Maxwell  = 0; dx_dash_dt = 0;
end





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

% sum of all probabilities
p1_0 = dS*sum(p1); 
p2_0 = dS*sum(p2); % p2_1 = dS*sum((sign(s+params.dr).*abs(s+params.dr).^params.estiff).*p2);
p3_0 = dS*sum(p3);
p3_1 = dS*sum((s+params.drp3).*p3);

if params.UseNegativeKstiff
    p1_1_pos = dS*sum(s.*p1.*(s>=0));
    p1_1_neg = dS*sum(s.*p1.*(s<=0));
    p2_1_pos = dS*sum((sign(s+params.dr).*abs(s+params.dr).^params.estiff).*p2.*(s >= -params.dr));
    p2_1_neg = dS*sum((sign(s+params.dr).*abs(s+params.dr).^params.estiff).*p2.*(s < -params.dr));
    F_active = params.kstiff3*(p3_1) + params.kstiff2*(p2_1_pos) + params.kstiff2_n*(p2_1_neg) + params.kstiff1*(p1_1_pos) + params.kstiff1_n*(p1_1_neg);     
    p1_1 = p1_1_neg + p1_1_pos;
    p2_1 = p2_1_pos + p2_1_neg;
else
    p1_1 = dS*sum(s.*p1);
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
if params.UsePassive && ~params.UseTitinIdentifiedPassive
    Lsc0    = 1.51;
    % gamma = 7.5;
    F_passive = F_passive + params.k_pas*max(SL-LSE - Lsc0, 0)^params.gamma; 
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

% default
F_total = F_active + F_passive;

% F_total = max(0, F_total);


if params.FudgeVmax && t > -1
    if params.kSE*LSE >= 0
        Force = max(params.MaxSlackNegativeForce, params.kSE*LSE);
        velHS = (Force - F_total)/params.mu;
    else
        Force = params.MaxSlackNegativeForce;
        % velHS = -(818*exp(-17*(SL - LSE)) + 22);
        velHS = -(params.FudgeA*SL^2 +params.FudgeB*SL + params.FudgeC);
        % velHS = -params.vmax;
    end
    dLSEdt = vel - velHS;
else
    Force = sign(LSE)*abs(max(params.MaxSlackNegativeForce, params.kSE*LSE)).^params.ekSE;
    % Force = params.kSE*LSE*(LSE >= 0) + params.kSEn*LSE*(LSE < 0);
    if (Force - F_total) > 0
        velHS = (Force - F_total)/params.mu;    
    else
        velHS = (Force - F_total)/params.mu_neg;
    end
    dLSEdt = vel - velHS;
end

Force = Force + params.mu2*vel + F_Maxwell;

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

%%
% quasi-equilibrium binding factor functions
% TODO move to evalModel for optim
% MgATP = params.MgATP;
% Pi = params.Pi;
% MgADP = params.MgADP;

g1 = 1; g2 = 1; f1 = 0; f2 = 1;

RTD = g2*params.kah*PT;
RDT = g2*params.kamh*PD;

RD1 = params.ka*PD*N_overlap; % to loosely attachemnt state

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
    
    % % Ensure monotonic increasing y (clip or enforce)
    % for i = 2:length(args)
    %     if args(i) <= args(i-1)
    %         args(i) = args(i-1) + 1e-3;
    %     end
    % end

    % Monotonic cubic interpolation    
    % if params.UseStrainDep4R1D
    %     R1D = params.kd*p1.*exp(-(s.^2) / (2*params.kd_sigma^2));
    % elseif isfield(params, 'PieceWiseStrainDepR1D') && ~isempty(params.PieceWiseStrainDepR1D)
        R1D = params.kd*p1.*ppval(params.PieceWiseStrainDepR1D, s);
    % else        
    %     R1D = params.kd*p1;
    % end

    R12 = params.k1*p1.*ppval(params.PieceWiseStrainDep, s);
    
    % if isfield(params, 'PieceWiseStrainDepR21') && ~isempty(params.PieceWiseStrainDepR21)    
        R21 = params.k_1.*p2.*ppval(params.PieceWiseStrainDepR21, s);
    % else
    %     R21 = p2*0;    
    % end
    
    % if isfield(params, 'A2_PieceWiseStrainDepX')
    %     error('Not implemented atm!');
    %     f2 = @(x)pchip(params.A2_PieceWiseStrainDepX, params.A2_PieceWiseStrainDepParams, x);
    % else
        % f2 = f;
    % elseif isfield(params, 'PieceWiseStrainDep2')
        R2 = params.k2*p2.*ppval(params.PieceWiseStrainDep2, s + params.dr2 - params.dr);
    % else
    %     R2 = params.k2*p2.*ppval(params.PieceWiseStrainDep, s + params.dr2 - params.dr);
    % end
    % plot(params.PieceWiseStrainDepX,params.PieceWiseStrainDepParams, 'o', s, f(s), 'x-', s, f(s+ params.dr2 - params.dr), '--')

elseif params.UseUniformTransitionFunc
    % the cycle goes: PT (ATP bound) <-> PD(ready) <-> P1 <-> P2 -> P3 -> PT
    % dPUdT_TransitionRates;
    sd = @(kx, alphaL, alphaR, dr,eL, eR) min(1e4, kx*(exp((alphaL*(s-dr)).^eL).*(s<dr) + exp((alphaR*(s-dr)).^eR).*(s>=dr)));
    R1D = p1.*sd(params.kd, params.alpha0_L, params.alpha0_R, params.dr0, 2, 2);
    
    R12 = p1.*sd(params.k1, params.alpha1, 0, params.dr1, 2, 2); % P1 to P2
    R21 = f1*p2.*sd(params.k_1, params.alpha_1, params.alpha_1, params.dr_1, 2, 2); % p2 to p1
    
    R2 = p2.*sd(params.k2, params.alpha2_L, params.alpha2_R, params.dr2, params.e2L, params.e2R);
else 
    % s_b = s;
    % s = s*( max(-1, -Force))
    strainDep = @(alpha, dr) min(1e4, exp((alpha*(s+dr)).^params.StrainExp));																		

    R12 = params.k1*p1.*exp(-params.alpha1*s); % P1 to P2
    R21 = f1*params.k_1*p2.*strainDep(params.alpha_1, params.dr_1); % p2 to p1

    if params.UseWDetachment
        R2 = params.k2_L*exp(-params.alpha2_L*(s - (params.dr2_L-1))) + params.k2*exp(-params.alpha2*(s - (params.dr2-1)).^2)  + params.k2_R*exp(params.alpha2_R*(s - (params.dr2_R-1)));
        R2 = min(1e4, R2.*p2);
    else
        if params.UseNegativeForceRip
            f = @(x,A,s) (x <= 1).* (1 + (A - 1).*(1 - x).^s) + (x > 1);
            R1D = params.kd*p1.*(strainDep(params.alpha0_L, params.dr0).*(s<= 0) ...
        + strainDep(params.alpha0_R, params.dr0).*(s> 0))*f(F_SR, 2, 3); %(exp(-params.alpha1*s)) + params.TK*(s>params.TK0).*s.*p1; % p1 to PU - detachment rate + tearing constant           
            kL = min(1e4, params.k2_L*((s+params.dr2_L)<=0).*(1 - exp(-(s+params.dr2_L)*params.alpha2_L)).^2)*f(F_SR, 2, 3);
         else
           R1D = params.kd*p1.*(strainDep(params.alpha0_L, params.dr0).*(s<= 0) ...
        + strainDep(params.alpha0_R, params.dr0).*(s> 0)); %(exp(-params.alpha1*s)) + params.TK*(s>params.TK0).*s.*p1; % p1 to PU - detachment rate + tearing constant            
            kL = min(1e4, params.k2_L*((s+params.dr2_L)<=0).*(1 - exp(-(s+params.dr2_L)*params.alpha2_L)).^2);           
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

%% Super relaxed transitions
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
    
    % AS_src = zeros(size(p2));
    % AS_trg = zeros(size(p2));
    % mask = (s < 0 & s > -params.dr);
    % AS_src(mask) = abs(s(mask)) * (params.a2RAS / params.dr);
    % trg = ((params.a2RAS / params.dr) - abs(s(mask)));
    % AS_trg(mask) = trg/sum(trg)*sum(AS_src);


    %% --- Parameters ---
    sp2 = s;
    % params.slope = 1000; % parametrized rate (s^-1/nm)
    % params.d_actin = 5.5e-3; % 5.5 nm
    % params.s_threshold = 5.5e-3;    
    % p2: N x 1 probability array
    % p0: scalar (or array) detached probability
    % s: N x 1 strain vector (must be monotonically increasing)
    
    % 1. Calculate Hopping Rate and Mass Leaving
    RA_K = max(0, params.slope/params.dr * (s - params.s_threshold_R)) + max(0, -params.slope/params.dr * (s + params.s_threshold_L));
    % RA_K = params.slope/params.dr * max(0, abs(sp2) - params.s_threshold);
    dp2_RAm = p2 .* RA_K;    
    
    % 2. Calculate Target Positions (Shifted toward center)
    s_target = (s > params.s_threshold_R).*(s - params.d_actin) + (s < - params.s_threshold_L).*(s+params.d_actin);  
    % s_target = sp2 - sign(sp2) .* params.d_actin;
    
    % 3. Find Indices using 'discretize'
    % This tells us which bin each target_s falls into
    L = discretize(s_target, s);
    R = L + 1;
    
    % 4. Bound check (Essential for safety)
    L = max(1, min(length(s)-1, L));
    R = L + 1;
            
    % 5. Linear Interpolation Weights
    dist = s(R) - s(L);
    w_R = (s_target - sp2(L)) ./ dist;
    w_L = 1 - w_R;
    % 6. Accumulate into p2 (Vectorized)
    % accumarray sums all mass falling into the same bin index
    dp2_RAL = accumarray(L, dp2_RAm .* w_L, [length(s), 1]);
    dp2_RAR = accumarray(R, dp2_RAm .* w_R, [length(s), 1]); 
    % plot(s, p2*100, s, dp2_RAm, '--',s, dp2_RAR,  ':', s, dp2_RAL,  ':', s, dp2_RAR + dp2_RAL -dp2_RAm, '--', LineWidth=2)
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
    if params.UseA1AttachmentKernel
        % s_0 = 1 + round(-s(1)/dS, 6); % strain position in space at 0
        % s_i0 = round(s_0);
        % if s_i0 > 0 && s_i0 <= length(s)
        %     halfSpan = ceil(Ax / dS);    % in bins
        % 
        %     x = (-halfSpan:halfSpan) * dS + s(s_i0);
        %     % jj = floor(max(1,s_i0-halfSpan)):ceil(min(length(dp1), s_i0+halfSpan));
        %     jj = (-halfSpan:halfSpan);
        %     x = jj*dS + s(s_i0);
        %     sigma = Ax/3;
        %     K = exp(-(x.^2)/(2*sigma^2));
        %     K = K/sum(K);
        %     K_vp = jj + s_i0 >= 1 & jj + s_i0 <= length(dp1); % Kernel valid positions in bounds
        %     s_pos = jj(K_vp) + s_i0;
        %     dp1(s_pos) = dp1(s_pos) + (RD1/dS)*K(K_vp)';
        % else
        %     warning([num2str(t) ' : Strain array is not around 0!']);
        % end
        % Find nearest grid index of zero
        s_i0 = 1 + floor(-s(1)/dS);   % much cheaper than round(x,6)    
        if s_i0 > params.cached.halfSpan && s_i0 <= length(dp1)-params.cached.halfSpan    
            % Valid positions (no clipping)
            s_pos = s_i0 +  params.cached.jj;    
            % Apply kernel mass - dS already divided in K0
            dp1(s_pos) = dp1(s_pos) + (RD1) * params.cached.K0;    
        else
            % Need boundary handling
            s_pos = s_i0 +  params.cached.jj;
            % Logical mask for valid indices
            vp = s_pos >= 1 & s_pos <= length(dp1);
            % Apply kernel mass - dS already divided in K0    
            dp1(s_pos(vp)) = dp1(s_pos(vp)) + (RD1) * params.cached.K0(vp);
        end

    elseif Ax == 0    
            [s_i0, s_i1, s_i0k] = attachmentPoint(s(1), dS, ss);
            dp1(s_i0) = dp1(s_i0) + s_i0k*(RD1/dS); % attachment
            dp1(s_i1) = dp1(s_i1) + (1-s_i0k)*(RD1/dS); % attachment
            % dp1(s_i0) = dp1(s_i0) + s_i0k*(RD1/dS)*(F_passive/10); % attachment
            % dp1(s_i1) = dp1(s_i1) + (1-s_i0k)*(RD1/dS)*(F_passive/10); % attachment
    else
        % inputs: dp1, RD1, dS, Ax, x0
        % grid indices: assume s_j = j*dS
        % [s_i0, s_i1, s_i0k] = attachmentPoint(s(1), params.dS, params.ss);% nearest left grid index
        % x_center = 0; % strain position in space at 0;
        
        % lets get that earlier in the code
        s_0 = 1 + round(-s(1)/dS, 6); % strain position in space at 0
        
        %% choose kernel support in indices
        halfspan = ceil(Ax/dS); 
        att = zeros(size(dp1));
        for jj = floor(max(1,s_0-halfspan)):ceil(min(length(dp1), s_0+halfspan))
            % jj = max(1, min(length(dp1), jj));
            xj = s(jj);
            dist = abs(xj - 0.0); % distance from 0
            % if dist <= Ax
            w = max(0,1 - dist/Ax);   % triangular weight
                % w = 1;
                % dp1(jj) = dp1(jj) + (RD1/dS) * w*dS/Ax;        
            att(jj) = w;
        end
        % if isempty(jj)
        %     % warning('%f: Strain array is not around 0!', t);
        %     disp([num2str(t) ' : Strain array is not around 0!']);
        % end
        %%
        dp1 = dp1 + att/sum(att)*RD1/dS;
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
        dp2(s_pop) = -p2(s_pop)*1e6;
        % dp1(s_pop) = -p1(s_pop)*1e6;
    end
end


if Ns == 2 
    f = [dp1; dp2; dU_SR; dNP; dSL;dLSEdt;dPD;dU_SRD;dx_dash_dt];
    outputs = [Force, F_active, F_passive, N_overlap, p1_0, p2_0, p1_1, p2_1, PT, F_Maxwell];    
elseif Ns == 3
    f = [dp1; dp2; dp3; dU_SR; dNP; dSL;dLSEdt;dPD;dU_SRD; dx_dash_dt];
    outputs = [Force, F_active, F_passive, N_overlap, p1_0, p2_0, p3_0, p1_1, p2_1, p3_1, PT, F_Maxwell];
end

rates = [RTD, RDT, RD1, sum([R1D, R12,R21,R2, XB_Ripped], 1)*dS, RSR2PT, RPT2SR, RSRD2PD, RPD2SRD, RSR2SRD, RSRD2SR, RT2, sum(R2D)*dS];

if params.DryRun
    f = zeros(size(PU));
    outputs = [0, 0, 0, 1, 0, 0, 0, 0, 1];
    rates = zeros(1, 13);
    return;
end

%% breakpints
if t > 2.765 % && any(PD > 0)
    % P_SR < 0 % && dU_SR < 0 
    % || any(~isreal(f)) || t > 0.012 % || t > 0 && (p1_0 + p2_0 + PD + P_SR) > 1
    numberofthebeast = 6678;
    
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