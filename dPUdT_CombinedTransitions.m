function [f, outputs, rates] = dPUdT_CombinedTransitions(t,PU,params)
% ODE function for the d/dt operator for the cross-bridge model of half-sarcomere.
%  first 2N-1 entries of PU represent p1(s,t)
%  second 2N-1 entries represent p2(s,t)
%  third-last entry is U_NR, the fraction of myosin heads in not-relaxed state
% then second last NP and last SL 

vel = params.Vums;

% Decompose State Variables from PU vector
ss = params.ss; % space size (length of the s for each of p1-p3)
p1 = PU(1:ss);
p2 = PU(ss+1:2*ss);
if params.UseSuperRelaxed
    P_SR = PU(2*ss+1);
else
    P_SR = 0;
end

% if ~params.UseCa
    NP = 0;

if ~params.UseSLInput
    SL = PU(2*ss + 3);
    dSL = vel;
else
    if t >= params.datatable(end-1, 1)
        % if the sim time is over the datatable length, hold the SL
        SL = params.datatable(end, 2);
        dSL = 0;
    elseif t <= params.datatable(1, 1)
        SL = params.datatable(1, 2);
        dSL = 0;
    else
        % TODO make it faster in sorted list
        % https://stackoverflow.com/questions/20166847/faster-version-of-find-for-sorted-vectors-matlab
        i = find(params.datatable(:, 1) >= t,1,'First');    
    %     i = min(length(params.datatable(:, 1))-1, i);
        SL = params.datatable(i, 2);
        if i == 1
            dSL = 0;
        else
            dSL = (params.datatable(i, 2) - params.datatable(i-1, 2))/((params.datatable(i, 1) - params.datatable(i-1, 1)));
        end
    end
    vel = dSL;
end    

    
% P_SRs = 1 - P_SR;
LSE = PU(2*ss + 4); % length of the serial stiffness
PD = PU(2*ss + 5);


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
s = params.s' + (-(SL - LSE) + params.LXBpivot)/2;
% s = params.s' - (-(SL - LSE) + params.LXBpivot)/2;
s = flipud(-s);
dS = params.dS;

% estimate the position of the actual index of zero strain
% IMPORTANT: s MUST be around 0 somewhere!
s_p0 = 1 + round(-s(1)/params.dS, 6); % strain position in space at 0

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
elseif ceil(s_p0) == params.ss
    s_i0 = params.ss - 1;
    s_i1 = params.ss;
    s_i0k = 0;
elseif ceil(s_p0) > params.ss
    s_i0 = params.ss - 1;
    s_i1 = params.ss;
    s_i0k = 0;
    error(sprintf('Out of bounds! At %0.6fs and SL %0.2fum, the s(end) was %0.2f', t, SL, s(end)));
else
    % if params.UseSpaceInterpolation
        s_i0 = floor(s_p0);
        s_i1 = ceil(s_p0);
        s_i0k = s_i1 - s_p0;
end

% sum of all probabilities
p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum((s+params.dr).*p2);

% non-hydrolized ATP in non-super relaxed state
PT = max(0, 1*(1.0 - NP) - (p1_0 + p2_0 + PD + P_SR)); % unattached permissive fraction - 
    F_active = params.kstiff2*(p2_1) + params.kstiff1*(p1_1);     

if params.UseOverlapFactor
    N_overlap = (N_overlap - (p1_0 + p2_0))/(P_SR + PT + PD); % overlap factor
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
    Force = max(params.MaxSlackNegativeForce, params.kSE*LSE);
    % Force = params.kSE*LSE*(LSE >= 0) + params.kSEn*LSE*(LSE < 0);
    velHS = (Force - F_total)/params.mu;
    dLSEdt = vel - velHS;
end

%% TRANSITIONS
% plotStateTransitionsFlag = true;
if params.justPlotStateTransitionsFlag
    % s = s - (s(end) - s(1))/2; 
    s = -0.05:0.001:0.05;
    F_total = -5:0.1:80;
    F_passive = -5:0.1:10;
    p1 = ones(size(p1));
    p2 = ones(size(p2));
    
    PD = 1;PT = 1;P_SR = 1;
end

%%
% quasi-equilibrium binding factor functions
% TODO move to evalModel for optim
MgATP = params.MgATP;
Pi = params.Pi;
MgADP = params.MgADP;

g1 = 1; g2 = 1; f1 = 0; f2 = 1;

RTD = g2*params.kah*PT;
RD1 = params.ka*PD*N_overlap; % to loosely attachemnt state

if params.UseUniformTransitionFunc
    % the cycle goes: PT (ATP bound) <-> PD(ready) <-> P1 <-> P2 -> P3 -> PT
    % dPUdT_TransitionRates;
    sd = @(kx, alphaL, alphaR, dr,eL, eR) min(1e4, kx*(exp((alphaL*(s-dr)).^eL).*(s<dr) + exp((alphaR*(s-dr)).^eR).*(s>=dr)));
    R1D = p1.*sd(params.kd, params.alpha0, params.alpha0, params.dr0, 2, 2);
    
    R12 = p1.*sd(params.k1, params.alpha1, 0, params.dr1, 2, 2); % P1 to P2
    R21 = f1*p2.*sd(params.k_1, params.alpha_1, params.alpha_1, params.dr_1, 2, 2); % p2 to p1
    
    R2T = p2.*sd(params.k2, params.alpha2_L, params.alpha2_R, params.dr2, params.e2L, params.e2R);
else 
	strainDep = @(alpha, dr) min(1e4, exp((alpha*(s+dr)).^2));																		
    R1D = params.kd*p1.*strainDep(params.alpha0, params.dr0); %(exp(-params.alpha1*s)) + params.TK*(s>params.TK0).*s.*p1; % p1 to PU - detachment rate + tearing constant

    R12 = params.k1*p1.*exp(-params.alpha1*s); % P1 to P2
    R21 = f1*params.k_1*p2.*strainDep(params.alpha_1, params.dr_1); % p2 to p1

    kL = min(1e4, params.k2_L*((s+0)<=0).*(1 - exp(-(s+0)*params.alpha2_L)).^2);
    kR = max(0, params.k2_R*(s-params.dr2_R)).^params.alpha2_R; %.*(s>0.002);
    R2T = p2.*(params.k2 + kL + kR);
end

% to PT state directly
XB_Ripped = params.k2rip*p2.*min(1e9, max(0, exp(params.alphaRip*(s+params.dr3))));

if params.UsePassiveForSR
    F_SR = F_passive;
else
    F_SR = F_total;
end

if params.UseSuperRelaxed && params.UseDirectSRXTransition
    % A2 to SR directly: A2 => SR <-> UD
    RSR2PT = params.kmsr*exp(F_total/params.sigma1)*P_SR;    
    RPT2SR = params.ksr0*exp(-F_total/params.sigma2)*PT;
    dU_SR = -RSR2PT  + RPT2SR + sum(R2T)*dS;    
elseif params.UseSuperRelaxed
    % A2 -> UD <-> SR
    RSR2PT = params.kmsr*exp(F_total/params.sigma1)*P_SR;
    RPT2SR = params.ksr0*exp(-F_total/params.sigma2)*PT;
    dU_SR = -RSR2PT  + RPT2SR;
else 
    dU_SR = 0;
    RSR2PT = 0;
    RPT2SR = 0;    
end

if params.justPlotStateTransitionsFlag    
    plotStateTransitions;
    error('Quitting after plotting states');
end
%% governing flows
% PT - calculated as complement of sum of all probabilities
% state 0: unattached, ATP-cocked
dPD = + RTD - RD1 + sum(R1D)*dS;
dp1 = - R1D -  R12 + R21; % state 1: loosely attached, just sitting&waiting
dp2 = + R12 - R21  - R2T - XB_Ripped; % strongly attached, post-ratcheted: hydrolyzed ATP to ADP, producing Pi - ready to ratchet


% if ~params.UseMutualPairingAttachment && params.UseSpaceInterpolation
    dp1(s_i0) = dp1(s_i0) + s_i0k*(RD1/dS); % attachment
    dp1(s_i1) = dp1(s_i1) + (1-s_i0k)*(RD1/dS); % attachment

% if ~params.UseCa
    dNP = 0;

f = [dp1; dp2; dU_SR; dNP; dSL;dLSEdt;dPD];
outputs = [Force, F_active, F_passive, N_overlap, R2T', p1_0, p2_0, p1_1, p2_1, PT];
rates = [RTD, RD1, sum([R1D, R12,R21,XB_Ripped], 1)*dS, RSR2PT, RPT2SR];

%% breakpints
if t >= 1.20162642e+00 % || t > 0 && (p1_0 + p2_0 + PD + P_SR) > 1
    numberofthebeast = 667;
end

% disp('oj')
