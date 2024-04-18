function [f, outputs] = dPUdTCa(t,PU,params)
% ODE function for the d/dt operator for the cross-bridge model of half-sarcomere.
%  first 2N-1 entries of PU represent p1(s,t)
%  second 2N-1 entries represent p2(s,t)
%  third 2N-1 entries represent p3(s,t)
%  third-last entry is U_NR, the fraction of myosin heads in not-relaxed state
% then second last NP and last SL 

vel = params.Vums;

% Decompose State Variables from PU vector
ss = params.ss; % space size (length of the s for each of p1-p3)
p1 = PU(1:ss);
p2 = PU(ss+1:2*ss);
p3 = PU(2*ss +1:3*ss);
if params.UseSuperRelaxed
    U_NSR = PU(3*ss+1);
else
    U_NSR = 1;
end
% if ~params.UseCa
    NP = 0;

% if ~params.UseSLInput
    SL = PU(3*ss + 3);
    dSL = vel;

    
U_SR = 1 - U_NSR;
LSE = PU(3*ss + 4); % length of the serial stiffness
PuR = PU(3*ss + 5);


% Sarcomere geometry
if params.UseOverlap
    L_thick = 1.67; % Length of thick filament, um
    L_hbare = 0.10; % Length of bare region of thick filament, um
    L_thin  = 1.20; % Length of thin filament, um
    % deltaR  = 0.010; % um    
    L_T_HS1 = min(L_thick*0.5, SL*0.5);
    L_T_HS2 = max(SL*0.5 - (SL-L_thin),L_hbare*0.5);
    L_ov = L_T_HS1 - L_T_HS2; % Length of single overlap region
    N_overlap = L_ov*2/(L_thick - L_hbare);
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
    s_i1 = 1;
    s_i0k = 0;
    % warning(sprintf('Out of bounds! At %0.6fs and SL %0.2fum, the s(1) was %0.2f', t, SL, s(1)));
elseif floor(s_p0) == 0
    s_i0 = 1;
    s_i1 = 1;
    s_i0k = 0;
elseif ceil(s_p0) == params.ss
    s_i0 = params.ss;
    s_i1 = 1;
    s_i0k = 0;
elseif ceil(s_p0) > params.ss
    error(sprintf('Out of bounds! At %0.6fs and SL %0.2fum, the s(end) was %0.2f', t, SL, s(end)));
else
    if params.UseSpaceInterpolation
        s_i0 = floor(s_p0);
        s_i1 = ceil(s_p0);
        s_i0k = s_i1 - s_p0;
    else
        s_i0 = round(s_p0);
        s_i1 = 0;
    end
end

% sum of all probabilities
p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum((s+params.dr).*p2);
p3_0 = dS*sum(p3); p3_1 = dS*sum((s+params.dr).*p3);

% non-hydrolized ATP in non-super relaxed state
PuATP = max(0, N_overlap*(1.0 - NP) - (p1_0 + p2_0 + p3_0 + PuR)); % unattached permissive fraction - 
PuATP_NSR = PuATP*U_NSR; 

% if ~params.F_act_UseP31
    % the dr shift is used here instead of UseP31Shift
    % p3_1_stroke = dS*sum((s+params.dr).*p3);   
    % F_active = params.kstiff2*p3_0*params.dr + params.kstiff1*( p2_1 + p3_1_stroke);     
    % F_active = params.kstiff2*p3_0*params.dr + params.kstiff1*(p3_1_stroke);     
    F_active = params.kstiff2*(p3_1 + p2_1) + params.kstiff1*(p1_1);     


if params.UsePassive
    Lsc0    = 1.51;
    gamma = 7.5;
    F_passive = params.k_pas*max(SL - Lsc0, 0)^gamma; 
else
    F_passive = 0;
end

F_total = F_active + F_passive;

% if params.UseSerialStiffness && ~params.UseSlack    
    if LSE >= 0
        Force = params.kSE*LSE;
    else
        % Slack, no force
        Force = 0;
    end
    velHS = (Force - F_total)/params.mu;
    dLSEdt = vel - velHS;

%% TRANSITIONS

% quasi-equilibrium binding factor functions
% TODO move to evalModel for optim
MgATP = params.MgATP;
Pi = params.Pi;
MgADP = params.MgADP;

% g1 = (MgADP/params.K_D)/(MgADP/params.K_D + MgATP/params.K_T1);
% g2 = (MgATP/params.K_T1)/(MgADP/params.K_D + MgATP/params.K_T1);
% % g4 = MgATP/(MgATP + params.K_T3);
% f1 = (Pi/params.K_Pi)/(1 + Pi/params.K_Pi); f2 = 1/(1 + Pi/params.K_Pi); 

g1 = 0; g2 = 1; f1 = 0; f2 = 1;

% the cycle goes: PuATP (ATP bound) <-> PuR(ready) <-> P1 <-> P2 -> P3 -> PuATP

PuATP2PuR = g2*params.kah*PuATP_NSR;
PuATP2PuR_r = params.kadh*PuR;
PU2p1 = params.ka*PuR; % to loosely attachemnt state
PU2p1_r = params.kd*p1.*(exp(-params.alpha1*s)) + params.TK*(s>params.TK0).*s.*p1; % p1 to PU - detachment rate + tearing constant
p12p2 = params.k1*p1; % P1 to P2
p12p2_r = f1*params.k_1*(exp(+params.alpha1*s).*p2); % backward flow from p2 to p1
p22p3 = f2*params.k2*p2.*min(params.alpha2, (exp(abs(params.alpha2*(s+params.dr).^params.alpha3)))); % P2 to P3
p22p3_r = g1*params.k_2*p3; % reverse flow from p3 to p2

% if params.UseTORNegShift
    % creates instable oscillations without suppressing the force enough.
    % XB_TOR = g2*params.k3*(exp(params.alpha3*(s-params.s3).^2).*p3) + p3.*min((s>params.TK0).*s*[params.TK], params.TK);
    % XB_TOR = g2*params.k3*(exp(abs(params.alpha3*(s))).*p3) + p3.*min((s>params.TK0).*s*[params.TK], params.TK);
    XB_TOR = g2*params.k3*p3 + p3.*min((s>params.TK0).*s*[params.TK], params.TK);

% if ~params.UseAtpOnUNR
    dU_NSR = params.ksr0*(exp(F_active/params.sigma0))*U_SR - params.kmsr*U_NSR*PuATP;

%% governing flows
% state 0: unattached, ATP-cocked
dPuR = PuATP2PuR - PuATP2PuR_r - PU2p1 + sum(PU2p1_r)*dS;
dp1   = - PU2p1_r -  p12p2 + p12p2_r; % state 1: loosely attached, just sitting&waiting
dp2   = + p12p2 - p12p2_r  - p22p3 +  p22p3_r; % strongly attached: hydrolyzed ATP to ADP, producing Pi - ready to ratchet
dp3   = + p22p3 - p22p3_r - XB_TOR; % post-ratcheted: ADP bound, still attached

% if ~params.UseMutualPairingAttachment && params.UseSpaceInterpolation
    dp1(s_i0) = dp1(s_i0) + s_i0k*(PU2p1/dS); % attachment
    dp1(s_i1) = dp1(s_i1) + (1-s_i0k)*(PU2p1/dS); % attachment

% if ~params.UseCa
    dNP = 0;

f = [dp1; dp2; dp3; dU_NSR; dNP; dSL;dLSEdt;dPuR];
outputs = [Force, F_active, F_passive, N_overlap, XB_TOR', p1_0, p2_0, p3_0, p2_1, p3_1, PuATP];
%% breakpints
if t > 2.765
    numberofthebeast = 666;
end
