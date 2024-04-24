function [f, outputs] = dPUdTCaSimpleAlternative2State(t,PU,params)
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
    % if t < 2.77
    %     dU_NSR = (0.5 - U_NSR)*1e3;
    %     % U_NSR = 1;
    % else
    %     dU_NSR = (0.001 - U_NSR)*1e3;
    %     % U_NSR = 0.001;    
    % end
    
else
    P_SR = 1;
end


% if ~params.UseCa
    NP = 0;

% if ~params.UseSLInput
    SL = PU(2*ss + 3);
    dSL = vel;

    
% P_SRs = 1 - P_SR;
LSE = PU(2*ss + 4); % length of the serial stiffness
PuR = PU(2*ss + 5);


% Sarcomere geometry
if params.UseOverlap
    L_thick = 1.67; % Length of thick filament, um
    L_hbare = 0.10; % Length of bare region of thick filament, um
    L_thin  = 1.20; % Length of thin filament, um
    % deltaR  = 0.010; % um    
    L_T_HS1 = min(L_thick*0.5, SL*0.5);
    L_T_HS2 = max((SL-LSE)*0.5 - ((SL-LSE)-L_thin),L_hbare*0.5);
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
    % if params.UseSpaceInterpolation
        s_i0 = floor(s_p0);
        s_i1 = ceil(s_p0);
        s_i0k = s_i1 - s_p0;
end

% sum of all probabilities
p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum((s+params.dr).*p2);
% p2_1_overstroke = dS*sum((s+params.dr).*p2.*(s>0));
% p2_1_understroke = dS*sum((s+params.dr).*p2.*(s<0));

% non-hydrolized ATP in non-super relaxed state
PuATP = max(0, N_overlap*(1.0 - NP) - (p1_0 + p2_0 + PuR + P_SR)); % unattached permissive fraction - 
% PuATP_NSR = PuATP*P_SR; 

% if ~params.F_act_UseP31
    % the dr shift is used here instead of UseP31Shift
    % p3_1_stroke = dS*sum((s+params.dr).*p3);   
    % F_active = params.kstiff2*p3_0*params.dr + params.kstiff1*( p2_1 + p3_1_stroke);     
    % F_active = params.kstiff2*p3_0*params.dr + params.kstiff1*(p3_1_stroke);     
    F_active = params.kstiff2*(p2_1) + params.kstiff1*(p1_1);     
    % F_active = params.kstiff2*(p3_1 + p2_1_understroke) + params.kstiff3*(p2_1_overstroke) + params.kstiff1*(p1_1);
% if t > 2.77
%     F_active = 0;
% end

if params.UsePassive
    Lsc0    = 1.51;
    % gamma = 7.5;
    F_passive = params.k_pas*max(SL-LSE - Lsc0, 0)^params.gamma; 
else
    F_passive = 0;
end

F_total = F_active + F_passive;

% if params.UseSerialStiffness && ~params.UseSlack    
    % if LSE >= 0
    %     Force = params.kSE*LSE;
    % else
    %     % Slack, no force
    %     Force = 0;
    % end
    Force = max(-0.1, params.kSE*LSE);
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

g1 = 1; g2 = 1; f1 = 0; f2 = 1;

% the cycle goes: PuATP (ATP bound) <-> PuR(ready) <-> P1 <-> P2 -> P3 -> PuATP
strainDep = @(alpha, dr) exp((alpha*(s+dr)).^2);
PuATP2PuD = g2*params.kah*PuATP;
PuATP2PuD_r = params.kadh*PuR;
PU2p1 = params.ka*PuR; % to loosely attachemnt state
PU2p1_r = params.kd*p1.*strainDep(params.alpha0, params.dr0); %(exp(-params.alpha1*s)) + params.TK*(s>params.TK0).*s.*p1; % p1 to PU - detachment rate + tearing constant
% p12p2 = params.k1*p1.*strainDep(params.alpha1, params.dr1); % P1 to P2
p12p2 = params.k1*p1.*exp(-params.alpha1*s); % P1 to P2
p12p2_r = f1*params.k_1*p2.*strainDep(params.alpha_1, params.dr_1); % backward flow from p2 to p1

% if params.UseTORNegShift
    % creates instable oscillations without suppressing the force enough.
    % XB_TOR = g2*params.k3*(exp(params.alpha3*(s-params.s3).^2).*p3) + p3.*min((s>params.TK0).*s*[params.TK], params.TK);
    % XB_TOR = g2*params.k3*(exp(abs(params.alpha3*(s))).*p3) + p3.*min((s>params.TK0).*s*[params.TK], params.TK);
    XB_TOR = g2*params.k2*p2.*min(1e9, max(1, strainDep(params.alpha2, params.dr2))) + p2.*min((s>params.TK0).*s*[params.TK], params.TK);

    % XB_Ripped = params.k2rip*p2.*min(1e9, max(1, strainDep(params.alpha2, params.dr2)));

% if ~params.UseAtpOnUNR
    % dU_NSR = params.ksr0*(exp(F_total/params.sigma0))*U_SR - params.kmsr*U_NSR*PuATP;
    % dU_NSR = params.ksr0*(F_total/params.sigma0)*U_SR - params.kmsr*U_NSR*PuATP;
    % dU_NSR = 10*(F_total/(F_total + 5))*U_SR - 100*U_NSR*PuATP;
    % dU_NSR = params.ksr0*(F_total/(F_total + params.sigma0))*U_SR - 100*U_NSR*PuATP*(exp(-F_total/params.sigma2));
    % dU_NSR = params.ksr0*(exp(F_total/params.sigma1))*U_SR - params.kmsr*U_NSR*PuATP*(exp(-F_total/params.sigma2));
% model 01
    % dU_NSR = params.ksr0*(F_total/params.sigma1)*U_SR - params.kmsr*U_NSR*PuATP*(exp(-F_total/params.sigma2));
    %model 11
    if params.UseSuperRelaxed
        dU_SR = - params.ksr0*exp(F_total/params.sigma1)*P_SR + params.kmsr*PuATP*(exp(-F_total/params.sigma2));
    else 
        dU_SR = 0;
    end

%% governing flows
% state 0: unattached, ATP-cocked
dPuR = PuATP2PuD - PuATP2PuD_r - PU2p1 + sum(PU2p1_r)*dS;
dp1   = - PU2p1_r -  p12p2 + p12p2_r; % state 1: loosely attached, just sitting&waiting
dp2   = + p12p2 - p12p2_r  - XB_TOR; % strongly attached, post-ratcheted: hydrolyzed ATP to ADP, producing Pi - ready to ratchet

% if ~params.UseMutualPairingAttachment && params.UseSpaceInterpolation
    dp1(s_i0) = dp1(s_i0) + s_i0k*(PU2p1/dS); % attachment
    dp1(s_i1) = dp1(s_i1) + (1-s_i0k)*(PU2p1/dS); % attachment

% if ~params.UseCa
    dNP = 0;

f = [dp1; dp2; dU_SR; dNP; dSL;dLSEdt;dPuR];
outputs = [Force, F_active, F_passive, N_overlap, XB_TOR', p1_0, p2_0, p1_1, p2_1, PuATP];
%% breakpints
if t > 2.91
    numberofthebeast = 666;
end
% disp('oj')
