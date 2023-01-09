function [f, outputs] = dPUdTCa(t,PU,params)
% ODE function for the d/dt operator for the cross-bridge mode.
%  first 2N-1 entries of PU represent p1(s,t)
%  second 2N-1 entries represent p2(s,t)
%  third 2N-1 entries represent p3(s,t)
%  third-last entry is U_NR, the fraction of myosin heads in not-relaxed state
% then second last NP and last SL 

MgATP = params.MgATP;
Pi = params.Pi;
MgADP = params.MgADP;
vel = params.Vums;
Ca_i = params.Ca;
ss = params.ss; % space size (length of the s for each of p1-p3)



freq = 1;
T = 1/freq;

% Decompose State Variables from PU vector
p1 = PU(1:ss);
p2 = PU(ss+1:2*ss);
p3 = PU(2*ss +1:3*ss);
U_NR = PU(3*ss+1);
if params.UseCa
    NP = PU(3*ss + 2);
else
    NP = 0;
end
if ~params.UseSLInput
    SL = PU(3*ss + 3);
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
    
U_SR = 1 - U_NR;
LSE = PU(3*ss + 4);


% Sarcomere geometry
if params.UseOverlap
    L_thick = 1.67; % Length of thick filament, um
    L_hbare = 0.10; % Length of bare region of thick filament, um
    L_thin  = 1.20; % Length of thin filament, um
    deltaR  = 0.010; % um    
    L_T_HS1 = min(L_thick*0.5, SL*0.5);
    L_T_HS2 = max(SL*0.5 - (SL-L_thin),L_hbare*0.5);
    L_ov = L_T_HS1 - L_T_HS2; % Length of single overlap region
    N_overlap = L_ov*2/(L_thick - L_hbare);
else
    N_overlap = 1;
end

% calculation of moments of strain distributions
s = params.s';%(-N:1:0)'*dS;
dS = params.dS;
p1_0 = dS*sum(p1);% p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
p3_0 = dS*sum(p3); p3_1 = dS*sum((s+params.dr).*p3);

Pu = N_overlap*(1.0 - NP) - (p1_0 + p2_0 + p3_0); % unattached permissive fraction

% quasi-equilibrium binding factor functions
% TODO move to evalModel for optim
% g1 = (MgADP/params.K_D)/(1 + MgADP/params.K_D + MgATP/params.K_T1);
% g2 = (MgATP/params.K_T1)/(1 + MgADP/params.K_D + MgATP/params.K_T1);
% % g3 = MgATP/(MgATP + K_T2);
% g4 = MgATP/(MgATP + params.K_T3);
% f1 = (Pi/params.K_Pi)/(1 + Pi/params.K_Pi); f2 = 1/(1 + Pi/params.K_Pi); 

g1 = 0;
g2 = 1;%0.95;
g4 = 1;%0.604;
f1 = 0;f2 = 1;
% Force model
kstiff1 = params.kstiff1; 
kstiff2 = params.kstiff2;
% F_active = kstiff2*p3_0/100 - max(-kstiff1*(p2_1 + p3_1 ), 0);
F_active = kstiff1*p2_1 + kstiff2*p3_1;

if params.UsePassive
    Lsc0    = 1.51;
    gamma = 7.5;
    F_passive = params.k_pas*max(SL - Lsc0, 0)^gamma; 
else
    F_passive = 0;
end

F_total = F_active + F_passive;

% we do nont know the velocity here, so we do that up a level
% Force = kstiff2*p3_0 + kstiff1*(( p2_1 + p3_1 )^g(20)) + mu*vel;
% muscle model
if params.UseSerialStiffness
    
    if ~params.UseSlack
        Force = params.kSE*LSE;
    elseif LSE >= 0
        Force = params.kSE*LSE;
    else
        % soft slack spring
%         Force = LSE/params.kSE;
        Force = 0;
    end
        
    velHS = (Force - F_total)/params.mu;
    dLSEdt = vel - velHS;
elseif params.UseSlack
    vmax = params.vmax;
    if vel < -vmax
        % slacking - lengthtening
        Force = F_passive;
        velHS = -vmax;
        dLSEdt = vel + vmax;
    elseif vel > -vmax && LSE < 0
        % slacking / shortening
        Force = F_passive;
        velHS = -vmax;
        dLSEdt = vel + vmax;
    else % vel > -vmax && LSE >= 0
        Force = F_total;
        velHS = vel;
        dLSEdt = 0;
    end
        
        
    
%     if LSE > 0
%     Force = params.kSE*LSE;
%     velHS = (Force - F_total)/params.mu;
%     dLSEdt = vel - velHS;
%     else
%     Force = params.kSE/1000*LSE;
%     velHS = (Force - F_total)/params.mu;
%     dLSEdt = vel - velHS;
%     end
    
%     if F_active > 0 && LSE > -1e-3
%         % normal
%         Force = F_total;
%         velHS = vel;
%         % return to norm
%         dLSEdt = -LSE*1000;
%     else
%         % slack - use a compliant spring
%         Force = params.kSE/10*LSE;
%         velHS = (Force - F_total)/params.mu;
%         dLSEdt = vel - velHS;
%     end
else
    % like 10x faster, does not cause any oscillations
    Force = F_total;
    velHS = vel;
    dLSEdt = 0;
end
    

% Estimating space derivatives, upwind differencing
if velHS > 0
    dp1ds = [(p1(1) - 0); p1(2:end) - p1(1:end-1)]/dS;
    dp2ds = [(p2(1) - 0); p2(2:end) - p2(1:end-1)]/dS;
    dp3ds = [(p3(1) - 0); p3(2:end) - p3(1:end-1)]/dS;
elseif velHS < 0
    dp1ds = [p1(2:end) - p1(1:end-1); (0 - p1(end))]/dS;
    dp2ds = [p2(2:end) - p2(1:end-1); (0 - p2(end))]/dS;    
    dp3ds = [p3(2:end) - p3(1:end-1); (0 - p3(end))]/dS;
else 
    % just optim, because its multiplied by 0 anyway
    dp1ds = 0;dp2ds = 0;dp3ds = 0;
end

% dU_NR = + ksr0*U_SR - kmsr*U_NR*Pu  ; 
dU_NR = + g4*params.ksr0*exp(F_total/params.sigma0)*U_SR - params.kmsr*U_NR*Pu; 
% dU_NR = ksr0*exp(F_active/sigma0)*U_SR - kmsr*U_NR*Pu;
% dU_NR = + ksr0*exp(F_active/sigma0)*U_SR*(1 + 3*U_NR) - kmsr*U_NR*(1 + 3*U_SR)*Pu  ; 
% dU_NR = + ksr*(1/(1.0 - MgATP/10))*(exp(F_active/sigma0))*U_SR - 50*kmsr*(1.0 - g3)*U_NR*Pu  ; 
% dU_NR = + ksr0*(1 + F_active/sigma0 )*U_SR - kmsr*U_NR*Pu  ; 
% dU_NR = + ksr0*U_SR - kmsr*exp(-F_active/sigma0)*U_NR*Pu  ; 
dp1   = -velHS/2*dp1ds - f1*params.kd*p1 - f2*params.k1*(exp(-params.alpha1*s).*p1) ...
    + params.k_1*(exp(+params.alpha1*s).*p2);
dp2   = -velHS/2*dp2ds + f2*params.k1*(exp(-params.alpha1*s).*p1) ...
    - params.k_1*(exp(+params.alpha1*s).*p2) - params.k2*(exp(-params.alpha2*s).*p2) ...
    + g1*params.k_2*p3  ;

         
% XB_TOR = max(-1, g2*params.k3*(exp(params.alpha3*(s-params.s3).^2).*p3));
XB_TOR = g2*params.k3*(exp(params.alpha3*(s-params.s3).^2).*p3);
if any(XB_TOR < -1) 
    a = 1;
end
dp3   = -velHS/2*dp3ds + params.k2*(exp(-params.alpha2*s).*p2) ...
    - g1*params.k_2*p3 - XB_TOR;
% dp1(N+1) = dp1(N+1) + ka*Pu*U_NR/dS; % attachment
% dp1(N+1) = dp1(N+1) + ka*Pu*(1.0 - (p1_0 + p2_0 + p3_0))*U_NR/dS; % attachment
dp1(params.s_i0) = dp1(params.s_i0) + params.ka*Pu*(params.Amax - (p1_0 + p2_0 + p3_0))*U_NR/dS; % attachment
% dp1(params.s_i0) = dp1(params.s_i0) + ka*Pu*U_NR/dS; % attachment

if params.UseCa
    Jon  = params.k_on*Ca_i*NP*(1 + params.K_coop*(1 - NP));
    Joff = params.k_off*(Pu/N_overlap)*(1 + params.K_coop*NP);
    dNP = - Jon + Joff; % dN_LV / dt 
else
    dNP = 0;
end

% dLse = Kse*Lse

f = [dp1; dp2; dp3; dU_NR; dNP; dSL;dLSEdt];
if t > 0.1
    a = 1;
end
outputs = [Force, F_active, F_passive, N_overlap, XB_TOR'];