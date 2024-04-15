function [f, outputs] = dPUdTCa(t,PU,params)
% ODE function for the d/dt operator for the cross-bridge model of half-sarcomere.
%  first 2N-1 entries of PU represent p1(s,t)
%  second 2N-1 entries represent p2(s,t)
%  third 2N-1 entries represent p3(s,t)
%  third-last entry is U_NR, the fraction of myosin heads in not-relaxed state
% then second last NP and last SL 

vel = params.Vums;
Ca_i = params.Ca;

freq = 1;
T = 1/freq;

% Decompose State Variables from PU vector
ss = params.ss; % space size (length of the s for each of p1-p3)
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
%     if t > 2.76
%         a = 1;
%     end
    vel = dSL;
end
    
U_SR = 1 - U_NR;
LSE = PU(3*ss + 4); % length of the serial stiffness


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

if isnan(s_p0)
    s_i0 = 1; % just hold on..
    s_i1 = 1;
    s_i0k = 0;
elseif floor(s_p0) < 0
    error(sprintf('Out of bounds! At %0.6fs and SL %0.2fum, the s(1) was %0.2f', t, SL, s(1)));
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

if params.UseSpaceDiscretization
    try
    s = s - s(s_i0); % subtrck the space difference
    catch e
        disp('wat')
    end
    if abs(s(s_i0)) > 1e-4
        disp('')
    end
end


% calculation of moments of strain distributions
try
% sum of all probabilities
p1_0 = dS*sum(p1);% p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
p3_0 = dS*sum(p3); 
catch e
    disp(e.message);
end

if params.UseP31Shift
    % mean of p3 probability with shifted centre by the contraction dr
    if params.UseKstiff3
        p3_1_stroke = dS*sum((s+params.dr).*p3.*(s+params.dr >= 0));
        p3_1_overstroke = dS*sum((s+params.dr).*p3.*(s+params.dr < 0));
    else
        p3_1_stroke = dS*sum((s+params.dr).*p3);   
        p3_1_overstroke = 0;
    end
else
    p3_1_stroke = dS*sum(s.*p3);
end

Pu = N_overlap*(1.0 - NP) - (p1_0 + p2_0 + p3_0); % unattached permissive fraction

% quasi-equilibrium binding factor functions
% TODO move to evalModel for optim
MgATP = params.MgATP;
Pi = params.Pi;
MgADP = params.MgADP;

% g1 = (MgADP/params.K_D)/(1 + MgADP/params.K_D + MgATP/params.K_T1);
% removed + 1 in the denominator, because there is no unbound fraction
g1 = (MgADP/params.K_D)/(MgADP/params.K_D + MgATP/params.K_T1);
% g2 = (MgATP/params.K_T1)/(1 + MgADP/params.K_D + MgATP/params.K_T1);
g2 = (MgATP/params.K_T1)/(MgADP/params.K_D + MgATP/params.K_T1);

% g3 = MgATP/(MgATP + K_T2);
g4 = MgATP/(MgATP + params.K_T3);
f1 = (Pi/params.K_Pi)/(1 + Pi/params.K_Pi); f2 = 1/(1 + Pi/params.K_Pi); 

% g1 = 0;
% g2 = 1;%0.95;
% g4 = 1;%0.604;
% f1 = 0;f2 = 1;
% Force model

% F_active = kstiff2*p3_0/100 - max(-kstiff1*(p2_1 + p3_1 ), 0);
if params.F_act_UseP31
    % together with UseP31Shift
    F_active = params.kstiff1*p2_1 + params.kstiff2*p3_1_stroke + params.kstiff3*p3_1_overstroke;
else
    % the dr shift is used here instead of UseP31Shift
    F_active = params.kstiff2*p3_0*params.dr + params.kstiff1*( p2_1 + p3_1_stroke); 
end

if params.UsePassive
    Lsc0    = 1.51;
    gamma = 7.5;
    F_passive = params.k_pas*max(SL - Lsc0, 0)^gamma; 
else
    F_passive = 0;
end

F_total = F_active + F_passive;

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
else
    % like 10x faster, does not cause any oscillations
    Force = F_total;
    velHS = vel;
    dLSEdt = 0;
end

% dU_NR = + ksr0*U_SR - kmsr*U_NR*Pu  ; 
if params.UseAtpOnUNR
    dU_NR = + g4*params.ksr0*exp(F_total/params.sigma0)*U_SR - params.kmsr*U_NR*Pu; 
else
    dU_NR = params.ksr0*(exp(F_active/params.sigma0))*U_SR - params.kmsr*U_NR*Pu;
end
% dU_NR = ksr0*exp(F_active/sigma0)*U_SR - kmsr*U_NR*Pu;
% dU_NR = + ksr0*exp(F_active/sigma0)*U_SR*(1 + 3*U_NR) - kmsr*U_NR*(1 + 3*U_SR)*Pu  ; 
% dU_NR = + ksr*(1/(1.0 - MgATP/10))*(exp(F_active/sigma0))*U_SR - 50*kmsr*(1.0 - g3)*U_NR*Pu  ; 
% dU_NR = + ksr0*(1 + F_active/sigma0 )*U_SR - kmsr*U_NR*Pu  ; 
% dU_NR = + ksr0*U_SR - kmsr*exp(-F_active/sigma0)*U_NR*Pu  ; 

%% TRANSITIONS
if exist('plotTransitions', 'var')
    % debug the transitions
    s = -0.2:0.02:0.2; p1 = ones(size(s)); p2 = p1;p3 = p1;
end

p12PU = f1*params.kd*p1; % p1 to PU
p12p2 = f2*params.k1*(exp(-params.alpha1*s).*p1); % P1 to P2
p12p2_r = params.k_1*(exp(+params.alpha1*s).*p2); % backward flow from p2 to p1
p22p3 = params.k2*(exp(-params.alpha2*s).*p2); % P2 to P3
p22p3_r = g1*params.k_2*p3; % reverse flow from p3 to p2


% XB_TOR = max(-1, g2*params.k3*(exp(params.alpha3*(s-params.s3).^2).*p3));
if params.UseTORNegShift
    % creates instable oscillations without suppressing the force enough.
    % params.TK ~ 1e3
    % XB_TOR = g2*params.k3*(exp(params.alpha3*(s-params.s3).^2).*p3) + p3.*(s>0).*min(exp(s*params.TK)/(1+params.TK), params.TK);

    XB_TOR = g2*params.k3*(exp(params.alpha3*(s-params.s3).^2).*p3) + p3.*min((s>params.TK0).*s*[params.TK], params.TK);
else
%     % test using interp
%     try
%     p3i = interp1(s, p3, s+params.s3, 'linear');
%     p3i(isnan(p3i)) = 0;
    XB_TOR = g2*params.k3*(exp(params.alpha3*(s+params.s3).^2).*p3);
%     catch e
%         disp('vile');
%     end
end

if exist('plotTransitions', 'var')
    % debug transitions only
    figure(87);
    plot(s, p12PU, s, p12p2, s, p12p2_r,'--', s, p22p3, s,  p22p3_r, '--', s, XB_TOR, 'Linewidth', 2);
    legend(...
    'p12PU',... = f1*params.kd*p1; % p1 to PU
    'p12p2',... = f2*params.k1*(exp(-params.alpha1*s).*p1); % P1 to P2
    'p12p2_r',..._r = params.k_1*(exp(+params.alpha1*s).*p2); % backward flow from p2 to p1
    'p22p3',... = params.k2*(exp(-params.alpha2*s).*p2); % P2 to P3
    'p22p3_r',... = g1*params.k_2*p3; % reverse flow from p3 to p2
    'TOR'...
    );
% ylim([0, 5000]);
end
%%
% governing flows
dp1   = - p12PU -  p12p2 + p12p2_r;
dp2   = + p12p2 - p12p2_r  - p22p3 +  p22p3_r;
dp3   = + p22p3 - p22p3_r - XB_TOR;

try
if params.UseMutualPairingAttachment
    dp1(s_i0) = dp1(s_i0) + params.ka*Pu*(params.Amax - (p1_0 + p2_0 + p3_0))*U_NR/dS; % attachment
else
    if params.UseSpaceInterpolation
        dp1(s_i0) = dp1(s_i0) + s_i0k*(params.ka*Pu*U_NR/dS); % attachment
        dp1(s_i1) = dp1(s_i1) + (1-s_i0k)*(params.ka*Pu*U_NR/dS); % attachment
    else
        dp1(s_i0) = dp1(s_i0) + (params.ka*Pu*U_NR/dS); % attachment
    end
end
catch e
    disp('vole')
end

if params.UseCa
    Jon  = params.k_on*Ca_i*NP*(1 + params.K_coop*(1 - NP));
    Joff = params.k_off*(Pu/N_overlap)*(1 + params.K_coop*NP);
    dNP = - Jon + Joff; % dN_LV / dt 
else
    dNP = 0;
end

% dLse = Kse*Lse

f = [dp1; dp2; dp3; dU_NR; dNP; dSL;dLSEdt];
if t  > 0.18 | any(isnan(PU)) | any(isnan(f))
    % just for placing a breakpoint here
    a = 1;
end
if SL > 2.01 & t > 0.05
    a = 3;
end
outputs = [Force, F_active, F_passive, N_overlap, XB_TOR'];