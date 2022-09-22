function f = dPUdTCa(~,PU,params, g)
% ODE function for the d/dt operator for the cross-bridge mode.
%  first 2N-1 entries of PU represent p1(s,t)
%  second 2N-1 entries represent p2(s,t)
%  third 2N-1 entries represent p3(s,t)
%  third-last entry is U_NR, the fraction of myosin heads in not-relaxed state
% then second last NP and last SL 

MgATP = params.MgATP;
Pi = params.Pi;
MgADP = params.MgADP;
vel = params.Velocity;
Ca_i = params.Ca;
N = params.N;


% Dan's parameters:
% from cross bridge model identrification
g0 = [ 1.5*0.3977    2.0478    1.4903    0.3765    0.5219    0.2726    1.25  1.0471    0.2382    0.9342];

K_coop = 5.7;
k_on   = g0(1)*100;
k_off  = g0(2)*1.5*100;

freq = 1;
T = 1/freq;

% State Variables
p1 = PU(1:1*N+1);
p2 = PU(1*N+2:2*N+2);
p3 = PU(2*N+3:3*N+3);
U_NR = PU(3*N+4);
if isfield(params, 'UseCa') && params.UseCa
    NP = PU(3*N + 5);
else
    NP = 0;
end

SL = PU(3*N + 6);

U_SR = 1 - U_NR;


% Sarcomere geometry
if isfield(params, 'UseOverlap') && params.UseOverlap
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

dr = +g(12)*0.01; % Power-stroke Size; Units: um
% calculation of moments of strain distributions
s = params.s';%(-N:1:0)'*dS;
dS = params.dS;
p1_0 = dS*sum(p1);% p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
Pu = N_overlap*(1.0 - NP) - (p1_0 + p2_0 + p3_0); % unattached permissive fraction

% strain-associated parameters
alpha1 = g(16)*50;
alpha2 = g(17)*50;
alpha3 = g(9)*5000;
s3     = dr;

% dissociation constants
K_Pi = 15;
K_T1 = g(11)*1; % (mM) ATP binding for detachment
% K_T2 = 0.05; % (mM) ATP binding to P state
K_T3 = g(15)*4; % (mM)
K_D = 0.194; % MgADP dissociation constant from Yamashita etal (Circ Res. 1994; 74:1027-33).

% quasi-equilibrium binding factors
g1 = (MgADP/K_D)/(1 + MgADP/K_D + MgATP/K_T1);
g2 = (MgATP/K_T1)/(1 + MgADP/K_D + MgATP/K_T1);
% g3 = MgATP/(MgATP + K_T2);
g4 = MgATP/(MgATP + K_T3);
f1 = (Pi/K_Pi)/(1 + Pi/K_Pi); f2 = 1/(1 + Pi/K_Pi); 

% rate constants
ka  = g(1)*100 ;
kd  = f1*g(2)*5; 
k1  = g(3)*f2*500;%
k_1 = g(4)*50;%
k2  = g(5)*500;
k_2 = g1*10; % not identified
k3  = g2*g(10)*25;%;


% Force model
kstiff1 = g(13)*2500; 
kstiff2 = g(14)*200;
F_active = kstiff2*p3_0 - max(-kstiff1*(p2_1 + p3_1 ), 0).^g(20);

% we do nont know the velocity here, so we do that up a level
% Force = kstiff2*p3_0 + kstiff1*(( p2_1 + p3_1 )^g(20)) + mu*vel;

% transitions between super relaxed state and non relaxed state
ksr0   = g(6)*10 ; % 
sigma0 = g(7)*20;
kmsr   = g(8)*10; % 
% kmsr   = g(8)*250*(1-g3); % 

% phi = mod(t,T)/T;
% Ts = (0.36/0.3197)*(T^2.2)/(0.39^2.2 + T^2.2); % ratio of relaxation time to 1-Hz relaxation time
% Ca0 = 0.050; % micro M
% a = 0.558*(1 + g0(6)*(freq-0.5)); % micro M
% % a = 0.550*(1 + (freq-0.5)/(freq) ); % micro M
% b = 7.6516/Ts;
% c = 0.2893*Ts;
% d = 162.07/Ts;
% Ca_i = a*( 0.5*(1-tanh(b*(T*phi-c))) - exp(-(d*T*phi))  ) + Ca0;


Amax = g(18)*0.5;
% dU_NR = + ksr0*U_SR - kmsr*U_NR*Pu  ; 
% dU_NR = + ksr0*g4*exp(F_active/sigma0)*U_SR - kmsr*U_NR*Pu; 
dU_NR = ksr0*exp(F_active/sigma0)*U_SR - kmsr*U_NR*Pu;
% dU_NR = + ksr0*exp(F_active/sigma0)*U_SR*(1 + 3*U_NR) - kmsr*U_NR*(1 + 3*U_SR)*Pu  ; 
% dU_NR = + ksr*(1/(1.0 - MgATP/10))*(exp(F_active/sigma0))*U_SR - 50*kmsr*(1.0 - g3)*U_NR*Pu  ; 
% dU_NR = + ksr0*(1 + F_active/sigma0 )*U_SR - kmsr*U_NR*Pu  ; 
% dU_NR = + ksr0*U_SR - kmsr*exp(-F_active/sigma0)*U_NR*Pu  ; 
dp1   = - kd*p1 - k1*(exp(-alpha1*s).*p1) + k_1*(exp(+alpha1*s).*p2);
dp2   = + k1*(exp(-alpha1*s).*p1) - k_1*(exp(+alpha1*s).*p2) - k2*(exp(-alpha2*s).*p2) + k_2*p3  ;
dp3   = + k2*(exp(-alpha2*s).*p2) - k_2*p3 - k3*(exp(alpha3*(s+s3).^2).*p3);
% dp1(N+1) = dp1(N+1) + ka*Pu*U_NR/dS; % attachment
% dp1(N+1) = dp1(N+1) + ka*Pu*(1.0 - (p1_0 + p2_0 + p3_0))*U_NR/dS; % attachment
dp1(params.s_i0) = dp1(params.s_i0) + ka*Pu*(Amax - (p1_0 + p2_0 + p3_0))*U_NR/dS; % attachment
% dp1(params.s_i0) = dp1(params.s_i0) + ka*Pu*U_NR/dS; % attachment

Jon  = k_on*Ca_i*NP*(1 + K_coop*(1 - NP));
Joff = k_off*(Pu/N_overlap)*(1 + K_coop*NP);
dNP = - Jon + Joff; % dN_LV / dt
dSL = params.Vums;

f = [dp1; dp2; dp3; dU_NR; dNP; dSL];

if ~params.PlotProbsOnStep
    return;
end

% figure(params.PlotProbsOnStep);hold on;
% cla;hold on;
plot(s,p1,s,p2,s,p3,'x-', 'linewidth',1.5);
ylabel('Probability density ($\mu$m$^{-1}$)','interpreter','latex','fontsize',16);
xlabel('strain, $s$ ($\mu$m)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14);
set(gca,'xlim',[s(1) s(N)]);
legend('$p_1(s)$','$p_2(s)$','$p_3(s)$','interpreter','latex','fontsize',16,'location','northwest');
drawnow;