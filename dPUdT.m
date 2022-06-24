function f = dPUdT(~,PU,N,dS,MgATP,Pi,MgADP,g,SL)
% ODE function for the d/dt operator for the cross-bridge mode.
%  first 2N-1 entries of PU represent p1(s,t)
%  second 2N-1 entries represent p2(s,t)
%  third 2N-1 entries represent p3(s,t)
%  last entry is U_NR, the fraction of myosin heads in not-relaxed state

% State Variables
p1 = PU(1:1*N+1);
p2 = PU(1*N+2:2*N+2);
p3 = PU(2*N+3:3*N+3);
U_NR = PU(3*N+4);
U_SR = 1 - U_NR;

% calculation of moments of strain distributions
s = (-N:1:0)'*dS;
p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
Pu = 1.0 - (p1_0 + p2_0 + p3_0); % unattached permissive fraction

% definition of parameters
alpha1 = 25;
alpha2 = 25;
alpha3 = 10*g(9)*276;
s3     = -g(12)*0.05;
K_Pi = 15;
K_T1 = g(11)*0.50; % (mM) 
K_T2 = 0.1; % (mM)
K_T3 = g(15)*4.0; % (mM)
K_D = 0.194; % MgADP dissociation constant from Yamashita etal (Circ Res. 1994; 74:1027-33).
g1 = (MgADP/K_D)/(1 + MgADP/K_D + MgATP/K_T1);
g2 = (MgATP/K_T1)/(1 + MgADP/K_D + MgATP/K_T1);
g3 = (MgATP/K_T2)/(1 + MgATP/K_T2);
f1 = (Pi/K_Pi)/(1 + Pi/K_Pi); f2 = 1/(1 + Pi/K_Pi); 
ka  = 0.1*g3*g(1)*373*1 ;
kd  = f1*g(2)*103;
k1  = 10*g(3)*f2*40*5;%
k_1 = 10*g(4)*17;%
k2  = 100*g(5)*420/5;
k_2 = 100*g1*2.8; % not identified
k3  = 0.010*g2*g(10)*44;%;

% Force model
dr = 0.01; % Power-stroke Size; Units: um
kstiff1 = g(13)*1500; 
kstiff2 = g(14)*10000;
F_active = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 );
% F_total = F_active + 2;

% transitions between super relaxed state and non relaxed state
ksr    = g(6)*9; % 
sigma0 = g(7)*33;
kmsr   = g(8)*250; % 

% dU_NR = + ksr*(1/(1.0 - MgATP/10))*(exp(F_active/sigma0))*U_SR - 50*kmsr*(1.0 - g3)*U_NR*Pu  ; 
dU_NR = + ksr*exp(F_active/sigma0)*U_SR - kmsr*(K_T3/(MgATP + K_T3))*U_NR*Pu  ; 
dp1   = - kd*p1 - k1*(exp(-alpha1*s).*p1) + k_1*(exp(+alpha1*s).*p2);
dp2   = + k1*(exp(-alpha1*s).*p1) - k_1*(exp(+alpha1*s).*p2) - k2*(exp(-alpha2*s).*p2) + k_2*p3  ;
dp3   = + k2*(exp(-alpha2*s).*p2) - k_2*p3 - k3*(exp(alpha3*(s+s3).^2).*p3);
dp1(N+1) = dp1(N+1) + ka*Pu*U_NR/dS; % attachment

f = [dp1; dp2; dp3; dU_NR];
