function f = dPUdT(~,PU,N,dS,MgATP,Pi,MgADP,g,vel)
% ODE function for the d/dt operator for the cross-bridge mode.
%  first 2N-1 entries of PU represent p1(s,t)
%  second 2N-1 entries represent p2(s,t)
%  third 2N-1 entries represent p3(s,t)
%  last entry is U_NR, the fraction of myosin heads in not-relaxed state

% State Variables
p1 = PU(1:2*N+1);
p2 = PU(2*N+2:4*N+2);
p3 = PU(4*N+3:6*N+3);
U_NR = PU(6*N+4);
U_SR = 1 - U_NR;
Ls = PU(6*N+5); % length of series element
Lh = PU(6*N+6); % 1/2-sarcomere length

% Sarcomere and Xbridge geometry parameters
dr = 0.01; % Power-stroke Size; Units: um
SL = 2*Lh;
L_thick = 1.67; % Length of thick filament, um
L_bare = 0.10; % Length of bare region of thick filament, um
L_thin  = 1.20; % Length of thin filament, um
LthickHS1 = min(L_thick*0.5, SL*0.5);
LthickHS2 = max(SL*0.5 - (SL-L_thin),L_bare*0.5);
L_OV = LthickHS1-LthickHS2;
OV = 2*L_OV/(L_thick - L_bare);
% OV = 1;

% calculation of moments of strain distributions
s = (-N:1:N)'*dS;
p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
Pu = OV - (p1_0 + p2_0 + p3_0); % unattached permissive fraction, assuming N = 0

% strain-associated parameters
alpha1 = g(16)*50;
alpha2 = g(17)*50;
alpha3 = g(9)*10000;
s3     = 0.0025;

% muscle model
mu = 0.1;
kSE = 25000; 

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
ka  = g(1)*50 ;
kd  = f1*g(2)*5; 
k1  = g(3)*f2*1000;%
k_1 = g(4)*10;%
k2  = g(5)*1000;
k_2 = g1*10; % not identified
k3  = g2*g(10)*100;%;

% Force model
kstiff2 = g(13)*5000; 
kstiff3 = g(14)*27500;
F_active = kstiff3*p3_1 + kstiff2*p2_1;

% muscle model
velHS = (kSE*Ls - F_active)/mu;% velocity of half-sarcomere
dLSEdt = vel - velHS;

% Estimating space derivatives, upwind differencing
dp1ds = zeros(N+1,1);
dp2ds = zeros(N+1,1);
dp3ds = zeros(N+1,1);
if velHS > 0
  dp1ds(2:2*N+1) = ( p1(2:2*N+1) - p1(1:2*N) )/dS;
  dp1ds(1)       = ( p1(1) - 0 )/dS;
  dp2ds(2:2*N+1) = ( p2(2:2*N+1) - p2(1:2*N) )/dS;
  dp2ds(1)       = ( p2(1) - 0 )/dS;
  dp3ds(2:2*N+1) = ( p3(2:2*N+1) - p3(1:2*N) )/dS;
  dp3ds(1)       = ( p3(1) - 0 )/dS;
else % velHS <= 0
  dp1ds(1:2*N) = ( p1(2:2*N+1) - p1(1:2*N) )/dS;
  dp1ds(2*N+1) = ( 0 - p1(2*N+1) )/dS;
  dp2ds(1:2*N) = ( p2(2:2*N+1) - p2(1:2*N) )/dS;
  dp2ds(2*N+1) = ( 0 - p2(2*N+1) )/dS;
  dp3ds(1:2*N) = ( p3(2:2*N+1) - p3(1:2*N) )/dS;
  dp3ds(2*N+1) = ( 0 - p3(2*N+1) )/dS;
end

% transitions between super relaxed state and non relaxed state
ksr0   = g(6)*10*g4 ; % 
sigma0 = g(7)*40;
kmsr   = g(8)*20; % 
% kmsr   = g(8)*250*(1-g3); % 

Amax = g(18)*1.0;
% dU_NR = + ksr0*U_SR - kmsr*U_NR*Pu  ; 
dU_NR = + ksr0*exp(F_active/sigma0)*U_SR - kmsr*U_NR*Pu  ; 
% dU_NR = + ksr0*exp(F_active/sigma0)*U_SR*(1 + 3*U_NR) - kmsr*U_NR*(1 + 3*U_SR)*Pu  ; 
% dU_NR = + ksr*(1/(1.0 - MgATP/10))*(exp(F_active/sigma0))*U_SR - 50*kmsr*(1.0 - g3)*U_NR*Pu  ; 
% dU_NR = + ksr0*(1 + F_active/sigma0 )*U_SR - kmsr*U_NR*Pu  ; 
% dU_NR = + ksr0*U_SR - kmsr*exp(-F_active/sigma0)*U_NR*Pu  ; 
dp1   = -velHS*dp1ds - kd*p1 - k1*(exp(-alpha1*s).*p1) + k_1*(exp(+alpha1*s).*p2);
dp2   = -velHS*dp2ds + k1*(exp(-alpha1*s).*p1) - k_1*(exp(+alpha1*s).*p2) - k2*(exp(-alpha2*s).*p2) + k_2*p3  ;
dp3   = -velHS*dp3ds + k2*(exp(-alpha2*s).*p2) - k_2*p3 - k3*(exp(alpha3*(s-s3).^2).*p3);
% dp1(N+1) = dp1(N+1) + ka*Pu*U_NR/dS; % attachment
% dp1(N+1) = dp1(N+1) + ka*Pu*(1.0 - (p1_0 + p2_0 + p3_0))*U_NR/dS; % attachment
dp1(N+1) = dp1(N+1) + ka*Pu*(Amax - (p1_0 + p2_0 + p3_0))*U_NR/dS; % attachment

f = [dp1; dp2; dp3; dU_NR; dLSEdt; velHS];

