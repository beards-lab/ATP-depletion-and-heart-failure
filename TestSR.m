% Script to set up and run simulation:
 
clear
 
%% Setting up problem
N = 40; % space (strain) discretization--number of grid points in half domain
Slim = 0.20;
dS = Slim/N;
s = (-N:1:0)*dS; % strain
p1 = zeros(N+1,1);
p2 = zeros(N+1,1);
pT = 0;
pD = 0;
pS = 1;
 
% State variable vector concatenates p1, p2, p2, and U_NR
PU0 = [p1; p2; pT; pD;pS];
% PU0 = [p1; p2; pT; pD];
 
 
%% Simulating steady state SL-dependence
SL = 1.8:0.05:2.2;

% SL = [1.9 2.0  2.1 2.2 2.3];
 
for i = 1:length(SL)
 
  % run model to steady state
  [t,PU] = ode15s(@dPUdT,[0 200],PU0,[],N,dS,SL(i), g);
 
  PUend = PU(end,:);
  p1 = PUend(1:1*N+1);  % probabilities of attatched state 1 on grid
  p2 = PUend(1*N+2:2*N+2); % probabilities of attatched state 2 on grid
  pT(i) = PUend(2*N+3); % probability of unattached state T
  pD(i) = PUend(2*N+4); % probability of unattached state D
  pS(i) = PUend(2*N+5); % probability of unattached state D
 
  p1_0(i) = dS*sum(p1); p1_1 = dS*sum(s.*p1(i));
  p2_0(i) = dS*sum(p2); p2_1 = dS*sum(s.*p2(i));
  pS(i) = 1 - p1_0(i) - p2_0(i) - pT(i) - pD(i); % SRX probability
  sop(i) = pS(i) + p1_0(i) + p2_0(i) + pT(i) + pD(i);
 
  dr = 0.01; % Power-stroke Size; Units: um
  kstiff1 = 1e4;
  F_active(i) = kstiff1*(p1_1 + p2_1 + dr*p2_0(i));

  kSE = 5000;
  LSE(i) = F_active(i)/kSE;
  % LSE =0;
    y = @(k_pas, x0, gamma, x) k_pas.*(x-x0).^gamma - 4*0 - x0*0 + 0.5e9.*(x-0.95).^13;    
    F_passive(i) = y(0.4, -0.4, 7.9, (SL(i)-LSE(i))/2); % from the ramp-ups

    % Overlap function
    OV(i) = min(1, 0.6 + (SL(i)-1.8)); % this is just a function that I made up    
    OF(i)     = (OV(i) - (p1_0 + p2_0))/(pS + pT + pD); % overlap factor
    % OF = OF;

  % F_passive(i) = 0*((SL(i)-1.6)/0.6)^1;
 
end
 
  F_total = F_active + F_passive;

figure(1); clf;
nexttile;
hold on;
plot(SLd, df, 'o-', LineWidth=2)
plot(SL,F_total,'o-', SL,F_passive,'x--', LineWidth=2)
% plot(SL, pS*100)
% nexttile;
% plot(SL, p1_0, SL, p2_0, SL, pT, SL, pD, SL, pS, LineWidth=2);legend('p1', 'p2', 'T', 'D', 'SR', Location='best')

%%
g = [  1.6536    0.1325    1.0119    1.0660]
g = [ 0.7591    1.0094    0.5783    1.3588];
g = [ 0.6791    0.9122    0.6484    1.2397];
g = [0.8770    0.9150    0.8324    1.2412];
g = [2.2429    0.9147    2.1301    1.2390];
resimulate(g)
% g = ones(4, 1);
% g = [x;1]
%%
resimulate(g, true)
%%
g = [4 3 2 1]
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.1, 'PlotFcns', @optimplotfval, 'MaxIter', 1500);
% g = [1, 1, 1, 1, 1, 1, 1, 1];
% g = [1.2539    0.4422];

x = fminsearch(@resimulate, g, options)
g = x';

%%
function Es = resimulate(g, plotThat)

if nargin < 2
    plotThat = false;
end

SL = [2.2 2.0400    2.0000    1.9600    1.9200  1.8800];
df = [76.5 68.3878   65.5438   59.2362   52.5507   43.9655];
N = 10; % space (strain) discretization--number of grid points in half domain
Slim = 0.20; dS = Slim/N; s = (-N:1:0)*dS; % strain
p1 = zeros(N+1,1); p2 = zeros(N+1,1); pT = 0; pD = 0; pS = 1;
 
% State variable vector concatenates p1, p2, p2, and U_NR
PU0 = [p1; p2; pT; pD;pS];

for i = 1:length(SL)
 
  % run model to steady state
  [t,PU] = ode15s(@dPUdT,[0 100],PU0,[],N,dS,SL(i), g);
 
  PUend = PU(end,:);
  p1 = PUend(1:1*N+1);  % probabilities of attatched state 1 on grid
  p2 = PUend(1*N+2:2*N+2); % probabilities of attatched state 2 on grid
  pT(i) = PUend(2*N+3); % probability of unattached state T
  pD(i) = PUend(2*N+4); % probability of unattached state D
  pS(i) = PUend(2*N+5); % probability of unattached state D
 
  p1_0(i) = dS*sum(p1); p1_1 = dS*sum(s.*p1);
  p2_0(i) = dS*sum(p2); p2_1 = dS*sum(s.*p2);
  pS(i) = 1 - p1_0(i) - p2_0(i) - pT(i) - pD(i); % SRX probability
  
  dr = 0.01; % Power-stroke Size; Units: um
  kstiff1 = 1e4*g(4);
  F_active(i) = kstiff1*(p1_1 + p2_1 + dr*p2_0(i));

  kSE = 5000;
  LSE(i) = F_active(i)/kSE;
  % LSE =0;
y = @(k_pas, x0, gamma, x) k_pas.*(x-x0).^gamma - 4*0 - x0*0 + 0.5e9.*(x-0.95).^13;    
F_passive(i) = y(0.4, -0.4, 7.9, (SL(i)-LSE(i))/2); % from the ramp-ups
 
end
 
  F_total = F_active + F_passive;

if plotThat
    %%
  figure(1); clf;
hold on;
plot(SL, df, 'o-', LineWidth=2)
plot(SL,F_total,'o-', SL,F_passive,'x--', LineWidth=2)
end

E = (F_total - df).^2;
Es = sum(E);

end


function f = dPUdTOld(t,PU,N,dS,SL)
% ODE function for the d/dt operator for the cross-bridge mode.
%  first 2N-1 entries of PU represent p1(s,t)
%  second 2N-1 entries represent p2(s,t)
%  The remaining entries are pT and pD, the proabilities of the T and D
%  states
%
%  SL = sarcomere length

% Overlap function
OV = min(1, 0.6 + (SL-1.8)); % this is just a function that I made up

% State Variables
p1 = PU(1:1*N+1);  % probabilities of attatched state 1 on grid
p2 = PU(1*N+2:2*N+2); % probabilities of attatched state 2 on grid
pT = PU(2*N+3); % probability of unattached state T
pD = PU(2*N+4); % probability of unattached state D
pS = PU(2*N+5); % probability of unattached state D

% calculation of moments of strain distributions
s = (-N:1:0)'*dS;
p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
% pS = 1 - p1_0 - p2_0 - pT - pD; % SRX probability

% definition of parameters
kTD = 25; % (1/sec) transition from T to D
kD1 = 25; % (1/sec) transition from D to 1
k12 = 25; % (1/sec) transition from 1 to 2
k2T = 2; % (1/sec) transition from 2 to T (detachment)

% Force model
dr = 0.01; % Power-stroke Size; Units: um
kstiff1 = 1e4;
F_active = kstiff1*(p1_1 + p2_1 + dr*p2_0);

% transitions between super relaxed state and non relaxed state
kST    = 10; % transition from S to T
sigma0 = 25;
kTS    = 600; % transition from T to SRX

OF     = (OV - (p1_0 + p2_0))/(pS + pT + pD); % overlap factor
OF = OV;

dp1      = - k12*p1 ; % transition from 1 to 2
dp1(N+1) = dp1(N+1) + kD1*OF*pD/dS; % attachment
dp2      = + k12*p1 - k2T*p2; %
dpT      = dS*sum(k2T*p2) - kTD*pT - kTS*pT + kST*exp(F_active/sigma0)*pS;
dpS = + kTS*pT - kST*exp(F_active/sigma0)*pS;
% dpT      = dS*sum(k2T*p2) - kTD*pT - kTS*pT + kST*(F_active/sigma0 + 0.01)*pS;
dpD      = + kTD*pT - kD1*pD*OF;
if t > 90
    breakpoints = true;
end
f = [dp1; dp2; dpT; dpD;dpS];
end

function f = dPUdT(~,PU,N,dS,SL, g)
% ODE function for the d/dt operator for the cross-bridge mode.
%  first 2N-1 entries of PU represent p1(s,t)
%  second 2N-1 entries represent p2(s,t)
%  The remaining entries are pT and pD, the proabilities of the T and D
%  states
%
%  SL = sarcomere length
 
% Overlap function
OV = min(1, 0.6 + (SL-1.8)); % this is just a function that I made up
 
% State Variables
p1 = PU(1:1*N+1);  % probabilities of attatched state 1 on grid
p2 = PU(1*N+2:2*N+2); % probabilities of attatched state 2 on grid
pT = PU(2*N+3); % probability of unattached state T
pD = PU(2*N+4); % probability of unattached state D
pS = PU(2*N+5); % probability of unattached state D

% calculation of moments of strain distributions
s = (-N:1:0)'*dS;
p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
pS = 1 - p1_0 - p2_0 - pT - pD; % SRX probability
 
% definition of parameters
kTD = 100; % (1/sec) transition from T to D
kD1 = 100; % (1/sec) transition from D to 1
k12 = 25; % (1/sec) transition from 1 to 2
k2T = 10; % (1/sec) transition from 2 to T (detachment)
kSE = 5000;

% Force model
dr = 0.01; % Power-stroke Size; Units: um
kstiff1 = 1e4;
F_active = kstiff1*(p1_1 + p2_1 + dr*p2_0);


LSE = F_active/kSE;
y = @(k_pas, x0, gamma, x) k_pas.*(x-x0).^gamma - 4*0 - x0*0 + 0.5e9.*(x-0.95).^13;    
F_passive = y(0.4, -0.4, 7.9, (SL-LSE)/2); % from the ramp-ups
% F_passive = 15*((SL-1.6)/0.6)^1; % Original dan's

F_total = F_active + F_passive;
 
% transitions between super relaxed state and non relaxed state
kST    = 5*g(1); % transition from S to T
sigma0 = 15*g(2);
kTS    = 1000*g(3); % transition from T to SRX
 
OF     = (OV - (p1_0 + p2_0))/(pS + pT + pD); % overlap factor
% OF = OF;

dp1      = - k12*p1 ; % transition from 1 to 2
dp1(N+1) = dp1(N+1) + kD1*OF*pD/dS; % attachment
dp2      = + k12*p1 - k2T*p2; %
dpT      = dS*sum(k2T*p2) - kTD*pT - kTS*pT + kST*exp((F_total/sigma0)^1)*pS;
% dpT      = dS*sum(k2T*p2) - kTD*pT - kTS*pT + kST*(F_total/sigma0)*pS;
% dpS = + kTS*pT - kST*exp((F_total/sigma0)^1)*pS;
dpD      = + kTD*pT - kD1*OF*pD;
 
f = [dp1; dp2; dpT; dpD;0];
end