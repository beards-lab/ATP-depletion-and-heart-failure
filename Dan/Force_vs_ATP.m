clear
% g0 = ones(1,14);
load g0

%% Data 
% Fmax (normalized.) versus [MgATP] (mM) from Ebus et al.(2001)
data = ...
    [0.010     0.020      0.050     0.10      0.50        5.
     1.5925    1.6826     1.5898    1.4657    1.2884      0.99732];

%% Setting up problem
N = 50; % space (strain) discretization--number of grid points in half domain
Slim = 0.040; 
dS = Slim/N;
s = (-N:1:N)*dS; % strain 
p1 = zeros(2*N+1,1);
p2 = zeros(2*N+1,1);
p3 = zeros(2*N+1,1);
U_NR = 1;
LSE = 0; % length of series element
% State variable vector concatenates p1, p2, p2, and U_NR
PU0 = [p1; p2; p3; U_NR; LSE];

% Set metabolite concentrations, 
MgATP = logspace(-2,log10(8),75);
% MgATP = data(1,:);
MgADP = 0.0; 
Pi    = 0; 

%% Simulating sliding and kinetics via Strang operator splitting

% moments and force
dr = +g0(12)*0.01; % Power-stroke Size; Units: um
kstiff1 = g0(13)*2500; 
kstiff2 = g0(14)*20000;

for k = 1:length(MgATP)

  % Zero velocity:
  [t,PU] = ode15s(@dPUdT,[0 1],PU0,[],N,dS,MgATP(k),Pi,MgADP,g0,0);  
  PU = PU(end,:);
  p1 = PU(1:2*N+1);
  p2 = PU(2*N+2:4*N+2);
  p3 = PU(4*N+3:6*N+3);
  p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
  p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
  p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
  % F_active = kstiff2*p3_0 + kstiff1*(  p2_1 + p3_1 ) ; % initial steady-state force
  F_active (k)= kstiff2*p3_1 + kstiff1*p2_1;
  
end

%% plots

figure(5); clf; axes('position',[0.15 0.15 0.8 0.8]); 
semilogx(MgATP,F_active,'b-','linewidth',1.5); hold on;
semilogx(data(1,:),data(2,:).*57,'ko','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1]);
xlabel('[MgATP] (mM)','interpreter','latex','fontsize',16);
ylabel('Max. Force (kPa)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'ylim',[0 120]); 
box on;


% figure(5); clf; axes('position',[0.15 0.15 0.8 0.8]); 
% semilogx(MgATP,XBCR,'b-','linewidth',1.5);
% xlabel('[MgATP] (mM)','interpreter','latex','fontsize',16);
% ylabel('XBCR','interpreter','latex','fontsize',16);
% set(gca,'fontsize',14,'ylim',[0 25]); 
% box on;

