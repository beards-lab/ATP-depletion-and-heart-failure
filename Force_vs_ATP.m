clear
% g0 = ones(1,14);
load g0

%% Data 
% Fmax (normalized.) versus [MgATP] (mM) from Ebus et al.(2001)
data = ...
    [0.0098736     0.019874      0.04959     0.098478      0.49024        5.063
      1.5925       1.6826       1.5898       1.4657       1.2884      0.99732];

%% Setting up problem
N = 40; % space (strain) discretization--number of grid points in half domain
Slim = 0.075; 
dS = Slim/N;
s = (-N:1:0)*dS; % strain 
p1 = zeros(N+1,1);
p2 = zeros(N+1,1);
p3 = zeros(N+1,1);
U_NR = 1;
% State variable vector concatenates p1, p2, p2, and U_NR
PU0 = [p1; p2; p3; U_NR];

% Set metabolite concentrations, 
MgATP = logspace(-2,log10(8),50);
% MgATP = linspace(0.01,5,250);
MgADP = 0.0; 
Pi    = 0; 

%% Simulating sliding and kinetics via Strang operator splitting

% moments and force
dr = 0.01; % Power-stroke Size; Units: um
kstiff1 = g0(13)*1500; 
kstiff2 = g0(14)*10000;

for k = 1:length(MgATP)

  % Zero velocity:
  [t,PU] = ode15s(@dPUdT,[0 1],PU0,[],N,dS,MgATP(k),Pi,MgADP,g0);  
  PU = PU(end,:);
  p1 = PU(1:1*N+1);
  p2 = PU(1*N+2:2*N+2);
  p3 = PU(2*N+3:3*N+3);
  p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
  p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
  p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
  F_active(k) = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 ) ;
  
end

%% plots

figure(4); clf; axes('position',[0.15 0.15 0.8 0.8]); 
semilogx(MgATP,F_active,'b-','linewidth',1.5); hold on;
semilogx(data(1,:),data(2,:).*55,'ko','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1]);
xlabel('[MgATP] (mM)','interpreter','latex','fontsize',16);
ylabel('Max. Force (kPa)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'ylim',[0 100]); 
box on;


% figure(5); clf; axes('position',[0.15 0.15 0.8 0.8]); 
% semilogx(MgATP,XBCR,'b-','linewidth',1.5);
% xlabel('[MgATP] (mM)','interpreter','latex','fontsize',16);
% ylabel('XBCR','interpreter','latex','fontsize',16);
% set(gca,'fontsize',14,'ylim',[0 25]); 
% box on;

