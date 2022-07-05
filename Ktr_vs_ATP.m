 clear
% g0 = ones(1,12);
load g0

%% Data -- need to verify these data
Ktr_mean = [25.8033 29.0 37.7928];
Ktr_err  = [2.0167  1.30    1.9308];

%% Setting up problem
N = 50; % space (strain) discretization--number of grid points in half domain
Slim = 0.075; 
dS = Slim/N;
s = (-N:1:0)*dS; % strain 
p1 = zeros(N+1,1);
p2 = zeros(N+1,1);
p3 = zeros(N+1,1);
U_NR = 1; % unknown
% State variable vector concatenates p1, p2, p2, and U_NR
PU0 = [p1; p2; p3; U_NR];

% Set metabolite concentrations, 
% MgATP = logspace(-2,1,25);
MgATP = 1:9;
MgADP = 0; 
Pi    = 0; 

%% Simulating sliding and kinetics via Strang operator splitting

% moments and force
dr = 0.01; % Power-stroke Size; Units: um
kstiff1 = g0(13)*1500; 
kstiff2 = g0(14)*10000; 
alpha3 = g0(9)*125;
s3     = 0.010;
K_T1 = g0(11)*1.0; % (mM) 
K_D = 0.194; % MgADP dissociation constant from Yamashita etal (Circ Res. 1994; 74:1027-33).

Tspan = [0:0.001:0.12];

for k = 1:length(MgATP)

  % Zero velocity:
  [t,PU] = ode15s(@dPUdT,Tspan,PU0,[],N,dS,MgATP(k),Pi,MgADP,g0);  
%   PU = PU(end,:);
  p1 = PU(:,1:1*N+1);
  p2 = PU(:,1*N+2:2*N+2);
  p3 = PU(:,2*N+3:3*N+3);
  p1_0 = dS*sum(p1'); p1_1 = dS*sum(s'.*p1');
  p2_0 = dS*sum(p2'); p2_1 = dS*sum(s'.*p2');
  p3_0 = dS*sum(p3'); p3_1 = dS*sum(s'.*p3');
  F_active(k,:) = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 ) ;
  
  Frel = F_active(k,:)./F_active(k,end);
  Ktr(k) = 1/interp1(Frel,Tspan,1-exp(-1)); % time constant for Frel(1/Ktr) = 1-exp(-1)

end

%% plots

figure(6); clf; axes('position',[0.15 0.15 0.8 0.8]); hold on;
plot(Tspan,F_active(2,:)./F_active(2,end),'r-','linewidth',1.5);
plot(Tspan,F_active(4,:)./F_active(4,end),'g-','linewidth',1.5);
plot(Tspan,F_active(8,:)./F_active(8,end),'b-','linewidth',1.5);
% plot(Tspan,F_active,'linewidth',1.5);
xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
ylabel('Force (rel.)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'ylim',[0 1.1]);  box on;
ldg = legend('2','4','8 mM','location','northwest');
title(ldg,'[MgATP]');

axes('position',[0.5 0.25 0.4 0.4]); hold on;
plot(MgATP,Ktr,'k-','linewidth',1.5);
errorbar(MgATP([2 4 8]),Ktr_mean,Ktr_err,'ko','linewidth',1.5,'markersize',6);
xlabel('[MgATP] (mM)','interpreter','latex','fontsize',6);
ylabel('$K_{tr}$ (sec.$^{-1}$)','interpreter','latex','fontsize',6);
set(gca,'fontsize',9,'xlim',[0 10],'ylim',[0 45]); box on;

% figure(7); clf; axes('position',[0.15 0.15 0.8 0.8]); hold on;
% plot(MgATP,Ktr,'k-','linewidth',1.5);
% errorbar(MgATP([2 4 8]),Ktr_mean,Ktr_err,'ko','linewidth',1.5,'markersize',12);
% xlabel('[MgATP] (mM)','interpreter','latex','fontsize',16);
% ylabel('$K_{tr}$ (sec.$^{-1}$)','interpreter','latex','fontsize',16);
% set(gca,'fontsize',14,'xlim',[0 10],'ylim',[0 45]); box on;
