clear
load g0

%% Setting up problem
N = 50; % space (strain) discretization--number of grid points in half domain
Slim = 0.075; 
dS = Slim/N;
s = (-N:1:0)*dS; % strain 
p1 = zeros(1*N+1,1);
p2 = zeros(1*N+1,1);
p3 = zeros(1*N+1,1);
U_NR = 0;
% State variable vector concatenates p1, p2, p2, and U_NR
PU0 = [p1; p2; p3; U_NR];

% Set metabolite concentrations, 
MgATP = 8;
MgADP = 0; 
Pi    = 0; 

% moments and force parameters
mu = 0.1; % viscosity
dr = 0.01; % Power-stroke Size; Units: um
kstiff1 = g0(13)*1500; 
kstiff2 = g0(14)*10000;

%   Set the outer timestep based on space step:
SLi = 1.21; % half sarcomere length at which experiment starts (micron)
SLo = 1.10; % half sarcomere length at which experiment ends (micron)


%% Simulating zero-velocity initial state
[t,PU] = ode15s(@dPUdT,[0 1],PU0,[],N,dS,MgATP,Pi,MgADP,g0,SLi);
PU = PU(end,:);
p1 = PU(1:1*N+1);
p2 = PU(1*N+2:2*N+2);
p3 = PU(2*N+3:3*N+3);
p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
F_active = kstiff2*p3_0*dr + kstiff1*(  p2_1 + p3_1 ) ; % initial steady-state force

T_sim = 0;
F_sim = F_active;
SL_sim = SLi;

%% Simulating sliding and kinetics via Strang operator splitting

% (1.) first run unloaded until active force drops to zero
forcepositive = true;
while forcepositive
    
  vel = -F_active/mu;
  dt = dS/abs(vel);
  
%   SL(i+1) = SL(i) + dt*vel;
  % advection (sliding step)
  PU(1:1*N+0)     = PU(2:1*N+1); PU(N+1) = 0;
  PU(1*N+2:2*N+1) = PU(1*N+3:2*N+2); PU(2*N+2) = 0;
  PU(2*N+3:3*N+2) = PU(2*N+4:3*N+3); PU(3*N+3) = 0;
  % simulate kinetics for full step
  [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP,Pi,MgADP,g0,1.1);
  PU = PU(end,:); 
  p1 = PU(1:1*N+1);
  p2 = PU(1*N+2:2*N+2);
  p3 = PU(2*N+3:3*N+3);
  p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
  p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
  p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
  F_active = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1) ;
  forcepositive = F_active > 25;
  F_sim = [F_sim, F_active];
  T_sim = [T_sim, T_sim(end)+dt ];
  SL_sim = [SL_sim, SL_sim(end)+dt*vel];
  
end

%% (2.) next run at constant velocity until SL = SLo
vel = -5.0*1.1; % micron per sec
dt = dS/abs(vel);
tend = 0.1/abs(vel); % ending time of simulation
Nstep = round(tend/dt);

% while SL_sim(end) > SLo
for i = 1:75
  
%   SL(i+1) = SL(i) + dt*vel;
  % advection (sliding step)
  PU(1:1*N+0)     = PU(2:1*N+1); PU(N+1) = 0;
  PU(1*N+2:2*N+1) = PU(1*N+3:2*N+2); PU(2*N+2) = 0;
  PU(2*N+3:3*N+2) = PU(2*N+4:3*N+3); PU(3*N+3) = 0;
  % simulate kinetics for full step
  [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP,Pi,MgADP,g0,1.1);
  PU = PU(end,:); 
  p1 = PU(1:1*N+1);
  p2 = PU(1*N+2:2*N+2);
  p3 = PU(2*N+3:3*N+3);
  p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
  p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
  p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
  F_active = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1) ;
  F_sim = [F_sim, F_active];
  T_sim = [T_sim, T_sim(end)+dt ];
  SL_sim = [SL_sim, SL_sim(end)+dt*vel];

end


%% plots

% strain distributions
p1 = PU(1:1*N+1);
p2 = PU(1*N+2:2*N+2);
p3 = PU(2*N+3:3*N+3);
figure(8);
plot(s,p1,s,p2,s,p3,'linewidth',1.5);
ylabel('Probability density ($\mu$m$^{-1}$)','interpreter','latex','fontsize',16);
xlabel('strain, $s$ ($\mu$m)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14);
set(gca,'xlim',[-Slim 0]);
legend('$p_1(s)$','$p_2(s)$','$p_3(s)$','interpreter','latex','fontsize',16,'location','northwest');

%
figure(9)
plot(T_sim*1e3,F_sim,'linewidth',1.5);
ylabel('Force (kPa)','interpreter','latex','fontsize',16);
xlabel('time (ms)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'Xlim',[0 120]);

figure(10); clf
plot(T_sim*1e3,SL_sim,'linewidth',1.5); hold on;
plot([0 120]*1e3,[SLo SLo],'k--');
ylabel('Length ($\mu$m)','interpreter','latex','fontsize',16);
xlabel('time (ms)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'Xlim',[0 120]);
