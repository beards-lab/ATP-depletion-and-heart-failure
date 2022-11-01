% clear
% load g0

%% Setting up problem
N = 50; % space (strain) discretization--number of grid points in half domain
Slim = 0.040; 
dS = Slim/N;
s = (-N:1:N)*dS; % strain 
p1 = zeros(2*N+1,1);
p2 = zeros(2*N+1,1);
p3 = zeros(2*N+1,1);
U_NR = 0;
Ls = 0.0; % initial length of series element (micron)
Lh = 1.0; % initial length of half sarcomere (micron)
% State variable vector concatenates p1, p2, p2, and U_NR
PU0 = [p1; p2; p3; U_NR; Ls; Lh];

% Set metabolite concentrations, 
MgATP = 5;
MgADP = 0; 
Pi    = 0; 

% moments and force parameters
mu = 0.1; % viscosity
dr = 0.01; % Power-stroke Size; Units: um
kstiff2 = g0(13)*5000; 
kstiff3 = g0(14)*27500;
kSE = 25000;

%% Simulating zero-velocity initial state
[t,PU] = ode15s(@dPUdT,[0 1],PU0,[],N,dS,MgATP,Pi,MgADP,g0,0);
PU = PU(end,:);
p1 = PU(1:2*N+1);
p2 = PU(2*N+2:4*N+2);
p3 = PU(4*N+3:6*N+3);
p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
% F_active = kstiff2*p3_0 + kstiff1*(  p2_1 + p3_1 ) ; % initial steady-state force
F_active = kstiff3*p3_1 + kstiff2*p2_1;
Ls = PU(6*N+5); % length series element
Lh = PU(6*N+6); % length series element
F_SE = kSE*Ls;  % force series element

T_sim = 0;
F_passive = 1;
F_sim = F_active + F_passive;
ML_sim = Lh; % initial muscle length
HS_sim = Lh; % initial half-sarcomere length

%% Simulating sliding and kinetics via Strang operator splitting
  
  
  % ramp for 40 ms
  RampLength = 0.01;
  vel = RampLength/0.020;
  Nstep = ceil(RampLength/dS);
  dt = RampLength/vel/Nstep;

 for i = 1:Nstep
    % simulate kinetics for dt
    [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP,Pi,MgADP,g0,0);
    PU = PU(end,:); 
    p1 = PU(1:2*N+1);
    p2 = PU(2*N+2:4*N+2);
    p3 = PU(4*N+3:6*N+3);
    p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
    p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
    F_active = kstiff3*p3_1 + kstiff2*p2_1;
    Ls = PU(6*N+5); % length series element
    Lh = PU(6*N+6); % 1/2 sarcomere length
    F_SE = kSE*Ls;  % force series element
    F_passive = 1 + 40*(Lh-1);
    F_sim = [F_sim, F_SE+F_passive];
    T_sim = [T_sim, T_sim(end)+dt ];
    ML_sim = [ML_sim, ML_sim(end)+dt*0];
    HS_sim = [HS_sim, Lh];
  end
  
  for j = 1:7 % number of ramps
    for i = 1:Nstep
      % simulate kinetics for dt
      [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP,Pi,MgADP,g0,vel);
      PU = PU(end,:); 
      p1 = PU(1:2*N+1);
      p2 = PU(2*N+2:4*N+2);
      p3 = PU(4*N+3:6*N+3);
      p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
      p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
      p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
      F_active = kstiff3*p3_1 + kstiff2*p2_1;
      Ls = PU(6*N+5); % length series element
      Lh = PU(6*N+6); % 1/2 sarcomere length
      F_SE = kSE*Ls;  % force series element
      F_passive = 1 + 40*(Lh-1);
      F_sim = [F_sim, F_SE+F_passive];
      T_sim = [T_sim, T_sim(end)+dt ];
      ML_sim = [ML_sim, ML_sim(end)+dt*vel];
      HS_sim = [HS_sim, Lh];
    end
  
    % steady for 20 ms
    for i = 1:Nstep
      % simulate kinetics for dt
      [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP,Pi,MgADP,g0,0);
      PU = PU(end,:); 
      p1 = PU(1:2*N+1);
      p2 = PU(2*N+2:4*N+2);
      p3 = PU(4*N+3:6*N+3);
      p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
      p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
      p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
      F_active = kstiff3*p3_1 + kstiff2*p2_1;
      Ls = PU(6*N+5); % length series element
      Lh = PU(6*N+6); % 1/2 sarcomere length
      F_SE = kSE*Ls;  % force series element
      F_passive = 1 + 40*(Lh-1);
      F_sim = [F_sim, F_SE+F_passive];
      T_sim = [T_sim, T_sim(end)+dt ];
      ML_sim = [ML_sim, ML_sim(end)+dt*0];
      HS_sim = [HS_sim, Lh];
    end
  end
  
  % last ramp at 3x speed
  RampLength = 0.03;
  vel = RampLength/0.020;
  Nstep = ceil(RampLength/dS);
  dt = RampLength/vel/Nstep;
  
  for i = 1:Nstep
    % simulate kinetics for dt
    [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP,Pi,MgADP,g0,vel);
    PU = PU(end,:); 
    p1 = PU(1:2*N+1);
    p2 = PU(2*N+2:4*N+2);
    p3 = PU(4*N+3:6*N+3);
    p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
    p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
    F_active = kstiff3*p3_1 + kstiff2*p2_1;
    Ls = PU(6*N+5); % length series element
    Lh = PU(6*N+6); % 1/2 sarcomere length
    F_SE = kSE*Ls;  % force series element
    F_passive = 1 + 40*(Lh-1);
    F_sim = [F_sim, F_SE+F_passive];
    T_sim = [T_sim, T_sim(end)+dt ];
    ML_sim = [ML_sim, ML_sim(end)+dt*vel];
    HS_sim = [HS_sim, Lh];
  end
  
  % steady for 200 ms
  for i = 1:25*Nstep
    % simulate kinetics for dt
    [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP,Pi,MgADP,g0,0);
    PU = PU(end,:); 
    p1 = PU(1:2*N+1);
    p2 = PU(2*N+2:4*N+2);
    p3 = PU(4*N+3:6*N+3);
    p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
    p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
    F_active = kstiff3*p3_1 + kstiff2*p2_1;
    Ls = PU(6*N+5); % length series element
    Lh = PU(6*N+6); % 1/2 sarcomere length
    F_SE = kSE*Ls;  % force series element
    F_passive = 1 + 40*(Lh-1);
    F_sim = [F_sim, F_SE+F_passive];
    T_sim = [T_sim, T_sim(end)+dt ];
    ML_sim = [ML_sim, ML_sim(end)+dt*0];
    HS_sim = [HS_sim, Lh];
  end
  
%% plots

% Overlap function
SL = 2*HS_sim;
L_thick = 1.67; % Length of thick filament, um
L_bare = 0.10; % Length of bare region of thick filament, um
L_thin  = 1.20; % Length of thin filament, um
LthickHS1 = min(L_thick*0.5, SL*0.5);
LthickHS2 = max(SL*0.5 - (SL-L_thin),L_bare*0.5);
L_OV = LthickHS1-LthickHS2;
OV = 2*L_OV/(L_thick - L_bare);
% OV = 1;

% strain distributions
p1 = PU(1:2*N+1);
p2 = PU(2*N+2:4*N+2);
p3 = PU(4*N+3:6*N+3);

figure(1);
plot(s,p1,s,p2,s,p3,'linewidth',1.5);
ylabel('Probability density ($\mu$m$^{-1}$)','interpreter','latex','fontsize',16);
xlabel('strain, $s$ ($\mu$m)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14);
set(gca,'xlim',[-Slim +Slim],'xtick',-Slim:0.02:Slim);
legend('$p_1(s)$','$p_2(s)$','$p_3(s)$','interpreter','latex','fontsize',16,'location','northwest');

figure(2); clf
plot(td,ld,'linewidth',1.5); hold on;
plot(T_sim*1e3,ML_sim,'linewidth',1.5); hold on;
ylabel('Length ($\mu$m)','interpreter','latex','fontsize',16);
xlabel('time (ms)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'Xlim',[0 450]);

figure(3);clf;hold on;
plot(td,fd,'linewidth',1.5);
plot(T_sim*1e3,F_sim,'linewidth',1.5);
ylabel('Force (kPa)','interpreter','latex','fontsize',16);
xlabel('time (ms)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'Xlim',[0 450]);

figure(4)
plot(T_sim*1e3,HS_sim,'linewidth',1.5);
ylabel('$L_h$, Half S Length ($\mu$m)','interpreter','latex','fontsize',16);
xlabel('time (ms)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'Xlim',[0 450]);

figure(5);
plot(T_sim*1e3,OV,'linewidth',1.5);
ylabel('OV','interpreter','latex','fontsize',16);
xlabel('time (ms)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'Xlim',[0 450]);

