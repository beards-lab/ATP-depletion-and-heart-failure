% clear
% load g0

%% Setting up problem
N = 40; % space (strain) discretization--number of grid points in half domain
Slim = 0.040; 
dS = Slim/N;
s = (-N:1:N)*dS; % strain 
p1 = zeros(2*N+1,1);
p2 = zeros(2*N+1,1);
p3 = zeros(2*N+1,1);
U_NR = 0;
LSE = 0; % length of series element
% State variable vector concatenates p1, p2, p2, and U_NR
PU0 = [p1; p2; p3; U_NR; LSE];

% Set metabolite concentrations, 
MgATP = 5;
MgADP = 0; 
Pi    = 0; 

% moments and force parameters
mu = 0.01; % viscosity
dr = 0.01; % Power-stroke Size; Units: um
kstiff1 = g0(13)*2500; 
kstiff2 = g0(14)*20000;
kSE = 10000;


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
F_active = kstiff2*p3_1 + kstiff1*p2_1;
LSE = PU(6*N+5); % length series element
F_SE = kSE*LSE;  % force series element

T_sim = 0;
F_sim = F_active;
SL_sim = 1.0; % initial SL

%% Simulating sliding and kinetics via Strang operator splitting
  
  
  % ramp for 40 ms
  RampLength = 0.01;
  vel = RampLength/0.020;
  Nstep = ceil(RampLength/dS);
  dt = RampLength/vel/Nstep;
  
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
      F_active = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1) ;
      LSE = PU(6*N+5); % length series element
      F_SE = kSE*LSE;  % force series element
      F_sim = [F_sim, F_SE];
      T_sim = [T_sim, T_sim(end)+dt ];
      SL_sim = [SL_sim, SL_sim(end)+dt*vel];
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
      F_active = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1) ;
      LSE = PU(6*N+5); % length series element
      F_SE = kSE*LSE;  % force series element
      F_sim = [F_sim, F_SE];
      T_sim = [T_sim, T_sim(end)+dt ];
      SL_sim = [SL_sim, SL_sim(end)+dt*0];
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
    F_active = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1) ;
    LSE = PU(6*N+5); % length series element
    F_SE = kSE*LSE;  % force series element
    F_sim = [F_sim, F_SE];
    T_sim = [T_sim, T_sim(end)+dt ];
    SL_sim = [SL_sim, SL_sim(end)+dt*vel];
  end
  
  % steady for 200 ms
  for i = 1:5*Nstep
    % simulate kinetics for dt
    [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP,Pi,MgADP,g0,0);
    PU = PU(end,:); 
    p1 = PU(1:2*N+1);
    p2 = PU(2*N+2:4*N+2);
    p3 = PU(4*N+3:6*N+3);
    p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
    p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
    F_active = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1) ;
    LSE = PU(6*N+5); % length series element
    F_SE = kSE*LSE;  % force series element
    F_sim = [F_sim, F_SE];
    T_sim = [T_sim, T_sim(end)+dt ];
    SL_sim = [SL_sim, SL_sim(end)+dt*0];
  end
  
%% plots

% strain distributions
p1 = PU(1:2*N+1);
p2 = PU(2*N+2:4*N+2);
p3 = PU(4*N+3:6*N+3);
figure(8);
plot(s,p1,s,p2,s,p3,'linewidth',1.5);
ylabel('Probability density ($\mu$m$^{-1}$)','interpreter','latex','fontsize',16);
xlabel('strain, $s$ ($\mu$m)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14);
set(gca,'xlim',[-Slim +Slim],'xtick',-Slim:0.02:Slim);
legend('$p_1(s)$','$p_2(s)$','$p_3(s)$','interpreter','latex','fontsize',16,'location','northwest');

%
figure(9)
plot(T_sim*1e3,F_sim,'linewidth',1.5);
ylabel('Force (kPa)','interpreter','latex','fontsize',16);
xlabel('time (ms)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'Xlim',[0 400]);

figure(10); clf
plot(T_sim*1e3,SL_sim,'linewidth',1.5); hold on;
% plot([0 120]*1e3,[SLo SLo],'k--');
ylabel('Length ($\mu$m)','interpreter','latex','fontsize',16);
xlabel('time (ms)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'Xlim',[0 400]);
