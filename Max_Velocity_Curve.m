clear
% g0 = ones(1,12);
load g0


%% Setting up problem
N = 40; % space (strain) discretization--number of grid points in half domain
Slim = 0.080; 
dS = Slim/N;
s = (-N:1:0)*dS; % strain 
p1 = zeros(N+1,1);
p2 = zeros(N+1,1);
p3 = zeros(N+1,1);
U_NR = 0;
% State variable vector concatenates p1, p2, p2, and U_NR
PU0 = [p1; p2; p3; U_NR];

% Set metabolite concentrations, 
MgATP = [ 0.02:0.02:0.1, 0.15:0.05:0.45, 0.5:0.5:2.5];
MgADP = 0; 
Pi    = 0; 

%% Simulating sliding and kinetics via Strang operator splitting

% moments and force
dr = 0.01; % Power-stroke Size; Units: um
kstiff1 = g0(13)*2500; 
kstiff2 = g0(14)*200;
mu = 0.50;
F_load = 0; % afterload force

for j = 1:length(MgATP)

  Force = 1; % initial value (arbitrary)
  vel = -6.35/(1 + (0.40/MgATP(j))); % initial guess
  % vel = 2*1.1;

  while abs(F_load-Force) > 0.0001
    
    dt = dS/abs(vel);
    tend = 0.40/abs(vel); % ending time of simulation
    Nstep = round(tend/dt);
    % simulate kinetics for 1/2 timestep
    [t,PU] = ode15s(@dPUdT,[0 dt/2],PU0,[],N,dS,MgATP(j),Pi,MgADP,g0);
    PU = PU(end,:); 
    for i = 1:(Nstep-1)
      % advection (sliding step)
      PU(1:1*N+0)     = PU(2:1*N+1); PU(N+1) = 0;
      PU(1*N+2:2*N+1) = PU(1*N+3:2*N+2); PU(2*N+2) = 0;
      PU(2*N+3:3*N+2) = PU(2*N+4:3*N+3); PU(3*N+3) = 0;
      % simulate kinetics for full step
      [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP(j),Pi,MgADP,g0);
      PU = PU(end,:); 
    end
    % final advection (sliding step)
    PU(1:1*N+0)     = 0.5*(PU(2:1*N+1) + PU(1:1*N+0));         PU(N+1) = 0.5*(0 + PU(N+1));
    PU(1*N+2:2*N+1) = 0.5*(PU(1*N+3:2*N+2) + PU(1*N+2:2*N+1)); PU(2*N+2) = 0.5*(0 + PU(2*N+2));
    PU(2*N+3:3*N+2) = 0.5*(PU(2*N+4:3*N+3) + PU(2*N+3:3*N+2)); PU(3*N+3) = 0.5*(0 + PU(3*N+3));
    % final 1/2 timestep for kinetics
    [t,PU] = ode15s(@dPUdT,[0 dt/2],PU,[],N,dS,MgATP(j),Pi,MgADP,g0);
    PU = PU(end,:);
  
    p1 = PU(1:1*N+1);
    p2 = PU(1*N+2:2*N+2);
    p3 = PU(2*N+3:3*N+3);
    p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
    p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
    Force = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 ) + mu*vel;
  
    vel = vel + 1*(F_load - Force);
    [Force vel]
  
  end
  vel_max(j) = vel;
  [MgATP(j) vel]
end

MgATP = [0 MgATP];
vel_max = [0 vel_max];

%%

figure(6); clf; axes('position',[0.15 0.15 0.8 0.8]); 
plot(MgATP,-vel_max.*1.1,'bo-','linewidth',1.5); hold on;
xlabel('[MgATP] (mM)','interpreter','latex','fontsize',16);
ylabel('Max. Velocity (ML/sec.)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'ylim',[0 8]); 
box on;





