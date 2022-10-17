clear

%% Data ([ATP] = 8, 4, 6 mM)

Data = [0		56.4048	63.6074	61.3192
        0.5		51.812	51.8626	47.4794
        1		37.4459	35.9182	31.387
        2		17.8025	13.5516	10.2112
        3		11.443	8.34	6.3895
        4		6.2643	2.8669	2.5781
        5		3.2759	1.6526	1.594
        6		2.212	1.1823	1.2117];

% Data = [0		56.4048	63.6074	61.3192
%         1		37.4459	35.9182	31.387
%         2		17.8025	13.5516	10.2112
%         3		11.443	8.34	6.3895
%         6		2.212	1.1823	1.2117];
    
ML = 1.1; % half sarcomere length (microns)

%% Setting up problem
N = 20; % space (strain) discretization--number of grid points in half domain
Slim = 0.060; 
dS = Slim/N;
s = (-N:1:0)*dS; % strain 

% Initial variables for Force-velocity experiment
p1 = zeros(N+1,1);
p2 = zeros(N+1,1);
p3 = zeros(N+1,1);
U_NR = 1;
% State variable vector concatenates p1, p2, p2, and U_NR
PU0 = [p1; p2; p3; U_NR];

% Set metabolite concentrations, 
MgATP = [8 4 2];
MgADP = 0; 
Pi    = 0; 

%% Parameters and computing initial error
% g0 = ones(1,12);
load g0

% moments and force
dr = +g0(12)*0.01; % Power-stroke Size; Units: um
kstiff1 = g0(13)*2500; 
kstiff2 = g0(14)*200; 
mu = 0.50; % viscosity

% Non-zero velocities
vel = (-Data(:,1)).*ML; % micron per sec

for k = [1 3]

  % Zero velocity:
  [t,PU] = ode15s(@dPUdT,[0 0.11],PU0,[],N,dS,MgATP(k),Pi,MgADP,g0);
  PU = PU(end,:);
  p1 = PU(1:1*N+1);
  p2 = PU(1*N+2:2*N+2);
  p3 = PU(2*N+3:3*N+3);
  p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
  p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
  p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
  Force(1,k) = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 );

  for j = 2:length(vel)
    j;
  %   Set the outer timestep based on space step:
    dt = dS/abs(vel(j));
    tend = 0.11/abs(vel(j)); % ending time of simulation
    Nstep = round(tend/dt);
    % simulate kinetics for 1/2 timestep
    [t,PU] = ode15s(@dPUdT,[0 dt/2],PU0,[],N,dS,MgATP(k),Pi,MgADP,g0);
    PU = PU(end,:); 
    for i = 1:(Nstep-1)
      % advection (sliding step)
      PU(1:1*N+0)     = PU(2:1*N+1); PU(N+1) = 0;
      PU(1*N+2:2*N+1) = PU(1*N+3:2*N+2); PU(2*N+2) = 0;
      PU(2*N+3:3*N+2) = PU(2*N+4:3*N+3); PU(3*N+3) = 0;
      % simulate kinetics for full step
      [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP(k),Pi,MgADP,g0);
      PU = PU(end,:); 
    end
    % final advection (sliding step)
    PU(1:1*N+0)     = 0.5*(PU(2:1*N+1) + PU(1:1*N+0));         PU(N+1) = 0.5*(0 + PU(N+1));
    PU(1*N+2:2*N+1) = 0.5*(PU(1*N+3:2*N+2) + PU(1*N+2:2*N+1)); PU(2*N+2) = 0.5*(0 + PU(2*N+2));
    PU(2*N+3:3*N+2) = 0.5*(PU(2*N+4:3*N+3) + PU(2*N+3:3*N+2)); PU(3*N+3) = 0.5*(0 + PU(3*N+3));
    % final 1/2 timestep for kinetics
    [t,PU] = ode15s(@dPUdT,[0 dt/2],PU,[],N,dS,MgATP(k),Pi,MgADP,g0);
    PU = PU(end,:);
  
    p1 = PU(1:1*N+1);
    p2 = PU(1*N+2:2*N+2);
    p3 = PU(2*N+3:3*N+3);
    p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
    p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
    Force(j,k) = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 ) + mu*vel(j);
    
  end

end

E0 = sum((Force(:,1)-Data(:,2)).^2) + ...
     sum((Force(:,3)-Data(:,4)).^2)
 
% Fmax (normalized.) versus [MgATP] (mM) from Ebus et al.(2001)
iso_data = ...
    [0.010     0.020      0.050     0.10      0.50        5.
     1.5925    1.6826     1.5898    1.4657    1.2884      0.99732];
MgATP_iso = iso_data(1,:);

for k = 1:length(MgATP_iso)
  % Zero velocity:
  [t,PU] = ode15s(@dPUdT,[0 1],PU0,[],N,dS,MgATP_iso(k),Pi,MgADP,g0);  
  PU = PU(end,:);
  p1 = PU(1:1*N+1);
  p2 = PU(1*N+2:2*N+2);
  p3 = PU(2*N+3:3*N+3);
  p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
  p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
  p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
  F_iso(k) = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 ) ;
end

F_data = iso_data(2,:).*57;
E0 = E0 + sum((F_iso-F_data).^2)

% Ktr data from Beard et al.
Ktr_mean = [37.7928 29.0 25.8033];

% initial state vector for Ktr exp
p1 = zeros(N+1,1);
p2 = zeros(N+1,1);
p3 = zeros(N+1,1);
U_NR = 1;
PU0 = [p1; p2; p3; U_NR];
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
  p3_0 = dS*sum(p3'); p3_1 = dS*sum((s+dr)'.*p3');
  F_active_ktr = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 ) ;
  
  Frel = F_active_ktr./F_active_ktr(end);
  Fexp = 1 - exp( -Tspan.*Ktr_mean(k) );
  E0 = E0 + 10*sum((Frel-Fexp).^2);
 
  Ktr(k) = 1/interp1(Frel,Tspan,1-exp(-1)); % time constant for Frel(1/Ktr) = 1-exp(-1)

end

E0 = E0 + 2*sum((Ktr-Ktr_mean).^2)

%% Search parameters to reduce error

gN = g0;


ii = [3 16; 5 17; 13 14; 10 11; 1 2; 6 15; 7 8; 9 12; 4 18;
      1 2; 3 4; 5 6; 7 8; 9 10; 11 12; 13 14; 15 18; 16 17; ];
  
% ii = [3 16; 5 17; 1 2; 6 15; 7 8; 9 12; 4 18;
%       1 2; 3 4; 5 6; 7 8; 11 12; 15 18; 16 17; ];

for n = 1:14
n

for l = 1:5
  gN(ii(n,:)) = g0(ii(n,:)).*(1 + 0.0025*randn(1,2));
%   gN([6 7 8]) = g0([6 7 8]).*(1 + 0.01*randn(1,3));
%   gN([3 16]) = g0([3 16]).*(1 + 0.05*randn(1,2));
%   gN([5 17]) = g0([5 17]).*(1 + 0.05*randn(1,2));
%   gN([13 14]) = g0([13 14]).*(1 + 0.02*randn(1,2));

  % moments and force
  dr = +gN(12)*0.01; % Power-stroke Size; Units: um
  kstiff1 = gN(13)*2500; 
  kstiff2 = gN(14)*200; 

  % initial state vector for Force-velocity exp
  p1 = zeros(N+1,1);
  p2 = zeros(N+1,1);
  p3 = zeros(N+1,1);
  U_NR = 1;
  PU0 = [p1; p2; p3; U_NR];

  for k = [1 3]
    % Zero velocity:
    [t,PU] = ode15s(@dPUdT,[0 0.11],PU0,[],N,dS,MgATP(k),Pi,MgADP,gN);
    PU = PU(end,:);
    p1 = PU(1:1*N+1);
    p2 = PU(1*N+2:2*N+2);
    p3 = PU(2*N+3:3*N+3);
    p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
    p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
    Force(1,k) = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 ) ;

    % Non-zero velocities
    for j = 2:length(vel)
    %   Set the outer timestep based on space step:
      dt = dS/abs(vel(j));
      tend = 0.11/abs(vel(j)); % ending time of simulation
      Nstep = round(tend/dt);
      % simulate kinetics for 1/2 timestep
      [t,PU] = ode15s(@dPUdT,[0 dt/2],PU0,[],N,dS,MgATP(k),Pi,MgADP,gN);
      PU = PU(end,:); 
      for i = 1:(Nstep-1)
        % advection (sliding step)
        PU(1:1*N+0)     = PU(2:1*N+1); PU(N+1) = 0;
        PU(1*N+2:2*N+1) = PU(1*N+3:2*N+2); PU(2*N+2) = 0;
        PU(2*N+3:3*N+2) = PU(2*N+4:3*N+3); PU(3*N+3) = 0;
        % simulate kinetics for full step
        [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP(k),Pi,MgADP,gN);
        PU = PU(end,:); 
      end
      % final advection (sliding step)
      PU(1:1*N+0)     = 0.5*(PU(2:1*N+1) + PU(1:1*N+0));         PU(N+1) = 0.5*(0 + PU(N+1));
      PU(1*N+2:2*N+1) = 0.5*(PU(1*N+3:2*N+2) + PU(1*N+2:2*N+1)); PU(2*N+2) = 0.5*(0 + PU(2*N+2));
      PU(2*N+3:3*N+2) = 0.5*(PU(2*N+4:3*N+3) + PU(2*N+3:3*N+2)); PU(3*N+3) = 0.5*(0 + PU(3*N+3));
      % final 1/2 timestep for kinetics
      [t,PU] = ode15s(@dPUdT,[0 dt/2],PU,[],N,dS,MgATP(k),Pi,MgADP,gN);
      PU = PU(end,:);
  
      p1 = PU(1:1*N+1);
      p2 = PU(1*N+2:2*N+2);
      p3 = PU(2*N+3:3*N+3);
      p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
      p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
      p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
      Force(j,k) = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 ) + mu*vel(j) ;
      
      if j==3
        p3_0;
      end
  
    end
  end
  
  EN = sum((Force(:,1)-Data(:,2)).^2) + ...
       sum((Force(:,3)-Data(:,4)).^2);

  for k = 1:length(MgATP_iso)
    % Zero velocity:
    [t,PU] = ode15s(@dPUdT,[0 1],PU0,[],N,dS,MgATP_iso(k),Pi,MgADP,gN);  
    PU = PU(end,:);
    p1 = PU(1:1*N+1);
    p2 = PU(1*N+2:2*N+2);
    p3 = PU(2*N+3:3*N+3);
    p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
    p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);
    F_iso(k) = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 ) ;
  end

  F_data = iso_data(2,:).*57;
  EN = EN + sum((F_iso-F_data).^2);
  
  % initial state vector for Ktr exp
  p1 = zeros(N+1,1);
  p2 = zeros(N+1,1);
  p3 = zeros(N+1,1);
  U_NR = 1;
  PU0 = [p1; p2; p3; U_NR];
  Tspan = [0:0.001:0.12];

  for k = 1:length(MgATP)

    % Zero velocity:
    [t,PU] = ode15s(@dPUdT,Tspan,PU0,[],N,dS,MgATP(k),Pi,MgADP,gN);  
    p1 = PU(:,1:1*N+1);
    p2 = PU(:,1*N+2:2*N+2);
    p3 = PU(:,2*N+3:3*N+3);
    p1_0 = dS*sum(p1'); p1_1 = dS*sum(s'.*p1');
    p2_0 = dS*sum(p2'); p2_1 = dS*sum(s'.*p2');
    p3_0 = dS*sum(p3'); p3_1 = dS*sum((s+dr)'.*p3');
    F_active_ktr = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 ) ;
  
    Frel = F_active_ktr./F_active_ktr(end);
    Fexp = 1 - exp( -Tspan.*Ktr_mean(k) );
    EN = EN + 10*sum((Frel-Fexp).^2);
    
    Ktr(k) = 1/interp1(Frel,Tspan,1-exp(-1)); % time constant for Frel(1/Ktr) = 1-exp(-1)

  end

  EN = EN + 2*sum((Ktr-Ktr_mean).^2);
 
  if EN < E0
    E0 = EN;
    g0 = gN;
    [ii(n,:) E0]
    save g0 g0
  end

end
end
