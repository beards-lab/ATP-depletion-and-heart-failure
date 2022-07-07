clear

LoadData;


    

%% Setting up problem

N = 20; % space (strain) discretization--number of grid points in half domain
Slim = 0.075; 
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

ML = 1.1; % half sarcomere length (microns)

%% Parameters and computing initial error
% g0 = ones(1,12);
load g0
% 
    % moments and force
    dr = 0.01; % Power-stroke Size; Units: um
    kstiff1 = g0(13)*1500; 
    kstiff2 = g0(14)*10000;     
    
vel = (-Data_ATP(:,1)).*ML; % micron per sec    

for k = [1 3]

  % Zero velocity:
%   [t,PU] = ode15s(@dPUdT,[0 1],PU0,[],N,dS,MgATP(k),Pi,MgADP,g0);
%   PU = PU(end,:);
%   p1 = PU(1:1*N+1);
%   p2 = PU(1*N+2:2*N+2);
%   p3 = PU(2*N+3:3*N+3);
%   p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
%   p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
%   p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
%   F_active(1,k) = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 );
%   
%   F_active2 = evaluateModel(0, 1, PU0,N,dS,MgATP(k),Pi,MgADP,g0);
%   disp(F_active(1,k) - F_active2)
% F_active(1,k) = evaluateModel([0],MgATP(k),Pi,MgADP,g0);
% 
%   for j = 2:length(vel)
%     j
%   %   Set the outer timestep based on space step:
%     dt = dS/abs(vel(j));
%     tend = 0.20/abs(vel(j)); % ending time of simulation
%     Nstep = round(tend/dt);
%     % simulate kinetics for 1/2 timestep
%     [t,PU] = ode15s(@dPUdT,[0 dt/2],PU0,[],N,dS,MgATP(k),Pi,MgADP,g0);
%     PU = PU(end,:); 
%     for i = 1:(Nstep-1)
%       % advection (sliding step)
%       PU(1:1*N+0)     = PU(2:1*N+1); PU(N+1) = 0;
%       PU(1*N+2:2*N+1) = PU(1*N+3:2*N+2); PU(2*N+2) = 0;
%       PU(2*N+3:3*N+2) = PU(2*N+4:3*N+3); PU(3*N+3) = 0;
%       % simulate kinetics for full step
%       [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP(k),Pi,MgADP,g0);
%       PU = PU(end,:); 
%     end
%     % final advection (sliding step)
%     PU(1:1*N+0)     = 0.5*(PU(2:1*N+1) + PU(1:1*N+0));         PU(N+1) = 0.5*(0 + PU(N+1));
%     PU(1*N+2:2*N+1) = 0.5*(PU(1*N+3:2*N+2) + PU(1*N+2:2*N+1)); PU(2*N+2) = 0.5*(0 + PU(2*N+2));
%     PU(2*N+3:3*N+2) = 0.5*(PU(2*N+4:3*N+3) + PU(2*N+3:3*N+2)); PU(3*N+3) = 0.5*(0 + PU(3*N+3));
%     % final 1/2 timestep for kinetics
%     [t,PU] = ode15s(@dPUdT,[0 dt/2],PU,[],N,dS,MgATP(k),Pi,MgADP,g0);
%     PU = PU(end,:);
%   
%     p1 = PU(1:1*N+1);
%     p2 = PU(1*N+2:2*N+2);
%     p3 = PU(2*N+3:3*N+3);
%     p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
%     p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
%     p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
%     F_active(j,k) = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 ) ;
%     
%   end
    for j = 1:length(vel)
        F_active(j,k) = evaluateModel(vel(j), 1, MgATP(k),Pi,MgADP,g0);
    end
end



E0 = sum(abs(F_active(:,1)-Data_ATP(:,2)).^2) + ...
     sum(abs(F_active(:,3)-Data_ATP(:,4)).^2);
 
%%
for k = 1:length(MgATP_iso)
  % Zero velocity:
%   [t,PU] = ode15s(@dPUdT,[0 1],PU0,[],N,dS,MgATP_iso(k),Pi,MgADP,g0);  
%   PU = PU(end,:);
%   p1 = PU(1:1*N+1);
%   p2 = PU(1*N+2:2*N+2);
%   p3 = PU(2*N+3:3*N+3);
%   p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
%   p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
%   p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
%   F_iso(k) = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 ) ;
  F_iso(k) = evaluateModel(0, 1, MgATP_iso(k),Pi,MgADP,g0);
end


E0 = E0 + sum(abs(F_iso-F_data).^2);
%%


% initial state vector for Ktr exp
p1 = zeros(N+1,1);
p2 = zeros(N+1,1);
p3 = zeros(N+1,1);
U_NR = 1;
PU0 = [p1; p2; p3; U_NR];
Tspan = [0:0.001:0.12];

for k = 1:length(MgATP)

%   % Zero velocity:
%   [t,PU] = ode15s(@dPUdT,Tspan,PU0,[],N,dS,MgATP(k),Pi,MgADP,g0);  
% %   PU = PU(end,:);
%   p1 = PU(:,1:1*N+1);
%   p2 = PU(:,1*N+2:2*N+2);
%   p3 = PU(:,2*N+3:3*N+3);
%   p1_0 = dS*sum(p1'); p1_1 = dS*sum(s'.*p1');
%   p2_0 = dS*sum(p2'); p2_1 = dS*sum(s'.*p2');
%   p3_0 = dS*sum(p3'); p3_1 = dS*sum(s'.*p3');
%   F_active_ktr2 = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 ) ;
%   
  % Tspan array returns F active array
  F_active_ktr = evaluateModel(0, Tspan, MgATP(k),Pi,MgADP,g0);
  % cehck
%   sum(F_active_ktr2 - F_active_ktr)
  
  Frel = F_active_ktr./F_active_ktr(end);
  Ktr(k) = 1/interp1(Frel,Tspan,1-exp(-1)); % time constant for Frel(1/Ktr) = 1-exp(-1)

end

E0 = E0 + sum(abs(Ktr-Ktr_mean).^2)

%% Search parameters to reduce error

gN = g0;
ii = [1 2 3;
	  4 5 6;
      7 8 9;
      10 11 12;
      13 14 15;
      2 3 4;
      5 6 7;
      8 9 10;
      11 12 13;
      14 15 1;
      3 4 5;
      6 7 8;
      9 10 11;
      12 13 14;
      15 1 2;
 	  ];
  
% ii = [1 2 3;
%       4 5 6;
%       7 8 10;
%       11 12 13;
%       14 15 1;
%       2 3 4;
%       5 6 7;
%       8 10 11;
%       12 13 14;
%       15 1 2;
%       3 4 5;
%       6 7 8;
%       10 11 12;
%       13 14 15;
%       ];

for n = 1:15
n

for l = 1:5
  gN(ii(n,:)) = g0(ii(n,:)).*(1 + 0.05*randn(1,3));
%   gN([11 12]) = g0([11 12]).*(1 + 0.02*randn(1,2));

  % moments and force
  kstiff1 = gN(13)*1500; 
  kstiff2 = gN(14)*10000; 

  % initial state vector for Force-velocity exp
  p1 = zeros(N+1,1);
  p2 = zeros(N+1,1);
  p3 = zeros(N+1,1);
  U_NR = 1;
  PU0 = [p1; p2; p3; U_NR];

  for k = [1 3]
    % Zero velocity:
    [t,PU] = ode15s(@dPUdT,[0 1],PU0,[],N,dS,MgATP(k),Pi,MgADP,gN);
    PU = PU(end,:);
    p1 = PU(1:1*N+1);
    p2 = PU(1*N+2:2*N+2);
    p3 = PU(2*N+3:3*N+3);
    p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
    p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
    F_active(1,k) = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 ) ;

    % Non-zero velocities
    for j = 2:length(vel)
    %   Set the outer timestep based on space step:
      dt = dS/abs(vel(j));
      tend = 0.20/abs(vel(j)); % ending time of simulation
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
      p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
      F_active(j,k) = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 ) ;
      
      if j==3
        p3_0;
      end
  
    end
  end
  
  EN = sum(abs(F_active(:,1)-Data_ATP(:,2)).^2) + ...
       sum(abs(F_active(:,3)-Data_ATP(:,4)).^2);
%   EN = sum(abs(F_active(:,1)-Data(:,2)).^2);

  for k = 1:length(MgATP_iso)
    % Zero velocity:
    [t,PU] = ode15s(@dPUdT,[0 1],PU0,[],N,dS,MgATP_iso(k),Pi,MgADP,gN);  
    PU = PU(end,:);
    p1 = PU(1:1*N+1);
    p2 = PU(1*N+2:2*N+2);
    p3 = PU(2*N+3:3*N+3);
    p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
    p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
    F_iso(k) = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 ) ;
  end

  F_data = iso_data(2,:).*57;
  EN = EN + sum(abs(F_iso-F_data).^2);
  
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
    p3_0 = dS*sum(p3'); p3_1 = dS*sum(s'.*p3');
    F_active_ktr = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 ) ;
  
    Frel = F_active_ktr./F_active_ktr(end);
    Ktr(k) = 1/interp1(Frel,Tspan,1-exp(-1)); % time constant for Frel(1/Ktr) = 1-exp(-1)

  end

  EN = EN + sum(abs(Ktr-Ktr_mean).^2);
  %%
  fcn = @dPUdT;
  EN1 = evaluateProblem(fcn, gN)
 %%
  if EN < E0
    E0 = EN;
    g0 = gN;
    [n l E0]
    save g0 g0
  end

end
end
