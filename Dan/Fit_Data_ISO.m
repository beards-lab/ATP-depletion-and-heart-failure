clear

%% Setting up problem
N = 20; % space (strain) discretization--number of grid points in half domain
Slim = 0.050; 
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
  p3_0 = dS*sum(p3); p3_1 = dS*sum((s-dr).*p3);
  F_iso(k) = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 ) ;
end

F_data = iso_data(2,:).*57;
E0 = sum((F_iso-F_data).^2)

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
  p3_0 = dS*sum(p3'); p3_1 = dS*sum((s-dr)'.*p3');
  F_active_ktr = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 ) ;
  
  Frel = F_active_ktr./F_active_ktr(end);
  Fexp = 1 - exp( -Tspan.*Ktr_mean(k) );
  E0 = E0 + 10*sum((Frel-Fexp).^2);
 
  Ktr(k) = 1/interp1(Frel,Tspan,1-exp(-1)); % time constant for Frel(1/Ktr) = 1-exp(-1)

end

E0 = E0 + 2*sum(abs(Ktr-Ktr_mean).^2)

%% Search parameters to reduce error

gN = g0;
% 1 2 3 4 5 10 11 15
  
% ii = [1 2;
%       3 4;
%       5 6;
%       7 8;
%       10 11;
%       15 1;
%       2 3;
%       4 5;
%       6 7;
%       8 10;
%       11 15;
%  	  ];
  
ii = [1 2;
      3 4;
      5 10;
      11 15;
      14 18;
 	  ];
  
for n = 1:4
n

for l = 1:5
  gN(ii(n,:)) = g0(ii(n,:)).*(1 + 0.05*randn(1,2));
%   gN([6 7 8]) = g0([6 7 8]).*(1 + 0.05*randn(1,3));
%   gN([9 12 13]) = g0([9 12 13]).*(1 + 0.10*randn(1,3));

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

  for k = 1:length(MgATP_iso)
    % Zero velocity:
    [t,PU] = ode15s(@dPUdT,[0 1],PU0,[],N,dS,MgATP_iso(k),Pi,MgADP,gN);  
    PU = PU(end,:);
    p1 = PU(1:1*N+1);
    p2 = PU(1*N+2:2*N+2);
    p3 = PU(2*N+3:3*N+3);
    p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
    p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    p3_0 = dS*sum(p3); p3_1 = dS*sum((s-dr).*p3);
    F_iso(k) = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 ) ;
  end

  F_data = iso_data(2,:).*57;
  EN = sum(abs(F_iso-F_data).^2);
  
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
    p3_0 = dS*sum(p3'); p3_1 = dS*sum((s-dr)'.*p3');
    F_active_ktr = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 ) ;
  
    Frel = F_active_ktr./F_active_ktr(end);
    Fexp = 1 - exp( -Tspan.*Ktr_mean(k) );
    EN = EN + 10*sum((Frel-Fexp).^2);
    
    Ktr(k) = 1/interp1(Frel,Tspan,1-exp(-1)); % time constant for Frel(1/Ktr) = 1-exp(-1)

  end

  EN = EN + 2*sum(abs(Ktr-Ktr_mean).^2);
 
  if EN < E0
    E0 = EN;
    g0 = gN;
    [ii(n,:) E0]
    save g0 g0
  end

end
end
