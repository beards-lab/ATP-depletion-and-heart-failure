% clear;

% setting up problem
N = 50; % number of space step
L = 0.60; % domain size (micros)
ds = L / N;
s  = 0:ds:(L-ds);
a = zeros(1,N); % initial p.d.f of attached positions
X = [0 0]; % initial length and L1


% rampup duration (s)
% rd = 1.00;
dl = (2.4-1.6)/2; % delta L of half-sarcomere (um)
V = dl/rd; % highest half-sarcomere velocity
% V = 40/2; % highest half-sarcomere velocity
% k = 10; % series spring constant
% Lo = 0.10;

if ~exist('opt_mods', 'var')
    opt_mods = ones(1, 10);
end
% Model parameters
r_a = opt_mods(1);
r_d = opt_mods(2)*(1).*ones(1,N) ;
% r_d = (1/50).*(s./0.2).^4;
% r_d = (1/50).ones(1,N) + (1/25)*ones(1,N).(s>0.25);
% plot(s,r_d)
mu = opt_mods(3)*10;
ks = opt_mods(4)*4;
k1 = opt_mods(5)*0.2;
c=opt_mods(6)*13.1; % titin linear koefficient
gamma=opt_mods(7)*4.7; % titin exponent
alpha1 = opt_mods(8)*10;
e = opt_mods(9)*2;
  
Tend = 60; % length of steady-state simulation - get rid of all transients
[t,x] = ode15s(@dadt,[0 Tend],a,[],N,ds,r_d,r_a, e);
a = x(end,:);

% figure(1); plot(s,a); pause

% ramp sretch time course
Tend_ramp = rd; % length of ramp
a1 = ds*sum(s.*a);
Fsim = 100*a1;
Ls = 0;
Fvs = 0; % viscous force set
Ftits = 0; % TITIN passive force set
Tsim = 0;
dt = ds/V; % numerical time step
aN = zeros(1,N); % new updated a vector for upwind diff.

for i = 1:Tend_ramp/dt 
%   a0 = ds*sum(a); % 0th moment
  [t,x] = ode15s(@dadt,[0 dt],a,[],N,ds,r_d,r_a, e);
  a = x(end,:);

  a1 = ds*sum(s.*a);
  F1 = k1*a1;
%   V1 = V - F1/mu;
  % UPWIND differencing for sliding
  aN(1) = a(1) - (dt*V/ds)*a(1);
  for j = 2:length(s)
     aN(j) = a(j) - (dt*V/ds)*a(j) + (dt*V/ds)*a(j-1);
  end
  a = aN;
  L = i*dt*V;
%   plot(s,a); pause

   % Fv = 0;
  [t,X] = ode15s(@dL1dT,[0 dt],X,[],V, ks/mu);
  X = X(end,:);
  L = X(1);
  L1 = X(2);
  Fv = ks*(L - L1); % dashpot viscous force
  
  Ftit = c*(L*2)^gamma; % nonlinear titin stiffness

  a1 = ds*sum((exp(alpha1*s)-1).*a); 
  Tsim = [Tsim, i*dt];
  Fsim = [Fsim, k1*a1];
  Ftits = [Ftits, Ftit];
  Fvs = [Fvs, Fv];
  Ls = [Ls, L];

end
% clf;hold on;
% plot(Tsim, Fsim + Ftits + Fvs, 'Linewidth', 2)
% plot(Tsim, Ls, Tsim, Fsim, Tsim, Ftits, Tsim, Fvs);
%
% plot(s,a); pause
% steady time course after ramp
Tend_relax = rd*100; % length of relaxation time
dt = rd/10;
for i = 1:Tend_relax/dt 
  [t,x] = ode15s(@dadt,[0 dt],a,[],N,ds,r_d,r_a, e);
  a = x(end,:);
  a1 = ds*sum((exp(alpha1*s)-1).*a); 

  [t,X] = ode15s(@dL1dT,[0 dt],X,[],0, ks/mu);
  X = X(end,:);
  L = X(1);
  L1 = X(2);
  Fv = ks*(L - L1); % dashpot viscous force
  
  Tsim = [Tsim, Tend_ramp + i*dt];
  Fsim = [Fsim, k1*a1];
  % the length stays constant
  Ftits = [Ftits, Ftits(end)];
  % TODO fix the viscous dynamics
  Fvs = [Fvs, Fv];
  Ls = [Ls, Ls(end)];
%   plot(s,a); pause 
end
%% calc error and plot

% interp to data points
Tsims = Tsim + 2; % time shifted
Ftot = Fsim + Ftits + Fvs;
Ftot_int = interp1(Tsims, Ftot, datatable(:,1));    
E = (Ftot_int - datatable(:,3)).^2;
E(isnan(E)) = 0; % zero outside bounds
Es = sum(E);


if plotEach
    % figure(2); clf; hold on;    
    cla;hold on;
    plot(datatable(:, 1), datatable(:, 3), 'o-', 'Linewidth', 1.5);
    % plot(Tsims,Ftot, '-', 'Linewidth', 1)
    plot(Tsims, Ls, Tsims, Fsim, Tsims, Ftits, Tsims, Fvs);

    % plot interpolated to data
    plot(datatable(:, 1),Ftot_int, 'x-', 'Linewidth', 1.5)
    plot(datatable(:, 1), E, '--');
    xlim([0, datatable(end, 1)])
    legend('Data', 'Ls', 'Fbinded', 'Ftitin passive', 'Viscous', 'Ftot passive', 'Error', 'Location', 'NorthWest') 
end