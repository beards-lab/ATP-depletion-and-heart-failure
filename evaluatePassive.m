% setting up problem
N = 50; % number of space step
L = 0.60; % domain size (micros)
ds = L / N;
s  = 0:ds:(L-ds);
a = zeros(1,N); % initial p.d.f of attached positions

s_rel0=1.62;
gamma=4.7;
c=13.1;
    
% rampup duration (s)
rd = 0.100;
dl = (2.4-1.6)/2; % delta L
V = dl/rd; % highest half-sarcomere velocity

r_a = 1.0;
% tau_d = 50 ./ (1 + (s./0.10).^4); % detatchment time constant
% tau_d = 250.*exp(-20*s);
% r_d   = 500*(exp( 10*(s).^6 ) - 1);
% r_d = (1/50).*ones(1,N);
% r_d = (1/50).*ones(1,N).*s.^6 ./ (s.^6 + 0.25.^6);
k_rd = 1/500;
r_d = k_rd.*ones(1,N).*(s./0.2).^4;
% plot(s,r_d)
% mu = 100;
k1 = 20;

Tend = 60; % length of steady-state simulation
[t,x] = ode15s(@dadt,[0 Tend],a,[],N,ds,r_d,r_a);
a = x(end,:);

% ramp sretch time course
Tend_ramp = rd; % length of ramp
a1 = ds*sum(s.*a);
Fsim = 100*a1;
Tsim = 0;
dt = ds/V; % numerical time step
aN = zeros(1,N); % new updated a vector for upwind diff.
for i = 1:Tend_ramp/dt 
%   a0 = ds*sum(a); % 0th moment
  [t,x] = ode15s(@dadt,[0 dt],a,[],N,ds,r_d,r_a);
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
%   plot(s,a); pause
    
  a1 = ds*sum(s.*a); 
  Ft = c * a1^gamma;
  Tsim = [Tsim, i*dt];
  Fsim = [Fsim, Ft];

end

% plot(s,a); pause
% steady time course after ramp
Tend_relax = 200; % length of relaxation time
dt = 1;
for i = 1:Tend_relax/dt 
  [t,x] = ode15s(@dadt,[0 dt],a,[],N,ds,r_d,r_a);
  a = x(end,:);
%   a1 = ds*sum((exp(4*s/0.4)-1).*a); 
  a1 = ds*sum(s.*a); 
  Ft = c * a1^gamma;
  Tsim = [Tsim, Tend_ramp + i*dt];
  Fsim = [Fsim, Ft];
%   plot(s,a); pause
  
end
Tsim = Tsim + 2;
figure(3); plot(Tsim*1000,Fsim)