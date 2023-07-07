% clear

% setting up problem
N = 50; % number of space step
L = 0.60; % domain size (micros)
ds = L / N;
s  = 0:ds:(L-ds);
a = zeros(1,N); % initial p.d.f of attached positions
X = [0 0]; % initial length and L1

V = 0.008/2; % slowest half-sarcomere velocity - 100s ramp
V = 0.8/2; % mid half-sarcomere velocity - 1s ramp
% V = 40/2; % highest half-sarcomere velocity - 20ms ramp
% k = 10; % series spring constant
% Lo = 0.10;

r_a = (1/250);
r_d = (1/25).*ones(1,N) ;

% plot(s,r_d)
mu = 3;
ks = 6*0;
k1 = 2;

Tend = 6000; % length of steady-state simulation
[t,x] = ode15s(@dadt,[0 Tend],a,[],N,ds,r_d,r_a);
a = x(end,:);

figure(1); plot(s,a);% pause

% ramp sretch time course
Tend_ramp = 0.4/V; % length of ramp
a1 = ds*sum(s.*a);
Fsim = 100*a1;
Tsim = 0;
dt = ds/V; % numerical time step
aN = zeros(1,N); % new updated a vector for upwind diff.
p_a = 0;
c_a = 0;
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

  a1 = ds*sum((exp(4*s/0.4)-1).*a); 
  Ftitin = k1*a1;
  [t,X] = ode15s(@dL1dT,[0 dt],X,[],V);
  X = X(end,:);
  L = X(1);
  L1 = X(2);
  Fs = ks*(L - L1);

  Tsim = [Tsim, i*dt];
  Fsim = [Fsim, Ftitin + Fs];
  p_a = [p_a ds*sum(a)];
  c_a = [c_a ds*sum(a.*s)];

end

% plot(s,a); pause
% steady time course after ramp
Tend_relax = 100; % length of relaxation time
dt = 0.1;
for i = 1:Tend_relax/dt 
  [t,x] = ode15s(@dadt,[0 dt],a,[],N,ds,r_d,r_a);
  a = x(end,:);
%   plot(s,a); pause

  a1 = ds*sum((exp(4*s/0.4)-1).*a); 
  Ftitin = k1*a1;
  [t,X] = ode15s(@dL1dT,[0 dt],X,[],0);
  X = X(end,:);
  L = X(1);
  L1 = X(2);
  Fs = ks*(L - L1);
  
  Tsim = [Tsim, Tend_ramp + i*dt];
  Fsim = [Fsim, Ftitin + Fs];
  p_a = [p_a ds*sum(a)];
  c_a = [c_a ds*sum(a.*s)];  

end
%%
figure(1); clf;plot(s,a);title('Profile of attached');xlabel('distance from attached');ylabel('Attached fraction');
figure(2); plot(Tsim,Fsim, Tsim, p_a, Tsim, c_a);xlabel('time'); legend('Force', 'Fracrtion attached', 'center of attached');
figure(3); semilogy(Tsim,Fsim, Tsim, p_a, Tsim, c_a);xlabel('time'); legend('Force', 'Fracrtion attached', 'center of attached');
