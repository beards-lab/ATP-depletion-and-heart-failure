clear

% setting up problem
N = 100; % number of space step
L = 0.60; % domain size (micros)
ds = L / N;
s  = 0:ds:(L-ds);
a = zeros(1,N); % initial p.d.f of attached positions
X = [0 0]; % initial length and L1

% V = 1*0.40/100; % slowest half-sarcomere velocity
% V = 1000*0.40/100; % slowest half-sarcomere velocity

V = 1000*0.40/100; % slowest half-sarcomere velocity

r_a = 100;
r_d = 0.01 ;

k1 = 0.1; % titin force constant
kp = 200; % parallel nonlinear force constant

% Obtain initial steady state
Tend = 5./r_a; % length of steady-state simulation
[t,x] = ode15s(@dadt,[0 Tend],a,[],N,ds,r_d,r_a);
a = x(end,:);

% ramp sretch time course
Tend_ramp = 0.4/V; % length of ramp
a1 = ds*sum(s.*a);
Fsim = k1*a1;
Lsim = 0.8; % length of 1/2 sarcomere
Tsim = 2;
dt = ds/V; % numerical time step
for i = 1:Tend_ramp/dt 
%   a0 = ds*sum(a); % 0th moment
  [t,x] = ode15s(@dadt,[0 dt],a,[],N,ds,r_d,r_a);
  a = x(end,:);

  % UPWIND differencing for sliding
  aN(1) = a(1) - (dt*V/ds)*a(1);
  for j = 2:length(s)
     aN(j) = a(j) - (dt*V/ds)*a(j) + (dt*V/ds)*a(j-1);
  end
  a = aN;

  a1 = ds*sum((exp(s/0.10)-1).*a); 
  Ftitin = k1*a1;
  L = 0.8+i*dt*V;
  Fp     = kp*(L-0.8)^4;

  Tsim = [Tsim, 2+i*dt];
  Fsim = [Fsim, Ftitin+Fp];
  Lsim = [Lsim, L];

end
figure(2); clf; hold on; plot(s,a); % plotting the a distribution after ramp

% steady time course after ramp
Tend_relax = 100; % length of relaxation time
dt = 0.1;
for i = 1:Tend_relax/dt 
  [t,x] = ode15s(@dadt,[0 dt],a,[],N,ds,r_d,r_a);
  a = x(end,:);

  a1 = ds*sum((exp(s/0.10)-1).*a); 
  Ftitin = k1*a1;
  L = 1.2;
  Fp     = kp*(L-0.8)^4;

  Tsim = [Tsim, 2+ Tend_ramp + i*dt];
  Fsim = [Fsim, Ftitin+Fp ];
  Lsim = [Lsim, L];

end

figure(1); plot(Tsim, Fsim,'g','LineWidth',2)
figure(2); plot(s,a)
figure(3); semilogx(Tsim,Fsim,'r','LineWidth',2); 
figure(4); semilogy(Tsim,Fsim);

