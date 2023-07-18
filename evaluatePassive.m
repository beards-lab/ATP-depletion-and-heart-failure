%% Evaluates passive model - script
% Expects input parameters in the workspace

% use simInit, simRamp, simRecover to run just parts of the code for faster
% optim

if ~(exist('simInit', 'var') && ~simInit)
    % setting up problem
    N = 100; % number of space step
    L = 0.60; % domain size (micros)
    ds = L / N;
    s  = 0:ds:(L-ds);
    a = zeros(1,N); % initial p.d.f of attached positions    
    
    % rampup duration (s)
    % rd = 1.00;
    dl = (2.4-1.6)/2; % delta L of half-sarcomere (um)
    % V = 40/2; % highest half-sarcomere velocity
    % k = 10; % series spring constant
    % Lo = 0.10;
    
    if ~exist('opt_mods', 'var')
        opt_mods = ones(1, 10);
    end
    % Model parameters
    r_a = opt_mods(1)*(100);
    r_d = opt_mods(2)*(0.01);
    % r_d = (1/50).*(s./0.2).^4;
    % r_d = (1/50).ones(1,N) + (1/25)*ones(1,N).(s>0.25);
    % plot(s,r_d)
    k1 = opt_mods(3)*0.1; % titin force constant
    kp = opt_mods(4)*200; % parallel nonlinear force constant

    
    gamma=opt_mods(5)*4; % titin exponent
    s0 = opt_mods(6)*4;
    beta = opt_mods(7)*2;
    
    % phi = opt_mods(10)*1;
    L0 = opt_mods(8)*0.8;
    s1 = opt_mods(9)*10;
    
    FtitFun = @(p_a, L)kp*(L-L0)^gamma; % nonlinear titin stiffness of unattached
    FattFun = @(p_a, a)k1*ds*sum((exp(s1*s)-1).*a); % Force of attached
      
    Tend = 5/r_a; % length of steady-state simulation - get rid of all transients
    [t,x] = ode15s(@dadt,[0 Tend],a,[],N,s, ds,r_a,r_d, beta, s0);
    a = x(end,:);
    
    
    % figure(1); plot(s,a); pause    
    % Tsim = 0;Fatts = 0;Ftits = 0; Fvs = 0;
    Fvs = 0; % viscous force set
    Ftits = 0; % TITIN passive force set
    Tsim = -2;

end % end init phase

if ~(exist('simRamp', 'var') && ~simRamp)
    X = [0 0]; % initial length and L1
    V = dl/rd; % highest half-sarcomere velocity
    p_a = ds*sum(a); % probability (fraction) of attached
    p_u = 1 - p_a; % probability (fraction) of unattached
    Ftit = FtitFun(p_a, L); % nonlinear titin stiffness of unattached
    Fatt = FattFun(p_a, a); % Force of attached
    p_as = [p_a];
    catts = [0];

    % ramp sretch time course
    Tend_ramp = rd; % length of ramp
    Fatts = [Fatt];
    Ls = 0;
    dt = ds/V; % numerical time step
    aN = zeros(1,N); % new updated a vector for upwind diff.

    for i = 1:Tend_ramp/dt 
    %   a0 = ds*sum(a); % 0th moment
    [t,x] = ode15s(@dadt,[0 dt],a,[],N,s, ds,r_a,r_d, beta, s0);
      a = x(end,:);
    
      % Fatt = ds*sum(s.*a);
      % F1 = k1*Fatt;
    %   V1 = V - F1/mu;
      % UPWIND differencing for sliding
      aN(1) = a(1) - (dt*V/ds)*a(1);
      for j = 2:length(s)
         aN(j) = a(j) - (dt*V/ds)*a(j) + (dt*V/ds)*a(j-1);
      end
      a = aN;
      % L = i*dt*V;
    %   plot(s,a); pause
    
      %  % Fv = 0;
      % [t,X] = ode15s(@dL1dT,[0 dt],X,[],V, ks/mu);
      % X = X(end,:);
      % L = X(1);
      % L1 = X(2);
      % Fv = ks*(L - L1); % dashpot viscous force
      Fv = 0;      
      L = L0 + i*dt*V;

      p_a = ds*sum(a); % probability (fraction) of attached
      p_u = 1 - p_a; % probability (fraction) of unattached
    Ftit = FtitFun(p_a, L); % nonlinear titin stiffness of unattached
    Fatt = FattFun(p_a, a); % Force of attached
    
    if Tsim(end)  > 0.8
      breakpointhere = true;
    end
    
      Tsim = [Tsim, i*dt];
      Fatts = [Fatts, Fatt]; % attached
      Ftits = [Ftits, Ftit]; % titin alon
      Fvs = [Fvs, Fv];
      Ls = [Ls, L];
      p_as = [p_as, p_a];
      catts = [catts, ds*sum(s.*a)]; % center of attached

    end
end % end ramp

% figure(0405);hold on;
% plot(s, a*ds);
% title('Attached distribution at the end of the ramp');

% steady time course after ramp
if rd ~= 0
    dt = rd/10;
else
    % do not simulate here at all, just defensive
    dt = 1;
end

if ~(exist('simRecover', 'var') && ~simRecover)
    Tend_relaxDt = min(rd*100, 300)/dt; % length of relaxation time
else
    % if not simRecover, simulate just a minor peak
    Tend_relaxDt = rd/dt; % length of relaxation time
end

    for i = 1:Tend_relaxDt
      [t,x] = ode15s(@dadt,[0 dt],a,[],N,s, ds,r_a,r_d, beta, s0);
      a = x(end,:);
    
      % [t,X] = ode15s(@dL1dT,[0 dt],X,[],0, ks/mu);
      % if Tsim(end) + dt > 4
      %     breakpointhere = true;
      % end
      % X = X(end,:);
      % L = X(1);
      % L1 = X(2);
      % Fv = ks*(L - L1); % dashpot viscous force
      Fv = 0;
    
      p_a = ds*sum(a); % probability (fraction) of attached
      p_u = 1 - p_a; % probability (fraction) of unattached
    Ftit = FtitFun(p_a, L); % nonlinear titin stiffness of unattached
    Fatt = FattFun(p_a, a); % Force of attached
      
      
      Tsim = [Tsim, Tend_ramp + i*dt];
      Fatts = [Fatts, Fatt];
      % the length stays constant
      Ftits = [Ftits, Ftit];
      % TODO fix the viscous dynamics
      Fvs = [Fvs, Fv];
      Ls = [Ls, Ls(end)];
      p_as = [p_as, p_a];      
      catts = [catts, ds*sum(s.*a)]; % center of attached
      
    %   plot(s,a); pause 
    end
%% calc error and plot

% interp to data points
Tsims = Tsim + 2; % time shifted
Ftot = Fatts + Ftits + Fvs;
if ~(exist('calcInterpE', 'var') && ~calcInterpE)
    Ftot_int = interp1(Tsims, Ftot, datatable(:,1));    
    E = (Ftot_int - datatable(:,3)).^2;
    E(isnan(E)) = 0; % zero outside bounds
    Es = 1e3*sum(E)/length(E);
end


if plotEach
    % figure(2); clf; hold on;    
    % cla;
    hold on;
    plot(datatable(:, 1), datatable(:, 3), 'o-', 'Linewidth', 1.5);
    % plot(Tsims,Ftot, '-', 'Linewidth', 1)
    plot(Tsims, Ls, Tsims, Fatts, Tsims, Ftits, Tsims, Fvs);

    % plot interpolated to data
    plot(datatable(:, 1),Ftot_int, 'x-', 'Linewidth', 1.5)
    % yl = ylim;
    % plot(datatable(:, 1), E, '--');
    plot(Tsims, p_as, '--');
    plot(Tsims, catts, ':', 'Linewidth', 1.5);
    xlim([0, Tsims(end)])
    % ylim(yl);
    ylim('auto');
    legend('Data', 'Ls', 'Fbinded', 'Ftitin passive', 'Viscous', 'Ftot passive', 'Attached (1)', 'Center of attached', 'Location', 'NorthEast') 
end