function [Etot, E] = evaluateProblem(fcn, g, drawPlots, evalParts)


if nargin < 4
    evalParts = [1, 1, 1, 1];
end
E = zeros(4, 1);

LoadData;

%% Set up environment
MgADP = 0; 
Pi    = 0;
F_active_0 = zeros(length(vel), length(MgATP));
t_sl0 = 0.11; % time at which the SL = 2.2 - time to stop the experiment
t_ss = 0.3; %% steady state time


%% force x velocity
F_active = F_active_0;
if evalParts(1)
% save from previous run
% FAb - F_active(j,k)    
opts = struct('N', 20, 'Slim', 0.06, 'PlotProbsOnFig', 0);
for k = [1 2 3]
    for j = 1:length(vel)
        F_active(j,k) = evaluateModel(fcn, vel(j), t_sl0, MgATP(k),Pi,MgADP,g, opts);
    end
end

% normalize by number of data points
E(1) = (sum(abs(F_active(:,1)-Data_ATP(:,2)).^2) + ...
       sum(abs(F_active(:,2)-Data_ATP(:,3)).^3) + ...
       sum(abs(F_active(:,3)-Data_ATP(:,4)).^2))/size(Data_ATP, 1)/(size(Data_ATP, 2) -1);
   
end

   %% linearized method for vmax
   k = 3; % atp of 2
   if F_active == F_active_0
       % the first part has not been run, we have to evaluate additionally
       opts = struct('N', 20, 'Slim', 0.06, 'PlotProbsOnFig', 0);
       F_active(end-1, k) = evaluateModel(fcn, vel(end-1), t_sl0, MgATP(k),Pi,MgADP,g, opts);
       F_active(end, k) = evaluateModel(fcn, vel(end), t_sl0, MgATP(k),Pi,MgADP,g, opts);
   end
   % take the last two
   
   % dy/dx = y/x, thus y = x*dy/dx
   dx = F_active(end-1, k) - F_active(end, k);
   dy = abs(vel(end)) - abs(vel(end-1));
   v_f0 = F_active(end-1, k)*dy/dx + abs(vel(end-1));

   dx_data = Data_ATP(end-1, k+1) - Data_ATP(end, k+1);
   dy_data = abs(vel(end)) - abs(vel(end-1));
   v_f0_data = Data_ATP(end-1, k+1)*dy/dx + abs(vel(end-1));
%    [v_f0 v_f0_data (v_f0-v_f0_data).^2]


%% force x iso MgATP

F_iso = zeros(1, length(MgATP_iso));
if evalParts(2)
    opts = struct('N', 10, 'Slim', 0.04, 'PlotProbsOnFig', 0);

    for k = 1:length(MgATP_iso)
      % Zero velocity:
      F_iso(k) = evaluateModel(fcn, 0, t_ss, MgATP_iso(k),Pi,MgADP,g, opts);
    end
    
    % weight per data point
    E(2) = (sum(abs(F_iso-F_data).^2))/length(MgATP_iso);
end
%% time constant of MgATP
Tspan = [0:0.01:0.2];
Frel = zeros(length(MgATP), length(Tspan));
Ktr = zeros(1, length(MgATP));

if evalParts(3)
opts = struct('N', 20, 'Slim', 0.02, 'PlotProbsOnFig', 0);
for k = 1:length(MgATP)

  % Tspan array returns F active array
  F_active_ktr = evaluateModel(fcn,0, Tspan, MgATP(k),Pi,MgADP,g);
  
  if F_active_ktr(end) == 0
      Ktr(k) = 0;
      continue;
  end
  
  Frel(k, :) = F_active_ktr./F_active_ktr(end);
  
  % get the time constant
  Ktr(k) = 1/interp1(Frel(k, :),Tspan,1-exp(-1)); % time constant for Frel(1/Ktr) = 1-exp(-1)
end

E(3) = (sum(abs(Ktr-Ktr_mean).^2))/length(MgATP);
end
%% Zero-force velocity
K_m = 0;

if evalParts(4)
%%
% have to search for velocity, which gives us a zero force
% currently the target is about 6 um/s
MgATP_vmax = 2;
targetVel = 6;

tol = 1e-2;

% scaling of the error function
errScale = 1e2;
% search for maximal velocity, that gives us fa = 0
i = 0;
e = Inf;
step = targetVel*2;
% the loop will search up to step*2
v_max = step;
maxIter = 5;
% real tolerance achievable:
realTol = targetVel*2*(1/2)^maxIter;

estimateVmax = 0;
if ~estimateVmax
    
    
    % debug only
    range = [v_max*(1/2)^(maxIter), v_max + step*(1 -(1/2)^(maxIter))];
    opts = struct('N', 30, 'Slim', 0.06, 'PlotProbsOnFig', 0);
    while abs(e) > tol && i < maxIter
        i = i + 1;
        e = evaluateModel(fcn, -v_max, t_ss, MgATP_vmax,Pi,MgADP,g ,opts);
%         e = 10;
    %     if i == 0 && e > 0 
    %         % not found in current range, cancelling the search
    %         break;
    %     end
        step = abs(step)/2;
        if e > 0
            % we havent reached the pivot
            v_max = v_max + step;
        else
            % we went too far, reduce the speed
            v_max = v_max - step;
        end
    end

else
    v_max = v_f0;
    e = 0;
    targetVel = v_f0_data;
end
% v_max
% e
% v_max - e

errScale = 1;
% disp("v_max - v_f0: " + num2str(v_max - v_f0))

%%
% if abs(e) > tol 
    % we have not found a zero force up till max speed of step*2
%     E(4) = (abs(e) > tol)*(exp(errScale*abs(e)-1e-4) - 1);
% else
    % we have found a zero force
    E(4) = 0*(v_max - targetVel).^2;
%%
    % search for MgATP, that gives us fa = 0 at 1/2 of the max velocity
    i = 0;
    e = Inf;
    step = MgATP_vmax/4;
    K_m = step;
    maxIter = 7;
    % real tolerance
    tol_min = step*step*2*(1/2)^maxIter;
    tol = 0.01;
    % debug only
    range = [K_m*(1/2)^(maxIter), K_m + step*(1 -(1/2)^(maxIter))];

    while abs(e) > tol && i < maxIter
        i = i + 1;
        opts = struct('N', 20, 'Slim', 0.06, 'PlotProbsOnFig', 0);
        e = evaluateModel(fcn, -v_max/2, t_ss, K_m,Pi,MgADP,g, opts);
% e = 11;
        step = abs(step)/2;
        if e > 0
            % Km too high, reduce
            K_m = K_m - step;
        else
            % we went too far, increase again
            K_m = K_m + step;
        end
    end


    % f = @(MgATP_halfv)(k_m_ADP0(1) > MgATP_halfv)*(k_m_ADP0(1)/k_m_ADP0(1)*scale - MgATP_halfv/k_m_ADP0(1)*scale)^2 ...
    %     + (k_m_ADP0(2) < MgATP_halfv)*(k_m_ADP0(2)/k_m_ADP0(2)*scale - MgATP_halfv/k_m_ADP0(2)*scale)^2;

    % divide the errscale by 2 so we are on right track
%     E(5) = (abs(e) > tol)*(exp(errScale/2*abs(e)-tol) - 1) ...
%         + (k_m_ADP0(1) > K_m)*(k_m_ADP0(1)/k_m_ADP0(1)*errScale - K_m/k_m_ADP0(1)*errScale)^2 ...
%         + (k_m_ADP0(2) < K_m)*(k_m_ADP0(2)/k_m_ADP0(2)*errScale - K_m/k_m_ADP0(2)*errScale)^2;

    % take the 0.07 as an average of 0.04 and 0.1
    E(5) = 100*(K_m/0.07 - 1).^2;

% end

end
%% Return

penalty = sum(max(0, -g))*1000;
E(6) = penalty;
Etot = sum(E);



%% plot?

if ~drawPlots
    return;
end

disp("With Km of " + num2str(K_m) + " at max velosity of " + v_max)
%% Force velocity
figure(101); clf; axes('position',[0.1 0.6 0.35 0.35]); hold on;
% figure(102);hold on;
plot(F_active(:,1), -vel./ML,'b-','linewidth',1.5);
plot(F_active(:,2), -vel./ML,'g-','linewidth',1.5);
plot(F_active(:,3), -vel./ML,'r-','linewidth',1.5);
ylabel('Velocity (ML/s)','interpreter','latex','fontsize',16);
xlabel('Force (kPa)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14); 
axis([0 65 0 6]);
title("With Km of " + num2str(K_m) )
box on;
plot(Data_ATP(:,2),Data_ATP(:,1),'bo','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
plot(Data_ATP(:,3),Data_ATP(:,1),'go','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
plot(Data_ATP(:,4),Data_ATP(:,1),'ro','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
ldg = legend('8','4','2 mM');
title(ldg,'[MgATP]');

% MgATP
% figure(104); clf; 
axes('position',[0.6 0.6 0.35 0.35]); 
title("and V_max of " + num2str(v_max));
semilogx(MgATP_iso,F_iso,'b-','linewidth',1.5); hold on;
semilogx(MgATP_iso,F_data,'ko','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1]);
xlabel('[MgATP] (mM)','interpreter','latex','fontsize',16);
ylabel('Max. Force (kPa)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'ylim',[0 100]); 
box on;

% KTR
% figure(106); clf; 
axes('position',[0.1 0.13 0.85 0.33]); hold on;
plot(Tspan,Frel(3, :),'r-','linewidth',1.5);
plot(Tspan,Frel(2, :),'g-','linewidth',1.5);
plot(Tspan,Frel(1, :),'b-','linewidth',1.5);
% plot(Tspan,F_active,'linewidth',1.5);
xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
ylabel('Force (rel.)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'ylim',[0 1.1]);  box on;
ldg = legend('2','4','8 mM','location','northwest');
title(ldg,'[MgATP]');
axes('position',[0.6 0.2 0.3 0.2]); hold on;
plot(MgATP,Ktr,'k.-','linewidth',1.5);
errorbar(MgATP,Ktr_mean,Ktr_err,'ko','linewidth',1.5,'markersize',6);
xlabel('[MgATP] (mM)','interpreter','latex','fontsize',6);
ylabel('$K_{tr}$ (sec.$^{-1}$)','interpreter','latex','fontsize',6);
set(gca,'fontsize',9,'xlim',[0 10],'ylim',[0 45]); box on;

return
%% Eval the N

% save from previous run
FAb = F_active;

tic 
t = []
EE = []
Nspan = 10:5:51;
for N = Nspan
    opts = struct('N', N, 'Slim', 0.06, 'PlotProbsOnFig', 1);

% insides
%     k = 3
%     for j = 1:length(vel)
%         F_active(j,k) = evaluateModel(fcn, vel(j), 0.11, MgATP(k),Pi,MgADP,g, opts);
%     end

N
t = [t toc];
EE = [EE sum(abs(FAb(:, k) - F_active(:, k)))];
end
figure();clf;subplot(121);
plot(Nspan, EE, 'x-',Nspan, t, '+-');xlabel('N');legend('Error', 'Time');
subplot(122);plot(EE, t, '*-');xlabel('Error');ylabel('time');

%% show in time
opts = struct('N', 50, 'Slim', 0.06, 'PlotProbsOnFig', 0);
tspan = 0:0.02:t_ss;
clear f;
for i = 1:length(vel)
    for j = 2:length(tspan)
        f(i,j) = evaluateModel(fcn, vel(i), tspan(j), 2,Pi,MgADP,g, opts);
    end
end   
%
figure(104);clf;hold on;
for i = 1:size(f, 1)
    plot(tspan, f(i, :), '*-')
    l{i} = "v = " + num2str(vel(i));
end
l{length(l)+1} = "zero force";
plot([0, tspan(end)], [0 0], 'r--');
legend(l)
    