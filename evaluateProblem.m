function [Etot, E] = evaluateProblem(fcn, g, drawPlots, evalParts)


if nargin < 4
    % force X velocity, force X concetration, ktr, Km, plot zeroforce
    % velocity to ATP conc
    evalParts = [1, 1, 1, 1, 0];
end
E = zeros(4, 1);

LoadData;

%% Limit the ATP to single val
Data_ATP = Data_ATP(:, 1:2);
MgATP = MgATP(1);

%% Set up environment
MgADP = 0; 
Pi    = 0;
F_active_0 = zeros(length(vel), max(length(MgATP), 3));
t_sl0 = [0 0.1]; % time at which the SL = 2.2 at velocity of 1 ML/s - time to stop the experiment
t_ss = [0 0.3]; %%  steady state time
params0 = getParams([], g);
params0.UseSerialStiffness = false;
% params0.UsePassive
PU0 = [];

%% force x velocity
F_active = F_active_0;
if evalParts(1)
%%
params = params0;
params.OutputAtSL = Inf; % this is solved by limited sim time
params.Slim = 0.2;
params.N = 120;
% update the dS
params = getParams(params, g);
for k = length(Data_ATP(1, 2:end))
    params.MgATP = Data_ATP(1, k);
    for j = 1:length(vel)
        params.Velocity = vel(j);
        if vel(j) == 0 
%             if isfield(params, 'PU0'), params = rmfield(params, 'PU0');end
            [F_active(j,k) out] = evaluateModel(fcn, t_ss, params);
            params.PU0 = out.PU(end, :);
        else
            [F_active(j,k) out] = evaluateModel(fcn, t_sl0/abs(vel(j)), params);
        end        
    end
end
% save the init PU to save some time
PU0 = params.PU0;

% sum the error per ATP level
E(1) = 0;
for i = 1:length(Data_ATP(1, 1:end-1))
    E(1) = E(1) + sum(abs(F_active(:,i)-Data_ATP(:,i+1)).^2);
end

% normalize by number of data points
E(1) = E(1)/size(Data_ATP, 1)/(size(Data_ATP, 2) -1);
   
end


%% linearized method for vmax
if evalParts(5)
   k = 1; % atp of 8
   if F_active == F_active_0
       % the first part has not been run, we have to evaluate additionally
       params = params0;
       params.MgATP = MgATP(k);
       params.Velocity = vel(end-1);
       F_active(end-1, k) = evaluateModel(fcn, t_ss, params);
       params.Velocity = vel(end);
       F_active(end, k) = evaluateModel(fcn, t_ss, params);
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

end
%% force x iso MgATP

F_iso = zeros(1, length(MgATP_iso));
if evalParts(2)
    %%
    params = params0;
%     params.PU0 = PU0;
%     params.N = 5;
%     params.Slim = 0.03;
%     params.OutputAtSL = Inf;
    % udpate for dS
%     params = getParams(params, g);
    for k = 1:length(MgATP_iso)
      % Zero velocity:
      params.MgATP = MgATP_iso(k);
      [F_iso(k) out] = evaluateModel(fcn, t_ss, params);
      params.PU0 = out.PU(end, :); % speeds the next initialization a bit
%       F_iso(k) - out.Force(end)
    end
    
    % weight per data point
    E(2) = (sum(abs(F_iso-F_data_r*F_iso(end)).^2))/length(MgATP_iso);
end
%% time constant of MgATP
% Tspan = [0:0.01:0.2];
% Frel = zeros(length(MgATP), length(Tspan));
Ktr = zeros(1, length(MgATP));
KtrOuts = cell(length(MgATP), 2);
if evalParts(3)
%%     
params = params0;
% if ~isempty(PU0)
% different slim, cant use now
%     params.PU0 = PU0;
% end

% cut all attached
% TODO elaborate on that, inlcude length etc...
% read about the protocol
% params.PU0(1:params.ss*3 + 4) = 0;
% params.N = 20;
params.ValuesInTime = true;
% params.OutputAtSL = Inf;
% params.Velocity = [-20, 0, 0];
% times = [-0.2, 0, 0.5];
% update dS
% params = getParams(params);

for k = 1:length(MgATP)
  params.MgATP = MgATP(k);
  
  if ~params.UseKtrProtocol
      % just an approximation
      % Tspan array returns F active array
      [~, out] = evaluateModel(fcn, t_ss, params);
  else
      v = 150; % ML/s
%       times  = [-1e3, 0,  2, 20, 22.5, 25 , 25.5, 1e3]/1000 - 25.5e-3;
      pos_ML = [1   , 1,0.8,0.8, 1.05,1.05,     1, 1 ];
    
      % putting the numbers as a difference
      times = cumsum([0   , 1,0.2/v,0.018,0.25/v, 0.005, 0.05/v, 1]);
      params.Velocity = diff(pos_ML)./diff(times);
      [~, out] = evaluateModel(fcn, times - times(end-1), params);
%       figure;hold on;
%       plot(out.t, out.SL)
%       plot(out.t, out.FXB/out.FXB(end))
  end
  
  if out.FXB(end) == 0
      Ktr(k) = 0;
      continue;
  end
  
%   Frel(k, :) = out.Force./out.Force(end);
  Frel = out.FXB./out.FXB(end);
  
  % get the time constant
%   Ktr(k) = 1/interp1(Frel,out.t,1-exp(-1)); % time constant for Frel(1/Ktr) = 1-exp(-1)
%   find val at time 1-exp(-1)
  i = find(Frel >= 1-exp(-1), 1);
  Ktr(k) = 1/out.t(i);
  KtrOuts{k, 1} = out.t;
  KtrOuts{k, 2} = Frel;
  
end

E(3) = (sum(abs(Ktr-Ktr_mean).^2))/length(MgATP);
end
%% Zero-force velocity
K_m = 0;
v_max = 0;

if evalParts(4)
%%
% have to search for velocity, which gives us a zero force
% currently the target is about 6 um/s
MgATP_vmax = 2;
targetVel = 6;

tol = 1e-3;

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
    params = struct('N', 30, 'Slim', 0.06, 'PlotProbsOnFig', 0, 'ValuesInTime', 0);
%     params.PU0 = PU0;
    params.MgATP = MgATP_vmax;
    params = getParams(params, g);
    
    while abs(e) > tol && i < maxIter
        i = i + 1;
        params.Velocity = -v_max;
        
        [e out] = evaluateModel(fcn, t_ss, params);
        params.PU0 = out.PU(end, :);
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
    tol = 0.00001;
    % debug only
    range = [K_m*(1/2)^(maxIter), K_m + step*(1 -(1/2)^(maxIter))];
    
    params.Velocity = -v_max/2;
    while abs(e) > tol && i < maxIter
        i = i + 1;
        params.MgATP = K_m;
        e = evaluateModel(fcn, t_ss, params, g, opts);
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

if evalParts(6)
%%
%     load and compare with bakers step-up experiment
% generated by LoadBakesExp.m
load data/bakers_rampup8.mat
params = params0;
params.Slim = 0.06;
params.N = 40;

params.OutputAtSL = Inf;
params.ValuesInTime = true;

params.Velocity = velocitytable(1:end-1, 2);
params.MgATP = 8;
params.SL0 = 2.0;
params.ML = 2.0;

% if ~isempty(PU0)
%     params.PU0 = PU0;
% end
% update dS and PU0
params = getParams(params, g);

try
[F out] = evaluateModel(fcn, velocitytable(:, 1), params);
t_exp = datatable(:, 1);

% Li = interp1(out.t, out.SL, t_exp);
Fi = interp1(out.t, out.Force, t_exp);
e = (datatable(:, 3) - Fi).^2;
se = mean(e);
catch e
    se = NaN;
end
E(6) = se;
%%
% figure(1);clf;
% subplot(211);plot(out.t, out.SL, t_exp, datatable(:, 2), t_exp, Li, 'x-');xlim([0 t_exp(end)]);
% subplot(212);plot(out.t, out.Force, t_exp, datatable(:, 3), t_exp, Fi, 'x', 'Linewidth', 2, 'MarkerSize', 10);
% yyaxis right;plot(t_exp, e,[0 t_exp(end)],[se se],'Linewidth', 1);xlim([0 t_exp(end)]);



end
%% Return

penalty = sum(max(0, -g))*1000;
if any(g < 1e-3) || any(g > 10)
    penalty = Inf;
end
E(end+1) = penalty;
Etot = sum(E);

% %% save best so far?
% load e_best;
% if Etot < e_best
%     save e_best Etot;
%     writematrix([g;E], [datestr(now(), 'DD_MM_YYYY HH_mm_ss') num2str(Etot) '.csv']);
% end

%% plot?

if ~drawPlots
    return;
end

disp("With Km of " + num2str(K_m) + " at max velosity of " + v_max)
%% Force velocity
figure(101); clf; axes('position',[0.1 0.6 0.35 0.35]); hold on;
% figure(); clf; axes('position',[0.1 0.6 0.35 0.35]); hold on;
% figure(102);hold on;
plot(F_active(:,1), -vel,'b-','linewidth',1.5);
if size(Data_ATP, 2) > 2
    plot(F_active(:,2), -vel,'g-','linewidth',1.5);
    plot(F_active(:,3), -vel,'r-','linewidth',1.5);
end
ylabel('Velocity (ML/s)','interpreter','latex','fontsize',16);
xlabel('Force (kPa)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14); 
% axis([0 65 0 6]);
axis([-10 65 0 6]);
title('Force-velocity')
box on;grid on;
plot(Data_ATP(:,2),Data_ATP(:,1),'bo','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
if size(Data_ATP, 2) > 2
    plot(Data_ATP(:,3),Data_ATP(:,1),'go','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
    plot(Data_ATP(:,4),Data_ATP(:,1),'ro','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
    ldg = legend('8','4','2 mM');
    title(ldg,'[MgATP]');
end


% MgATP
% figure(104); clf; 
axes('position',[0.6 0.6 0.35 0.35]); 
semilogx(MgATP_iso,F_iso,'b-','linewidth',1.5); hold on;
semilogx(MgATP_iso,F_data_r*F_iso(end),'ko','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1]);
xlabel('[MgATP] (mM)','interpreter','latex','fontsize',16);
ylabel('Max. Force (kPa)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14); 
box on;
title('Max-force to ATP');

% KTR
% figure(106); clf; 
axes('position',[0.1 0.13 0.35 0.35]); hold on;
if size(KtrOuts, 1) > 1 
    plot(KtrOuts{3, 1},KtrOuts{3, 2},'r-','linewidth',1.5);
    plot(KtrOuts{2, 1},KtrOuts{2, 2},'g-','linewidth',1.5);
end
plot(KtrOuts{1, 1},KtrOuts{1, 2},'b-','linewidth',1.5);
% plot(Tspan,F_active,'linewidth',1.5);
xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
ylabel('Force (rel.)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'ylim',[0 1.1]);  box on;
title('Speed of the transient');
if size(KtrOuts, 1) > 1 
    ldg = legend('2','4','8 mM','location','northwest');
    title(ldg,'[MgATP]');
end

axes('position',[0.25 0.15 0.15 0.2]); hold on;
plot(MgATP,Ktr,'k.-','linewidth',1.5);
errorbar([8 4 2],Ktr_mean,Ktr_err,'ko','linewidth',1.5,'markersize',6);
xlabel('[MgATP] (mM)','interpreter','latex','fontsize',6);
ylabel('$K_{tr}$ (sec.$^{-1}$)','interpreter','latex','fontsize',6);
% set(gca,'fontsize',9,'xlim',[0 10],'ylim',[0 45]); box on;
set(gca,'fontsize',9,'xlim',[0 10]); box on;

% axes('position',[0.1 0.13 0.45 0.33]); hold on;
% figure;
%%
if evalParts(5)
    %%
    v_atp = [];
    params = struct('N', 50, 'Slim', 0.06, 'PlotProbsOnFig', 0);
    atp_range = 0:0.2:2;
    % more datapoints for lower bounds, relax near SS
    atp = (exp(atp_range)-1)/exp(atp_range(end))*atp_range(end)*1.5;
%     plot(atp, ones(length(atp)), 'bx--');
    v_atp = zeros(size(atp));
    maxIter = 8;
    % real tolerance achievable:
    % realTol = v_max/2*2*(1/2)^maxIter
    params = params0;
    

    for a = 1:length(atp)
        disp("Running K_m " + num2str(round(a/length(atp)*100)) + "%...");
        % the loop will search up to step*2
        v = v_max/2; step = v;
        
        params.MgATP = atp(a);
        for i = 1:maxIter
            params.Velocity = -v;
            e = evaluateModel(fcn, t_ss, params, g,opts);
            step = abs(step)/2;
            if e > 0 v = v + step; else v = v - step;end;
        end
        v_atp(a) = v;
    end
else
    atp = [0, K_m, 2];
    v_atp = [0, v_max/2, v_max];
end
%%
% figure(102);clf;hold on;
if evalParts(5)
%     axes('position',[0.6 0.13 0.10 0.35]);cla;hold on;
    axes('position',[0.6 0.13 0.35 0.35]);cla;hold on;
    set(gca,'fontsize',14);box on;
    title("Maximum sliding velocity, Km " + num2str(K_m));
    fill([0.04 0.1 0.1 0.04], [0 0 4 4], [0.8 0.98 0.8], 'LineStyle', 'None');
    fill([1.96 2.0 2.0 1.96], [2 2 8 8], [0.8 0.8 0.98], 'LineStyle', 'None');
    plot([0, K_m, 2], [0, v_max/2, v_max], 's--', 'Color', [0 0 0.3], 'MarkerSize', 8, 'MarkerFaceColor', [0 0 0.3]);
    plot(atp, v_atp, 'bo-', 'MarkerSize', 4, 'MarkerFaceColor', [0 0 1]);

    plot([0, 2], [v_max/2, v_max/2], '--', 'Color', [1 0 0]);
    text(1.5, v_max/2 - 0.3, 'half-max velocity', 'Color', [1 0 0])
    xlim([0, 2]);
    xlabel('MgATP');
    ylabel('Zero force velocity ML/s');
    legend('K_m area', 'V_{max} area', 'v_{max} simplified', 'v_{max}', 'Location', 'NorthWest', 'fontsize',9)
end
if evalParts(6)
    % Bakers step-up
%     axes('position',[0.75 0.13 0.2 0.35]);cla;hold on;
    axes('position',[0.6 0.13 0.35 0.35]);cla;hold on;
%     axes('position',[0.75 0.13 0.2 0.35]);cla;hold on;
    cla;plot(out.t, out.Force, t_exp, datatable(:, 3), t_exp, Fi, '|', 'Linewidth', 2, 'MarkerSize', 8);
    axis([0, 0.5, 50, 90]);
%     yyaxis right;plot(t_exp, e,[0 t_exp(end)],[se se],'Linewidth', 1);xlim([0 t_exp(end)]);
end
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

%% show in time

opts = struct('N', 50, 'Slim', 0.06, 'PlotProbsOnFig', 0);
opts.ValuesInTime = 1;

for i = 2:length(vel)
    opts.PlotProbsOnFig = 1;
    evaluateModel(fcn, vel(i), abs(2*vel(i)),2, Pi,MgADP,g, opts);
end
    