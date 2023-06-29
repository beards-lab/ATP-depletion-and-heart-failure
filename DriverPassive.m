% driver passive

% set the params

% set the modifiers
% let it optim
% clear;

% Documentation use
% mods = {'r_a', 'r_d', 'mu', 'ks', 'k1', 'c', 'gamma', 'alpha1', 'e'};


% optimized for 1s
opt_mods = [  1.1797,     0.6314,    1.1477,    0.5181,    0.5833,    1.9550,    1.6055, 1];
% optimized for 100s
opt_mods =   [  0.0298    0.0290    0.3837    1.2774    1.0610    2.0034     1.7758, 1, 1];
% reoptimized for 100s with faster sampling around the peak
opt_mods =   [   0.0051    0.0379    0        1         2.3957    2.1260    1.8250    1.3943   0.4839];

% opt tuned for 1s
opt_mods =   [  0.1032    0.0339    1         0         1.0787    1.0346    1.2138    1.2506    1.1316];

% reopt for 100s, incl the exponent
opt_mods =   [0.0108    0.0338    0.1         10    2.2214    1.8489    1.7039    1.1048    0.5377];
% opt_mods =   [  0.0298    0.10290    0.006837    8.2774    1.0610    2.0034     1.7758, 1];
% opt_mods =   [  0.06298    0.060    0.3837    1.2774    4.0610    2.0034     1.7758, 1];

% optimized for 0.1s
% opt_mods = [0.0298    0.0290    200   0.1   1.0610    2.0034    1.7758, 1];

% optimized for 100ms and 100s
% opt_mods = [ 0.0022    0.0403    0.5058 0.0003    3.7887    1.5490 3.1429    1.0322    0.6723];

% optimized for 100ms and 100s
opt_mods = [1.1568    0.7497    2.0208    0.2414  0.5852    1.0600    1.1421    1.6024    1.0790];
% handtuned
opt_mods = [1.1568    0.7497    .20208    2.414/5  0.5852    1.0600    1.1421    1.6024    1.0790];

% optimized for peaks
% opt_mods = [1.1021    2.0296    1.3551    0.2645    0.2657   1.0506    0.9401    2.3257    0.7581];

% optimized for peaks with some overshoot
% opt_mods = [0.9081    3.2012   19.9710    1.4491    0.3239    0.3479    0.0401    2.3739    0.7189];

% optimized for 
plotEach = true;
figure(101);clf;

% ramp duration
rd = 1; 

% mods = {'r_a', 'r_d', 'mu', 'ks', 'k1', 'c', 'gamma', 'alpha1', 'e', 'phi', L0};
% opt_mods = [0.05081    15.2012   19.9710    0.14491    0.013239    0.3479    0.401    2.3739    0.7189];
opt_mods = [.001    .002        1    1      1  0.01    2.8        1        1   20     5.5];

% optim result
opt_mods = [0.0007    0.0013    0.8387    2.1727    2.0435    0.0133 4.1448    0.0010    0.0043   23.8029    5.9188];

% optim result handtuned
opt_mods = [0.0007    0.0013    0.018387    0.21727    50.435    0.0133 4.1448    0.0010    0.0043   4.8029    5.9188];

% optim result again, handtuned again
% mods = {'r_a',      'r_d',     'mu',     'ks',     'k1',     'c', 'gamma', 'alpha1',    'beta',   'phi',    L0};
opt_mods = [0.0007    0.0015    0.0184    0.2173   48.1429    0.010019 4.2373    0.0010    0.0046    3.5433    5.9446]

r = 1;
% second round of handtuning
opt_mods = [0.2*r    0.0514*r    0.0184    0.2173   120.1429*r    0.010019 4.2373    0.0010    2.6    3.5433    5.9446]

% enlarge the struct
on = ones(1, 11);
on(1:length(opt_mods)) = opt_mods;
opt_mods = on;

%
figure(10101)
tic
datastruct = load(['data/bakers_passiveStretch_' num2str(rd*1000) 'ms.mat']);
datatable = datastruct.datatable;
time_end = datatable(end, 1);
toc 
tic
evaluatePassive;
Es
toc
xlim([0, 2+min(rd*3, 200)])
% xlim([2 2.2])
% ylim([0, 0.02])
%%
figure(10102);semilogy(Tsim, Fatts);
%% compare peaks and steady state to data
peaks_sim = [];ss_sim = []; % sim peaks and sim steady state
peaks_data = [];ss_data = []; % data peaks and steady state

rds = [0.02, 0.1, 1, 10, 100];
rd_i = 1;
ft_sim = cell(length(rds), 1);ft_data = cell(length(rds), 1);
for rd = rds
    disp(['Processing ' num2str(rd*1000) 'ms...'])
    datastruct = load(['data/bakers_passiveStretch_' num2str(rd*1000) 'ms.mat']);
    datatable = datastruct.datatable; time_end = 200;
    figure(rd_i + 40);
    evaluatePassive;
    xlim([0, 2+min(rd*3, 200)])
    peaks_sim = [peaks_sim, max(Ftot)]; ss_sim = [ss_sim, Ftot(end)];
    peaks_data = [peaks_data, max(datatable(:, 3))]; ss_data = [ss_data, datatable(end, 3)];
    
    ft_sim{rd_i} = Ftot_int;
    ft_data{rd_i} = datatable;
    rd_i = rd_i + 1;
end

% save('passiveData.mat', 'peaks_data', 'ss_data');
%
figure(14);clf;
Np = length(peaks_sim); % number of peaks
x = 1:Np; % xaxis
set(gca,'ColorOrderIndex',1)
semilogx(rds, peaks_data, 'o-', rds, peaks_sim, '+-', 'LineWidth',2, 'MarkerSize',8);
hold on;set(gca,'ColorOrderIndex',1)
semilogx(rds, ss_data, 's-', rds, ss_sim, 'x-', 'LineWidth',2, 'MarkerSize',8);
legend('Peak (data)', 'Peak (model)', 'End (data)', 'End (model)');
xlabel('Ramp-up time (s)');ylabel('Passive tension (kPa)');
title({"Peak and steady state","during passive stretch"})
set(gca,'fontsize',16);
%%
figure(15);clf;
ldat = [];lsim = [];
for rd_i = 1:length(rds)
    set(gca,'ColorOrderIndex',rd_i)
    ldat = [ldat semilogx(ft_data{rd_i}(:, 1), ft_data{rd_i}(:, 3), '-', 'Linewidth', 2)];
    hold on;
%     set(gca,'ColorOrderIndex',rd_i)
%     semilogx(ft_data{rd_i}(:, 1), ft_sim{rd_i}(:, 1), ':', 'Linewidth', 2);
end

for rd_i = 1:length(rds)
%     set(gca,'ColorOrderIndex',rd_i)
%     l = [l semilogx(ft_data{rd_i}(:, 1), ft_data{rd_i}(:, 3), '-', 'Linewidth', 1)];
%     hold on;
    set(gca,'ColorOrderIndex',rd_i)
    lsim = [lsim semilogx(ft_data{rd_i}(:, 1), ft_sim{rd_i}(:, 1), ':', 'Linewidth', 2)];
end

title({" ","Tension response to ramp-up"});
xlabel('time (s)')
ylabel('Passive tension (kPa)')
ylim([0 15])
xlim([1.9 200])
legend([lsim(1) ldat], 'Ramp-up 20ms model', 'Ramp-up 20ms data', 'Ramp-up 100ms', 'Ramp-up 1s', 'Ramp-up 10s', 'Ramp-up 100s');
set(gca,'fontsize',16);
%% gradient search
figure(1010);clf;
options = optimset('Display','iter', 'TolFun', 1e-6, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 500);
% g = [1, 1, 1, 1, 1, 1, 1, 1];
% g = [1.2539    0.4422];
opt_mods =    [0.0022    0.0403    0.5058     0.0003    3.7887    1.5490     3.1429    1.0322    2]
opt_mods = ones(1, 9);
%%
datatables = [];rds = [];
% ident all params
x0 = opt_mods;
% optimfun = @(g)evalPassiveCost(g, datatables, rds);

% pick some vars only
mods = {'r_a', 'r_d', 'mu', 'ks', 'k1', 'c', 'gamma', 'alpha1', 'beta', 's0', 'L0'};
sel = [1 2 5 8 9 10 11]
x0 = opt_mods(sel);

% use only selected subset
optimfun = @(g)evalPassiveCost(insertAt(opt_mods, g, sel) , datatables, rds);
%%
% optimfun = @(g)evalPassiveCost(g, datatables, rds);

optimfun(x0)
%%
% ident just a subset
% ps = [1:2, 5, 6, 7, 8, 9];
% disp('Optimizing for ') 
% disp(mods(ps))
% x0 = opt_mods(ps);
% optimfun = @(x0)evalPassiveCost([x0(1, [1, 2]), 1, 0, x0(1, [3, 4, 5, 6, 7])], datatable, rd); % optim just a subset
%
% x0 = [0.0008    0.0017    1.0000    1.0000    1.0095    0.0119 3.6611    0.9297    1.0620   17.3403    5.7232];
x = fminsearch(optimfun, x0, options);
x

%%
opt_mods(sel) = x

%% EvalPassive test
evalPassiveCost(ones(1, 9), [], [])

%%

evalPassivePeaks(opt_mods, rds)

optimfun = @(g)evalPassivePeaks(g, rds);
x = fminsearch(optimfun, opt_mods, options);


%% GA
ga_Opts = optimoptions('ga', ...
    'PopulationSize',180, ...            % 250
    'Display','iter', ...
    'MaxStallGenerations',5, ...  % 10
    'UseParallel',true);

mods = {'r_a', 'r_d', 'mu', 'ks', 'k1', 'c', 'gamma', 'alpha1', 'e'};

Ng = length(sel);

% optimfun = @(g)evalPassiveCost(g, datatable, rd);

% default bounds struct - based on the mods
upp = 5; lwr = 0.05;
lb = lwr*ones(1, Ng);
ub = upp*ones(1, Ng);

[p_OptimGA,Res_OptimGA,~,~,FinPopGA,FinScoreGA] = ...
    ga(optimfun,Ng, ...
    [],[],[],[],lb,ub,[],ga_Opts);

save env;
%% PSO
pso_opts = optimoptions('particleswarm', ...    
    'Display','iter', ...
    'UseParallel',true, 'InitialSwarmMatrix', opt_mods);
x = particleswarm(optimfun, Ng, lb, ub, pso_opts)
%% optimfun(p_OptimGA)
options = optimset('Display','iter', 'TolFun', 1e-6, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 500);
x2 = fminsearch(optimfun, p_OptimGA, options)
save x2
optimfun(x2)


%%
function Ep = evalPassivePeaks(opt_mods, rds)
    if any(opt_mods<0)
        Ep = Inf;
        return;
    end
%%    
    load('passiveData.mat'); % 'peaks_data', 'ss_data'
    plotEach = false;
    calcInterpE = false;

    Es = zeros(1, length(rds)); % error set
    peaks_sim = zeros(1, length(rds)); % peak values
    ss_sim = 0; %% steady state val
    
    % init sim
    % figure(1);clf;hold on;
    tic
    simInit = true;simRamp = false;simRecover = false;
    rd = 0;
    evaluatePassive;
    aInit = a;
    % toc

    for i_rd = 1:length(rds)-1
        % ramp 10 - take peak
        a = aInit;
        rd = rds(i_rd);
        simInit = false;simRamp = true;simRecover = false;
        % figure;
        % plotEach = true;
        evaluatePassive;
        % plot(Tsim,Ftot)
        peaks_sim(i_rd) = max(Ftot);
        % toc
    end

    % longest ramp - take peak and steady state
    a = aInit;
    rd = rds(end);
    simInit = false;simRamp = true;simRecover = true;
    % figure;
    evaluatePassive;
    peaks_sim(end) = max(Ftot);
    ss_sim = Ftot(end); % steady state value
    Es(end) = max(Ftot); % TODO error error peak
    toc

    Ep = sum((peaks_sim - peaks_data).^2) + length(rds)*(ss_sim - mean(ss_data))^2
    
%%
    figure(15);cla;hold on;
    plot(1:length(rds), peaks_data, 'o-', 'LineWidth',2)
    plot(1:length(rds), peaks_sim, 'x--', 'LineWidth',2)
    plot([1 length(rds)], mean(ss_data)*ones(2, 1), '-', 'LineWidth',2);
    plot([1 length(rds)], ss_sim*ones(2, 1), '--', 'LineWidth',2);
end

function Es = evalPassiveCost(opt_mods, datatable, rd)
if any(opt_mods<0)
    Es = Inf;
    return;
end
% normalize to ones as input
opt_mods = opt_mods.*[ 0.0007    0.0014    0.0184    0.2173   45.3794    0.0100 4.2373    0.0009    0.0048    3.4744    5.8306];
plotEach = true;
figure(1);

% datastruct = load('data/bakers_passiveStretch_1000ms.mat');
% datatable = datastruct.datatable;
% evaluatePassive;

% ramp duration
subplot(211);
rd = 0.1; 
datastruct = load(['data/bakers_passiveStretch_' num2str(rd*1000) 'ms.mat']);
datatable = datastruct.datatable;
time_end = datatable(end, 1);

evaluatePassive;
Es100 = Es;

% ramp duration
subplot(212);
rd = 100; 
datastruct = load(['data/bakers_passiveStretch_' num2str(rd*1000) 'ms.mat']);
datatable = datastruct.datatable;
time_end = datatable(end, 1);

evaluatePassive;
Es100000 = Es;

Es = Es100000 + Es100;
end


