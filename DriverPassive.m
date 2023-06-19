% driver passive

% set the params

% set the modifiers
% let it optim
clear;

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

% optimized for 
plotEach = true;
figure(101);clf;

% ramp duration
rd = 1; 

tic
datastruct = load(['data/bakers_passiveStretch_' num2str(rd*1000) 'ms.mat']);
datatable = datastruct.datatable;
time_end = datatable(end, 1);
toc 
tic

evaluatePassive;
Es
toc
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
    evaluatePassive;
    peaks_sim = [peaks_sim, max(Ftot)]; ss_sim = [ss_sim, Ftot(end)];
    peaks_data = [peaks_data, max(datatable(:, 3))]; ss_data = [ss_data, datatable(end, 3)];
    
    ft_sim{rd_i} = Ftot_int;
    ft_data{rd_i} = datatable;
    rd_i = rd_i + 1;
end
%%
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
datatables = [];rds = [];
% ident all params
x0 = opt_mods;
optimfun = @(g)evalPassiveCost(g, datatables, rds);


% ident just a subset
% ps = [1:2, 5, 6, 7, 8, 9];
% disp('Optimizing for ') 
% disp(mods(ps))
% x0 = opt_mods(ps);
% optimfun = @(x0)evalPassiveCost([x0(1, [1, 2]), 1, 0, x0(1, [3, 4, 5, 6, 7])], datatable, rd); % optim just a subset
%
x = fminsearch(optimfun, x0, options);
x

%% EvalPassive test
evalPassiveCost(ones(1, 9), [], [])

%% GA
ga_Opts = optimoptions('ga', ...
    'PopulationSize',100, ...            % 250
    'Display','iter', ...
    'MaxStallGenerations',4, ...  % 10
    'UseParallel',true);

mods = {'r_a', 'r_d', 'mu', 'ks', 'k1', 'c', 'gamma', 'alpha1', 'e'};

Ng = length(mods);

optimfun = @(g)evalPassiveCost(g, datatable, rd);

% default bounds struct - based on the mods
upp = 4; lwr = 0.05;
lb = upp*ones(1, Ng);
ub = upp*ones(1, Ng);

[p_OptimGA,Res_OptimGA,~,~,FinPopGA,FinScoreGA] = ...
    ga(optimfun,Ng, ...
    [],[],[],[],lb,ub,[],ga_Opts);

save env;
%% PSO
pso_opts = optimoptions('particleswarm', ...    
    'Display','iter', ...
    'UseParallel',true);
x = particleswarm(optimfun, Ng, lb, ub, pso_opts)
%% optimfun(p_OptimGA)
options = optimset('Display','iter', 'TolFun', 1e-6, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 500);
x2 = fminsearch(optimfun, p_OptimGA, options)
save x2
optimfun(x2)


%%
function Es = evalPassiveCost(opt_mods, datatable, rd)
if any(opt_mods<0)
    Es = Inf;
    return;
end
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


