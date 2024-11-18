% test sensitivity 
params0 = getParams();
clf
% ModelParamsInitNiceSlack;
% ModelParamsInit2;


cost_sap = []; % SA plus
cost_sam = []; % SA minus

% ModelParamsInitOptim_slack4
% ModelParamsOptim_tmp
% ModelParamsInitOptim_slackAll
% ModelParamsOptim_tf2_slackLast
% ModelParamsOptim_tf2_slackOnsetAll_LeftOnly
% ModelParamsOptim_DtKtr;
% ModelParamsOptim_DtKtr_OV
% ModelParamsOptim_DtKtr_OV2
% ModelParamsOptim_fudgeSlackVelocities
% ModelParamsInitOptim_slack4
ModelParams_SR

params0.UseOverlapFactor = false;
% params0.ksr = params0.ksr0;

params0.RunSlackSegments = 'First';
% params0.RunSlackSegments = 'FirstAndLastExtended';

% need to run the sim first to get the velocity table in the workspace. 
% velocitytable = load('data/bakers_slack8mM_all.mat').velocitytable;
% params0.ResetSRat = [velocitytable(4, 1), 0.0;...
%     velocitytable(8, 1), 0.2;...
%     velocitytable(12, 1), 0.3;...
%     velocitytable(16, 1), 0.4;...
%     velocitytable(20, 1), 0.7 ];
params0.ResetSRat = [];

% params0.RunSlackSegments = 'stairs-up';

% params0.ka = 1e-3;
% params0.kd = 1e-2;
% params0.kadh = 0;

% params0.xrate = 1.6;
% params0.ksr0 = 5;
% params0.sigma1 = 1e6;
% % to SR
% params0.kmsr = 1.5;
% params0.sigma2 = 1e6;

params0.UseDirectSRXTransition = false;

% params0.Lsc0 = 1.51;
% params0.e2R = 1;
% ModelParamsOptim_tf2_slackFirst
% ModelParamsOptim_tf2_slackFirstLast
% ModelParamsOptim_tmp.m

% params0.L_thick = 1.67; % Length of thick filament, um
% params0.L_hbare = 0.10; % Length of bare region of thick filament, um
% params0.L_thin  = 1.20; % Length of thin filament, um

params0.drawPlots = true;
params0.drawForceOnset = true;
params0.PlotEachSeparately = true;
params0.drawForceOnset = true;
params0.ShowResidualPlots = false;
params0.justPlotStateTransitionsFlag = false;

params0.UsePassiveForSR = false;
% RSR2PT = params.ksr0*exp(F_total/)*P_SR;
% RPT2SR = params.kmsr*exp(-F_total/params.sigma2)*PT;

% params0.ksr0 = 1.78;
% params0.sigma1 = 18.7;
% params0.kmsr = 51.9;
% params0.sigma2 = 1e6;

% params0.dr = 0.01;
% params0.kstiff2 = 1.8e4;
% params0.kstiff1 = 10.2e3;

params0.MaxSlackNegativeForce = 0;
params0.FudgeVmax = true;
params0.FudgeVmax = 1;
params0.vmax1 = 20;
params0.FudgeA = 0;
params0.FudgeB = 121.922;
params0.FudgeC = -209.18;

% params0.vmax = 18;

% params0.kSE = 1.3e4;
% params0.kSEn = 1.3e1;
% params0.mu = 1e-1;


params0.ghostLoad = 'DtKtr_OV';
tic

RunBakersExp;
toc
sum(E)
%% experiment with SRD
figure(1);clf;
params0.UseSuperRelaxedADP = false;
params0.UseSuperRelaxed = true;
params0.UsePassiveForSR = false;

params0.kmsr = 20;
params0.ksr0 = 1;
params0.sigma1 = 1e6;
params0.sigma2 = 1e6;

params0.ksrd0 = params0.ksr0;
params0.kmsrd = params0.kmsr;
params0.ksr2srd = .1;
% params0.RunForceVelocity = 0;
params0.RunStairs = false;
params0.RunSlack = true;
params0.RunForceVelocity = false;

params0.RunSlackSegments = 'ramp-down';
% params0.RunSlackSegments = 'ramp-up';
params0.RunSlackSegments = 'stairs-down';

RunBakersExp;
%%
ModelParamsInit_TF2_slack4;
ModelParamsOptim_tf2_slackLast;
ModelParamsOptim_tmp;
%%
figure(1001); clf; hold on;
params0.RunSlackSegments = 'All';
params0.PlotEachSeparately = true;
params0.ShowStatePlots = false;
params0.ShowResidualPlots = false;
c0 = isolateRunBakersExp(params0);
%%
% all parameters
params0.mods = {'MaxSlackNegativeForce', ... 
'L_thick', ... 
'L_hbare', ... 
'L_thin', ... 
'e2L', ... 
'e2R', ... 
'K_coop', ... 
'k_on', ... 
'k_off', ... 
'vmax', ... 
'kah', ... 
'kadh', ... 
'ka', ... 
'kd', ... 
'k1', ... 
'k_1', ... 
'k2', ... 
'k_2', ... 
'k3', ... 
'ksr0', ... 
'sigma1', ... 
'sigma2', ... 
'kmsr', ... 
'kmsrd', ... 
'ksrd', ... 
'ksr2srd', ... 
'sigma_srd1', ... 
'sigma_srd2', ... 
'K_Pi', ... 
'K_T1', ... 
'dr', ... 
'kstiff1', ... 
'kstiff2', ... 
'kstiff3', ... 
'K_T3', ... 
'K_D', ... 
'alpha0', ... 
'alpha1', ... 
'alpha_1', ... 
'alpha2', ... 
'alpha2_L', ... 
'alpha2_R', ... 
'alpha3', ... 
's3', ... 
'alphaRip', ... 
'k2rip', ... 
'dr0', ... 
'dr1', ... 
'dr_1', ... 
'dr2', ... 
'dr3', ... 
'Amax', ... 
'mu', ... 
'kSE', ... 
'k_pas', ... 
'k2_L', ... 
'sigma0', ... 
'TK', ... 
'TK0', ... 
'gamma', ... 
'k2_R', ... 
'dr2_R', ... 
'Lsc0', ... 
'xrate'};
%%
params0.mods = {'k_on', ...
    'k_off', 'vmax', 'kah', 'kadh', 'ka', 'kd', 'k1', 'k_1', 'k2', 'k_2', 'k3', 'ksr0', 'sigma0', 'kmsr', 'K_Pi', 'K_T1', 'dr', 'kstiff1', 'kstiff2', 'kstiff3', 'K_T3', 'K_D', 'alpha0', 'alpha1', 'alpha_1', 'alpha2', 'alpha3', 's3', 'alphaRip', 'k2rip', 'dr0', 'dr_1', 'dr2', 'dr3', 'Amax', 'mu', 'kSE', 'k_pas', 'k2_L', 'ss', 'TK', 'TK0', 'dr1', 'sigma1', 'sigma2', 'gamma', 'alpha2_L', 'k2_R', 'dr2_R', 'alpha2_R', 'e2R', 'Lsc0', 'e2L', 'xrate'}
params0.mods = {'kstiff1', 'kstiff2', 'kmsr', 'ksr0', 'sigma1', 'sigma2', 'ka', 'kd', 'alpha2', 'k1', 'sigma1', 'kah', 'k2rip', 'alphaRip', 'alpha1', 'k2','k2_R','dr2_R', 'alpha2_R'};
params0.mods = {'kstiff1', 'kstiff2'};

% full
% params0.mods = {'kadh', 'ka', 'kd', 'k1', 'k2', 'ksr0', 'kmsr', 'dr', 'kstiff1', 'kstiff2', 'alpha0', 'alpha1', 'alpha2', 'alpha3', 'alphaRip', 'dr2', 'dr3', 'kSE', 'k_pas', 'TK', 'TK0', 'dr1', 'sigma1', 'gamma', 'alpha2_L', 'k2_R', 'k2_L', 'dr2_R', 'alpha2_R'};

%reduced
% params0.mods = {'kadh', 'ka', 'kd', 'k1', 'k2', 'ksr0', 'kmsr', 'dr', 'kstiff1', 'kstiff2', 'alpha0', 'alpha1', 'kSE', 'k_pas', 'sigma1', 'gamma', 'alpha2_L', 'k2_R', 'k2_L', 'dr2_R', 'alpha2_R'};

%reduced 2nd round
params0.mods = {'dr1', 'dr2','kd', 'k1', 'ksr0', 'kmsr', 'kstiff1', 'alpha1', 'alpha0', 'sigma1', 'k2_R', 'dr2_R', 'alpha2_R', 'alpha2_L'};
params0.g = ones(size(params0.mods));

% reduced
params0.mods = {'dr1', 'dr2','kd', 'k1', 'ksr0', 'kmsr', 'kstiff1', 'alpha1', 'alpha0', 'sigma1', 'alpha2_R', 'alpha2_L'};

% prep for the trans changes
params0.mods = {'dr1', 'alpha1', 'k1', 'alpha2_L', 'k2', 'dr2', 'alpha2_R', 'e2R', 'e2L'};

% extended
params0.mods = {'dr1', 'alpha1', 'k1', 'alpha2_L', 'k2', 'dr2', 'alpha2_R', 'e2R', 'e2L', 'kd', 'ksr0', 'kmsr', 'kstiff1', 'kstiff2', 'k_pas', 'gamma', 'Lsc0'};

% only left 
params0.mods = {'dr1', 'alpha1', 'k1', 'alpha2_L', 'k2', 'dr2', 'e2L', 'kd', 'ksr0', 'kmsr', 'sigma1', 'kstiff1', 'kstiff2', 'k_pas', 'gamma', 'Lsc0', 'kSE'};

params0.mods = {'dr1', 'alpha1', 'k1', 'alpha2_L', 'k2', 'dr2', 'e2L', 'kd', 'ksr0', 'kmsr', 'sigma1', 'kstiff1', 'kstiff2', 'kSE'};

% reduced
% params0.mods = {'dr1', 'alpha1', 'alpha2_L', 'k2', 'dr2', 'e2L', 'kd', 'ksr0', 'kmsr', 'kstiff1', 'kstiff2', 'k_pas', 'gamma', 'Lsc0', 'kSE'};

% params0.mods = {'k_pas', 'gamma', 'Lsc0', 'kSE', 'ksr0', 'kmsr', 'xrate'};

% params0.mods = {'xrate'};
% params0.mods = {'k_pas', 'gamma', 'Lsc0'};
params0.mods = {'dr1', 'alpha1', 'k1', 'alpha2_L', 'k2', 'dr2', 'e2L', 'kd', 'ksr0', 'kmsr', 'sigma1', 'kstiff1', 'kstiff2', 'kSE'};

% params0.mods = {'L_thick', 'L_hbare', 'L_thin'};
% 
% params0.mods = {'vmax', 'kSE'}
% 
% params0.mods = {'FudgeB','FudgeC'};

% new attempt
params0.mods = {'ksr', 'kmsr', 'sigma2', 'kstiff2', 'ka', 'kd', 'k1', 'k2', 'kah', 'alpha2_R', 'alpha2_L', 'e2R', 'e2L'};
params0.g = [0.7231    1.5436    0.7716    0.9494    0.8632    0.5957    1.2720    0.8273    0.2280    2.8035    0.2681    3.6707    0.5337];

% second round
params0.mods = {'ksr', 'kmsr', 'sigma2', 'kstiff1', 'kstiff2', 'ka', 'kd', 'k1', 'k2', 'kah', 'alpha2_R', 'e2R'};
params0.g = [0.7231    1.5436    0.7716    1    0.9494    0.8632    0.5957    1.2720    0.8273    0.2280    2.8035   3.6707];
params0.g = [0.7206    1.5516    0.7663    0.9437    0.9286    0.8858    0.5387    1.2715    0.8249    0.2298    2.4325    5.4854];

% last slack and params0.RunForceLengthEstim - cost 350
params0.mods = {'alpha0', 'alpha2_R','kmsr','gamma','kd','k1','ka','dr2_R','sigma1', 'k2_R','k2','alpha1', 'alpha2_L','kah','kstiff2','k_pas', 'kstiff1','xrate'};
params0.g = [0.9765    1.1966    1.6759    0.2419    1.1190    0.8689    1.3629    1.4356    0.9853    0.3401    0.9596    0.4014    1.0851    0.5159    1.0188    0.5793    2.0921    1.2786];
% only last slack
params0.g = [0.9578    1.4289    1.6017    0.2347    1.1474    0.8426    1.3470    1.3556    0.9765    0.3708    0.9696    0.4100    1.0063    0.5208    1.0063    0.5792    1.9811    1.3405];
%%
params0.g = ones(size(params0.mods));
saSet = 1:length(params0.mods);

params0.ghostLoad = '';
SAFact = 1.01;
c0 = isolateRunBakersExp(params0);
params0.PlotEachSeparately = false;
params0.ShowStatePlots = false;
params0.ShowResidualPlots = false;
p0 = params0;

for i_m = saSet
    % reset params
    params0 = p0;

    if params0.g(i_m) == 0
        % replace 0 with one - 
       params0.g(i_m) = 1;
    else
        params0.g(i_m) = params0.g(i_m)*SAFact;
    end

    fprintf('Mod %s is up to %g %%..', params0.mods{i_m}, params0.g(i_m)*100);
    % figure(i_m)
    cost = isolateRunBakersExp(params0);
    cost_sap(i_m) = cost;
    % 
    params0.g(i_m) = params0.g(i_m)/SAFact*(1+ 1-SAFact);
    fprintf('costing %1.4e€ and down to %g...', cost, params0.g(i_m)*100);
    cost = isolateRunBakersExp(params0);
    cost_sam(i_m) = cost;

    fprintf('costing %1.4e€. \n', cost);
end
params0 = p0;
%% visualize one way
figure(404)
err_1 = min(cost_sap)
figure;bar([cost_sap]');hold on;plot([0 length(saSet)+1], [c0 c0])    
xticks(1:length(params0.mods))
xticklabels(params0.mods)

%% visualize two way
err_1 = min(cost_sap, cost_sam)
figure;bar([cost_sap; cost_sam]');hold on;plot([0 length(saSet)+1], [c0 c0])    
xticks(1:length(params0.mods))
xticklabels(params0.mods)

%% only effective ones
params0.mods(cost_sap ~= c0)
%% sort by reduced and most increased
% c0 - min(cost_sam, cost_sap)
cost_cmb = cost_sap;
% cost_cmb = min(cost_sam, cost_sap);

cost_diff = c0 - cost_cmb
% better value
[val_bett cost_diff_better_i] = sort(cost_diff.*(cost_diff > 0), 'descend');
% most influential
[val_infl cost_diff_infl_i] = sort(cost_diff.*(cost_diff < 0), 'ascend');

param_ord = [cost_diff_better_i(val_bett > 0) cost_diff_infl_i(val_infl < 0)]
clf;
bar(cost_cmb(param_ord));xticks(1:length(params0.mods))
% hold on;plot([0 length(saSet)+1], [c0 c0]) 
xticklabels(params0.mods(param_ord))

%% select where it can be reduced OR most influential
mods = params0.mods(param_ord(1:20))
% mods = params0.mods(1:20)
g = ones(size(params0.mods));
% params0.mods = mods;
% saSet = 1:length(params0.mods);
%% sort by most influential
clf;
[val_infl cost_diff_infl_i] = sort(abs(cost_diff), 'descend')
sorted_infl = cost_diff_infl_i(val_infl>0);

bar(cost_cmb(sorted_infl));xticks(1:length(params0.mods))
xticklabels(params0.mods(sorted_infl))
hold on;plot([0 length(sorted_infl)+1], [c0 c0])    

% params0.mods = params0.mods(sorted_infl)
%% Optim fminsearch

options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.1, 'PlotFcns', @optimplotfval, 'MaxIter', 1500);
% g = [1, 1, 1, 1, 1, 1, 1, 1];
% g = [1.2539    0.4422];

params0.ShowResidualPlots = false;

g = params0.g;
% g = ones(length(params0.mods), 1);
optimfun = @(g)evaluateBakersExp(g, params0);
x = fminsearch(optimfun, params0.g, options)
params0.g = x;

params0.PlotEachSeparately = true;
RunBakersExp;
% params0.g = params0.g(1:26);
% save tmpFminSlackOpt;
% writeParamsToMFile('ModelParams_surroFmin.m', params0);
%% optim surrogateopt

options = optimoptions('surrogateopt','Display','iter', 'MaxTime', 60*60*6, 'UseParallel',true);

optparams = params0;
% optparams.mods = optparams.mods(sorted_infl(1:20))
optparams.BreakOnODEUnstable = true;
optparams.RunForceVelocity = false;
optparams.RunForceLengthEstim = false;
optparams.RunSlackSegments = 'All';

g = ones(length(optparams.mods), 1);
% retrieve from table
% T = load("optTable.mat").T;

lb = T.lb./T.val; ub = T.ub./T.val;
optparams.mods = cellstr(T.mods');

optimfun = @(g)evaluateBakersExp(g, optparams);

[x,fval,exitflag,output,trials] = surrogateopt(optimfun, lb,ub, options);
save tmpOpt3
%% pattern search
% Run the pattern search optimization
% Create optimization options with suitable settings for high-dimensional optimization
options = optimoptions('patternsearch', ...
                       'UseParallel', true, ...
                       'Display', 'iter', ...
                       'PollMethod', 'GSSPositiveBasis2N', ...  % Polling method
                       'SearchMethod', 'MADSPositiveBasis2N', ... % Search method
                       'MaxIterations', 50, ...
                       'MaxFunctionEvaluations', 2000, ...
                       'MeshTolerance', 1e-3);  % Tolerance for the mesh size

[x, fval] = patternsearch(optimfun, g, [], [], [], [], lb, ub, [], options);



%%
clf;

params0.g = ones(length(params0.mods), 1);
params0.mods = optparams.mods;
params0.PlotEachSeparately = true;
% params0.RunSlackSegments = 'All';
% params0.RunForceVelocity = false;
% params0.RunForceLengthEstim = false;
% params0.ShowStatePlots = true;
% params0.justPlotStateTransitionsFlag = false;
%%
clear;
clf;
load tmpOpt3.mat
% load tmpOpt2.mat
%% load tmpOpt3.mat
clf;
params0.BreakOnODEUnstable = false;
params0.useHalfActiveForSR = false;
params0.RunForceVelocity = true;
params0.RunForceLengthEstim = false;
params0.PlotEachSeparately = true;
params0.EvalFitSlackOnset = true;
params0.ShowStatePlots = true;
params0.justPlotStateTransitionsFlag = false;

 
params0.kd = 10;
params0.alpha0 =100;
params0.alpha0_L = 100;
params0.alpha0_R = 10;
% params0.FudgeVmax = true;
tic
RunBakersExp;
toc
E
sum(E)

% writeParamsToMFile('ModelParams_surro.m', params0);
%%
% T = table('VariableNames',["mod", "lb", "val", "ub"])
% T = table()
clear t_mods t_val;
for i = 1:length(optparams.mods)
    t_mods(i) = string(optparams.mods{i});
    t_val(i) = optparams.(optparams.mods{i});
end

T = table(t_mods', t_val'*0.1, t_val', t_val'*10);
T.Properties.VariableNames = ["mods", "lb", "val", "ub"];

% manually edit the values
save("optTable2.mat", "T");

%% Experiment run
% params0.g = p_OptimGA;
% params0.RunSlackSegments = 'Fourth-rampuponly';
params0.RunSlackSegments = 'All';
params0.EvalFitSlackOnset = true;
params0.PlotEachSeparately = true;
params0.ShowStatePlots = true;
params0.justPlotStateTransitionsFlag = false;

% params0.Lsc0 = 1.51;
% params0.xrate = 1;

% params0.EvalFitSlackOnset = true;
% params0.drawForceOnset = true;

params0.RunForceLengthEstim = false;
params0.RunForceVelocity = true;
params0.RunForceLengthEstim = false;
% experimenting
% params0.FudgeVmax = true;
% params0.FudgeA = 0;
% params0.FudgeB = 121.922;
% params0.FudgeC = -209.18;

tic
clf;
LoadData;
RunBakersExp
toc
sum(E)

% writeParamsToMFile('ModelParamsInit_FudgedSlack.m', params0);

%%
%% Running GA

% parpool
ga_Opts = optimoptions('ga', ...
    'PopulationSize',25, ...            % 250
    'Display','iter', ...
    'MaxStallGenerations',6, ...  % 10
    'UseParallel',false, 'PlotFcn','gaplotbestf');

[p_OptimGA,Res_OptimGA,~,~,FinPopGA,FinScoreGA] = ...
    ga(optimfun,size(g, 2), ...
    [],[],[],[],ones(size(g))*0.01,ones(size(g))*20,[],ga_Opts);

%%
writeParamsToMFile('ModelParamsOptim_tf2_slackLast.m', params0);
writeParamsToMFile('ModelParamsOptim_tf2_slackFirst.m', params0);
writeParamsToMFile('ModelParamsOptim_tf2_slackFirstLast.m', params0);
writeParamsToMFile('ModelParamsOptim_tf2_slackFirstLast_LeftOnly.m', params0);
writeParamsToMFile('ModelParamsOptim_tf2_slackOnsetAll_LeftOnly.m', params0);
writeParamsToMFile('ModelParamsOptim_tmp.m', params0);
writeParamsToMFile('ModelParamsOptim_DtKtr.m', params0);
writeParamsToMFile('ModelParamsOptim_DtKtr_OV.m', params0);
writeParamsToMFile('ModelParamsOptim_DtKtr_OV2.m', params0);
writeParamsToMFile('ModelParamsOptim_fudgeSlack.m', params0);
writeParamsToMFile('ModelParamsOptim_fudgeSlackVelocities.m', params0);
writeParamsToMFile('ModelParamsOptim_FSV_ForceOnset.m', params0);
writeParamsToMFile('ModelParamsOptim_FSV_ForceOnset_tmp.m', params0);
writeParamsToMFile('ModelParamsOptim_DoublePeakOptim.m', params0);
clf;RunBakersExp;
%% show

figure(8340);clf;
LoadData;
params0.ksr = params0.ksr0;
params0.kmsr
% params0.L_thick = 1.67; % Length of thick filament, um
% params0.L_hbare = 0.10; % Length of bare region of thick filament, um
% params0.L_thin  = 1.20; % Length of thin filament, um

% params0.mods = {};
% params0.Lsc0    = 1.51;
% params0.RunForceVelocity = false;
% params0.xrate = 3;

% params0.kmsr = params0.kmsr*10;
% params0.ksr0 = params0.ksr0*2;
% params0.gamma = 4;
% params0.Lsc0 = 1.51;
% params0.k_pas = 50;
params0.RunSlack = true;
% params0.RunForceVelocity = false;
% params0.RunForceVelocityTime = false;
% params0.PlotEachSeparately = true;
% params0.justPlotStateTransitionsFlag = false;
% params0.RunSlackSegments = 'FirstAndLast';
params0.RunSlackSegments = 'FirstAndLastExtended';
params0.RunSlackSegments = 'All';
params0.UsePassiveForSR = false;

% params0.ShowStatePlots = true;
params0.drawForceOnset = true;

% params0.UseDirectSRXTransition = false;true

% params0.kSE = 5.245151e3;
% params0.SL0 = 2.2*rsl0;

% params0.mu = 1e-2;
params0.UseDirectSRXTransition = false;
 
% rather extreme parametrization - max rate in, max dependence out
% a bit weird shape - the force depends on force
% from SR
params0.ksr0 = 100/20;
params0.sigma1 = 10e0;
% to SR
params0.kmsr = 6000/20*2;
params0.sigma2 = 1e6;

% Other way around
% from SR
params0.ksr0 = 100;
params0.sigma1 = 1e6;
% to SR
params0.kmsr = 5000;
params0.sigma2 = 20e0;
% 
% Other way around - less steep
% from SR
% params0.ksr0 = 100;
% params0.sigma1 = 1e6;
% to SR
% params0.kmsr = 600;
% params0.sigma2 = 40e0;
% 
% serial setup
% from SR
% params0.ksr0 = 100;
% params0.sigma1 = 40e0;
% to SR - backflow
% params0.kmsr = 500;
% params0.sigma2 = 1e6;
% params0.UseDirectSRXTransition = true;

% Dans tests
% out of SR
params0.ksr0 = 10;
params0.sigma1 = 50e6;
% to SR
params0.kmsr = 30;
params0.sigma2 = 1e6;
params0.UseDirectSRXTransition = false;

params0.ka = 25;
params0.kd = 1;

%
F_SR = 80;PT = 0.15;
RSR2PT = params0.ksr0*exp(F_SR/params0.sigma1);
RPT2SR = params0.kmsr*exp(-F_SR/params0.sigma2);
r_SR = RPT2SR/RSR2PT*PT
%
params0.justPlotStateTransitionsFlag = false;
% params0.UseTitinInterpolation = false;
params0.ghostSave = '';
% params0.ghostSave = 'tmp';
params0.ghostLoad = 'tmp';
% params0.alpha2_L = -100;
% params0.mu = 1e-1;
% params0.kstiff1 =7.2e3; 
% params0.kstiff2  = 1.5e4;

% params0.FudgeVmax = true;
% params0.FudgeA = 0;% 2.1750e+03;
% %
% params0.FudgeB = 120;
% params0.FudgeC = -210;

% SL = 1.7:0.1:2.2;
% velHS = -(params0.FudgeB*SL + params0.FudgeC);
% plot(SL, velHS);

%
% params0.g = [1 1];
% for fudgeB and C x =    0.8414    1.0641
%
% params0.g = [1.4747    0.9379]
params0.FudgeVmax = true;
% params0.justPlotStateTransitionsFlag = false;
tic
RunBakersExp;
toc
sum(E)
%% plot overlay
figure(222);clf;
dropstart = velocitytable([3, 7, 11, 15, 19], 1);

% sr
% plot(out.t-dropstart, out.SR, '-', LineWidth=1);    hold on; set(gca,'ColorOrderIndex',1);

% p1
plot(out.t-dropstart, out.p2_0, '-', LineWidth=1);     hold on; set(gca,'ColorOrderIndex',1);

xlim([0 0.3])
%%
LXB = 1.6:0.1:3.8;
    L_thick = params.L_thick;% = 1.67; % Length of thick filament, um
    L_hbare = params.L_hbare;% = 0.10; % Length of bare region of thick filament, um
    L_thin = params.L_thin;  %= 1.20; % Length of thin filament, um

    % deltaR  = 0.010; % um    
    L_T_HS1 = min(L_thick*0.5, LXB*0.5);
    L_T_HS2 = max(LXB*0.5 - (LXB-L_thin),L_hbare*0.5);
    L_ov = L_T_HS1 - L_T_HS2; % Length of single overlap region
    N_overlap = L_ov*2/(L_thick - L_hbare);

hold on;plot(LXB, N_overlap.^2, '-', LineWidth=2)
xlabel('SL');ylabel('OV');legend('OV', 'OV^2')



function cost = isolateRunBakersExp(params0)
LoadData;
    RunBakersExp;
    cost = sum(E);
end