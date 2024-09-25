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
ModelParamsOptim_DtKtr_OV
% ModelParamsOptim_DtKtr_OV2

params0.RunSlackSegments = 'All';
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
params0.ShowResidualPlots = false;
params0.justPlotStateTransitionsFlag = false;


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
params0.vmax = 18;
params0.kSE = 1.3e4;
params0.kSEn = 1.3e1;
params0.mu = 1e-1;


params0.ghostLoad = 'DtKtr_OV';
tic
RunBakersExp;
toc
sum(E)
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
params0.mods = {'k_on', 'k_off', 'vmax', 'kah', 'kadh', 'ka', 'kd', 'k1', 'k_1', 'k2', 'k_2', 'k3', 'ksr0', 'sigma0', 'kmsr', 'K_Pi', 'K_T1', 'dr', 'kstiff1', 'kstiff2', 'kstiff3', 'K_T3', 'K_D', 'alpha0', 'alpha1', 'alpha_1', 'alpha2', 'alpha3', 's3', 'alphaRip', 'k2rip', 'dr0', 'dr_1', 'dr2', 'dr3', 'Amax', 'mu', 'kSE', 'k_pas', 'k2_L', 'ss', 'TK', 'TK0', 'dr1', 'sigma1', 'sigma2', 'gamma', 'alpha2_L', 'k2_R', 'dr2_R', 'alpha2_R', 'e2R', 'Lsc0', 'e2L', 'xrate'}
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
params0.mods = {'dr1', 'alpha1', 'k1', 'alpha2_L', 'k2', 'dr2', 'e2L', 'kd', 'ksr0', 'kmsr', 'sigma1', 'kstiff1', 'kstiff2', 'kSE', 'L_thick', 'L_hbare', 'L_thin'};

params0.mods = {'L_thick', 'L_hbare', 'L_thin'};

params0.mods = {'vmax', 'kSE'}

params0.mods = {'FudgeB','FudgeC'};

params0.g = ones(size(params0.mods));
saSet = 1:length(params0.mods);

%%
params0.ghostLoad = '';
p0 = params0;
SAFact = 1.01;
c0 = isolateRunBakersExp(params0);
params0.PlotEachSeparately = false;
params0.ShowStatePlots = false;
params0.ShowResidualPlots = false;
for i_m = saSet
    % reset params
    params0 = p0;
    params0.g(i_m) = params0.g(i_m)*SAFact;
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

%%
% c0 - min(cost_sam, cost_sap)
cost_cmb = min(cost_sam, cost_sap);
% cost_cmb = cost_sap;
cost_diff = c0 - cost_cmb
% better value
[val_bett cost_diff_better_i] = sort(cost_diff.*(cost_diff > 0), 'descend');
% most influential
[val_infl cost_diff_infl_i] = sort(cost_diff.*(cost_diff < 0), 'ascend');

param_ord = [cost_diff_better_i(val_bett > 0) cost_diff_infl_i(val_infl < 0)]
clf;
bar(cost_cmb(param_ord));xticks(1:length(params0.mods))
xticklabels(params0.mods(param_ord))

%% select where it can be reduced OR most influential
params0.mods = params0.mods(param_ord);
params0.mods = params0.mods(1:8);
params0.g = ones(size(params0.mods));
% saSet = 1:length(params0.mods);
%% Optim

options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.1, 'PlotFcns', @optimplotfval, 'MaxIter', 1500);
% g = [1, 1, 1, 1, 1, 1, 1, 1];
% g = [1.2539    0.4422];

params0.ShowResidualPlots = false;

g = params0.g;
optimfun = @(g)evaluateBakersExp(g, params0);
x = fminsearch(optimfun, params0.g, options)
params0.g = x;

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
%% show

figure(8340);clf;
LoadData;

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
params0.RunSlackSegments = 'All';
% params0.ShowStatePlots = true;
% params0.drawForceOnset = true;
% params0.UseDirectSRXTransition = false;

params0.kSE = 5.245151e3;
% params0.SL0 = 2.2*rsl0;
params0.vmax = 10;
params0.vmax1 = 20;
params0.mu = 1e-1;
% params0.sigma1 = 1e6;
% params0.sigma2 = 3e1;
% params0.ksr0 = 50;
% params0.kmsr = 50.78;
% params0.mu = 1e-1;
% params0.kstiff1 =7.2e3; 
% params0.kstiff2  = 1.5e4;
% 
params0.FudgeVmax = true;
params0.FudgeA = 0;% 2.1750e+03;
%
params0.FudgeB = 120;
params0.FudgeC = -210;

SL = 1.7:0.1:2.2;
velHS = -(params0.FudgeB*SL + params0.FudgeC);
plot(SL, velHS);
%
% params0.g = [1 1];
% for fudgeB and C x =    0.8414    1.0641
%
% params0.g = [1.4747    0.9379]
tic
RunBakersExp;
toc
sum(E)



function cost = isolateRunBakersExp(params0)

    RunBakersExp;
    cost = sum(E);
end