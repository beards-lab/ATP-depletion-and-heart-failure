% test sensitivity 
params0 = getParams();
% ModelParamsInitNiceSlack;
% ModelParamsInit2;


cost_sap = []; % SA plus
cost_sam = []; % SA minus

% ModelParamsInitOptim_slack4
ModelParamsOptim_tmp
% ModelParamsInitOptim_slackAll
% ModelParamsOptim_tf2_slackLast
params0.RunSlackSegments = 'All';
% params0.Lsc0 = 1.51;
% params0.e2R = 1;
% ModelParamsOptim_tf2_slackFirst
% ModelParamsOptim_tf2_slackFirstLast
% ModelParamsOptim_tmp.m

params0.drawPlots = true;
params0.drawForceOnset = true;
params0.PlotEachSeparately = true;
params0.ShowResidualPlots = false;
params0.justPlotStateTransitionsFlag = false;


% params0.ghostLoad = 'NiceFit_slack4';

RunBakersExp;
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
params0.ShowResidualPlots = true;
c0 = isolateRunBakersExp(params0);
%%
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

% reduced
params0.mods = {'dr1', 'alpha1', 'alpha2_L', 'k2', 'dr2', 'e2L', 'kd', 'ksr0', 'kmsr', 'kstiff1', 'kstiff2', 'k_pas', 'gamma', 'Lsc0', 'kSE'};

params0.mods = {'k_pas', 'gamma', 'Lsc0', 'kSE'};

% params0.mods = {'k_pas', 'gamma', 'Lsc0'};

params0.g = ones(size(params0.mods));
saSet = 1:length(params0.mods);

%%
params0.ghostLoad = '';
p0 = params0;
SAFact = 1.05;
for i_m = saSet
    % reset params
    params0 = p0;
    params0.g(i_m) = params0.g(i_m)*SAFact;
    fprintf('Mod %s is up to %g %%..', params0.mods{i_m}, params0.g(i_m)*100);
    % figure(i_m)
    cost = isolateRunBakersExp(params0);
    cost_sap(i_m) = cost;
    % 
    % params0.g(i_m) = params0.g(i_m)/SAFact*(1+ 1-SAFact);
    % fprintf('costing %1.4e€ and down to %g...', cost, params0.g(i_m)*100);
    % cost = isolateRunBakersExp(params0);
    % cost_sam(i_m) = cost;

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
%% show
clf;
figure;
% params0.mods = {};
% params0.Lsc0    = 1.51;
% params0.RunForceVelocity = false;
params0.RunSlack = true;
params0.RunForceVelocity = false;
params0.RunForceVelocityTime = false;
params0.PlotEachSeparately = true;
params0.justPlotStateTransitionsFlag = false;
params0.RunSlackSegments = 'FirstAndLast';
params0.RunSlackSegments = 'All';
params0.ShowStatePlots = true;
RunBakersExp;
sum(E)



function cost = isolateRunBakersExp(params0)

    RunBakersExp;
    cost = sum(E);
end