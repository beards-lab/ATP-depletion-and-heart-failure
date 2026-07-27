root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root));

params0 = getParams();
% opt2state_opt;

%%
% optfull4_opt
% optfull_FourthSlackTimebase_opt
optfull6_opt
fn = params0.fn;
% fn = {'FV_fnorm|FV_v|10', 'ktr|8', 'A|50', 'ktr_rmse|0-.5|.1', 'XTOR[1]|3-15|15,XTOR_vmax[1]|5-25,SRX_ss[1]|0.01-0.40|5,attached_ss[1]|0.10-0.45|10,PT_ss[1]|0.01-0.70', 't0_crossing|SLdiff|2', 'restretchSlopeStart|1', 'peak1_y|10', 'peak1_dSL|1', 'vall_y|10', 'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|1.0', 'ovrsht_dy|1', 'k2|100-1500|0.01,ksrd|0.-20|0.001,kmsrd|0.0-20|0.001', 'doublePeak|10', 'coolDownLS|0.05'};
% passiveCoupling_opt
% ModelParamsWPassive_NewOpt3
% optfull_FourthSlackTimebase_opt
% optfull5_opt

params0.PlotEachSeparately = 1;
figure(231);clf;
params0.MaxRunTime = 320;
params0.FV_velocities = -[0, 0.5, 1, 2, 3, 4, 5, 6];
params0.RunForceLengthEstim = false;
params0.RunForceVelocityTime = false;
params0.RunForceVelocity = true;
params0.RunSlack = true;
params0.RunSlackSegments = 'AllPar';
params0.RunKtr = true;
params0.RunStairs = true;
params0.RunSlackPassive = true;
params0.passive_velocitytableonfile = 'protocol_03_27_2026_ActivePNBMava_slack.mat';
params0.EvalFeatures = true;

%%
% fn = {'FV_fnorm|FV_v|10', 'ktr|10', 'A|50', 'ktr_rmse|0-.5|.1', 'XTOR[1]|3-15|15,XTOR_vmax[1]|5-25,SRX_ss[1]|0.01-0.40|5,attached_ss[1]|0.10-0.45|10,PT_ss[1]|0.01-0.70', 't0_crossing|SLdiff|2', 'restretchSlopeStart|1', 'peak1_y|10', 'peak1_dSL|1', 'vall_y|10', 'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|1.0', 'ovrsht_dy|1', 'k2|100-1500|0.01,ksrd|0.-20|0.001,kmsrd|0.0-20|0.001', 'doublePeak|10', 'coolDownLS|0.05'};
% PS_fn = {'PS_restretchPeak', 'PS_rampupRMSE|0.005', 'PS_holdDecayRMSE|0.01', 'PS_steady22|0.5', 'PS_steady20|0.5'};
% params0.fn = [fn, PS_fn];
%%
% params0.velocitytableonfile = 'protocol_03_27_2026_ActivePNBMava_slack.mat';
% params0.ka = 0;
tic
RunBakersExp;
toc

% w      = 1 ./ max(evalFeatureCost(features_data, features_model, PS_fn), 1e-6);
% PS_fnw = arrayfun(@(i) sprintf('%s|%.4g', PS_fn{i}, w(i)), 1:numel(PS_fn), 'UniformOutput', false);
% params0.fn = [params0.fn, PS_fnw];   % each PS_ term now contributes ~1 at the seed
% write
%%
if ~exist('features_ghost', 'var')
    features_ghost = features_model;
end
sum(plotFeatures(features_data, features_model, features_ghost, params0.fn))

%%
params0.RunSlackSegments = 'AllPar';
params0.MaxRunTime = 120;
% opt2state_v2_opt
params_2state_a2hop

% ModelParams_coupled_lowATP
% params_3state_seedB
% params0.k2 = 95/2;
params0.PlotEachSeparately = 1;
params0.FV_velocities = -[0, 0.5, 1, 2, 4];
figure;

tic
params0.justPlotStateTransitionsFlag = true;
RunBakersExp;
toc

%%
writeParamsToMFile('ThursdayNightFever', params0)