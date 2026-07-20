root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root));

params0 = getParams();
% opt2state_opt;

%%
% optfull4_opt
optfull_FourthSlackTimebase_opt
ModelParamsWPassive_NewOpt3

params0.PlotEachSeparately = 1;
figure(231);clf;
params0.MaxRunTime = 320;
% params0.FV_velocities = -[0, 0.5, 1, 2, 3, 4, 5, 6];
params0.RunForceLengthEstim = false;
params0.RunForceVelocityTime = false;
params0.RunForceVelocity = false;
params0.RunSlack = true;
params0.RunSlackSegments = 'All';
params0.RunKtr = false;
params0.RunStairs = false;
% params0.velocitytableonfile = 'protocol_03_27_2026_ActivePNBMava_slack.mat';
% params0.ka = 0;
tic
RunBakersExp;
toc

%%
if ~exist('features_ghost', 'var')
    features_ghost = features_model;
end
plotFeatures(features_data, features_model, features_ghost, params0.fn)

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
% writeParamsToMFile('SaturdayNightFever', params0)