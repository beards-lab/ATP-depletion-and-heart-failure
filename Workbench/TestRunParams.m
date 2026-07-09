params0 = getParams();
opt2state_opt;
%%
params0.PlotEachSeparately = 1;
figure(231);clf;
params0.MaxRunTime = 120;
% params0.FV_velocities = -[0, 0.5, 1, 2, 3, 4, 5, 6];
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
RunBakersExp;
toc

%%
% writeParamsToMFile('SaturdayNightFever', params0)