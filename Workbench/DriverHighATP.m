% DriverHighATP.m
% =====================================================================
% Run the current BEST high-ATP (8 mM) parametrization and report its fit to the
% 8 mM slack data.
%
%   PARAM FILE : params/params_reseeded_regavail_opt2.m
%                (current best; registration-availability config; feature cost
%                ~8.86 on the legacy target. NumberOfStates = 2.)
%
% This is the baseline the low-ATP mechanism is built on (see DriverLowATP.m).
%
% Run:  cd(root); addpath(genpath('.')); DriverHighATP
% See also: DriverLowATP, Analyses/LowATP_Mechanisms/REPORT.md, Auxiliary/loadParams.m
% =====================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root));
refreshPool(5);                                   % workers pick up any edited Model code

PARAM_FILE = 'params_reseeded_regavail_opt2';     % <-- the high-ATP best

% loadParams sandboxes the file's internal `params0`, so the variable name here is safe.
params0 = getParams(loadParams(PARAM_FILE), [], true, false);
params0.RunForceVelocity=false; params0.RunKtr=false; params0.RunStairs=false;
params0.RunForceVelocityTime=false; params0.RunForceLengthEstim=false; params0.RunSlack=true;
params0.EvalFeatures=true; params0.recalculateDataFeats=false;
params0.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
params0.PlotProbsOnFig=0; params0.justPlotStateTransitionsFlag=false;
params0.BreakOnODEUnstable=false; params0.MaxRunTime=120;
params0.ghostLoad=''; params0.ghostSave=''; params0.EvalPeaks=false; params0.EvalFitSlackOnset=false;
params0.RunSlackSegments='All';                   % serial -> clean contiguous trace for the plot
params0.PlotEachSeparately=1; params0.ShowStatePlots=0; params0.ShowResidualPlots=0;
params0.drawForceOnset=0; params0.PlotFeatureFitting=0;

fprintf('High-ATP best = params/%s.m\n', PARAM_FILE);
figure(1); clf; tiledlayout('flow','TileSpacing','compact','Padding','compact');
% [~, out, fm, ~] = runSlackExperiment(params0);    % plots model (blue) vs data (black)
RunBakersExp;
sgtitle('High-ATP (8 mM) best: model vs data (slack-restretch)');

figure(232);
if ~exist(features_ghost)
    features_ghost = features_model;
end
cost = plotFeatures(features_data, features_model, features_ghost, params0.fn)


% feature comparison vs 8 mM data
f8 = load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
m  = @(x) mean(double(x),'omitnan');
fprintf('\n%-20s %10s %10s\n','feature','model','data');
for f = {'steady','A','Am','ktr','peak1_y','vall_y','peak2','restretchSlopeStart','t0'}
    fprintf('%-20s %10.3g %10.3g\n', f{1}, m(fm.(f{1})), m(f8.(f{1})));
end
fprintf(['\nKnown 8 mM residuals (see SYNTHESIS.md): ktr too fast (~66 vs 49),\n' ...
         'restretch peak overshoot (peak2 >> steady), restretchSlope low (kSE under-set).\n']);
