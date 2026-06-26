%% RunK2Frontier.m
% =====================================================================
% LOW-ATP FIRST FRONTIER: can reducing k2 (ADP-release / detachment = gapp)
% reproduce the measured 2 mM-vs-8 mM cross-bridge signature?
%
% Method (decided with the user, 2026-06-25):
%   * SCORING = relative ratio.  The 8 mM baseline fit is itself imperfect,
%     so we score the model's ATP RESPONSE  model(reduced k2)/model(base k2)
%     against the measured data ratio  features(2mM)/features(8mM).  Any
%     multiplicative baseline offset cancels.
%   * RUNDOWN = ignored (files used as-is) for this first pass.
%   * TARGET  = broad feature set (force gain + ktr + peak/onset/stiffness).
%
% Base parametrization: params/params_reseeded_regavail_opt2.m (current best,
% feature cost ~8.86). ATP lever = a single multiplicative scale on params.k2.
% The slack protocol / velocity table comes from the NEW 8 mM data file, so
% the model is driven by exactly the protocol the data were measured under.
%
% Run:  cd(root); addpath(genpath('.')); run Analyses/LowATP_k2Frontier/RunK2Frontier
% Outputs: results/k2_frontier.mat  +  results/k2_frontier.png
% =====================================================================
clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
addpath(genpath(root));
resdir = fullfile(here, 'results');
if ~exist(resdir, 'dir'); mkdir(resdir); end

%% -------- load current best, configure slack-only on the new 8 mM protocol
params0 = getParams();
run(fullfile(root, 'params', 'params_reseeded_regavail_opt2.m'));
params0 = getParams(params0, [], true, false);

params0.RunForceVelocity = false; params0.RunKtr = false; params0.RunStairs = false;
params0.RunForceVelocityTime = false; params0.RunForceLengthEstim = false;
params0.RunSlack = true; params0.EvalFeatures = true; params0.recalculateDataFeats = false;
params0.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
params0.PlotEachSeparately = 0; params0.PlotProbsOnFig = 0; params0.PlotFeatureFitting = false;
params0.justPlotStateTransitionsFlag = false; params0.BreakOnODEUnstable = false; params0.MaxRunTime = 120;
params0.ghostLoad = ''; params0.ghostSave = ''; params0.EvalPeaks = false; params0.EvalFitSlackOnset = false;
try
    if isempty(gcp('nocreate')); parpool('Threads', 5); end
    params0.RunSlackSegments = 'AllPar';
catch
    params0.RunSlackSegments = 'All';
end
baseK2 = params0.k2;

%% -------- experimental target ratios (2 mM / 8 mM), broad feature set
f2 = load(fullfile(root, 'data', 'protocol_03_27_2026_2mM_slack.mat')).features_data;
f8 = load(fullfile(root, 'data', 'protocol_03_27_2026_8mM_slack.mat')).features_data;
feats = {'steady','A','Am','peak2','peak1_y','vall_y','ktr','t0','restretchSlopeStart','ktr2_overshoot'};
nF = numel(feats);
dR = zeros(nF, 1);
for i = 1:nF
    dR(i) = mean(double(f2.(feats{i})), 'omitnan') / mean(double(f8.(feats{i})), 'omitnan');
end

%% -------- k2 sweep (the "ATP decrease" axis)
scales = [1.0 0.8 0.65 0.5 0.4 0.3];
nS = numel(scales);
Mabs = nan(nF, nS);
fprintf('k2 frontier sweep (baseK2=%.2f)...\n', baseK2); tic
for j = 1:nS
    p = params0; p.k2 = baseK2 * scales(j);
    p = getParams(p, [], true, false);
    [~, ~, fm, ~] = runSlackExperiment(p);
    for i = 1:nF; Mabs(i, j) = mean(double(fm.(feats{i})), 'omitnan'); end
    fprintf('  k2*%.2f=%6.1f | steady=%.1f A=%.1f ktr=%.1f peak2=%.1f\n', ...
        scales(j), p.k2, Mabs(strcmp(feats,'steady'),j), Mabs(strcmp(feats,'A'),j), ...
        Mabs(strcmp(feats,'ktr'),j), Mabs(strcmp(feats,'peak2'),j));
end
fprintf('done in %.0f s.\n', toc);

%% -------- relative-ratio analysis
mR = Mabs ./ Mabs(:, 1);                  % model ratio vs baseline (scale=1)
LR = log(mR) - log(dR);                   % log-ratio residual (symmetric in over/under)
core = {'steady','A','peak2','peak1_y','ktr','t0'}; ic = find(ismember(feats, core));
costAll  = sum(LR.^2, 1);
costCore = sum(LR(ic,:).^2, 1);
[~, jb] = min(costCore);

% per-feature crossing scale (where model ratio meets its data target)
sStar = nan(nF, 1);
for i = 1:nF
    [u, iu] = unique(mR(i,:)); sc = scales(iu);
    if (dR(i)-min(u))*(dR(i)-max(u)) <= 0; sStar(i) = interp1(u, sc, dR(i)); end
end

fprintf('\nscale:        %s\n', num2str(scales, '%6.2f '));
fprintf('cost(broad):  %s\n', num2str(costAll, '%6.3f '));
fprintf('cost(core6):  %s\n', num2str(costCore,'%6.3f '));
fprintf('\nbest single k2 (core6): scale=%.2f (k2=%.1f), cost=%.3f\n', scales(jb), baseK2*scales(jb), costCore(jb));
fprintf('%-20s %7s %7s %9s\n', 'feature', 'dataR', sprintf('s=%.2f',scales(jb)), 's*match');
for i = 1:nF
    ss = '   --'; if ~isnan(sStar(i)); ss = sprintf('%5.2f', sStar(i)); end
    fprintf('%-20s %7.3f %7.3f %9s\n', feats{i}, dR(i), mR(i,jb), ss);
end

%% -------- frontier figure
fig = figure('Position',[100 100 1100 650]); tl = tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Low-ATP k2 frontier: model ratio vs data target (2mM/8mM)   [best core k2=%.2f]', scales(jb)));
for i = 1:nF
    nexttile; hold on; box on;
    plot(scales, mR(i,:), 'o-', 'LineWidth', 1.8, 'Color', [0.1 0.3 0.8]);
    yline(dR(i), '--', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.5);
    if ~isnan(sStar(i)); xline(sStar(i), ':', 'Color', [0.3 0.6 0.3]); end
    xline(scales(jb), '-', 'Color', [0.6 0.6 0.6]);
    set(gca, 'XDir', 'reverse');  % decreasing ATP -> left to right
    xlabel('k2 scale'); title(feats{i}, 'Interpreter','none');
    if i == 1; legend({'model ratio','data target','crossing','best k2'}, 'FontSize', 6, 'Location','best'); end
    ylim padded;
end
exportgraphics(fig, fullfile(resdir, 'k2_frontier.png'), 'Resolution', 130);

%% -------- save
save(fullfile(resdir, 'k2_frontier.mat'), 'scales', 'Mabs', 'mR', 'dR', 'feats', ...
     'core', 'ic', 'LR', 'costAll', 'costCore', 'jb', 'sStar', 'baseK2');
fprintf('\nSaved results/k2_frontier.mat + k2_frontier.png\n');
