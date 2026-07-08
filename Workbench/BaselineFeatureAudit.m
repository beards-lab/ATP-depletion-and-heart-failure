% BaselineFeatureAudit.m
% =========================================================================
% One-shot baseline feature audit for params_reseeded_opt (the current best
% snapshot used by RunBoundedFit_Optim). NO optimization is run.
%
% Purpose: dump data-vs-model for EVERY feature the extractors produce — not
% just those currently in params0.fn — so we can decide which unused features
% to add to the optim, which to rule out, and how to reweight.
%
% Run:  cd(root); addpath(genpath('.')); BaselineFeatureAudit
% =========================================================================

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root));
LoadData;

%% Load the current best baked snapshot (self-contained: all flags + fn)
params0 = getParams();
run(fullfile(root, 'params', 'params_reseeded_opt.m'));
params0.mods = {}; params0.g = [];
params0 = getParams(params0, [], true, false);   % rebuild init state from SL0

% one-shot baseline: serial slack (no pool needed), no per-tile plots
params0.RunSlackSegments = 'All';
params0.PlotEachSeparately = false;
params0.PlotFeatureFitting = false;
params0.BreakOnODEUnstable = false;

%% Run once
tic; RunBakersExp; toc

%% (1) Cost breakdown + full data-vs-model dump for the ACTUAL fn
fprintf('\n############ BASELINE: actual params0.fn ############\n');
reportFeatureCost(features_data, features_model, params0.fn, 2);

%% (2) Augmented audit: score EVERY data-fit feature (weight 1) so we can see
%      what each currently-unused feature would contribute, data vs model.
fn_all = { ...
    'FV_fnorm|FV_v', 'FV_f|FV_v', ...                              % FV
    'ktr', 'A', 'Am', 'ktr_rmse|0-0.2', 't0', 't0_crossing|SLdiff', ... % recovery/onset
    'ktr2_k', 'ktr2_overshoot', 'ktr2_rmse|0-0.2', ...             % oscillatory recovery
    'peak1_y', 'peak1_t', 'peak1_dSL', 'peak1_SL', ...             % restretch peak
    'vall_y', 'vall_t', 'peak2', ...                               % valley after peak
    'restretchSlopeStart', 'restretchSlopeEnd', ...                % restretch stiffness
    'steady', 'vall2_dy', 'vall2_t', 'ovrsht_dy', 'ovrsht_t', ...  % steady + slow under/overshoot
    'SLslack', 'SLdiff' ...                                        % SL geometry
};
fprintf('\n############ AUDIT: all data-fit features (unit weight) ############\n');
reportFeatureCost(features_data, features_model, fn_all, 2);

%% (3) plotFeatures figure for the actual fn
fig = figure(777); clf;
set(fig, 'Visible', 'on', 'Position', [60 60 1500 850]);
plotFeatures(features_data, features_model, [], params0.fn);
drawnow;
figDir = fullfile(root, 'Docs', 'figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end
exportgraphics(fig, fullfile(figDir, 'baseline_feature_audit.png'), 'Resolution', 150);
fprintf('\nSaved Docs/figures/baseline_feature_audit.png\n');
