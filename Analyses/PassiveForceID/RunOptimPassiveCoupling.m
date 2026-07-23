% RunOptimPassiveCoupling.m
% =========================================================================
% Joint PASSIVE + ACTIVE(fourth-slack) feature optimizer. Same engine as
% RunOptimFull (optimizeFeatures: random-subset block-coordinate descent),
% but with the coupling objective:
%   - PASSIVE whole sim (ka=0 titin/SE mechanics): restretchPeak(5),
%     rampupRMSE(5), holdDecayRMSE(5), steady22, steady20
%   - ACTIVE fourth slack: slackLSQE  (per-slack time-based trace RMSE)
% All features go through evalFeatureCost (OptimizeOn='Feats'), so they are
% parametrised exactly like every other feature.
%
% Feature weights are set so every group contributes EQUALLY at the seed
% (weight = 1/cost_i(g=1)); no group dominates the descent at the start.
%
% Rationale: the passive fit is non-unique (sloppy k_pas/gamma/Lsc0/kSE/
% kSE_M/eta_M manifold). Scoring passive AND the active fourth slack in one
% cost lets the active data pin the passive soft directions.
% See Model/experiments/runPassiveExperiment.m and runSlackExperiment.m.
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));
%%
cfg = struct();
cfg.baseSnap = 'params/optfull_FourthSlackTimebase_opt.m';   % seed (active fit)
cfg.baseSnap = 'params/passiveCoupling_opt.m';   % seed (active fit)
cfg.tag      = 'passiveCoupling';                            % -> params/passiveCoupling_opt.m

% ---- base params (sanctioned isolated load) ----
params0 = getParams(loadParams(cfg.baseSnap), [], true, false);
params0.mods = {}; params0.g = [];

% ---- experiments: passive (whole) + active fourth slack, feature-scored ----
params0.RunForceVelocity = false; params0.RunForceVelocityTime = false;
params0.RunKtr = false; params0.RunStairs = false; params0.RunForceLengthEstim = false;
params0.RunSlack = true;  params0.RunSlackSegments = 'AllPar';
params0.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
params0.RunSlackPassive = true;
params0.passive_velocitytableonfile = 'protocol_03_27_2026_ActivePNBMava_slack.mat';
params0.EvalFeatures = false;   % 'Feats' mode reads params0.fn directly
params0.EvalPeaks    = false;
params0.OptimizeOn   = 'Feats';
params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0;
params0.BreakOnODEUnstable = false;
params0.MaxRunTime = 60;

% ---- feature groups (model<->data comparable after slack trimming) ----
% PS_* = passive-slack features (runPassiveExperiment); slackLSQE = active Fourth slack.
baseFn = {'PS_restretchPeak|2', 'PS_rampupRMSE|0.03', 'PS_holdDecayRMSE|0.04', 'PS_steady22|2', 'PS_steady20|1', 'slackLSQE|0.2'};

% ---- weight so each group contributes equally at the seed ----
params0.fn = baseFn;
tic
[fn, w0, c0] = iWeightEqual(params0, baseFn);
toc
% params0.fn = fn;
fprintf('[%s] initial per-feature cost / weight (equalised):\n', cfg.tag);
for i = 1:numel(baseFn)
    fprintf('   %-14s cost=%9.4f  weight=%.5g\n', baseFn{i}, c0(i), w0(i));
end

% ---- tunable pool: passive/mechanical + active levers ----
cand = { ...
    'k_pas','gamma','Lsc0','kSE','ekSE','kSE_M','eta_M','mu','mu_neg', ...   % passive / SE / coupling
    'ka','kd','k1','k2','kah','kamh','kstiff1','kstiff2','dr','dr1','dr2', ...% active core
    'ksr0','kmsrd','sigma1','sigma2', ...                                    % SRX manifold
    'A0','v_ref_reg','tau_reg'};                                            % registration availability
keep = true(1, numel(cand));
for i = 1:numel(cand); keep(i) = (iSeedVal(params0, cand{i}) ~= 0); end
cfg.pool = cand(keep);

% ---- compulsory: coupling-critical mechanics + force/rate rebalancers ----
cfg.compulsory = {'kSE_M','eta_M','k_pas','kstiff2','k2'};

% ---- engine budgets ----
cfg.fn            = fn;        % required by optimizeFeatures (params0.fn is what Feats uses)
cfg.N_DRAW        = 12;
cfg.SIMPLEX_EVALS = 90;
cfg.SURR_EVALS    = 120;         % pure fminsearch
cfg.MaxRunTime    = 60;
cfg.DEBUG           = false;   % set true for a ~1-round smoke test
cfg.TIME_BUDGET_HRS = 12;
cfg.RESUME          = false;
cfg.params0         = params0;
%%
fprintf('[%s] pool=%d | compulsory={%s}\n', cfg.tag, numel(cfg.pool), strjoin(cfg.compulsory, ','));
optimizeFeatures(cfg);

%% ---- local helpers -----------------------------------------------------
function [fn, w, c] = iWeightEqual(params0, baseFn)
%IWEIGHTEQUAL Run the seed once and set weights so each feature-group's
%   weighted cost == 1 at g=1 (comparable at initial). costExp=2 matches the
%   'Feats' objective (evaluateBakersExp -> evalFeatureCost, 3-arg default).
    params0.PlotEachSeparately = true;
    RunBakersExp;   % populates features_model / features_data at the seed
    c = evalFeatureCost(features_data, features_model, baseFn);
    w = 1 ./ max(c, 1e-6);
    fn = cell(1, numel(baseFn));
    for i = 1:numel(baseFn); fn{i} = sprintf('%s|%.6g', baseFn{i}, w(i)); end
    plotFeatures(features_data, features_model, [], baseFn);
end

function v = iSeedVal(p, name)
%ISEEDVAL Read a scalar/array-element param value by name (handles Field__N).
    us = strfind(name, '__');
    if isempty(us)
        if isfield(p, name); v = p.(name); else; v = 0; end
    else
        f = name(1:us(1)-1); idx = str2double(name(us(1)+2:end));
        if isfield(p, f); arr = p.(f); v = arr(idx); else; v = 0; end
    end
end
