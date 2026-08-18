% RunOptimRealign.m — the "slight reoptim" with compliant realignment switched on.
%
% Seed: params/rskR2_w025_opt.m + M1 at the setting that put rsK on the data
% (k_cr = 6, tau_cr = 20 ms, F_cr = 25) + its two hand-found compensators
% (kA2re = 30, eta_M x0.3). That point is rsK x1.12 / ktr x1.04 / steady 76.2 but
% L2 61.6 against a baseline 6.51, with the whole excess in restretch transient
% SHAPE. The question this run answers: how much of that 55 can a joint refit find,
% and does the optimiser KEEP the mechanism or walk k_cr back to zero?
%
% The objective is the UNCHANGED params0.fn from the baseline snapshot, so the
% resulting L2 is directly comparable to the 6.509 baseline. Nothing is reweighted
% to flatter the mechanism.
%
% Pool = the levers with a reason to be here, not the full 70-parameter set:
%   mechanism      k_cr tau_cr F_cr            (so it CAN be walked back)
%   series-visco   eta_M kSE_M mu_neg mu kstiff1_n kstiff2_n
%                  -- RestretchRecoveryFit 6 showed this class moves rsK and the
%                     restretch shape while leaving ktr and steady untouched, which
%                     is exactly the two things M1 already gets right
%   restretch      kA2re kA2hop sA2hop
%   SRX system     ksr0 ksr2srd kmsrd ksrd2sr sigma2 sigma_srd1
%                  -- included so the p2->SRXT-under-stretch route gets a fair
%                     joint test rather than a one-at-a-time one
%   routing        SRXFromR2HighStrain sSRXrip   (seeded NONZERO: a multiplicative
%                     modifier can never move a zero base)
%   core           k2 kstiff2 ka kah
%
% Usage:  matlab -batch "TIME_HRS=1; run('Analyses/RestretchSRXRecoil/RunOptimRealign.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('TIME_HRS', 'var'); TIME_HRS = 1.0; end
if ~exist('NWORKERS', 'var'); NWORKERS = 5;   end
if ~exist('TAG',      'var'); TAG      = 'realign1'; end
if ~exist('SEED',     'var'); SEED     = ''; end          % '' => build the seed below
% Features to remove from the objective, by primary field name.
if ~exist('DROP_FEATS','var'); DROP_FEATS = {}; end

params0 = getParams();
if isempty(SEED)
    run('params/rskR2_w025_opt.m');
    % ---- seed the mechanism at the point that reached the data -----------
    params0.UseLoadRealign = true;
    params0.k_cr   = 6;
    params0.tau_cr = 0.020;
    params0.F_cr   = 25;
    params0.Adef_max = 0.9;
    params0.kA2re  = 30;
    params0.eta_M  = 0.3 * params0.eta_M;
    % seeded small but NONZERO so the optimiser can explore it in both directions
    params0.SRXFromR2HighStrain = 0.2;
    params0.sSRXrip = 0.015;
    params0.dSRXrip = 0.005;
    params0.SRXRecoilToD = false;  % p2 -> SRXT under stretch (the slower route)
else
    run(SEED);                     % continue from a previous snapshot
end

params0.mods = {}; params0.g = [];

% ---- objective edits ------------------------------------------------------
% ovrsht_dy is dropped on request: SlackDataAnalysis 5 already lists it under
% "do not fit" (ratio CV 115 %, sign flips +1.35/+0.34/+0.37/-0.01) and
% recommendation 4 is to retire ovrsht_* from cost functions. The same window is
% still constrained from both ends -- coolDownLS (trace LSQE over the cool-down)
% and rsK (its first-order rate) -- so removing it does not leave the recovery
% unscored. NOTE this makes L2 NOT comparable to the 6.509 baseline number;
% the baseline must be re-scored under the same fn (done below).
if ~isempty(DROP_FEATS)
    prim = cellfun(@(s) regexprep(strtok(s,'|'), '\[.*$',''), params0.fn, 'uni', 0);
    drop = ismember(prim, DROP_FEATS);
    fprintf('dropping %d objective term(s): %s\n', nnz(drop), strjoin(params0.fn(drop), ', '));
    params0.fn = params0.fn(~drop);
end

% optimizeFeatures requires a seed snapshot on disk (provenance + resume), so
% freeze the hand-built seed rather than passing an in-memory struct only.
SEEDSNAP = fullfile('params', ['realign_seed_' TAG '.m']);
writeParamsToMFile(SEEDSNAP, params0, [], ...
    'seed for RunOptimRealign: rskR2_w025_opt + UseLoadRealign(k_cr=6,tau=20ms,F_cr=25) + kA2re=30 + eta_M x0.3 + SRXFromR2HighStrain=0.2');

cfg = struct();
cfg.baseSnap = SEEDSNAP;
cfg.tag = TAG;
cfg.fn  = params0.fn;
cfg.pool = { ...
    'k_cr','tau_cr','F_cr', ...
    'eta_M','kSE_M','mu_neg','mu','kstiff1_n','kstiff2_n', ...
    'kA2re','kA2hop','sA2hop', ...
    'ksr0','ksr2srd','kmsrd','ksrd2sr','sigma2','sigma_srd1', ...
    'SRXFromR2HighStrain','sSRXrip', ...
    'k2','kstiff2','ka','kah','kSE','estiff','k1','k_1','dr2'};
% every round must be able to trade the master rate, the master force scale and
% the mechanism's own gain against each other
cfg.compulsory = {'k_cr','eta_M','kstiff2'};

keep = true(1, numel(cfg.pool));
for i = 1:numel(cfg.pool)
    f = cfg.pool{i};
    keep(i) = isfield(params0, f) && params0.(f) ~= 0;
end
if any(~keep)
    fprintf('dropping zero-base params: %s\n', strjoin(cfg.pool(~keep), ', '));
    cfg.pool = cfg.pool(keep);
end

cfg.N_DRAW        = 8;
cfg.SIMPLEX_EVALS = 50;
cfg.MaxRunTime    = 300;
cfg.SURR_EVALS    = 0;
cfg.KICK_FRAC     = 0.1;
cfg.DEBUG         = false;
cfg.TIME_BUDGET_HRS = TIME_HRS;
cfg.RESUME = isfile(fullfile(root, 'params', [cfg.tag '_state.mat']));

params0.RunForceVelocity = true;  params0.RunKtr = false; params0.RunSlack = true;
params0.RunStairs = false; params0.RunForceVelocityTime = false;
params0.RunSlackPassive = true;
params0.EvalFeatures = true; params0.BreakOnODEUnstable = false;
params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0;
params0.RunSlackSegments = 'AllPar';
params0.MaxRunTime = cfg.MaxRunTime;
params0.OptimizeOn = 'Feats';
cfg.params0 = params0;

if isempty(gcp('nocreate')); parpool('local', NWORKERS); end

% Re-score the ORIGINAL baseline under the SAME (edited) objective, otherwise the
% 6.509 quoted elsewhere is not the number to beat.
pbase = getParams(loadParams('params/rskR2_w025_opt.m'), [], true, false);
pbase.fn = cfg.fn; pbase.mods = {}; pbase.g = ones(1,0);
pbase.RunForceVelocity = true; pbase.RunKtr = false; pbase.RunSlack = true;
pbase.RunStairs = false; pbase.RunForceVelocityTime = false; pbase.RunSlackPassive = true;
pbase.EvalFeatures = true; pbase.RunSlackSegments = 'AllPar';
pbase.MaxRunTime = cfg.MaxRunTime; pbase.OptimizeOn = 'Feats';
[cB, ~, ctB, nmB] = evaluateBakersExp(1, pbase, true, false);
fprintf('[%s] BASELINE re-scored under this objective = %.4f\n', cfg.tag, cB);

pchk = cfg.params0; pchk.mods = {}; pchk.g = ones(1,0);
[c0, ~, ct0, nm0] = evaluateBakersExp(1, pchk, true, false);
fprintf('[%s] seed cost = %.4f\n', cfg.tag, c0);

% Show the decomposition, biggest first, with the FV and passive terms proven
% present — the whole point of keeping RunForceVelocity/RunSlackPassive on.
fprintf('[%s] %-22s %10s %10s\n', cfg.tag, 'feature', 'baseline', 'seed');
[~, ordc] = sort(ct0, 'descend');
for q = ordc(1:min(14, numel(ordc)))
    fprintf('[%s] %-22s %10.3f %10.3f\n', cfg.tag, nm0{q}, ctB(q), ct0(q));
end
fvi = find(contains(nm0, 'FV'), 1);
if isempty(fvi)
    error('RunOptimRealign:noFV', ...
        'No FV term is being scored — RunForceVelocity or the FV data is missing.');
end
fprintf('[%s] FV IS SCORED: %s = %.3f (baseline %.3f)\n', cfg.tag, nm0{fvi}, ct0(fvi), ctB(fvi));
if ~isfinite(c0) || c0 > 1e4
    error('RunOptimRealign:badSeed', 'seed cost %.4g looks like a penalty', c0);
end
fprintf('[%s] pool=%d | N_DRAW=%d | budget %.2f h | resume=%d\n', ...
    cfg.tag, numel(cfg.pool), cfg.N_DRAW, cfg.TIME_BUDGET_HRS, cfg.RESUME);

optimizeFeatures(cfg);
