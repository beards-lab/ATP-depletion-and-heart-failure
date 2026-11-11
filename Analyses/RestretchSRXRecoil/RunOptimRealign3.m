% RunOptimRealign3.m — third refit: Maxwell tension-only ON, FV protected.
%
% Two changes from RunOptimRealign2:
%
% 1. `UseMaxwellTensionOnly = 1` (§10b). Worth ~-0.99 on its own — about what the
%    whole 4 h realign2 run bought — and it acts on rsK/coolDownLS, so it should
%    REDUCE the pressure to pay for rsK out of FV.
%    Note this makes `x_M_slack` a live parameter, so it joins the pool.
%
% 2. FV PROTECTED, two ways, because realign2 spent it monotonically
%    (FV_fnorm 0.518 -> 0.881 -> 0.922 across its last rounds, ending as the
%    LARGEST term in the objective while rsK fell 0.687 -> 0.633):
%
%    (a) WEIGHT: FV_fnorm 10 -> FV_WEIGHT (default 30). This changes the
%        objective, so the baseline is re-scored under it below and BOTH numbers
%        are printed — the previous 6.1025 is NOT comparable.
%
%    (b) LEVERS: realign2's pool contained no FV mechanism at all. `A0`,
%        `v_ref_reg` and `tau_reg` — the registration-availability parameters that
%        ARE the FV shoulder — were absent, as were the R2(s) shape knots. The
%        optimiser could therefore only ever pay FV, never repair it. That is the
%        more likely cause of the monotone decline than any property of the
%        mechanism, and it is fixed here.
%
% Usage:
%   matlab -batch "TIME_HRS=6; run('Analyses/RestretchSRXRecoil/RunOptimRealign3.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('TIME_HRS',  'var'); TIME_HRS  = 6.0; end
if ~exist('NWORKERS',  'var'); NWORKERS  = 5;   end
if ~exist('TAG',       'var'); TAG       = 'realign3'; end
if ~exist('FV_WEIGHT', 'var'); FV_WEIGHT = 30;  end
if ~exist('DROP_FEATS','var'); DROP_FEATS = {'ovrsht_dy'}; end

params0 = getParams(); run('params/realign2_opt.m');

% ---- change 1: the Maxwell fix -----------------------------------------
params0.UseMaxwellTensionOnly = 1;
params0.x_M_slack = 0.01;          % the value that tested best (0.003 was worse)

% ---- change 2a: reweight the objective ----------------------------------
fn   = params0.fn;
prim = cellfun(@(s) regexprep(strtok(s,'|'), '\[.*$',''), fn, 'uni', 0);
drop = ismember(prim, DROP_FEATS);
if any(drop)
    fprintf('dropping %d objective term(s): %s\n', nnz(drop), strjoin(fn(drop), ', '));
    fn = fn(~drop); prim = prim(~drop);
end
iFV = find(strcmp(prim, 'FV_fnorm'));
if isempty(iFV)
    error('RunOptimRealign3:noFV', 'FV_fnorm not found in params0.fn — cannot protect it.');
end
for k = iFV(:)'
    parts = strsplit(fn{k}, '|');
    old   = parts{end};
    parts{end} = num2str(FV_WEIGHT);
    fn{k} = strjoin(parts, '|');
    fprintf('FV term reweighted: %s  (was weight %s)\n', fn{k}, old);
end
params0.fn = fn;
params0.mods = {}; params0.g = [];

SEEDSNAP = fullfile('params', ['realign_seed_' TAG '.m']);
writeParamsToMFile(SEEDSNAP, params0, [], ...
    sprintf('seed for RunOptimRealign3: realign2_opt + UseMaxwellTensionOnly + FV weight %g', FV_WEIGHT));

cfg = struct();
cfg.baseSnap = SEEDSNAP;
cfg.tag = TAG;
cfg.fn  = fn;

% ---- change 2b: give the optimiser the FV levers it never had ------------
cfg.pool = { ...
    'k_cr','tau_cr','F_cr', ...                                  % the mechanism
    'eta_M','kSE_M','mu_neg','mu','kstiff1_n','kstiff2_n','x_M_slack', ... % series visco
    'kA2re','kA2hop','sA2hop', ...                               % restretch transient
    'ksr0','ksr2srd','kmsrd','ksrd2sr','sigma2','sigma_srd1', ... % SRX
    'SRXFromR2HighStrain','sSRXrip', ...                         % the routing (still untested)
    'A0','v_ref_reg','tau_reg', ...                              % <-- FV SHOULDER, was missing
    'P_bound_max', ...                                           % <-- attachment cap, was missing
    'PieceWiseStrainDep2X_logOffset','PieceWiseStrainDepX_logOffset', ... % <-- R2/R12 shape -> FV
    'PieceWiseStrainDep2Params__3','PieceWiseStrainDep2Params__4','PieceWiseStrainDep2Params__5', ...
    'k2','kstiff2','ka','kah','kSE','estiff','k1','k_1','dr2'};
% every round must be able to trade the mechanism, the series element, the master
% force scale AND force-velocity against each other
cfg.compulsory = {'k_cr','eta_M','kstiff2','v_ref_reg'};

keep = true(1, numel(cfg.pool));
for i = 1:numel(cfg.pool)
    f = cfg.pool{i}; us = strfind(f,'__');
    if isempty(us)
        keep(i) = isfield(params0, f) && isscalar(params0.(f)) && params0.(f) ~= 0;
    else
        b = f(1:us(1)-1); ix = str2double(f(us(1)+2:end));
        keep(i) = isfield(params0, b) && numel(params0.(b)) >= ix && params0.(b)(ix) ~= 0;
    end
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

%% ---- reference points under THIS objective ------------------------------
pchk = cfg.params0; pchk.mods = {}; pchk.g = ones(1,0);
[cS, ~, ctS, fnamesS] = evaluateBakersExp(1, pchk, true, false);
iF = find(contains(fnamesS, 'FV_fnorm'), 1);
fprintf('[%s] SEED (realign2 + MaxwellTensionOnly) = %.4f   FV term %.4f\n', ...
        cfg.tag, cS, ctS(iF));

pb = getParams(loadParams('params/rskR2_w025_opt.m'), [], true, false);
pb.fn = fn; pb.mods = {}; pb.g = ones(1,0);
pb.RunForceVelocity = 1; pb.RunSlack = 1; pb.RunSlackPassive = 1;
pb.RunKtr = 0; pb.RunStairs = 0; pb.RunForceVelocityTime = 0;
pb.EvalFeatures = 1; pb.PlotEachSeparately = 0; pb.PlotFeatureFitting = 0;
pb.RunSlackSegments = 'AllPar'; pb.MaxRunTime = cfg.MaxRunTime; pb.OptimizeOn = 'Feats';
[cB, ~, ctB] = evaluateBakersExp(1, pb, true, false);
fprintf('[%s] BASELINE re-scored under this objective = %.4f   FV term %.4f\n', ...
        cfg.tag, cB, ctB(iF));
fprintf('[%s] (the old 6.1025 / 5.0601 are NOT comparable — FV weight changed)\n', cfg.tag);

if ~isfinite(cS) || cS > 1e4
    error('RunOptimRealign3:badSeed', 'seed cost %.4g looks like a penalty', cS);
end
fprintf('[%s] pool=%d | N_DRAW=%d | budget %.2f h | resume=%d\n', ...
    cfg.tag, numel(cfg.pool), cfg.N_DRAW, cfg.TIME_BUDGET_HRS, cfg.RESUME);

optimizeFeatures(cfg);
