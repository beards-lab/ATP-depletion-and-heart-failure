% RunOptimRsK.m — refit with the post-restretch redevelopment RATE in the
% objective.
%
% WHY THIS EXISTS
% ---------------
% At params/optfull7_opt_mov.m the model recovers force after the re-stretch
% 2.70x too fast (rsK 118 s^-1 vs data 43.7), while the recovery AMPLITUDE
% (rsA 13.6 vs 12.3) and the post-slack rate (ktr x1.05) are already right.
% The old objective could not see this: the only term scoring that window by
% trace shape is coolDownLS, which contributes 0.32 of a 28.18 feature cost
% (1.1%) because the 280 ms window is dominated by the plateau the model
% already gets right. Replacing that least-squares view with a first-order
% RATE makes the defect a first-class term worth ~8.5.
%
% USAGE
%   matlab -batch "run('Analyses/RestretchRecoveryFit/RunOptimRsK.m')"
% Set WEIGHT_RSK / TAGSUF from the caller to run a second instance at a
% different weight (tags are the isolation key -- two instances with
% different tags never share a file).

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

% WEIGHT_RSK: the optimiser minimises the L2 feature cost (evaluateBakersExp
% calls evalFeatureCost without costExp, so it defaults to 2), and rsK's
% normalised errors are 1.3-2.1 -- large enough that squaring makes them
% dominate. Measured share of the objective at this seed:
%     weight 1.00 -> 78.6%     weight 0.35 -> 56.2%     weight 0.20 -> 42.3%
%     weight 0.50 -> 64.7%     weight 0.25 -> 47.8%     weight 0.15 -> 35.5%
% 0.25 makes rsK by far the largest single term while leaving over half the
% weight to the rest of the fit, so a round cannot be accepted purely for
% chasing rsK into a corner.
if ~exist('WEIGHT_RSK', 'var'); WEIGHT_RSK = 0.25; end
if ~exist('TAGSUF',     'var'); TAGSUF     = 'w025'; end
if ~exist('TIME_HRS',   'var'); TIME_HRS   = 18; end
if ~exist('NWORKERS',   'var'); NWORKERS   = 6; end

cfg = struct();
cfg.baseSnap = 'params/optfull7_opt_mov.m';
cfg.tag      = ['rsk_' TAGSUF];

params0 = getParams(); run(cfg.baseSnap);
NoS = params0.NumberOfStates;

%% ---- objective: the previous feature list + the recovery RATE -----------
fn = params0.fn;

% Drop any stale post-restretch entry so re-seeding from a snapshot this
% driver itself wrote cannot double-count it. Must match on the PRIMARY FIELD
% NAME, not on a prefix of the whole string: startsWith(fn,'rsK') misses
% 'rsR2|0.4-1.0|0.5', which is then appended a second time below.
RS_TERMS = {'rsK', 'rsR2'};
primary  = cellfun(@(s) regexprep(strtok(s, '|'), '\[.*$', ''), fn, ...
                   'UniformOutput', false);
fn = fn(~ismember(primary, RS_TERMS));

% rsK: normalised |data-model| per cycle, summed over the 5 cycles.
% At the seed this scores ~8.5 at weight 1 -- deliberately the largest single
% term, because the run's question is "CAN the model do this, and what does
% it cost?". A small weight would answer neither.
fn{end+1} = sprintf('rsK|%g', WEIGHT_RSK);

% Guardrail: keep the fitted rise recognisably first-order. Without this the
% optimiser can win on rsK by making the window non-exponential (e.g. a slow
% creep with a fast foot) rather than by genuinely slowing the recovery.
fn{end+1} = 'rsR2|0.4-1.0|0.5';

params0.fn   = fn;
params0.mods = {};
params0.g    = [];

fprintf('[%s] objective (%d terms), rsK weight %.2f:\n', cfg.tag, numel(fn), WEIGHT_RSK);
for i = 1:numel(fn); fprintf('    %s\n', fn{i}); end

%% ---- candidate pool ------------------------------------------------------
% Same comprehensive pool as Analyses/RestretchFeatureFit/RunOptimFull.m.
% The levers most likely to matter here are the strain-dependent detachment
% knots (R2 reaches 2100-4250 /s at re-stretch strains, so every stretched
% bridge is gone in <0.5 ms and the recovery is pure fast re-attachment) and
% the attachment side (ka / A0 / v_ref_reg / tau_reg).
cand = { ...
    'ka','kd','k1','k_1','k2','k_2','kah','kamh', ...
    'kstiff1','kstiff2','kstiff1_n','kstiff2_n','kSE','estiff', ...
    'kSE_M','eta_M','mu_neg','mu', ...
    'dr','dr1','dr2', ...
    'ksr0','kmsrd','ksr2srd','ksrd2sr','sigma1','sigma2','sigma_srd1','sigma_srd2', ...
    'A0','v_ref_reg','tau_reg', ...
    'P_bound_max', ...
    'kA2hop','sA2hop', ...
    'L_thick','L_hbare','L_thin','LXBpivot', ...
    'd10_ref','SL_ref_lattice','R_thick','R_thin','d_optimal','sigma_lattice', ...
    'PieceWiseStrainDep2X_logOffset', 'PieceWiseStrainDepR1DX_logOffset', 'PieceWiseStrainDepR21X_logOffset', ...
    'PieceWiseStrainDepX__2','PieceWiseStrainDepX__3','PieceWiseStrainDepX__4', ...
    'PieceWiseStrainDep2X__2','PieceWiseStrainDep2X__3','PieceWiseStrainDep2X__4','PieceWiseStrainDep2X__5', ...
    'PieceWiseStrainDepR1DX__1','PieceWiseStrainDepR1DX__2','PieceWiseStrainDepR1DX__3','PieceWiseStrainDepR1DX__4', ...
    'PieceWiseStrainDepR21X__2','PieceWiseStrainDepR21X__3','PieceWiseStrainDepR21X__4',...
    'PieceWiseStrainDepParams__2','PieceWiseStrainDepParams__3','PieceWiseStrainDepParams__4', ...
    'PieceWiseStrainDep2Params__2','PieceWiseStrainDep2Params__3','PieceWiseStrainDep2Params__4','PieceWiseStrainDep2Params__5', ...
    'PieceWiseStrainDep2Params__6','PieceWiseStrainDep2Params__7', ...
    'PieceWiseStrainDepR1DParams__1','PieceWiseStrainDepR1DParams__2','PieceWiseStrainDepR1DParams__3','PieceWiseStrainDepR1DParams__4', ...
    'PieceWiseStrainDepR21Params__2','PieceWiseStrainDepR21Params__3','PieceWiseStrainDepR21Params__4'};

if NoS == 3
    cand = [cand, {'k3','k3m','kstiff3','drp3','alpha3','s3'}];
end

% a multiplicative modifier g can never move a zero base
keep = true(1, numel(cand));
for i = 1:numel(cand); keep(i) = (seedval(params0, cand{i}) ~= 0); end
dropped = cand(~keep);
cfg.pool = cand(keep);
if ~isempty(dropped)
    fprintf('[%s] dropped %d zero-base params: %s\n', cfg.tag, numel(dropped), strjoin(dropped, ', '));
end

% Compulsory = the levers every round has to be able to trade against.
% k2 (master rate), kstiff2 (master force scale), ka (peak-neutral force
% compensation) are the standard three. eta_M is added here because the lever
% scan found the Maxwell dashpot to be the ONLY handle on rsK that leaves ktr
% and isometric force untouched (eta_M x0.6 -> d ln rsK = -0.337 at dFeat
% +5.2, ktr and steady unchanged) — and rsK is now ~48% of the objective, so
% every round needs access to it or it cannot trade.
cfg.compulsory = {'k2','kstiff2','ka','eta_M'};
if NoS == 3; cfg.compulsory = [cfg.compulsory, {'k3','kstiff3'}]; end

cfg.fn            = fn;
cfg.N_DRAW        = 10;
cfg.SIMPLEX_EVALS = 60;
% MaxRunTime guards against wandering into pathologically slow parameter
% regions -- it is NOT a budget for a healthy evaluation. At this snapshot the
% slack init run alone is ~22 s (TimeBattery.m), and the first launch of this
% driver died with a 1016441 baseline (the 1e6 ODE-timeout penalty) because
% MaxRunTime=60 left no margin once the lever scan was competing for cores.
if ~exist('MAXRUNTIME', 'var'); MAXRUNTIME = 600; end
cfg.MaxRunTime    = MAXRUNTIME;
cfg.SURR_EVALS    = 0;
cfg.KICK_FRAC     = 0.1;
cfg.DEBUG         = false;
cfg.TIME_BUDGET_HRS = TIME_HRS;
cfg.RESUME        = isfile(fullfile(root, 'params', [cfg.tag '_state.mat']));

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

fprintf('[%s] pool=%d | N_DRAW=%d | budget %.1f h | resume=%d\n', ...
    cfg.tag, numel(cfg.pool), cfg.N_DRAW, cfg.TIME_BUDGET_HRS, cfg.RESUME);

%% ---- BASELINE SANITY GUARD ----------------------------------------------
% An optimiser seeded from a broken evaluation optimises the penalty, not the
% model, and will happily run all night doing it. Two distinct failure modes,
% so two distinct checks.

% (1) COVERAGE: can every fn entry actually be scored? A feature absent from
% features_data or features_model scores MISSING_FEATURE_COST = 100, which
% after weighting looks like an ordinary bad-fit number rather than an error
% ('rsK|0.25' with no data target scores exactly 25). This bites whenever the
% untracked data/*.mat files are older than Model/extractSlackAttributes.m --
% e.g. after moving the run to another machine.
pchk = cfg.params0; pchk.mods = {}; pchk.g = ones(1, 0);   % = optimizeFeatures' setMods
fprintf('[%s] checking every fn entry can be scored...\n', cfg.tag);
rep = checkFeatureCoverage(pchk);
if ~rep.ok
    if ~isempty(rep.missingModel)
        fprintf(2, '  ABSENT from features_model (stale Model/ code?):\n');
        fprintf(2, '    %s\n', rep.missingModel{:});
    end
    if ~isempty(rep.missingData)
        fprintf(2, '  ABSENT from features_data (stale data/*.mat?):\n');
        fprintf(2, '    %s\n', rep.missingData{:});
    end
    error('RunOptimRsK:missingFeatures', ...
        ['%d feature(s) in params0.fn cannot be scored; each contributes a ' ...
         'flat MISSING_FEATURE_COST=100 x weight and is INVISIBLE to the ' ...
         'optimiser. For rsK/rsR2: copy data/%s from the machine that has ' ...
         'them, or run Analyses/RestretchRecoveryFit/AddRestretchFeatsToData.m ' ...
         '(DRYRUN=false) to add them in place.'], ...
        numel(rep.missingModel) + numel(rep.missingData), params0.velocitytableonfile);
end
fprintf('[%s]   coverage OK (%d fn entries)\n', cfg.tag, numel(pchk.fn));

% (2) EVALUABILITY: does the seed integrate inside MaxRunTime?
fprintf('[%s] checking baseline evaluates...\n', cfg.tag);
c0 = evaluateBakersExp(1, pchk, true, true);
fprintf('[%s] baseline cost = %.4f\n', cfg.tag, c0);
if ~isfinite(c0) || c0 > 1e4
    error('RunOptimRsK:badBaseline', ...
        ['Baseline cost %.4g looks like an ODE-failure penalty (>=1e6) or a ' ...
         'broken seed. Raise cfg.MaxRunTime (now %d s) or check the snapshot. ' ...
         'Refusing to start a long run from a bad baseline.'], c0, cfg.MaxRunTime);
end

optimizeFeatures(cfg);

function v = seedval(p, name)
    us = strfind(name, '__');
    if isempty(us)
        if isfield(p, name); v = p.(name); else; v = 0; end
    else
        f = name(1:us(1)-1); idx = str2double(name(us(1)+2:end));
        if isfield(p, f) && numel(p.(f)) >= idx; v = p.(f)(idx); else; v = 0; end
    end
end
