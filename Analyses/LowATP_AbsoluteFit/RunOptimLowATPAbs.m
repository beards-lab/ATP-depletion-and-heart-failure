% RunOptimLowATPAbs.m — fit the LOW-ATP data on an ABSOLUTE target, jointly with high ATP.
%
% WHY THIS EXISTS
% ---------------
% Every low-ATP fit so far scored the 2 mM/8 mM RATIO, which cancels the high-ATP
% baseline's own error and therefore cannot make the model reproduce the 2 mM data
% itself. This driver scores ABSOLUTE features at BOTH conditions:
%
%     Et = w8 * E(8 mM, absolute)  +  w2 * E(2 mM, absolute)
%
% Joint rather than 2 mM alone, deliberately: the shared parameters are the only way
% to reach most of the 2 mM residual, and moving them against a 2 mM-only cost would
% silently destroy the 8 mM fit. The joint cost makes that trade explicit and priced.
%
% Implemented by params0.UseAtpJointFit (Model/evaluateBakersExp.m) -- a guarded
% recursive call, inert when off, so optimizeFeatures/featCost/iter-gallery are reused.
%
% WHERE THE 2 mM RESIDUAL LIVES (baseline rskR2_w025_opt + K_T2=1.97, K_srx=13):
%     total 2 mM absolute 24.96, of which
%       rsK 9.72 | ovrsht_dy 5.81 | coolDownLS 3.80 | vall2_dy 1.63 | vall_t 1.06
%       ...and the amplitude/kinetic features are ALREADY SOLVED by the two ATP gates:
%       steady 0.11 | peak1_y 0.15 | A 0.29 | ktr 0.43
%   => ~88% of the low-ATP residual is the POST-RESTRETCH RECOVERY cluster.
%   The data's rsK ATP effect is x0.36 (2 mM 14.8-18.9 vs 8 mM 37.8-50.0) -- the
%   largest ATP effect in the dataset, larger than ktr (x0.55) or force (x1.36).
%
% LEVER STRUCTURE (from ../RestretchRecoveryFit/conclusions.md §6)
%   k2    sets the rsK LEVEL           -- and K_T2 (UseAtpDetachR2) is an [ATP] gate ON k2,
%                                         so it is a low-ATP-only rsK lever, free at 8 mM.
%   eta_M removes the AMPLITUDE dependence -- shared, so it is priced by the 8 mM leg.
%   K_srx gates kmsrd (ktr), also low-ATP-only.
% Both gates are MgATP_ref-normalised => exactly inert at 8 mM for any Km, so the
% optimiser can spend them on the 2 mM leg at zero cost to the 8 mM leg.
%
% USAGE
%   matlab -batch "run('Analyses/LowATP_AbsoluteFit/RunOptimLowATPAbs.m')"
% Override before running: W2, TAGSUF, TIME_HRS, NWORKERS, ALLOW_CONTESTED.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('W2',       'var'); W2       = 1;      end   % weight on the 2 mM leg
if ~exist('W8',       'var'); W8       = 1;      end   % weight on the 8 mM leg
if ~exist('TAGSUF',   'var'); TAGSUF   = 'w1';   end
if ~exist('TIME_HRS', 'var'); TIME_HRS = 18;     end
if ~exist('NWORKERS', 'var'); NWORKERS = 6;      end
% The R2 ceiling / D0 apo pool are FALSIFIED for high ATP (../ApoPoolDetachment/
% conclusions.md: the ktr overstretch peak pins R2(s) and every capped variant is
% 1.4-3x too high). Off by default; opt in only to re-open that question knowingly.
if ~exist('ALLOW_CONTESTED', 'var'); ALLOW_CONTESTED = false; end

cfg = struct();
cfg.baseSnap = 'params/rskR2_w025_opt.m';
cfg.tag      = ['lowatpabs_' TAGSUF];

params0 = getParams(); run(cfg.baseSnap);
params0 = getParams(params0, [], true, false);
NoS = params0.NumberOfStates;

%% ---- objective ----------------------------------------------------------
fn = params0.fn;
% The snapshot carries a duplicated 'rsR2|0.4-1.0|0.5' (RunOptimRsK appends it every
% run and only de-duplicates rsK), which double-counts that guardrail. Drop dupes.
[~, iu] = unique(fn, 'stable'); fn = fn(iu);

% Low-ATP leg: the 2 mM day has SLACK ONLY -- no force-velocity, no passive protocol.
% Their features would otherwise register as missing (100 each) on every evaluation.
fn2 = fn(~(startsWith(fn, 'FV') | startsWith(fn, 'PS_')));

params0.fn   = fn;
params0.mods = {};
params0.g    = [];

%% ---- joint dual-condition absolute objective ----------------------------
params0.UseAtpJointFit  = true;
params0.AtpJoint_w8     = W8;
params0.AtpJoint_w2     = W2;
params0.AtpJoint_fn2    = fn2;
params0.AtpJoint_vt2    = 'protocol_03_27_2026_2mM_slack.mat';
params0.AtpJoint_MgATP2 = 2;
params0.AtpJoint_Pi2    = 0.8;
params0.AtpJoint_RunFV2 = false;
params0.AtpJoint_RunPS2 = false;

% ATP gates ON, seeded at the current best (ratio-fit optimum on the TNF baseline).
params0.UseAtpDetachR2 = true;  params0.K_T2  = 1.97;
params0.UseAtpKmsrd    = true;  params0.K_srx = 13;
params0.MgATP_ref      = 8;
% high-ATP leg config (unchanged from the rsK campaign)
params0.MgATP = 8; params0.MgADP = 0; params0.Pi = 0;
params0.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';

fprintf('[%s] joint objective: %d hi-terms + %d lo-terms, w8=%g w2=%g\n', ...
    cfg.tag, numel(fn), numel(fn2), W8, W2);

%% ---- candidate pool -----------------------------------------------------
% Shared params (priced by BOTH legs) -- same comprehensive pool as the rsK campaign.
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

% LOW-ATP-ONLY levers: both gates are normalised to 1 at MgATP_ref, so they move the
% 2 mM leg and leave the 8 mM leg bit-for-bit identical. These are what make an
% absolute low-ATP fit reachable without paying at high ATP.
cand = [cand, {'K_T2','K_srx'}];

if ALLOW_CONTESTED
    warning('RunOptimLowATPAbs:contested', ...
        ['Adding the R2 ceiling / D0 apo pool to the pool. These are FALSIFIED at ' ...
         'high ATP by the ktr overstretch peak (Analyses/ApoPoolDetachment). ' ...
         'The 8 mM leg should punish them -- that is the point of running joint.']);
    params0.UseR2Ceiling = true; params0.UseD0State = true;
    cand = [cand, {'k2max','K_T2c','k_D0T','K_D0T','k2rip'}];
end

if NoS == 3; cand = [cand, {'k3','k3m','kstiff3','drp3','alpha3','s3'}]; end

% a multiplicative modifier can never move a zero base
keep = true(1, numel(cand));
for i = 1:numel(cand); keep(i) = (seedval(params0, cand{i}) ~= 0); end
dropped = cand(~keep);
cfg.pool = cand(keep);
if ~isempty(dropped)
    fprintf('[%s] dropped %d zero-base params: %s\n', cfg.tag, numel(dropped), strjoin(dropped, ', '));
end

% Compulsory: the standard three (k2 master rate, kstiff2 master force, ka
% peak-neutral force compensation) + eta_M (the only rsK amplitude-slope handle)
% + BOTH ATP gates, so every round can trade on the low-ATP leg. Without the gates
% in every round the optimiser can only improve 2 mM by damaging 8 mM.
cfg.compulsory = {'k2','kstiff2','ka','eta_M','K_T2','K_srx'};
if NoS == 3; cfg.compulsory = [cfg.compulsory, {'k3','kstiff3'}]; end

cfg.fn              = fn;
cfg.N_DRAW          = 10;
cfg.SIMPLEX_EVALS   = 60;
if ~exist('MAXRUNTIME','var'); MAXRUNTIME = 600; end
cfg.MaxRunTime      = MAXRUNTIME;
cfg.SURR_EVALS      = 0;
cfg.KICK_FRAC       = 0.1;
cfg.DEBUG           = false;
cfg.TIME_BUDGET_HRS = TIME_HRS;
cfg.RESUME          = isfile(fullfile(root, 'params', [cfg.tag '_state.mat']));

params0.RunForceVelocity=true; params0.RunKtr=false; params0.RunSlack=true;
params0.RunStairs=false; params0.RunForceVelocityTime=false; params0.RunSlackPassive=true;
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false;
params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0;
params0.RunSlackSegments='AllPar';
params0.MaxRunTime=cfg.MaxRunTime;
params0.OptimizeOn='Feats';

cfg.params0 = params0;

if isempty(gcp('nocreate')); parpool('local', NWORKERS); end
fprintf('[%s] pool=%d | N_DRAW=%d | budget %.1f h | resume=%d\n', ...
    cfg.tag, numel(cfg.pool), cfg.N_DRAW, cfg.TIME_BUDGET_HRS, cfg.RESUME);

%% ---- BASELINE SANITY GUARD ---------------------------------------------
% A joint evaluation is ~2x a single one; if the seed is already a failure penalty
% the run would spend the night optimising the penalty. Refuse to start.
fprintf('[%s] checking baseline evaluates (joint)...\n', cfg.tag);
pchk = cfg.params0; pchk.mods = {}; pchk.g = ones(1,0);
[c0, ~, fc0, fn0] = evaluateBakersExp(1, pchk, true, false);
fprintf('[%s] joint baseline cost = %.4f\n', cfg.tag, c0);
if ~isfinite(c0) || c0 > 1e4
    error('RunOptimLowATPAbs:badBaseline', ...
        ['Joint baseline %.4g looks like an ODE-failure penalty or a broken seed. ' ...
         'Raise cfg.MaxRunTime (now %d s) or check the snapshot.'], c0, cfg.MaxRunTime);
end
[~, o] = sort(fc0, 'descend');
fprintf('[%s] top baseline terms:\n', cfg.tag);
for i = 1:min(10, numel(o)); fprintf('    %-26s %8.3f\n', fn0{o(i)}, fc0(o(i))); end

% Validate the whole setup (pool, objective, seed, joint wiring) without spending
% the time budget:  DRYRUN_SETUP = true; run(...)
if exist('DRYRUN_SETUP','var') && DRYRUN_SETUP
    fprintf('[%s] DRY RUN -- setup validated, not starting optimizeFeatures.\n', cfg.tag);
    return;
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
