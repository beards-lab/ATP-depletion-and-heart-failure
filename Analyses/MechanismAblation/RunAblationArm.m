% RunAblationArm.m -- stage 1: one budget-matched minor re-optimisation.
%
% Every arm and every control replicate gets EXACTLY the same search budget:
% one round, N_DRAW = 10, SIMPLEX_EVALS = 60, KICK 0.10, no surrogate. ~38 min.
%
% WHY BUDGET-MATCHED, AND WHY CONTROLS. The incumbent sits at a local minimum
% and every flip starts off one. Ranking a re-optimised flip against the
% un-re-optimised 0.9081 charges the flip for the restart and makes every
% mechanism look load-bearing. The control arms (no flip, same budget, different
% RNG) measure what one round buys for nothing; sigma over those replicates is
% the only honest significance threshold. This is the StairsPrediction lesson:
% the incumbent was where the search STOPPED, not a minimum.
%
% WHAT IS "RELEVANT PARAMETERS" HERE. cfg.compulsory = the arm's own mechanism
% parameters (`add`) plus its named compensators (`comp`), so the mechanism is
% allowed to be absorbed rather than merely removed. Everything else is drawn at
% random from the pool, 10 - numel(compulsory) per round -- see plan.md trap T3.
%
% Usage (one arm):
%   matlab -batch "ARM='C1'; run('Analyses/MechanismAblation/RunAblationArm.m')"
% Control replicate:
%   matlab -batch "ARM='control'; REP=1; run('Analyses/MechanismAblation/RunAblationArm.m')"
%
% Writes: params/abl_<name>_opt.m, params/abl_<name>_state.mat, params/abl_<name>_iter/
%
% See also: ablationSpec, RunAblationPreflight, RunAblationReport, plan.md

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('ARM','var'); error('RunAblationArm:noArm','set ARM (an id from ablationSpec, or ''control'')'); end
if ~exist('REP','var');       REP       = 1;              end
if ~exist('SEEDSNAP','var');  SEEDSNAP  = 'floored_seed'; end
if ~exist('NWORKERS','var');  NWORKERS  = 5;              end

% The budget. Identical for every arm; change it here or nowhere.
%
% SIZED 2026-09-11 AFTER THE FIRST ATTEMPT RAN 8 h AND FINISHED 9 OF 36 ARMS.
% What went wrong: MaxRunTime is a watchdog on ONE SIMULATION, not on one
% objective evaluation (evaluateModel.m:102 puts it in the ODE event vector).
% One evaluation runs ~18 simulations -- 8 force-velocity points, 5 slack
% segments, 5 passive -- so MAX_RUN_TIME = 300 permitted a single evaluation to
% take 90 minutes, and a 60-evaluation round to take days. A healthy evaluation
% is 29-104 s (measured, median ~35), i.e. ~2-6 s per simulation.
%
% MAX_RUN_TIME = 60 was then tried and was WORSE: it killed every arm, the
% unmodified control included, with "ODEslower took longer than MaxRunTime".
% Two things that watchdog is, which the average per-simulation time hid:
%   (1) it is WALL-CLOCK, so it trips on machine contention, not just on a
%       stiff parameter set. 60 s was survivable with two instances and fatal
%       with three competing for 12 cores and 31 GB.
%   (2) it THROWS (evaluateModel.m:253). optimizeFeatures catches it inside the
%       search and turns it into a high finite cost, but an uncaught one
%       anywhere else takes the whole arm down -- which is what happened.
% It also aborts the evaluation at the FIRST simulation to trip, so an
% evaluation costs at most (healthy prefix) + MAX_RUN_TIME, not 18x the cap.
%
% So: a cap set above the legitimate maximum (a full healthy evaluation is
% 29-104 s, so no single simulation exceeds ~104 s), plus two rounds so the
% wall-clock budget between them can actually fire, plus try/catch on every
% evaluation this script makes itself.
N_DRAW        = 10;
SIMPLEX_EVALS = 20;
SURR_EVALS    = 0;
KICK_FRAC     = 0.10;
MAX_RUN_TIME  = 150;    % above the largest healthy evaluation, with contention headroom
N_ROUNDS      = 2;
TIME_CAP_HRS  = 0.75;   % now reachable, because there is more than one round

%% ---- base ---------------------------------------------------------------
base = getParams();
snap = loadParams(SEEDSNAP);
f = fieldnames(snap); for i = 1:numel(f); base.(f{i}) = snap.(f{i}); end

base.RunForceVelocity = true;  base.RunKtr = false;   base.RunSlack = true;
base.RunStairs = false;        base.RunForceVelocityTime = false;
base.RunSlackPassive = true;   base.EvalFeatures = true;
base.BreakOnODEUnstable = false;
base.PlotEachSeparately = 0;   base.PlotFeatureFitting = 0;
base.RunSlackSegments = 'AllPar';
base.MaxRunTime = MAX_RUN_TIME; base.OptimizeOn = 'Feats';
base.mods = {}; base.g = [];

pool = iBasePool();

%% ---- the arm ------------------------------------------------------------
if strcmpi(ARM, 'control')
    name = sprintf('ctrl%d', REP);
    % `scale` matters here too: the control must spend the SAME 10-evaluation
    % line search as every arm or it is not budget-matched, and an
    % under-budgeted control makes every arm look better than it is.
    a = struct('id',sprintf('CTRL%d',REP), 'name',name, 'set',{{}}, ...
               'add',{{}}, 'drop',{{}}, ...
               'comp',{{'k2','kstiff2','ka','v_ref_reg','k_cr','eta_M'}}, ...
               'note','control: no flip, same budget, different RNG', ...
               'scale','k2');
    rng(1000 + REP);
else
    A = ablationSpec();
    k = find(strcmp({A.id}, ARM) | strcmp({A.name}, ARM), 1);
    if isempty(k); error('RunAblationArm:unknownArm','no arm "%s" in ablationSpec', ARM); end
    a = A(k); name = a.name;
    rng(12345);   % same draw seed for every flip, so arms differ by mechanism
end

params0 = base;
for j = 1:2:numel(a.set); params0.(a.set{j}) = a.set{j+1}; end

%% ---- pool surgery -------------------------------------------------------
pool = pool(~ismember(pool, a.drop));
pool = [pool, a.add(~ismember(a.add, pool))];

% The engine's own filter: a zero-valued entry is silently dropped, which is
% how an ON-side arm ends up on but undrawable (plan.md T1). Fail loudly here
% instead, because a silently-dropped mechanism parameter produces a confident
% null result.
keep = true(1, numel(pool));
for i = 1:numel(pool)
    v = iResolve(params0, pool{i});
    keep(i) = ~isnan(v) && v ~= 0;
end
droppedZero = pool(~keep);
pool = pool(keep);
bad = a.add(ismember(a.add, droppedZero));
if ~isempty(bad)
    error('RunAblationArm:inertMechanism', ...
        ['arm %s: %s is zero after the flip, so the mechanism is ON but inert ' ...
         'and undrawable. Give it a seed value in ablationSpec.'], a.id, strjoin(bad,', '));
end

compulsory = unique([a.add, a.comp], 'stable');
compulsory = compulsory(ismember(compulsory, pool));
if numel(compulsory) > N_DRAW - 2
    % leave at least two slots for random exploration
    compulsory = compulsory(1:N_DRAW-2);
end

TAG = ['abl_' name];
fprintf(['\n=== ablation arm %s (%s) ===\n  %s\n' ...
         '  set        : %s\n  add        : %s\n  drop       : %s\n' ...
         '  compulsory : %s\n  pool %d | N_DRAW %d | %d round x %d evals\n\n'], ...
    a.id, TAG, a.note, iKV(a.set), strjoin(a.add,', '), strjoin(a.drop,', '), ...
    strjoin(compulsory,', '), numel(pool), N_DRAW, N_ROUNDS, SIMPLEX_EVALS);
if ~isempty(droppedZero)
    fprintf('  (zero-valued, not drawn: %s)\n\n', strjoin(droppedZero,', '));
end

%% ---- pool ---------------------------------------------------------------
if isempty(gcp('nocreate'))
    try
        cl  = parcluster('local');
        jsl = fullfile(tempdir, sprintf('mjob_%s_%d', TAG, feature('getpid')));
        if ~isfolder(jsl); mkdir(jsl); end
        cl.JobStorageLocation = jsl;
        if cl.NumWorkers < NWORKERS; cl.NumWorkers = NWORKERS; end
        parpool(cl, NWORKERS);
    catch ME
        warning('RunAblationArm:pool','per-session pool failed (%s) -- falling back', ME.message);
        if isempty(gcp('nocreate')); parpool('local', NWORKERS); end
    end
end

%% ---- pre-flight ---------------------------------------------------------
pchk = params0; pchk.mods = {}; pchk.g = ones(1,0);
rep = checkFeatureCoverage(pchk);
if ~rep.ok
    fprintf(2, 'UNSCORABLE fn entries (each costs MISSING_FEATURE_COST=100):\n');
    if ~isempty(rep.missingModel); fprintf(2,'  model: %s\n', strjoin(rep.missingModel,', ')); end
    if ~isempty(rep.missingData);  fprintf(2,'  data : %s\n', strjoin(rep.missingData ,', ')); end
    error('RunAblationArm:coverage','objective has unscorable terms');
end
[c0, ct0] = iEval(pchk);
fprintf('[%s] raw flip = %.4f (frozen parameters)\n', TAG, c0);
% No hard guard here. The raw flip is EXPECTED to be terrible for a gate-like
% mechanism -- B1 lands at 46007, C9 at a NaN sentinel -- and the calibration
% below is what it is for. The guard is applied after calibration instead.

if exist('DRYRUN','var') && DRYRUN
    fprintf('[%s] DRYRUN -- not starting the search.\n', TAG); return;
end

%% ---- scalar pre-compensation -------------------------------------------
% One parameter, one line search, before the simplex. WHY: a flip that changes
% a gain lands the model where a 60-eval simplex cannot climb out of -- B1
% (drop the availability gate, A0 = 0.21) starts at 46007, of which 45000 is
% coolDownLS, an unbounded raw MSE. Charging an arm for a pure scale offset and
% calling the result "load-bearing" is the failure mode this whole design is
% meant to avoid. For ON-side arms the scaled parameter is the mechanism's own
% strength, so this also replaces the seed value guessed in ablationSpec with
% one the data picked.
%
% It is part of the budget and identical for every arm, controls included.
cal = struct('param', a.scale, 'grid', [], 'cost', [], 'best', 1, 'c', c0);
if ~isempty(a.scale) && ~isnan(iResolve(params0, a.scale))
    cal.grid = [0.1 0.2 0.35 0.5 0.7 1 1.4 2 3 5];
    v0 = iResolve(params0, a.scale);
    cal.cost = nan(size(cal.grid));
    fprintf('[%s] calibrating %s (base %.4g) over %d multipliers\n', TAG, a.scale, v0, numel(cal.grid));
    for i = 1:numel(cal.grid)
        if cal.grid(i) == 1; cal.cost(i) = c0; continue; end
        pc = iAssign(params0, a.scale, v0 * cal.grid(i));
        cal.cost(i) = iEval(iBare(pc));
        fprintf('        x%-5.2g -> %10.4f\n', cal.grid(i), cal.cost(i));
    end
    [cbest, ib] = min(cal.cost);
    if isfinite(cbest) && cbest < c0
        cal.best = cal.grid(ib); cal.c = cbest;
        params0 = iAssign(params0, a.scale, v0 * cal.best);
        fprintf('[%s] calibrated: %s x%.3g  -> %.4f (was %.4f)\n', TAG, a.scale, cal.best, cbest, c0);
    else
        fprintf('[%s] calibration found nothing better than x1 (%.4f)\n', TAG, c0);
    end
else
    fprintf('[%s] no scale parameter -- starting the search at the raw flip\n', TAG);
end
[ccal, ctcal] = iEval(iBare(params0));
if ~isfinite(ccal) || ccal > 1e4
    error('RunAblationArm:badStart', ...
        ['arm %s still sits at %.4g after calibrating %s -- no single scalar ' ...
         'rescues it, so a 60-eval simplex will not either. Re-seed the arm or ' ...
         'give it a companion revival.'], a.id, ccal, a.scale);
end

%% ---- search -------------------------------------------------------------
cfg = struct();
cfg.baseSnap = fullfile('params', [SEEDSNAP '.m']);   % provenance only; the
cfg.params0  = params0;                               % start point is params0
cfg.tag = TAG; cfg.fn = params0.fn;
cfg.pool = pool; cfg.compulsory = compulsory;
cfg.N_DRAW = N_DRAW; cfg.SIMPLEX_EVALS = SIMPLEX_EVALS; cfg.SURR_EVALS = SURR_EVALS;
cfg.KICK_FRAC = KICK_FRAC; cfg.MaxRunTime = MAX_RUN_TIME; cfg.N_ROUNDS = N_ROUNDS;
cfg.TIME_BUDGET_HRS = TIME_CAP_HRS; cfg.DEBUG = false;
cfg.RESUME = isfile(fullfile(root,'params',[TAG '_state.mat']));

state = optimizeFeatures(cfg);

%% ---- record -------------------------------------------------------------
OUTDIR = fullfile(here, 'results', 'arms');
if ~isfolder(OUTDIR); mkdir(OUTDIR); end
best = loadParams(fullfile('params', [TAG '_opt.m']));
pb = base; fb = fieldnames(best); for i = 1:numel(fb); pb.(fb{i}) = best.(fb{i}); end
[c1, ct1, nm1] = iEval(iBare(pb));
res = struct('id',a.id,'name',name,'tag',TAG,'note',a.note,'set',{a.set}, ...
             'add',{a.add},'drop',{a.drop},'compulsory',{compulsory}, ...
             'cost_raw',c0,'cost_start',ccal,'cost_end',c1,'cal',cal, ...
             'featCost_start',ctcal,'featCost_end',ct1,'featNames',{nm1}, ...
             'nPool',numel(pool),'dNParams',numel(a.add)-numel(a.drop));
save(fullfile(OUTDIR, [TAG '.mat']), 'res');
fprintf('\n[%s] raw %.4f -> calibrated %.4f -> reoptimised %.4f  (dParams %+d)\n', ...
        TAG, c0, ccal, c1, numel(a.add)-numel(a.drop));
fprintf('[%s] written: %s\n', TAG, fullfile(OUTDIR, [TAG '.mat']));

%% ---- helpers ------------------------------------------------------------
function s = iKV(kv)
p = cell(1, numel(kv)/2);
for i = 1:2:numel(kv)
    v = kv{i+1};
    if ischar(v); p{(i+1)/2} = sprintf('%s=%s', kv{i}, v);
    else;         p{(i+1)/2} = sprintf('%s=%g', kv{i}, v); end
end
s = strjoin(p, ', ');
if isempty(s); s = '(none)'; end
end

function p = iAssign(p, nm, v)
us = strfind(nm,'__');
if isempty(us); p.(nm) = v;
else; b = nm(1:us(1)-1); p.(b)(str2double(nm(us(1)+2:end))) = v; end
end

function p = iBare(p)
p.mods = {}; p.g = ones(1,0);
end

function [c, ct, nm] = iEval(p)
% One objective evaluation that cannot take an arm down with it. The
% MaxRunTime watchdog THROWS (evaluateModel.m:253) and a parfor turns that into
% "One or more futures resulted in an error"; optimizeFeatures already converts
% both to a high finite cost inside the search, and this does the same for the
% evaluations made outside it. A thrown watchdog means "this point is too slow
% to be worth costing", which is a bad point -- not a reason to lose the arm.
try
    [c, ~, ct, nm] = evaluateBakersExp(1, p, true, false);
catch ME
    fprintf(2, '        eval failed (%s) -> treated as a penalty\n', ME.message);
    c = 1e5; ct = []; nm = {};
end
end

function v = iResolve(p, nm)
us = strfind(nm,'__');
if isempty(us)
    if isfield(p,nm) && isscalar(p.(nm)) && isnumeric(p.(nm)); v = p.(nm); else; v = NaN; end
else
    b = nm(1:us(1)-1); ix = str2double(nm(us(1)+2:end));
    if isfield(p,b) && numel(p.(b)) >= ix; v = p.(b)(ix); else; v = NaN; end
end
end

function pool = iBasePool()
pool = { ...
    'ka','kd','k1','k_1','k2','k_2','kah','kamh', ...
    'kstiff1','kstiff2','kstiff1_n','kstiff2_n','kSE','estiff','ekSE', ...
    'kSE_M','eta_M','mu_neg','mu','x_M_slack', ...
    'dr','dr1','dr2', ...
    'ksr0','kmsrd','ksr2srd','ksrd2sr','sigma1','sigma2','sigma_srd1','sigma_srd2', ...
    'A0','v_ref_reg','tau_reg','P_bound_max', ...
    'k_cr','tau_cr','F_cr','SRXFromR2HighStrain','sSRXrip', ...
    'kA2re','kA2hop','sA2hop', ...
    'L_thick','L_hbare','L_thin','LXBpivot', ...
    'd10_ref','SL_ref_lattice','R_thick','R_thin','d_optimal','sigma_lattice', ...
    'Lsc0','gamma','k_pas', ...
    'PieceWiseStrainDep2X_logOffset','PieceWiseStrainDepR1DX_logOffset', ...
    'PieceWiseStrainDepR21X_logOffset','PieceWiseStrainDepX_logOffset', ...
    'PieceWiseStrainDepX__2','PieceWiseStrainDepX__3','PieceWiseStrainDepX__4', ...
    'PieceWiseStrainDep2X__2','PieceWiseStrainDep2X__3','PieceWiseStrainDep2X__4','PieceWiseStrainDep2X__5', ...
    'PieceWiseStrainDepR1DX__1','PieceWiseStrainDepR1DX__2','PieceWiseStrainDepR1DX__3','PieceWiseStrainDepR1DX__4', ...
    'PieceWiseStrainDepR21X__2','PieceWiseStrainDepR21X__3','PieceWiseStrainDepR21X__4', ...
    'PieceWiseStrainDepParams__2','PieceWiseStrainDepParams__3','PieceWiseStrainDepParams__4', ...
    'PieceWiseStrainDep2Params__2','PieceWiseStrainDep2Params__3','PieceWiseStrainDep2Params__4', ...
    'PieceWiseStrainDep2Params__5','PieceWiseStrainDep2Params__6','PieceWiseStrainDep2Params__7', ...
    'PieceWiseStrainDepR1DParams__1','PieceWiseStrainDepR1DParams__2', ...
    'PieceWiseStrainDepR1DParams__3','PieceWiseStrainDepR1DParams__4', ...
    'PieceWiseStrainDepR21Params__2','PieceWiseStrainDepR21Params__3','PieceWiseStrainDepR21Params__4'};
end
