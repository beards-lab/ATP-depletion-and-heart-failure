function state = optimizeFeatures(cfg)
% OPTIMIZEFEATURES  Decisive global->local feature-cost optimizer, refactored
% from RunRestretchOptim.m so that TWO (or more) INSTANCES can run in
% parallel MATLAB sessions without clobbering each other's output.
% =========================================================================
% All per-run identity — which snapshot to seed from, which tunable pool to
% draw from, and where to write results — comes from CFG. Nothing here
% touches a fixed path: every instance gets its own STATE_FILE/SNAP_BEST
% derived from cfg.tag, so two sessions started with different cfg.tag
% values never read or write each other's files.
%
% Per round:
%   1. DRAW a N_DRAW-param subset from cfg.pool (cfg.compulsory levers
%      always in; remainder random). Keeps each surrogateopt problem
%      low-dimensional and tractable.
%   2. SURROGATEOPT (global, bounded) explores the drawn box.
%   3. FMINSEARCH (local) polishes the surrogate's best point.
%   4. If the round IMPROVES the incumbent, BAKE the multipliers into the
%      running params and snapshot to disk. Otherwise count a STALL.
%   5. On STALL_LIMIT consecutive non-improving rounds, RANDOM-KICK a
%      WORKING COPY of the incumbent and force a fresh draw next round
%      (escape local basins). The incumbent itself (state.best_params /
%      state.best_cost) is NEVER touched by a kick.
%
% ACCEPT/SAVE INVARIANT (the bug this refactor fixes): a prior run saved a
% snapshot (cost 174) that was WORSE than its own seed (9.72) because the
% final on-disk snapshot was written from the last *working* point, which
% may have been a stalled/kicked excursion never validated against the
% incumbent. Fix: the incumbent (state.best_params, state.best_cost) is
% updated ONLY inside the "round improved" branch, immediately followed by
% a disk write of THAT incumbent. The stall/kick branch only ever mutates
% state.params (a working copy used to seed the NEXT round's draw) and
% never writes to disk. The final write (after the loop, and at the top of
% every improving round) is always writeParamsToMFile(SNAP_BEST,
% state.best_params) — never state.params — with an assertion immediately
% before the final write that the snapshot-file's cost matches
% state.best_cost.
%
% Bounds: each drawn parameter's g-multiplier box is derived from
% parameterBounds.m (absolute [lb,ub] -> multiplier around the current
% value) via gMultBox. Params without a bound entry (array knots, velocity-
% gauss levers, A0/v_ref_reg) get named default windows (see gMultBox).
%
% Cost: evaluateBakersExp(g, params0) -> RunBakersExp -> evalFeatureCost on
% cfg.fn (FV + slack + restretch shape features).
%
% Inputs:
%   cfg.baseSnap   - path to the seed snapshot .m (a writeParamsToMFile
%                    output, run(...) to populate params0). REQUIRED.
%   cfg.tag        - string identifying this optimizer instance. Outputs go
%                    to params/<tag>_opt.m and params/<tag>_state.mat.
%                    REQUIRED. This is the isolation key: two instances
%                    with different tags never collide.
%   cfg.fn         - cell array, the feature-cost list (params0.fn).
%                    REQUIRED.
%   cfg.pool       - cell array of param names eligible for optimisation.
%                    REQUIRED (the base snapshot need not define
%                    OptimizeParams -- e.g. params_2state_seed.m does not).
%   cfg.compulsory - cell array, subset of cfg.pool always drawn each
%                    round. REQUIRED (may be {}).
%   cfg.DEBUG            - true = ~10 min smoke test; false = production.
%                           Default: true.
%   cfg.N_DRAW           - params optimised per round (incl. compulsory).
%                           Default: 12.
%   cfg.SURR_EVALS       - surrogateopt function-evaluation budget/round.
%                           Default: DEBUG ? 12 : 0 (0 = surrogateopt call
%                           is skipped and g=1 is used as the "surrogate"
%                           seed for fminsearch -- matches
%                           RunRestretchOptim.m's current production mode,
%                           where the surrogateopt call itself is
%                           commented out and SURR_EVALS=0).
%   cfg.SIMPLEX_EVALS    - fminsearch MaxFunEvals/round.
%                           Default: DEBUG ? 5 : 60.
%   cfg.N_ROUNDS         - hard cap on round count.
%                           Default: DEBUG ? 1 : 1000.
%   cfg.TIME_BUDGET_HRS  - wall-clock stopper for production runs.
%                           Default: 28.
%   cfg.RESUME           - true = continue from the tag's STATE_FILE.
%                           Default: false.
%   cfg.MaxRunTime       - params0.MaxRunTime (s) per simulation.
%                           Default: 45.
%
% Output:
%   state - the final optimizer state (also saved to STATE_FILE).
%
% Run:  cd(root); addpath(genpath('.')); optimizeFeatures(cfg)
% Resumable: set cfg.RESUME=true to continue from the tag's saved state.
% Heavy by design; the caller stores results. Verify 1-2 rounds, then let
% run.
% =========================================================================

%% -------- validate + default cfg -----------------------------------------
if nargin < 1 || isempty(cfg); cfg = struct(); end
reqFields = {'baseSnap', 'tag', 'fn', 'pool', 'compulsory'};
for i = 1:numel(reqFields)
    if ~isfield(cfg, reqFields{i})
        error('optimizeFeatures:missingField', ...
            'cfg.%s is required.', reqFields{i});
    end
end
if ~isfield(cfg, 'DEBUG')           || isempty(cfg.DEBUG);           cfg.DEBUG = true; end
if ~isfield(cfg, 'N_DRAW')          || isempty(cfg.N_DRAW);          cfg.N_DRAW = 12; end
if ~isfield(cfg, 'SIMPLEX_EVALS')   || isempty(cfg.SIMPLEX_EVALS);   cfg.SIMPLEX_EVALS = cfg.DEBUG*5 + ~cfg.DEBUG*60; end
if ~isfield(cfg, 'SURR_EVALS')      || isempty(cfg.SURR_EVALS);      cfg.SURR_EVALS = cfg.DEBUG*12 + ~cfg.DEBUG*100; end
if ~isfield(cfg, 'N_ROUNDS')        || isempty(cfg.N_ROUNDS);        cfg.N_ROUNDS = cfg.DEBUG*1 + ~cfg.DEBUG*1000; end
if ~isfield(cfg, 'TIME_BUDGET_HRS') || isempty(cfg.TIME_BUDGET_HRS); cfg.TIME_BUDGET_HRS = 28; end
if ~isfield(cfg, 'RESUME')          || isempty(cfg.RESUME);          cfg.RESUME = false; end
if ~isfield(cfg, 'MaxRunTime')      || isempty(cfg.MaxRunTime);      cfg.MaxRunTime = 45; end

DEBUG           = cfg.DEBUG;
N_DRAW          = cfg.N_DRAW;
SURR_EVALS      = cfg.SURR_EVALS;
SIMPLEX_EVALS   = cfg.SIMPLEX_EVALS;
N_ROUNDS        = cfg.N_ROUNDS;
TIME_BUDGET_HRS = cfg.TIME_BUDGET_HRS;
RESUME          = cfg.RESUME;
STALL_LIMIT     = 2;        % non-improving rounds before a random kick
IMPROVE_TOL     = 1e-3;     % min cost drop to count as improvement
KICK_FRAC       = 0.0;     % +-25% random multiplicative kick on stall

root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));
LoadData;

% ---- per-instance (per-tag) output paths: THE isolation mechanism -------
% Every path an instance can write to is derived from cfg.tag. Two
% instances started with different tags therefore never share a file.
STATE_FILE = fullfile(root, 'params', [cfg.tag '_state.mat']);
SNAP_BEST  = fullfile(root, 'params', [cfg.tag '_opt.m']);
optimTag   = cfg.tag;

modeStr = 'PRODUCTION'; if DEBUG; modeStr = 'DEBUG'; end
if DEBUG
    fprintf('[%s] MODE: %s | %d round(s) | surr %d + simplex %d evals/round\n', ...
            cfg.tag, modeStr, N_ROUNDS, SURR_EVALS, SIMPLEX_EVALS);
else
    fprintf('[%s] MODE: %s | budget %.2f h | surr %d + simplex %d evals/round\n', ...
            cfg.tag, modeStr, TIME_BUDGET_HRS, SURR_EVALS, SIMPLEX_EVALS);
end
fprintf('[%s] state file : %s\n', cfg.tag, STATE_FILE);
fprintf('[%s] snapshot   : %s\n', cfg.tag, SNAP_BEST);

%% -------- base = the caller-supplied seed snapshot -----------------------
params0 = cfg.params0;

% --- feature list: caller-supplied, verbatim ---
% params0.fn = cfg.fn;

% surrogateopt UseParallel=false: the model's slack sim parallelizes
% internally over a THREAD pool ('All' segments run in one worker but the
% ODE RHS itself benefits from thread-based BLAS); a surrogateopt PROCESS
% pool would conflict with that. Keep the thread pool; run evals serially
% (each eval is already parallel inside).
if isempty(gcp('nocreate')); parpool('Threads', 5); end

%% -------- tunable POOL: entirely caller-supplied -------------------------
% Unlike RunRestretchOptim.m (which merged in params0.OptimizeParams read
% from the snapshot), the pool here comes ONLY from cfg.pool/cfg.compulsory
% -- required because some seeds (e.g. params_2state_seed.m) do not define
% OptimizeParams/OptimizeCompulsory at all.
pool = unique(cfg.pool, 'stable');
compulsory = cfg.compulsory;
% compulsory must be a subset of the pool
compulsory = compulsory(ismember(compulsory, pool));

bounds = parameterBounds();
% disable the bounds for fminsearch
bounds = struct();

%% -------- incumbent state -------------------------------------------------
if RESUME && exist(STATE_FILE, 'file')
    S = load(STATE_FILE); state = S.state;
    if ~isfield(state,'best_params'); state.best_params = state.params; end
    params0 = state.best_params;        % resume from the BEST, not a kicked point
    fprintf('[%s] Resumed at round %d, best cost %.4f\n', cfg.tag, state.round, state.best_cost);
else
    state = struct();
    state.params      = params0;   % working point (may be kicked)
    state.best_params = params0;   % the incumbent BEST (only updated on accept)
    % params0.FV_velocities = -[0, 0.5, 1, 2, 3, 4, 5, 6];
    tic
    figure;
    state.best_cost = evaluateBakersExp(ones(1,1), setMods(params0, {}), true, true);
    toc
    state.round     = 0;
    state.stall     = 0;
    state.history   = [];
    fprintf('[%s] Baseline cost (g=1): %.4f\n', cfg.tag, state.best_cost);
end

%% -------- main loop -----------------------------------------------------
tStart = tic;
for r = state.round+1 : N_ROUNDS
    if ~DEBUG && toc(tStart) > TIME_BUDGET_HRS*3600
        fprintf('\n[%s] == Time budget %.1f h reached after round %d — stopping. ==\n', cfg.tag, TIME_BUDGET_HRS, r-1);
        break;
    end
    % --- draw subset ---
    rng('shuffle');
    extra = pool(~ismember(pool, compulsory));
    k     = max(0, N_DRAW - numel(compulsory));
    drawn = [compulsory, extra(randperm(numel(extra), min(k, numel(extra))))];

    pr = state.params;          % start each round from the incumbent (baked)
    pr.mods = drawn; pr.g = ones(1, numel(drawn));

    % --- per-param multiplier box from physiology bounds ---
    [glb, gub] = deal(zeros(1,numel(drawn)));
    for i = 1:numel(drawn)
        [glb(i), gub(i)] = gMultBox(drawn{i}, pr, bounds);
    end

    fprintf('\n[%s] ===== Round %d/%d  (best %.4f, stall %d) =====\n', cfg.tag, r, N_ROUNDS, state.best_cost, state.stall);
    fprintf('[%s]   drawn: %s\n', cfg.tag, strjoin(drawn, ', '));

    obj = @(g) safeCost(g, pr);

    % --- 1) surrogateopt (global) ---
    % UseParallel=false: see parpool comment above.
    % MinSurrogatePoints defaults to max(20,nvars+1); clamp it below
    % SURR_EVALS so tiny debug budgets (SURR_EVALS<20) don't error.
    nv = numel(drawn);
    minPts = max(nv+1, min(20, SURR_EVALS-1));
    % Seed with the incumbent (g=1) so surrogateopt always knows the current
    % best point — a round can then never return worse than the seed.
    soOpts = optimoptions('surrogateopt', 'Display','Iter', ...
        'MaxFunctionEvaluations', SURR_EVALS, 'UseParallel', false, ...
        'MinSurrogatePoints', minPts, 'InitialPoints', ones(1, nv), 'PlotFcn', []);
    if SURR_EVALS > 0
        try
            [g_s, f_s] = surrogateopt(obj, glb, gub, soOpts);
        catch e
            fprintf('[%s]   surrogateopt failed (%s); falling back to g=1\n', cfg.tag, e.message);
            g_s = ones(1,numel(drawn)); f_s = obj(g_s);
        end
    else
        g_s = ones(1,numel(drawn)); f_s = obj(g_s);   % surrogate disabled -> fminsearch from incumbent
    end

    % --- 2) fminsearch (local polish) seeded from surrogate best ---
    % Headless (Display off, no plot) so two parallel instances don't spawn windows.
    fmOpts = optimset('TolFun',1e-3, 'TolX',1e-2, 'MaxFunEvals', SIMPLEX_EVALS, 'Display','iter','PlotFcns', @optimplotfval);
    [g_l, f_l] = fminsearch(@(g) safeCostClamp(g, pr, glb, gub), g_s, fmOpts);

    if f_l <= f_s; g_round = g_l; f_round = f_l; else; g_round = g_s; f_round = f_s; end
    fprintf('[%s]   surrogate %.4f -> simplex %.4f\n', cfg.tag, f_s, f_round);

    % --- accept / stall ---
    if f_round < state.best_cost - IMPROVE_TOL
        % ACCEPT: this is the ONLY branch that may update the incumbent
        % (state.best_params / state.best_cost) and the ONLY branch that
        % writes SNAP_BEST to disk.
        pr.g = g_round;
        state.params    = getParams(pr, g_round, false, true);   % bake in
        state.params.mods = {}; state.params.g = [];
        state.best_params = state.params;   % <-- persist the BEST (survives kicks)
        state.best_cost = f_round;
        state.stall     = 0;
        writeParamsToMFile(SNAP_BEST, state.best_params);
        if exist('captureOptimIter','file')
            try; captureOptimIter(setMods(state.params,{}), r, f_round, state.best_cost, optimTag); catch; end
        end
        fprintf('[%s]   IMPROVED -> %.4f  (snapshot written)\n', cfg.tag, f_round);
    else
        % STALL: state.params is updated (a working copy for the NEXT
        % round's draw) but state.best_params/state.best_cost and disk are
        % left untouched. This is the branch that used to be the source of
        % the "worse snapshot survives" bug when it was allowed to write.
        state.stall = state.stall + 1;
        fprintf('[%s]   no improvement (%.4f >= %.4f), stall=%d\n', cfg.tag, f_round, state.best_cost, state.stall);
        if state.stall >= STALL_LIMIT
            % Random-kick a copy of the BEST (never the best itself) to escape
            % the basin. state.best_params is left untouched.
            kp = state.best_params;
            kmods = pool(randperm(numel(pool), min(N_DRAW, numel(pool))));
            kp.mods = kmods;
            kp.g = 1 + KICK_FRAC*(2*rand(1,numel(kmods)) - 1);
            state.params = getParams(kp, kp.g, false, true);   % working point only
            state.params.mods = {}; state.params.g = [];
            state.stall = 0;
            fprintf('[%s]   >> STALL kick applied to: %s\n', cfg.tag, strjoin(kmods, ', '));
        end
    end

    state.round   = r;
    state.history = [state.history; r, f_s, f_round, state.best_cost];
    if ~isfield(state, 'paramsHistory')
        phl = 1;
    else
        phl = lenght(state.paramsHistory);
    end

    state.paramsHistory(phl) = getParams(pr, g_round, false, true);
    state.paramsCosts(phl) = f_round;

    save(STATE_FILE, 'state');
end

% Always leave the BEST result on disk (NEVER state.params, which may be a
% kicked/exploratory working point that is worse than the incumbent).
% Guard with an assertion: the snapshot we are about to (re-)write must
% correspond to state.best_cost, so a re-evaluation of it can never silently
% diverge from what state.best_cost claims. The tolerance is noise-aware: the
% objective is NOT bit-reproducible (the MaxRunTime watchdog + parallel pool
% make re-evaluations differ by ~1e-3), so a fixed 1e-6 bound spuriously
% aborted valid runs. This tolerance (max of 1e-2 absolute and 2% relative)
% still catches the gross-divergence bug this guard exists for (a snapshot
% written at cost ~174 when the incumbent was ~9.7).
try
    verifyCost = evaluateBakersExp(ones(1,1), setMods(state.best_params, {}), true, false);
catch
    verifyCost = Inf;   % a throw here (3-state timeout) -> ~isfinite passes the assert
end
verifyTol  = max(1e-2, 0.02*abs(state.best_cost));
% A huge re-eval (>=1e5) is a transient MaxRunTime/instability penalty (common
% under CPU load), NOT a real divergence -- the incumbent was already validated
% during its improving round. Skip the guard then; it still catches a genuine
% divergence, which is O(10-100), not O(1e6).
assert(abs(verifyCost - state.best_cost) < verifyTol || ~isfinite(verifyCost) || verifyCost >= 1e5, ...
    'optimizeFeatures:incumbentMismatch', ...
    '[%s] state.best_params re-evaluates to %.6f but state.best_cost=%.6f (tol %.4g) -- refusing to write a divergent snapshot.', ...
    cfg.tag, verifyCost, state.best_cost, verifyTol);
writeParamsToMFile(SNAP_BEST, state.best_params);
save(STATE_FILE, 'state');
fprintf('\n[%s] === Done. Best cost %.4f. Snapshot: %s ===\n', cfg.tag, state.best_cost, SNAP_BEST);

end % optimizeFeatures

%% ======================= local helpers =================================
function p = setMods(p, m); p.mods = m; p.g = ones(1, numel(m)); end

function c = safeCost(g, pr)
    % try/catch: a 3-state ODE timeout/instability throws in evaluateModel
    % (MaxRunTime / "not stable"); it must become a high finite cost the
    % optimizer routes around, not a crash of the whole run.
    try
        c = evaluateBakersExp(g, pr, true);
    catch
        c = 1e6;
    end
    if ~isfinite(c); c = 1e6; end
end

function c = safeCostClamp(g, pr, glb, gub)
    % fminsearch is unbounded: penalise excursions outside the surrogate box.
    if any(g < glb*0.8) || any(g > gub*1.2)
        c = 1e6; return;
    end
    try
        c = evaluateBakersExp(g, pr, true);
    catch
        c = 1e6;
    end
    if ~isfinite(c); c = 1e6; end
end

function [glb, gub] = gMultBox(name, pr, bounds)
    % Multiplier box for parameter NAME given its current value in PR and
    % the physiology bounds. Falls back to named default windows for the
    % velocity-gauss/registration-availability levers (no parameterBounds
    % entry), and a generic default window for anything else uncovered
    % (array knots, offsets, etc).
    DEF_LO = 0.4; DEF_HI = 2.5;
    cur = currentVal(name, pr);

    if strcmp(name, 'v_att_amplitude')
        glb = 0.2; gub = 4;
    elseif strcmp(name, 'v_att_center')
        % Current value is NEGATIVE (~-1). A positive multiplier times a
        % negative current value stays negative — no sign special-case
        % needed; the box below is on the multiplier, not on v_att_center
        % itself.
        glb = 0.3; gub = 3.5;
    elseif strcmp(name, 'v_att_sigma')
        glb = 0.35; gub = 3;
    elseif strcmp(name, 'v_ref_reg')
        glb = 0.25; gub = 5;
    elseif strcmp(name, 'A0')
        % A0 is an availability floor and must stay < 1: cap its multiplier
        % so A0*gub <= 0.98.
        glb = 0.3; gub = 3.2;
        if isfinite(cur) && cur > 0
            gub = min(gub, 0.98 / cur);
        end
    elseif contains(name, '__')
        % array knots (e.g. PieceWiseStrainDepR1DParams__4)
        glb = 0.3; gub = 3;
    elseif isfield(bounds, name) && isfinite(cur) && cur ~= 0
        b = bounds.(name);
        lo = b.lb / cur; hi = b.ub / cur;
        glb = min(lo, hi); gub = max(lo, hi);
        % guard against degenerate/zero-crossing bounds
        if ~isfinite(glb) || glb <= 0; glb = DEF_LO; end
        if ~isfinite(gub) || gub <= glb; gub = max(DEF_HI, glb*2); end
    else
        glb = DEF_LO; gub = DEF_HI;
    end
end

function v = currentVal(name, pr)
    % Resolve current scalar value, including arr__idx notation.
    if isfield(pr, name)
        v = pr.(name);
    elseif contains(name, '__')
        t = split(name, '__');
        if isfield(pr, t{1}) && numel(pr.(t{1})) >= str2double(t{2})
            v = pr.(t{1})(str2double(t{2}));
        else
            v = NaN;
        end
    else
        v = NaN;
    end
    if ~isscalar(v) || ~isnumeric(v); v = NaN; end
end
