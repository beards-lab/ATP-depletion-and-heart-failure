% RunSurrogateSimplex_Optim.m
% =========================================================================
% Decisive global->local optimizer for the reglu-availability 2-state fit.
%
% Per round:
%   1. DRAW a ~10-param subset from a larger mechanism-organised pool
%      (compulsory levers always in; remainder random) — keeps each
%      surrogateopt problem low-dimensional and tractable.
%   2. SURROGATEOPT (global, bounded, parallel) explores the drawn box.
%   3. FMINSEARCH (local) polishes the surrogate's best point.
%   4. If the round IMPROVES the incumbent, BAKE the multipliers into the
%      running params and snapshot to disk. Otherwise count a STALL.
%   5. On STALL_LIMIT consecutive non-improving rounds, RANDOM-KICK the
%      incumbent and force a fresh draw (escape local basins).
%
% Bounds: each drawn parameter's g-multiplier box is derived from
% parameterBounds.m (absolute [lb,ub] -> multiplier around the current
% value), so surrogateopt stays inside physiology. Params without a bound
% entry (array knots, offsets, A0/v_ref_reg/tau_reg) get a default window.
%
% Cost: evaluateBakersExp(g, params0) — feature cost (FV+ktr+slack) +
% literature OUTPUT bounds + soft AssertParams input guard.
%
% Run:  cd(root); addpath(genpath('.')); RunSurrogateSimplex_Optim
% Resumable: set RESUME=true to continue from the saved state .mat.
% NOTE: ~25-35 s/eval. A round = (surrogate budget + simplex budget) evals.
% Heavy by design; the user stores results. Verify 1-2 rounds, then let run.
% =========================================================================

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root));
LoadData;

STATE_FILE = fullfile(root, 'params', 'surrsimplex_state.mat');
SNAP_BEST  = fullfile(root, 'params', 'params_surrsimplex_best.m');
optimTag   = 'SurrSimplex';
SNAP_BEST  = fullfile(root, 'params', 'params_restretch_best.m');
optimTag   = 'restretch';

% ---- tunables of the OPTIMIZER itself ----
DEBUG = false;             % <<< true = ~10 min smoke test; false = grid run >>>
if DEBUG
    % Fast smoke test: proves the whole pipeline runs end-to-end in ~10 min.
    N_ROUNDS        = 1;
    TIME_BUDGET_HRS = 0.5;
    SURR_EVALS      = 12;  % < 20 (an 11-pt initial design + 1 adaptive eval)
    SIMPLEX_EVALS   = 5;
else
    % Production / grid: TIME_BUDGET_HRS is the real stopper — set it to the
    % grid wall-time. N_ROUNDS is just a high cap so the budget governs.
    N_ROUNDS        = 1000;
    TIME_BUDGET_HRS = 3;   % <<< set the grid wall-clock budget here (hours) >>>
    SURR_EVALS      = 120;
    SIMPLEX_EVALS   = 120;
end
N_DRAW         = 10;       % params optimised per round (incl. compulsory)
STALL_LIMIT    = 2;        % non-improving rounds before a random kick
IMPROVE_TOL    = 1e-3;     % min cost drop to count as improvement
KICK_FRAC      = 0.25;     % +-25% random multiplicative kick on stall
RESUME         = false;    % set true to continue from STATE_FILE
modeStr = 'PRODUCTION'; if DEBUG; modeStr = 'DEBUG'; end
fprintf('MODE: %s | budget %.2f h | surr %d + simplex %d evals/round\n', ...
        modeStr, TIME_BUDGET_HRS, SURR_EVALS, SIMPLEX_EVALS);

%% -------- base = the restretch-tuned snapshot (RegAvailShorteningOnly on) --
% From Analyses/RestretchFeatureFit: directional availability + restretch/FV
% rebalance took the feature cost 15.9 -> 6.2. Seed the grid run from that basin.
params0 = getParams();
run(fullfile(root, 'params', 'params_restretch_best.m'));
params0.RegAvailShorteningOnly = true;   % restretch/lengthening does not boost attachment
params0 = getParams(params0, [], true, false);

params0.RunForceVelocity = true; params0.RunKtr = false; params0.RunSlack = true;
params0.RunStairs = false; params0.RunForceVelocityTime = false;
params0.EvalFeatures = true; params0.BreakOnODEUnstable = false; params0.MaxRunTime = 120;
params0.PlotEachSeparately = 1; params0.PlotFeatureFitting = 1;
params0.RunSlackSegments = 'AllPar';
params0.FV_velocities = -[0 0.5 1 2 4];

% --- feature list: canonical bounds list, with two justified reweights ---
%   vall2_dy 0.1 -> 1.0  : model lost the ~13 kPa post-restretch sag entirely
%                          (regavail collapsed it); this miss must be visible.
%   ovrsht_dy 0.1 -> 0.01: data mean ~1.3 kPa -> normalisation is degenerate,
%                          and the model value is the 2-state oscillation artifact.
params0.fn = boundedOutputFn(0.001);
params0.fn = strrep(params0.fn, 'vall2_dy|0.1',  'vall2_dy|1.0');
params0.fn = strrep(params0.fn, 'ovrsht_dy|0.1', 'ovrsht_dy|1');
% ktr is iron-law-locked in 2-state (ktr ~= ka_eff + k2_eff): chasing the ktr
% VALUE just distorts force/duty, so drop the 'ktr|2' feature and the RunKtr
% protocol (disabled above). KEEP 'ktr_rmse|0-0.2' — it is a fit-quality bound
% (penalises non-exponential/oscillatory recovery), not the ktr value itself.
% isKtr = startsWith(params0.fn, 'ktr|');
% params0.fn = params0.fn(~isKtr);

if isempty(gcp('nocreate')); parpool('Threads', 5); end

%% -------- mechanism-organised tunable POOL ------------------------------
% Grouped by the hard-to-fit feature each lever targets (see header notes).
ktr_iron     = {'k2', 'ka', 'k_2', ...
                'PieceWiseStrainDep2Params__5', 'PieceWiseStrainDep2Params__6', ...
                'PieceWiseStrainDep2Params__7'};          % ktr / FV detachment shape
onset_delay  = {'kSE', 'ekSE', 'tau_reg'};                % t0 onset growth
slow_sag     = {'kSE_M', 'eta_M', 'mu', ...
                'PieceWiseStrainDepR1DParams__2', 'PieceWiseStrainDepR1DParams__3', ...
                'PieceWiseStrainDepR1DParams__4'}; % vall2/vall_y: __4 = positive-strain clamp (key restretch lever)
force_scale  = {'kstiff1', 'kstiff2', 'dr'};              % peak2 / A / steady
fv_shoulder  = {'A0', 'v_ref_reg'};                       % FV shape (regavail)
srx_pool     = {'ksr0', 'kmsr', 'SRXT_0', 'ksr2srd'};     % SRX_ss / PT_ss floors
OptimizeParams = {'k2', 'v_att_amplitude', 'A0', 'ka', 'kstiff2', ...
    'v_att_center', 'v_att_sigma', 'v_ref_reg', 'kstiff1', 'kSE', 'k_2', 'k1', ...
    'PieceWiseStrainDepR1DParams__4', 'PieceWiseStrainDepR1DParams__3', 'ksr0'};

pool = unique([ktr_iron, onset_delay, slow_sag, force_scale, fv_shoulder, srx_pool, OptimizeParams], 'stable');

% Always-in levers: the FV shoulder shape + the strongest force/duty knobs.
compulsory = {'A0', 'v_ref_reg', 'k2', 'kstiff2'};
compulsory = {'k2', 'v_att_amplitude', 'A0', 'ka', 'kstiff2'};

bounds = parameterBounds();

params0.fn = {'FV_fnorm|FV_v|10', 'ktr|2', 'A|50', 'ktr_rmse|0-.5|.1', ...
    'XTOR[1]|3-15|15,XTOR_vmax[1]|5-25,SRX_ss[1]|0.01-0.40|5,attached_ss[1]|0.10-0.45|10,PT_ss[1]|0.01-0.70',...
    't0_crossing|SLdiff|2', 'restretchSlopeStart|1', 'peak1_y|10', 'peak1_dSL|1', 'vall_y|10', ...
    'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|1.0', 'ovrsht_dy|1', ...
    'k2|100-1500|0.01,ksrd|0.-20|0.001,kmsrd|0.0-20|0.001'};
%% -------- incumbent state ----------------------------------------------
if RESUME && exist(STATE_FILE, 'file')
    S = load(STATE_FILE); state = S.state;
    params0 = state.params;            % already baked
    fprintf('Resumed at round %d, best cost %.4f\n', state.round, state.best_cost);
else
    state = struct();
    state.params    = params0;
    tic
    state.best_cost = evaluateBakersExp(ones(1,1), setMods(params0, {}), true, true);
    toc
    state.round     = 0;
    state.stall     = 0;
    state.history   = [];
    fprintf('Baseline cost (g=1): %.4f\n', state.best_cost);
end

%% -------- main loop -----------------------------------------------------
tStart = tic;
for r = state.round+1 : N_ROUNDS
    if toc(tStart) > TIME_BUDGET_HRS*3600
        fprintf('\n== Time budget %.1f h reached after round %d — stopping. ==\n', TIME_BUDGET_HRS, r-1);
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

    fprintf('\n===== Round %d/%d  (best %.4f, stall %d) =====\n', r, N_ROUNDS, state.best_cost, state.stall);
    fprintf('  drawn: %s\n', strjoin(drawn, ', '));

    obj = @(g) safeCost(g, pr);

    % --- 1) surrogateopt (global, parallel) ---
    % UseParallel=false: surrogateopt's own parallelism needs a PROCESS pool,
    % but the model's slack 'AllPar' uses the THREAD pool internally for its 5
    % chunks. Keep the thread pool; run surrogate evals serially (each eval is
    % already parallel inside). Process pools break the thread-based AllPar path.
    % MinSurrogatePoints defaults to max(20,nvars+1); clamp it below SURR_EVALS
    % so tiny debug budgets (SURR_EVALS<20) don't error.
    nv = numel(drawn);
    minPts = max(nv+1, min(20, SURR_EVALS-1));
    % Seed with the incumbent (g=1) so surrogateopt always knows the current
    % best point — a round can then never return worse than the seed.
    soOpts = optimoptions('surrogateopt', 'Display','iter', ...
        'MaxFunctionEvaluations', SURR_EVALS, 'UseParallel', false, ...
        'MinSurrogatePoints', minPts, 'InitialPoints', ones(1, nv), 'PlotFcn', 'surrogateoptplot');
    try
        [g_s, f_s] = surrogateopt(obj, glb, gub, soOpts);
    catch e
        fprintf('  surrogateopt failed (%s); falling back to g=1\n', e.message);
        g_s = ones(1,numel(drawn)); f_s = obj(g_s);
    end

    % --- 2) fminsearch (local polish) seeded from surrogate best ---
    fmOpts = optimset('Display','off', 'TolFun',1e-3, 'TolX',1e-2, 'MaxFunEvals', SIMPLEX_EVALS);
    [g_l, f_l] = fminsearch(@(g) safeCostClamp(g, pr, glb, gub), g_s, fmOpts);

    if f_l <= f_s; g_round = g_l; f_round = f_l; else; g_round = g_s; f_round = f_s; end
    fprintf('  surrogate %.4f -> simplex %.4f\n', f_s, f_round);

    % --- accept / stall ---
    if f_round < state.best_cost - IMPROVE_TOL
        pr.g = g_round;
        state.params    = getParams(pr, g_round, false, true);   % bake in
        state.params.mods = {}; state.params.g = [];
        state.best_cost = f_round;
        state.stall     = 0;
        writeParamsToMFile(SNAP_BEST, state.params);
        if exist('captureOptimIter','file')
            try; captureOptimIter(setMods(state.params,{}), r, f_round, state.best_cost, optimTag); catch; end
        end
        fprintf('  IMPROVED -> %.4f  (snapshot written)\n', f_round);
    else
        state.stall = state.stall + 1;
        fprintf('  no improvement (%.4f >= %.4f), stall=%d\n', f_round, state.best_cost, state.stall);
        if state.stall >= STALL_LIMIT
            % random kick the incumbent within bounds, then continue
            kp = state.params;
            kmods = pool(randperm(numel(pool), min(N_DRAW, numel(pool))));
            kp.mods = kmods;
            kp.g = 1 + KICK_FRAC*(2*rand(1,numel(kmods)) - 1);
            state.params = getParams(kp, kp.g, false, true);
            state.params.mods = {}; state.params.g = [];
            state.stall = 0;
            fprintf('  >> STALL kick applied to: %s\n', strjoin(kmods, ', '));
        end
    end

    state.round   = r;
    state.history = [state.history; r, f_s, f_round, state.best_cost];
    save(STATE_FILE, 'state');
end

% Always leave a usable result on disk (grid runs must produce a snapshot even
% if no round improved on the seed).
writeParamsToMFile(SNAP_BEST, state.params);
save(STATE_FILE, 'state');
fprintf('\n=== Done. Best cost %.4f. Snapshot: %s ===\n', state.best_cost, SNAP_BEST);

%% ======================= local helpers =================================
function p = setMods(p, m); p.mods = m; p.g = ones(1, numel(m)); end

function c = safeCost(g, pr)
    c = evaluateBakersExp(g, pr, true);
    if ~isfinite(c); c = 1e6; end
end

function c = safeCostClamp(g, pr, glb, gub)
    % fminsearch is unbounded: penalise excursions outside the surrogate box.
    if any(g < glb*0.8) || any(g > gub*1.2)
        c = 1e6; return;
    end
    c = evaluateBakersExp(g, pr, true);
    if ~isfinite(c); c = 1e6; end
end

function [glb, gub] = gMultBox(name, pr, bounds)
    % Multiplier box for parameter NAME given its current value in PR and
    % the physiology bounds. Falls back to a default window for params with
    % no bound entry (array knots, offsets, regavail levers).
    DEF_LO = 0.4; DEF_HI = 2.5;
    cur = currentVal(name, pr);
    base = strtok(name, '_');                 % crude; exact match preferred below
    if isfield(bounds, name) && isfinite(cur) && cur ~= 0
        b = bounds.(name);
        lo = b.lb / cur; hi = b.ub / cur;
        glb = min(lo, hi); gub = max(lo, hi);
        % guard against degenerate/zero-crossing bounds
        if ~isfinite(glb) || glb <= 0; glb = DEF_LO; end
        if ~isfinite(gub) || gub <= glb; gub = max(DEF_HI, glb*2); end
    else
        glb = DEF_LO; gub = DEF_HI;
    end
    % A0 is an availability floor and must stay < 1: cap its multiplier.
    if strcmp(name, 'A0') && isfinite(cur) && cur > 0
        gub = min(gub, 0.999 / cur);
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
