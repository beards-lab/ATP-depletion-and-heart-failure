% RunRestretchOptim.m
% =========================================================================
% Decisive global->local optimizer for the restretch/FV feature fit, seeded
% from params/params_restretch_best.m. Adapted from
% Workbench/RunSurrogateSimplex_Optim.m (that file is owned by the user and
% is NOT modified here) — same per-round subset draw -> surrogateopt
% (global) -> fminsearch (local polish) -> accept/stall -> random-kick
% architecture, retargeted at the restretch snapshot's own tunable pool.
%
% Per round:
%   1. DRAW a N_DRAW-param subset from params0.OptimizeParams (compulsory
%      levers from params0.OptimizeCompulsory always in; remainder random).
%      Keeps each surrogateopt problem low-dimensional and tractable.
%   2. SURROGATEOPT (global, bounded) explores the drawn box.
%   3. FMINSEARCH (local) polishes the surrogate's best point.
%   4. If the round IMPROVES the incumbent, BAKE the multipliers into the
%      running params and snapshot to disk. Otherwise count a STALL.
%   5. On STALL_LIMIT consecutive non-improving rounds, RANDOM-KICK the
%      incumbent and force a fresh draw (escape local basins).
%
% Bounds: each drawn parameter's g-multiplier box is derived from
% parameterBounds.m (absolute [lb,ub] -> multiplier around the current
% value) via gMultBox. Params without a bound entry (array knots, velocity-
% gauss levers, A0/v_ref_reg) get named default windows (see gMultBox).
%
% Cost: evaluateBakersExp(g, params0) -> RunBakersExp -> evalFeatureCost on
% the hardcoded params0.fn feature list set below (FV + slack + restretch
% shape features).
%
% Run:  cd(root); addpath(genpath('.')); RunRestretchOptim
% Resumable: set RESUME=true to continue from the saved state .mat.
% Heavy by design; the user stores results. Verify 1-2 rounds, then let run.
% =========================================================================

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));
LoadData;

STATE_FILE = fullfile(root, 'params', 'restretchopt2_state.mat');
SNAP_BEST  = fullfile(root, 'params', 'params_restretch_opt2.m');
optimTag   = 'restretchopt2';

% ---- tunables of the OPTIMIZER itself ----
DEBUG = false;              % <<< true = ~10 min smoke test; false = production run >>>
if DEBUG
    % Fast smoke test: proves the whole pipeline runs end-to-end in ~10 min.
    N_ROUNDS        = 1;
    SURR_EVALS      = 12;  % < 20 (an 11-pt initial design + 1 adaptive eval)
    SIMPLEX_EVALS   = 5;
else
    % Production: TIME_BUDGET_HRS is the real stopper — N_ROUNDS is just a
    % high cap so the wall-clock budget governs.
    N_ROUNDS        = 1000;
    TIME_BUDGET_HRS = 28;   % <<< set the wall-clock budget here (hours) >>>
    SURR_EVALS      = 0;
    SIMPLEX_EVALS   = 60;
end
N_DRAW         = 12;        % params optimised per round (incl. compulsory)
STALL_LIMIT    = 2;        % non-improving rounds before a random kick
IMPROVE_TOL    = 1e-3;     % min cost drop to count as improvement
KICK_FRAC      = 0.25;     % +-25% random multiplicative kick on stall
RESUME         = false;    % set true to continue from STATE_FILE
modeStr = 'PRODUCTION'; if DEBUG; modeStr = 'DEBUG'; end
if DEBUG
    fprintf('MODE: %s | %d round(s) | surr %d + simplex %d evals/round\n', ...
            modeStr, N_ROUNDS, SURR_EVALS, SIMPLEX_EVALS);
else
    fprintf('MODE: %s | budget %.2f h | surr %d + simplex %d evals/round\n', ...
            modeStr, TIME_BUDGET_HRS, SURR_EVALS, SIMPLEX_EVALS);
end

%% -------- base = the restretch-tuned snapshot ---------------------------
params0 = getParams();
run(fullfile(root, 'params', 'params_restretch_best.m'));
params0 = getParams(params0, [], true, false);

params0.RunForceVelocity = true; params0.RunKtr = false; params0.RunSlack = true;
params0.RunStairs = false; params0.RunForceVelocityTime = false;
params0.EvalFeatures = true; params0.BreakOnODEUnstable = false; params0.MaxRunTime = 45;
params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0;
params0.RunSlackSegments = 'AllPar';
params0.FV_velocities = -[0 0.5 1 2 4];
params0.MaxRunTime = 30;

% --- feature list: hardcoded per spec (verbatim from params_restretch_best) ---
params0.fn = {'FV_fnorm|FV_v|10', 'ktr|2', 'A|50', 'ktr_rmse|0-.5|.1', ...
    'XTOR[1]|3-15|15,XTOR_vmax[1]|5-25,SRX_ss[1]|0.01-0.40|5,attached_ss[1]|0.10-0.45|10,PT_ss[1]|0.01-0.70',...
    't0_crossing|SLdiff|2', 'restretchSlopeStart|1', 'peak1_y|10', 'peak1_dSL|1', 'vall_y|10', ...
    'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|1.0', 'ovrsht_dy|1', ...
    'k2|100-1500|0.01,ksrd|0.-20|0.001,kmsrd|0.0-20|0.001'};

% surrogateopt UseParallel=false: the model's slack sim parallelizes
% internally over a THREAD pool ('All' segments run in one worker but the
% ODE RHS itself benefits from thread-based BLAS); a surrogateopt PROCESS
% pool would conflict with that. Keep the thread pool; run evals serially
% (each eval is already parallel inside).
if isempty(gcp('nocreate')); parpool('Threads', 5); end

%% -------- tunable POOL: read from the snapshot, not hardcoded -----------
% params0.OptimizeParams / params0.OptimizeCompulsory live in
% params/params_restretch_best.m so the user can edit the snapshot's own
% tunable set without touching this driver.
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


if ~isfield(params0, 'OptimizeParams') || isempty(params0.OptimizeParams)
    error('RunRestretchOptim:noPool', ...
        'params0.OptimizeParams is missing/empty — define it in params_restretch_best.m');
end
pool = unique([ktr_iron, onset_delay, slow_sag, force_scale, fv_shoulder, srx_pool, OptimizeParams,params0.OptimizeParams], 'stable');

if isfield(params0, 'OptimizeCompulsory') && ~isempty(params0.OptimizeCompulsory)
    compulsory = params0.OptimizeCompulsory;
else
    compulsory = {};
end
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
    fprintf('Resumed at round %d, best cost %.4f\n', state.round, state.best_cost);
else
    state = struct();
    state.params      = params0;   % working point (may be kicked)
    state.best_params = params0;   % the incumbent BEST (only updated on accept)
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
    if ~DEBUG && toc(tStart) > TIME_BUDGET_HRS*3600
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

    % --- 1) surrogateopt (global) ---
    % UseParallel=false: see parpool comment above.
    % MinSurrogatePoints defaults to max(20,nvars+1); clamp it below
    % SURR_EVALS so tiny debug budgets (SURR_EVALS<20) don't error.
    nv = numel(drawn);
    minPts = max(nv+1, min(20, SURR_EVALS-1));
    % Seed with the incumbent (g=1) so surrogateopt always knows the current
    % best point — a round can then never return worse than the seed.
    soOpts = optimoptions('surrogateopt', 'Display','off', ...
        'MaxFunctionEvaluations', SURR_EVALS, 'UseParallel', false, ...
        'MinSurrogatePoints', minPts, 'InitialPoints', ones(1, nv), 'PlotFcn', []);
    % try
    %     [g_s, f_s] = surrogateopt(obj, glb, gub, soOpts);
    % catch e
    %     fprintf('  surrogateopt failed (%s); falling back to g=1\n', e.message);
        g_s = ones(1,numel(drawn)); f_s = obj(g_s);
    % end

    % --- 2) fminsearch (local polish) seeded from surrogate best ---
    fmOpts = optimset('Display','iter', 'TolFun',1e-3, 'TolX',1e-2, 'MaxFunEvals', SIMPLEX_EVALS, 'PlotFcns', @optimplotfval);
    [g_l, f_l] = fminsearch(@(g) safeCostClamp(g, pr, glb, gub), g_s, fmOpts);

    if f_l <= f_s; g_round = g_l; f_round = f_l; else; g_round = g_s; f_round = f_s; end
    fprintf('  surrogate %.4f -> simplex %.4f\n', f_s, f_round);

    % --- accept / stall ---
    if f_round < state.best_cost - IMPROVE_TOL
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
        fprintf('  IMPROVED -> %.4f  (snapshot written)\n', f_round);
    else
        state.stall = state.stall + 1;
        fprintf('  no improvement (%.4f >= %.4f), stall=%d\n', f_round, state.best_cost, state.stall);
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
            fprintf('  >> STALL kick applied to: %s\n', strjoin(kmods, ', '));
        end
    end

    state.round   = r;
    state.history = [state.history; r, f_s, f_round, state.best_cost];
    save(STATE_FILE, 'state');
end

% Always leave the BEST result on disk (NEVER state.params, which may be a
% kicked/exploratory working point that is worse than the incumbent).
writeParamsToMFile(SNAP_BEST, state.best_params);
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
