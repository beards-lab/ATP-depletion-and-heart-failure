% IdentifyLastSlack.m
%
% Two-phase parameter identification for the last slack-restretch cycle.
% Starting point: params/ModelParamsInitManualLastSlack.m
%
% Phase 1 — fminsearch on {kSE, mu, ka, kd, k1, k2, slope, k_pas} only.
%            Expected: cannot reproduce velocity-dependent double-peak.
% Phase 2 — same + k_catch_bond free (catch bond enabled).
%            Expected: double-peak reproduced.
%
% State is checkpointed to IdentifyLastSlack_state.mat so re-running
% the script continues from the last checkpoint.

STATE_FILE      = 'IdentifyLastSlack_state.mat';
MAX_EVALS       = 150;   % evals per call; re-run script to continue

fn_target   = {'ktr|SLslack','A|SLslack','peak1_y','peak1_dSL', ...
               'peak2','vall_y','vall2_dy','steady'};
mods_phase1 = {'kSE','mu','ka','kd','k1','k2','slope','k_pas'};

%% ── Setup ────────────────────────────────────────────────────────────────
clear params0;
params0 = getParams();
ModelParamsInitManualLastSlack;
LoadData;

params0.RunSlackSegments    = 'Last';
params0.RunForceVelocity    = false;
params0.RunForceLengthEstim = false;
params0.PlotEachSeparately  = false;
params0.ShowResidualPlots   = false;
params0.ghostLoad = ''; params0.ghostSave = '';
params0.EvalFeatures        = true;
params0.BreakOnODEUnstable  = false;
params0.fn   = fn_target;
params0.mods = mods_phase1;
params0.g    = ones(1, numel(mods_phase1));

opts = optimset('Display','iter','TolFun',1e-4, ...
                'MaxIter',MAX_EVALS,'MaxFunEvals',MAX_EVALS);

%% ── Load or init state ───────────────────────────────────────────────────
if exist(STATE_FILE,'file')
    load(STATE_FILE,'state');
    fprintf('Resuming: phase=%d, evals_so_far=%d, best_feat=%.4f\n', ...
        state.phase, state.total_evals, state.best_feat);
else
    state.phase       = 0;
    state.g_opt1      = ones(1, numel(mods_phase1));
    state.g_opt2      = [];
    state.best_feat   = Inf;
    state.total_evals = 0;
    state.costs_base  = [];
    state.costs1      = [];
    state.costs2      = [];
end

%% ── Phase 0: Baseline ────────────────────────────────────────────────────
if state.phase < 1
    fprintf('\n=== Phase 0: Baseline ===\n');
    [E_base, ~, fm0, fd0] = runSlackExperiment(params0);
    state.costs_base = evalFeatureCost(fd0, fm0, fn_target, 1);
    fprintf('E_slack=%.2f  FeatTot=%.3f\n', sum(E_base), sum(state.costs_base));
    printCosts(fn_target, state.costs_base);
    state.phase = 1;
    save(STATE_FILE,'state');
else
    fprintf('Phase 0 done (FeatTot=%.3f). Continuing.\n', sum(state.costs_base));
end

%% ── Phase 1: No catch bond ───────────────────────────────────────────────
if state.phase == 1
    fprintf('\n=== Phase 1: No catch bond (max %d evals) ===\n', MAX_EVALS);
    g_new = fminsearch(@(g) slack_cost(g, params0, fn_target), state.g_opt1, opts);

    p1 = params0; p1.g = g_new;
    [E1, ~, fm1, fd1] = runSlackExperiment(p1);
    c1 = evalFeatureCost(fd1, fm1, fn_target, 1);
    state.total_evals = state.total_evals + MAX_EVALS;

    fprintf('\nPhase 1 result: FeatTot=%.4f (prev best=%.4f, evals=%d)\n', ...
        sum(c1), state.best_feat, state.total_evals);
    printCosts(fn_target, c1);
    fprintf('g: %s\n', num2str(g_new,'%.3f '));

    improved = sum(c1) < state.best_feat - 0.01;
    state.g_opt1  = g_new;
    state.costs1  = c1;
    if improved
        state.best_feat = sum(c1);
    end

    % Advance to Phase 2 when gain per run is negligible or budget used
    if state.total_evals >= 4*MAX_EVALS || ~improved
        fprintf('→ Phase 1 done. Moving to Phase 2.\n');
        state.phase     = 2;
        state.best_feat = Inf;  % reset for phase 2 tracking
    else
        fprintf('→ Re-run script to continue Phase 1.\n');
    end
    save(STATE_FILE,'state');
end

%% ── Phase 2: ViscoSE (one-sided dashpot, vel>0) ──────────────────────────
% c_SE_visc adds viscous resistance to SE only during restretch (vel>0).
% This damps peak1 overshoot without slowing ktr (isometric, vel≈0).
mods_phase2 = [mods_phase1, {'c_SE_visc'}];

if state.phase == 2
    params0_v = params0;
    params0_v.UseViscoelasticSE = true;
    params0_v.c_SE_visc         = 30;     % base value; g multiplies this
    % SE time constant tau = (mu + c_SE_visc)/kSE ≈ (0.3+30)/3705 ≈ 8ms,
    % comparable to the ~20ms restretch duration → peak can be meaningfully delayed
    params0_v.mods = mods_phase2;

    if isempty(state.g_opt2)
        state.g_opt2 = [state.g_opt1, 1.0];
    end

    fprintf('\n=== Phase 2: +ViscoSE dashpot (max %d evals) ===\n', MAX_EVALS);
    g_new2 = fminsearch(@(g) slack_cost(g, params0_v, fn_target), state.g_opt2, opts);

    pv = params0_v; pv.g = g_new2;
    [~, ~, fm2, fd2] = runSlackExperiment(pv);
    c2 = evalFeatureCost(fd2, fm2, fn_target, 1);
    state.total_evals = state.total_evals + MAX_EVALS;

    improved2 = sum(c2) < state.best_feat - 0.01;
    state.g_opt2  = g_new2;
    state.costs2  = c2;
    if improved2; state.best_feat = sum(c2); end

    fprintf('\nPhase 2 result: FeatTot=%.4f (evals=%d)\n', sum(c2), state.total_evals);
    c_SE_opt = params0_v.c_SE_visc * g_new2(end);
    fprintf('c_SE_visc = %.4f\n', c_SE_opt);
    fprintf('g: %s\n', num2str(g_new2,'%.3f '));
    printCosts(fn_target, c2);

    if state.total_evals >= 8*MAX_EVALS || ~improved2
        state.phase = 3;
        fprintf('→ Phase 2 done. Moving to Phase 3 (Maxwell).\n');
    else
        fprintf('→ Re-run script to continue Phase 2.\n');
    end
    save(STATE_FILE,'state');
end

%% ── Phase 3: Maxwell dashpot comparison ─────────────────────────────────
% UseMaxwellDashpot adds a parallel spring-dashpot. Unlike ViscoSE it
% affects isometric force slightly. Compare to Phase 2 to see which wins.
mods_phase3 = [mods_phase1, {'kM', 'eta_M'}];

if state.phase == 3
    params0_m = params0;
    params0_m.UseMaxwellDashpot = true;
    params0_m.kM    = 20;    % base spring (kPa/µm); g multiplies
    params0_m.eta_M = 1.0;   % base dashpot; g multiplies
    params0_m.mods  = mods_phase3;

    if ~isfield(state,'g_opt3') || isempty(state.g_opt3)
        state.g_opt3 = [state.g_opt1, 1.0, 1.0];
    end

    fprintf('\n=== Phase 3: Maxwell dashpot comparison (max %d evals) ===\n', MAX_EVALS);
    g_new3 = fminsearch(@(g) slack_cost(g, params0_m, fn_target), state.g_opt3, opts);

    pm = params0_m; pm.g = g_new3;
    [~, ~, fm3, fd3] = runSlackExperiment(pm);
    c3 = evalFeatureCost(fd3, fm3, fn_target, 1);
    state.total_evals = state.total_evals + MAX_EVALS;
    state.g_opt3  = g_new3;
    state.costs3  = c3;

    fprintf('\nPhase 3 result: FeatTot=%.4f (evals=%d)\n', sum(c3), state.total_evals);
    fprintf('kM=%.3f  eta_M=%.3f\n', params0_m.kM*g_new3(end-1), params0_m.eta_M*g_new3(end));
    printCosts(fn_target, c3);

    state.phase = 4;
    save(STATE_FILE,'state');
end

%% ── Final summary ────────────────────────────────────────────────────────
if state.phase >= 3
    hdr = '%-20s  %8s  %8s  %8s';
    sep = repmat('-',1,54);
    has3 = isfield(state,'costs3') && ~isempty(state.costs3);
    if has3
        fprintf('\n%-20s  %8s  %8s  %8s  %8s\n','Feature','Baseline','NoCatch','ViscoSE','Maxwell');
        fprintf('%s\n', repmat('-',1,64));
    else
        fprintf(['\n' hdr '\n'],'Feature','Baseline','NoCatch','ViscoSE');
        fprintf('%s\n', sep);
    end
    for i = 1:numel(fn_target)
        c0 = getfield_safe(state.costs_base, i);
        c1 = getfield_safe(state.costs1,     i);
        c2 = getfield_safe(state.costs2,     i);
        if has3
            c3v = getfield_safe(state.costs3, i);
            fprintf('%-20s  %8.3f  %8.3f  %8.3f  %8.3f\n', fn_target{i}, c0,c1,c2,c3v);
        else
            fprintf('%-20s  %8.3f  %8.3f  %8.3f\n', fn_target{i}, c0,c1,c2);
        end
    end
    fprintf('%s\n', sep);
    if has3
        fprintf('%-20s  %8.3f  %8.3f  %8.3f  %8.3f\n','TOTAL',...
            sum(state.costs_base),sum(state.costs1),sum(state.costs2),sum(state.costs3));
    else
        fprintf('%-20s  %8.3f  %8.3f  %8.3f\n','TOTAL',...
            sum(state.costs_base),sum(state.costs1),sum(state.costs2));
    end
end

%% ── Local functions ──────────────────────────────────────────────────────
function E = slack_cost(g, params0, fn_target)
    p = params0; p.g = g;
    try
        [~, ~, fm, fd] = runSlackExperiment(p);
        E = sum(evalFeatureCost(fd, fm, fn_target, 1));
    catch
        E = 1e6;
    end
end

function printCosts(fn_list, costs)
    for i = 1:numel(fn_list)
        fprintf('  %-20s  %.3f\n', fn_list{i}, costs(i));
    end
    fprintf('  %-20s  %.3f\n','TOTAL', sum(costs));
end

function v = getfield_safe(arr, i)
    if numel(arr) >= i; v = arr(i); else; v = NaN; end
end
