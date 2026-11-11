% OptimLastSlackPieceWise.m
%
% Block-coordinate-descent optimizer for the last-slack restretch experiment.
% Starting point: params/ModelOptParamsLastSlack_PieceWiseFunc.m
%
% Strategy
% --------
%   Parameters are divided into functional groups. Each iteration:
%     1. Shuffle the group order.
%     2. Run fminsearch on each group (~50 evals).
%     3. Accept the new g_best if total cost improved.
%   Repeat until no group produces improvement > CONV_TOL in a full cycle.
%
% Key mechanism under test: UseLatticeSpacing — provides SL-dependent
% attachment (f_lattice) via inter-filament spacing. This is the primary
% candidate for fixing the "A|SLslack insensitive to SL" problem.
%
% Run interactively in MATLAB (not -batch) so you see live fminsearch output.

STATE_FILE    = 'OptimLastSlackPieceWise_state.mat';
STATE_FILE    = 'OptimLastSlackPieceWise_stateTimeCourse.mat';
STATE_FILE    = 'OptimLastSlackPieceWise_stateTimeCourseAll.mat';
EVALS_PER_BLOCK = 50;   % fminsearch evals per group per iteration
MAX_CYCLES    = 10;     % max full-group cycles before stopping
CONV_TOL      = 0.02;   % min cost improvement to continue a cycle

%% ── Feature targets ───────────────────────────────────────────────────────
fn_target = {'ktr|SLslack','A|SLslack','t0|SLslack', ...
             'peak1_y','peak1_dSL','peak2','steady', ...
             'XTOR|0.1','vall_y','restretchSlopeStart','vall2_dy'};

fn_target = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart'};

%% ── Parameter groups ──────────────────────────────────────────────────────
% Each group is a cell-array of param names (used via params0.mods / g).
% Use 'FieldName__k' for the k-th element of vector field FieldName.
%
% Knot positions (X vectors) are NOT optimized — only knot VALUES (Params).
% This keeps the search space lower-dimensional and physically interpretable.

groups = { ...
    % A: SE / mechanical — controls peak1 height and ktr time constant
    {'kSE', 'mu', 'c_SE_visc'}, ...

    % B: kinetic rates — controls ktr rate and steady-state bridge fraction
    {'ka', 'kd', 'k1', 'k2'}, ...

    % C: force scale + passive — steady force level, passive contribution
    {'kstiff2', 'k_pas', 'slope'}, ...

    % D: lattice spacing — SL-dependent activation (primary SL-sensitivity fix)
    %    d10_ref  : reference lattice spacing at SL_ref_lattice
    %    sigma_lattice : Gaussian width of attachment penalty vs excess spacing
    {'d10_ref', 'sigma_lattice', 'd_optimal'}, ...

    % E: R2 detachment knots in positive-strain (restretch) region
    %    Knots 5-9 of PieceWiseStrainDep2Params at strains [−0.004, 0, 0.006, 0.01, 0.0233]
    %    Controls the valley depth and peak2 amplitude
    {'PieceWiseStrainDep2Params__5', 'PieceWiseStrainDep2Params__6', ...
     'PieceWiseStrainDep2Params__7', 'PieceWiseStrainDep2Params__8', ...
     'PieceWiseStrainDep2Params__9'}, ...

    % F: R1D (p1 detachment) knots at positive strain
    %    Knots 5-7 of PieceWiseStrainDepR1DParams at strains [0.004, 0.01, 0.1]
    {'PieceWiseStrainDepR1DParams__5', 'PieceWiseStrainDepR1DParams__6', ...
     'PieceWiseStrainDepR1DParams__7'}, ...

    % G: R12 (p1→p2 / attachment) knots — all 5 knots of PieceWiseStrainDepParams
    %    Controls where and how fast bridges power-stroke
    {'PieceWiseStrainDepParams__1', 'PieceWiseStrainDepParams__2', ...
     'PieceWiseStrainDepParams__3', 'PieceWiseStrainDepParams__4', ...
     'PieceWiseStrainDepParams__5'}, ...
};

n_groups = numel(groups);

%% ── Shared model setup ────────────────────────────────────────────────────
clear params0;
params0 = getParams();
% ModelOptParamsLastSlack_PieceWiseFunc;   % loads fully-applied params (mods={})
% ModelOptParams_PerPartesAdAstra
ModelOptParams_SlackFeaturesMultiOpt
LoadData;

% Override for optimization run
params0.RunSlackSegments    = 'AllPar';
params0.RunForceVelocity    = true;

params0.RunForceLengthEstim = false;
params0.PlotEachSeparately  = false;
params0.ShowResidualPlots   = false;
params0.PlotFeatureFitting  = false;
params0.ghostLoad = ''; params0.ghostSave = '';
params0.EvalFeatures        = true;
params0.BreakOnODEUnstable  = false;
params0.fn                  = fn_target;
params0.BreakOnODEUnstable = true;


% Enable lattice spacing (key SL-dependence mechanism)
params0.UseLatticeSpacing   = true;
% d10_ref=0.037 µm and sigma_lattice=0.0006 µm are starting values from params file.
% Group D will tune these.
% RunBakersExp

params0.RunForceVelocity = false;
% params0.RunSlackSegments = 'Fourth';
params0.EvalFeatures = false;
%% ── State initialisation ──────────────────────────────────────────────────
if exist(STATE_FILE, 'file')
    load(STATE_FILE, 'state');
    fprintf('Resumed: cycle=%d, best_cost=%.4f\n', state.cycle, state.best_cost);
else
    % Build initial g=1 for every param across all groups
    all_params = [groups{:}];
    g0 = ones(1, numel(all_params));
    state.all_params  = all_params;
    state.g_best      = g0;
    state.best_cost   = Inf;
    state.cycle       = 0;
    state.costs_hist  = [];   % [cycle x n_features] history
    state.cycle_costs = [];   % [cycle x n_groups] per-group best cost
end

%% ── Evaluate baseline ─────────────────────────────────────────────────────
% if state.best_cost == Inf
    p = params0;
    params0.mods = state.all_params;
    params0.g    = state.g_best;
    % params0.RunForceVelocity = false;
    % params0.RunSlackSegments = 'Fourth';
    % params0.RunSlackSegments = 'AllPar';
    % [~, ~, features_model, features_data] = runSlackExperiment(params0);
    params0.PlotEachSeparately = true;
    params0.EvalFeatures = true;

    % params0.ghostSave = 'Fourth-preOp';
    params0.ghostLoad = 'Fourth-preOp';

    RunBakersExp;
    c0 = evalFeatureCost(features_data, features_model, fn_target);
    plotFeatures(features_data, features_model, [], fn_target);
    % state.best_cost  = sum(c0);
    % state.costs_hist = c0;

    % state.best_cost  = E(1);
    % state.costs_hist = E(1);
    E(1)
    
    fprintf('Baseline cost: %.4f\n', state.best_cost);
    printFeatures(fn_target, c0);
    params0 = p;
    % save(STATE_FILE, 'state');
% end

%% ── Main optimisation loop ────────────────────────────────────────────────
opts = optimset('Display','iter','TolFun',1e-4, ...
                'MaxIter', EVALS_PER_BLOCK, ...
                'MaxFunEvals', EVALS_PER_BLOCK,'PlotFcns','optimplotfval');

for cycle = state.cycle + 1 : MAX_CYCLES

    fprintf('\n====== Cycle %d / %d  (current best = %.4f) ======\n', ...
        cycle, MAX_CYCLES, state.best_cost);

    group_order = randperm(n_groups);   % shuffle group order each cycle
    cycle_improved = false;
    cycle_costs = zeros(1, n_groups);

    for gi = 1 : n_groups
        g_idx = group_order(gi);
        grp   = groups{g_idx};
        fprintf('\n-- Group %d (%s...) --\n', g_idx, grp{1});

        % Build mods/g slice for this group
        [grp_g0, grp_idx] = extractGroupG(state.g_best, state.all_params, grp);

        % Cost function: only the current group's g values move
        costfun = @(g_grp) block_cost(g_grp, grp_idx, state.g_best, ...
                                       state.all_params, params0, fn_target);

        g_grp_new = fminsearch(costfun, grp_g0, opts);

        % Evaluate result
        g_candidate = state.g_best;
        g_candidate(grp_idx) = g_grp_new;
        params0.mods = state.all_params;
        params0.g    = g_candidate;
        try
            params0.PlotEachSeparately = false;
            RunBakersExp;
            c = evalFeatureCost(features_data, features_model, fn_target);                        
            % [~, ~, fm, fd] = runSlackExperiment(params0);
            % c = evalFeatureCost(fd, fm, fn_target, 1);
            % new_cost = sum(c);
            new_cost = E(1);
        catch
            new_cost = state.best_cost;  % reject on ODE failure
            % c = state.costs_hist(end,:);
            E = state.costs_hist(end,:);
        end

        cycle_costs(gi) = new_cost;

        if new_cost < state.best_cost - 1e-6
            fprintf('  Improved: %.4f → %.4f\n', state.best_cost, new_cost);
            state.g_best   = g_candidate;
            state.best_cost = new_cost;
            cycle_improved  = true;
        else
            fprintf('  No improvement (%.4f >= %.4f)\n', new_cost, state.best_cost);
        end
    end

    % Record history
    state.cycle       = cycle;
    state.cycle_costs = [state.cycle_costs; cycle_costs];

    % Final eval for feature breakdown
    params0.mods = state.all_params;
    params0.g    = state.g_best;
    [~, ~, fmc, fdc] = runSlackExperiment(params0);
    cc = evalFeatureCost(fdc, fmc, fn_target, 1);
    % state.costs_hist = [state.costs_hist; cc];

    fprintf('\nEnd of cycle %d — best cost = %.4f\n', cycle, state.best_cost);
    printFeatures(fn_target, cc);
    printBestG(state.all_params, state.g_best, groups);

    save(STATE_FILE, 'state');

    if ~cycle_improved
        fprintf('No improvement in cycle %d — converged.\n', cycle);
        break;
    end
end

%% ── Final summary ─────────────────────────────────────────────────────────
fprintf('\n====== FINAL RESULT ======\n');
fprintf('Total cycles: %d,  Best cost: %.4f\n', state.cycle, state.best_cost);
printFeatures(fn_target, state.costs_hist(end,:));
fprintf('\nOptimised g values:\n');
printBestG(state.all_params, state.g_best, groups);

fprintf('\nWrite to params file with writeParamsToMFile.\n');
% writeParamsToMFile('ModelOptParams_SlackFeaturesMultiOpt', params0)
% writeParamsToMFile('ModelOptParams_SlackTimebaseMultiOpt', params0)
%% ── Local functions ───────────────────────────────────────────────────────
function E = block_cost(g_grp, grp_idx, g_full, all_params, params0, fn_target)
    g = g_full;
    g(grp_idx) = g_grp;
    params0.mods = all_params;
    params0.g    = g;
    try
        % [~, ~, fm, fd] = runSlackExperiment(p);
        LoadData;
        params0.PlotEachSeparately = false;
        RunBakersExp;
        % c0 = evalFeatureCost(features_data, features_model, fn_target);
        % E = sum(c0);
        E = E(1);
    catch ee
        E = 1e6;
    end
end

function [g_grp, grp_idx] = extractGroupG(g_full, all_params, grp)
    grp_idx = zeros(1, numel(grp));
    for k = 1:numel(grp)
        idx = find(strcmp(all_params, grp{k}), 1);
        if isempty(idx)
            error('Param "%s" not found in all_params.', grp{k});
        end
        grp_idx(k) = idx;
    end
    g_grp = g_full(grp_idx);
end

function printFeatures(fn_list, costs)
    fprintf('  %-26s  Cost\n', 'Feature');
    fprintf('  %s\n', repmat('-', 1, 36));
    for i = 1:numel(fn_list)
        fprintf('  %-26s  %.4f\n', fn_list{i}, costs(i));
    end
    fprintf('  %-26s  %.4f\n', 'TOTAL', sum(costs));
end

function printBestG(all_params, g_best, groups)
    for gi = 1:numel(groups)
        fprintf('  Group %d:\n', gi);
        for k = 1:numel(groups{gi})
            nm  = groups{gi}{k};
            idx = find(strcmp(all_params, nm), 1);
            fprintf('    %-40s  g=%.4f\n', nm, g_best(idx));
        end
    end
end
