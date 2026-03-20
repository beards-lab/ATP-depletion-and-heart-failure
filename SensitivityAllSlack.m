% SensitivityAllSlack.m
%
% OAT ±10% sensitivity: which params drive A|SLslack?
% Loads current best from OptimLastSlackPieceWise_state.mat.
% Tests all current optimised params + extended candidate list.
% Outputs: sorted keep/drop table + ready-to-paste group definitions.
%
% UseSuperRelaxed / UsePassiveForSR intentionally NOT enabled here —
% sensitivity runs with current mechanism set as-is.
%
% Results saved to SensitivityAllSlack_results.mat including full
% sens_c tensor (n_params × 2 × n_features) for Jacobian analysis.

STATE_FILE  = 'OptimLastSlackPieceWise_state.mat';
DELTA       = 0.10;
TARGET_FEAT = 'A|SLslack';
DROP_THRESH = 0.05;

%% ── Setup ────────────────────────────────────────────────────────────────
load(STATE_FILE, 'state');
clear params0; params0 = getParams();
ModelOptParamsLastSlack_PieceWiseFunc;
LoadData;

params0.RunSlackSegments    = 'AllPar';
params0.RunForceVelocity    = true;
params0.PlotEachSeparately  = false;
params0.EvalFeatures        = true;
params0.BreakOnODEUnstable  = true;
params0.UseLatticeSpacing = true;
% NOTE: UseSuperRelaxed / UsePassiveForSR intentionally NOT set here.
% Sensitivity runs with current state as-is to establish baseline
% sensitivities before enabling new mechanisms.
params0.mods = state.all_params;
params0.g    = state.g_best;

fn_target = {'FV_f|FV_v','ktr|SLslack','A|SLslack','t0|SLslack', ...
             'peak1_y','peak1_dSL','peak2','steady','XTOR|0.1', ...
             'vall_y','restretchSlopeStart'};

params0.fn = fn_target;

%% ── Baseline ─────────────────────────────────────────────────────────────
fprintf('Computing baseline...\n');
tic
RunBakersExp;
toc
c0       = evalFeatureCost(features_data, features_model, fn_target);
E0_total = sum(c0);
i_tgt    = find(strcmp(fn_target, TARGET_FEAT));
E0_tgt   = c0(i_tgt);
fprintf('Baseline: total=%.4f  %s=%.4f\n', E0_total, TARGET_FEAT, E0_tgt);

%% ── Candidate list (current params + extended) ───────────────────────────
extra_params = {'sigma1','sigma2','kmsr','ksr0', ...          % SRX regulation
                'k_pas','slope', ...                          % passive force
                'L_thin','L_thick','L_hbare', ...            % overlap geometry
                'd10_ref','sigma_lattice','d_optimal', ...   % lattice spacing
                'kstiff2','kstiff1', ...                     % XB force scale
                'ka','kd','k1','k2'};                        % kinetics
all_test = unique([state.all_params(:)', extra_params], 'stable');
fprintf('Total candidate params: %d (%d current + %d extra, after dedup)\n', ...
    numel(all_test), numel(state.all_params), numel(extra_params));

%% ── Detect zero-value params (cannot optimise by multiplication) ──────────
p_resolved = getParams(params0, state.g_best, false, true);
zero_params  = {};
nonzero_test = {};
for k = 1:numel(all_test)
    nm      = all_test{k};
    base_nm = regexprep(nm, '__\d+$', '');
    if isfield(p_resolved, base_nm)
        v   = p_resolved.(base_nm);
        tok = regexp(nm, '__(\d+)$', 'tokens');
        if ~isempty(tok)
            idx_arr = str2double(tok{1}{1});
            if idx_arr <= numel(v)
                v = v(idx_arr);
            else
                v = 0;
            end
        else
            v = v(1);   % take first element if array
        end
        if v == 0
            zero_params{end+1} = nm; %#ok<SAGROW>
            continue;
        end
    end
    nonzero_test{end+1} = nm; %#ok<SAGROW>
end
fprintf('\nZero-value params (excluded from sweep — set non-zero base first):\n');
for k = 1:numel(zero_params)
    fprintf('  %s\n', zero_params{k});
end
fprintf('Params to sweep: %d\n\n', numel(nonzero_test));

%% ── OAT sweep ────────────────────────────────────────────────────────────
n  = numel(nonzero_test);
nf = numel(fn_target);
% Use NaN as sentinel for failed evaluations (ODE unstable / error).
% NaN is excluded from impact calculation via omitnan — prevents false
% positives from error-inflated costs (e.g. 1e6 sentinel would make
% every broken param look maximally sensitive).
sens_total = nan(n, 2);          % total cost for each param × direction
sens_tgt   = nan(n, 2);          % target feature cost per param × direction
sens_c     = nan(n, 2, nf);      % full cost vector for Jacobian analysis
had_error  = false(n, 2);        % tracks which evals failed

for k = 1:n
    nm = nonzero_test{k};
    for s = [-1 1]
        g_test   = state.g_best;
        mods_use = state.all_params;

        % Find param in current mods or append as new
        idx = find(strcmp(state.all_params, nm));
        if ~isempty(idx)
            g_test(idx) = g_test(idx) * (1 + s*DELTA);
        else
            mods_use = [state.all_params(:)', {nm}];
            g_test   = [state.g_best,          1 + s*DELTA];
        end

        params0.mods = mods_use;
        params0.g    = g_test;
        col = (s > 0) + 1;   % col 1 = -10%, col 2 = +10%

        try
            RunBakersExp;
            c = evalFeatureCost(features_data, features_model, fn_target);
            sens_total(k, col)  = sum(c);
            sens_tgt(k, col)    = c(i_tgt);
            sens_c(k, col, :)   = c;
        catch ee
            % Leave as NaN — do NOT use 1e6, which would create false
            % sensitivity (|1e6 - E0| >> real cost differences).
            had_error(k, col) = true;
            fprintf('    ERROR [%d/%d %+d%%] %s: %s\n', k, n, s*100, nm, ee.message);
        end

        % Restore baseline mods
        params0.mods = state.all_params;
        params0.g    = state.g_best;
    end

    % Report: show NaN as '  err' for failed directions
    fmt_val = @(v) sprintf('%+.3f', v);
    v_neg = sens_tgt(k,1) - E0_tgt; v_pos = sens_tgt(k,2) - E0_tgt;
    s_neg = fmt_val(v_neg); if had_error(k,1); s_neg = '   err'; end
    s_pos = fmt_val(v_pos); if had_error(k,2); s_pos = '   err'; end
    fprintf('  [%d/%d] %-35s  %s(-10%%)=%s  %s(+10%%)=%s\n', ...
        k, n, nm, TARGET_FEAT, s_neg, TARGET_FEAT, s_pos);
end

%% ── Report ───────────────────────────────────────────────────────────────
% omitnan: if one direction errored but the other succeeded, use that one.
% If both errored, impact = NaN → listed in error section, not keep/drop.
impact = max(abs(sens_tgt - E0_tgt), [], 2, 'omitnan');
error_both = all(had_error, 2);   % params where both directions failed

[~, ord] = sort(impact, 'descend', 'MissingPlacement', 'last');

fprintf('\n%-35s  %8s  %8s  %8s  %8s  %8s\n', ...
    'Parameter','ΔTot(-)','ΔTot(+)','ΔTgt(-)','ΔTgt(+)','MaxImpact');
fprintf('%s\n', repmat('-', 1, 83));
keep_list  = {};
drop_list  = {};
error_list = {};
for k = 1:n
    i = ord(k);
    if error_both(i)
        error_list{end+1} = nonzero_test{i}; %#ok<SAGROW>
        tag = ' !! (both dirs errored)';
        fprintf('%-35s  %8s  %8s  %8s  %8s  %8s%s\n', nonzero_test{i}, ...
            'err','err','err','err','NaN', tag);
        continue;
    end
    % Format NaN entries (one direction errored) as 'err'
    fmt_tbl = @(v) sprintf('%8.3f', v);
    tot_neg = fmt_tbl(sens_total(i,1)-E0_total); if had_error(i,1); tot_neg = '     err'; end
    tot_pos = fmt_tbl(sens_total(i,2)-E0_total); if had_error(i,2); tot_pos = '     err'; end
    tgt_neg = fmt_tbl(sens_tgt(i,1)-E0_tgt);     if had_error(i,1); tgt_neg = '     err'; end
    tgt_pos = fmt_tbl(sens_tgt(i,2)-E0_tgt);     if had_error(i,2); tgt_pos = '     err'; end
    if impact(i) >= DROP_THRESH
        keep_list{end+1} = nonzero_test{i}; %#ok<SAGROW>
        tag = ' ✓';
    else
        drop_list{end+1} = nonzero_test{i}; %#ok<SAGROW>
        tag = ' ✗';
    end
    fprintf('%-35s  %s  %s  %s  %s  %8.3f%s\n', nonzero_test{i}, ...
        tot_neg, tot_pos, tgt_neg, tgt_pos, impact(i), tag);
end
fprintf('%s\n', repmat('-', 1, 83));
fprintf('(✓ keep | ✗ drop | !! both-dir error | threshold = %.2f)\n', DROP_THRESH);
if ~isempty(error_list)
    fprintf('Both-direction errors (exclude from groups): ');
    fprintf('%s  ', error_list{:}); fprintf('\n');
end
fprintf('\nZero-value (excluded): ');
fprintf('%s  ', zero_params{:});
fprintf('\n');

%% ── Ready-to-paste group definitions (≤6 params per group) ──────────────
MAX_PER_GROUP = 6;
fprintf('\n%% ── Paste into OptimLastSlackPieceWise.m (replace existing Group D onward) ──\n');
chunks = {};
for k = 1:MAX_PER_GROUP:numel(keep_list)
    chunks{end+1} = keep_list(k : min(k+MAX_PER_GROUP-1, numel(keep_list))); %#ok<SAGROW>
end
g_start = 4;   % first group index to replace (Group D)
for ci = 1:numel(chunks)
    fprintf('groups{%d} = {', g_start + ci - 1);
    chunk = chunks{ci};
    for j = 1:numel(chunk)
        if j < numel(chunk); sep = ', '; else; sep = ''; end
        fprintf('''%s''%s', chunk{j}, sep);
    end
    fprintf('};\n');
end

%% ── Save results ─────────────────────────────────────────────────────────
% sens_c layout: (n_params × 2 directions × n_features)
%   sens_c(:,1,:) = cost vectors for -10% perturbation
%   sens_c(:,2,:) = cost vectors for +10% perturbation
%   Jacobian approx: J(k,f) = (sens_c(k,2,f) - sens_c(k,1,f)) / (2*DELTA)
save('SensitivityAllSlack_results.mat', ...
    'nonzero_test', 'zero_params', 'fn_target', 'c0', ...
    'sens_total', 'sens_tgt', 'sens_c', 'had_error', ...
    'E0_total', 'E0_tgt', 'keep_list', 'drop_list', 'error_list', ...
    'impact', 'ord', 'DELTA', 'TARGET_FEAT', 'DROP_THRESH');
fprintf('\nSaved SensitivityAllSlack_results.mat\n');
