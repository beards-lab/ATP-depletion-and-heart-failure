% RunSensitivityAnalysis.m
% Calculates parameter sensitivities to identify the best candidates for optimization.
% Set ReRunJacobian = true to rerun the costly Jacobian computation, or
% false to load a previously computed result from 'ResidualsAndJacobian.mat'.

disp('--- Starting Sensitivity Analysis ---');
ReRunJacobian = false;

if ReRunJacobian
    % 1. Load the evaluated mechanism setup
    tic
    RunMechanismEvaluation;
    toc

    % Disable plotting for performance
    params0.PlotEachSeparately = 0;
    params0.justPlotStateTransitionsFlag = 0;
    params0.BreakOnODEUnstable = false;

    % 2. Define the exact features we are fitting against
    params0.fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', ...
                  'peak1_y', 'peak1_dSL', 'peak2', 'steady', ...
                  'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};

    % 3. Broad list of candidate parameters (~40)
    candidate_params = {
        'kd', 'ka', 'k1', 'k_1', 'k2', 'k_2', 'k3', 'k3m', 'kah', 'kamh', ...
        'kstiff1', 'kstiff2', 'kstiff3', 'k_pas', 'gamma', 'kSE', 'ekSE', 'mu', 'mu_neg', ...
        'dr', 'dr2', 'dr3', ...
        'ksr0', 'kmsr', 'sigma1', 'sigma2', 'ksrd', 'kmsrd', 'sigma_srd1', 'sigma_srd2', 'ksrd2sr', 'ksr2srd', ...
        'PieceWiseStrainDepParams__2', 'PieceWiseStrainDepParams__3', 'PieceWiseStrainDepParams__4', ...
        'PieceWiseStrainDepX__2', 'PieceWiseStrainDepX__3', 'PieceWiseStrainDepX__4', ...
        'PieceWiseStrainDepR1DParams__2', 'PieceWiseStrainDepR1DParams__3'
    };

    tp = tunableParams(params0);
    candidate_params = fieldnames(tp);

    params0.mods = candidate_params;
    params0.g    = ones(1, length(candidate_params));

    % 4. Compute Jacobian via finite differences
    disp('Computing Jacobian for candidate parameters (this may take a while)...');
    [Residuals, Jacobian] = ResidualAndJacobian(params0.g, params0, false);

    save('ResidualsAndJacobian.mat', 'Jacobian', 'Residuals', ...
         'candidate_params', 'params0', 'features_model', 'features_data');
    disp('Saved to ResidualsAndJacobian.mat');

else
    % FIX 1: actually unpack the struct fields
    disp('Loading precomputed Jacobian from ResidualsAndJacobian.mat ...');
    ResJac         = load('ResidualsAndJacobian.mat');
    Jacobian       = ResJac.Jacobian;
    Residuals      = ResJac.Residuals;         % FIX 2: was just bare 'Residuals'
    candidate_params = ResJac.candidate_params;
    params0        = ResJac.params0;
    features_model = ResJac.features_model;   % FIX 6: needed for analysis below
    features_data  = ResJac.features_data;
end

if isempty(Jacobian)
    error('Jacobian calculation failed or not loaded properly.');
end

%% --- Parameter nominal values (for elasticity scaling) ---
% FIX 4: use params0 not undefined 'params'
paramVals = zeros(1, length(candidate_params));
for i_canPan = 1:length(candidate_params)
    pname = candidate_params{i_canPan};
    if ~contains(pname, '__')
        paramVals(i_canPan) = params0.(pname);
    else
        arr = strsplit(pname, '__');
        if isfield(params0, arr{1}) && length(params0.(arr{1})) >= str2double(arr{2})
            paramVals(i_canPan) = params0.(arr{1})(str2double(arr{2}));
        else
            warning('Could not resolve candidate parameter "%s", defaulting to 1.', pname);
            paramVals(i_canPan) = 1;
        end
    end
end
paramVals(paramVals == 0) = eps; % guard against zero-valued params

%% --- Basic sensitivity: sum |J| per column ---
sensitivities = sum(abs(Jacobian), 1);
[sorted_sens, sort_idx] = sort(sensitivities, 'descend');
sorted_params = candidate_params(sort_idx);

disp('--- Raw Sensitivity Results ---');
for i_canPan = 1:length(sorted_params)
    fprintf('%2d. %-35s : %.6f\n', i_canPan, sorted_params{i_canPan}, sorted_sens(i_canPan));
end

top_K = 12;
fprintf('\nTop %d parameters to use in optimization:\n', top_K);
disp(sorted_params(1:top_K)');

%% --- Elasticity / Normalized Sensitivity Matrix ---
featureNames = cellfun(@(s) strtok(s, '|'), params0.fn, 'UniformOutput', false);
nFeatures    = length(featureNames);             % FIX 3: was length(features_model) -> 1
nParams      = length(paramVals);

% Model feature values at nominal point
features_model0Vals = cellfun(@(n) mean(abs(features_model.(n)) + eps), featureNames);

W_feat  = diag(1 ./ features_model0Vals);        % normalize rows by feature scale
% W_param removed: g multipliers are already dimensionless (all start at 1).
% J_{i,j} = dfeature_i/dg_j is already a fair cross-parameter comparison.
% Adding diag(|p_j0|) would incorrectly inflate large-magnitude params (e.g. kstiff2~13000 vs dr~0.01).
S_norm  = W_feat * Jacobian;                      % log-elasticity: (df/f) per unit multiplier change

%% --- SVD-based Identifiability Analysis ---
[U, Sigma, V] = svd(S_norm, 'econ');
singular_values = diag(Sigma);

% LOWERED threshold: 0.1% of max singular value (was 1%)
sv_threshold = 0.001 * max(singular_values);
k_sig = sum(singular_values > sv_threshold);

% Print all singular values with cumulative variance
total_variance = sum(singular_values.^2);
cum_var = cumsum(singular_values.^2) / total_variance * 100;
fprintf('\n--- Singular Values (threshold at 0.1%% of max = %.4f) ---\n', sv_threshold);
for sv_i = 1:length(singular_values)
    marker = '';
    if singular_values(sv_i) > sv_threshold; marker = ' <-- active'; end
    fprintf('  Mode %2d: %.4f   (cumul. var: %5.1f%%)%s\n', sv_i, singular_values(sv_i), cum_var(sv_i), marker);
end
fprintf('\nModel has %d parameters, %d effective degrees of freedom (SVD threshold = 1%%).\n', nParams, k_sig);

% Identifiability score per parameter: projection onto active subspace
Identifiability_Score = sum(V(:, 1:max(k_sig,1)).^2, 2);
[~, id_idx] = sort(Identifiability_Score, 'descend');
sorted_ident_params = candidate_params(id_idx);

fprintf('\nTop %d identifiable parameters (SVD-based):\n', top_K);
disp(sorted_ident_params(1:top_K)');

%% --- Missing Physics Detection (Orthogonal Projection) ---
ExpFeatures          = cellfun(@(n) mean(features_data.(n)), featureNames)';
features_model0Vals  = features_model0Vals';

Residual_raw       = ExpFeatures - features_model0Vals;
Residual_norm      = W_feat * Residual_raw;
Reachable_norm     = S_norm * (pinv(S_norm) * Residual_norm);
Missing_Physics_norm = Residual_norm - Reachable_norm;
Missing_Physics_Signal = diag(abs(features_model0Vals) + eps) * Missing_Physics_norm;

%% --- Visualizations ---
figure('Name', 'Singular Values');
bar(singular_values); set(gca, 'YScale', 'log');
xlabel('Mode'); ylabel('Singular value');
title('Information Content (Singular Values)');

figure('Name', 'Sensitivity Heatmap', 'Position', [100 100 1400 500]);
imagesc(S_norm);
colorbar; colormap(parula);
set(gca, 'XTick', 1:nParams, 'XTickLabel', candidate_params, 'TickLabelInterpreter', 'none', 'FontSize', 8);
set(gca, 'YTick', 1:nFeatures, 'YTickLabel', featureNames, 'TickLabelInterpreter', 'none');
xtickangle(45);
title('Normalized Feature–Parameter Sensitivity Matrix (S\_norm)');

figure('Name', 'SVD Analysis', 'Position', [100 100 1200 800]);

subplot(2,2,1);
bar(Identifiability_Score(id_idx));
set(gca, 'XTick', 1:nParams, 'XTickLabel', sorted_ident_params, 'TickLabelInterpreter', 'none', 'FontSize', 7);
xtickangle(45);
title('Parameter Identifiability Score'); ylabel('Score');

subplot(2,2,2);
scatter(V(:,1), V(:,2), 50, 'filled');
text(V(:,1)+0.02, V(:,2), candidate_params, 'Interpreter', 'none', 'FontSize', 7);
title('Parameter Collinearity (Mode 1 vs 2)');
xlabel('Mode 1'); ylabel('Mode 2'); grid on;

num_modes_to_plot = min(2, k_sig);
for m = 1:num_modes_to_plot
    subplot(2, 2, 2+m);
    barh(U(:, m));
    set(gca, 'YTick', 1:nFeatures, 'YTickLabel', featureNames, 'TickLabelInterpreter', 'none', 'FontSize', 8);
    title(sprintf('Feature Archetype %d', m));
    xlabel('Relative Change');
end

figure('Name', 'Missing Physics');
barh(Missing_Physics_Signal);
set(gca, 'YTick', 1:nFeatures, 'YTickLabel', featureNames, 'TickLabelInterpreter', 'none');
title('Structural Error (Cannot be fixed by tuning)');
xlabel('Unreachable Discrepancy Magnitude');

%% --- Save ---
save('SensitivityAnalysisResults.mat', 'candidate_params', 'Jacobian', 'Residuals', ...
     'sensitivities', 'sorted_params', 'sorted_sens', 'S_norm', ...
     'Identifiability_Score', 'sorted_ident_params', ...
     'Missing_Physics_Signal', 'singular_values');
disp('Results saved to SensitivityAnalysisResults.mat');
