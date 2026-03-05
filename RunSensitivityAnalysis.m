% RunSensitivityAnalysis.m
% Calculates parameter sensitivities to identify the best candidates for optimization.

disp('--- Starting Sensitivity Analysis ---');
ReRunJacobian = false;
if ReRunJacobian
    
    % 1. Load the evaluated mechanism setup
    RunMechanismEvaluation; 
    
    % Disable plotting for performance
    params0.PlotEachSeparately = 0;
    params0.justPlotStateTransitionsFlag = 0;
    
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
    
    params0.mods = candidate_params;
    params0.g = ones(1, length(params0.mods));
    
    % 4. Call ResidualAndJacobian to compute Jacobian using finite differences
    disp('Computing Jacobian for candidate parameters (this may take a while)...');
    % We pass 'false' to compute the Jacobian (S_Cache will be updated)
    [Residuals, Jacobian] = ResidualAndJacobian(params0.g, params0, false);
    
    save('ResidualsAndJacobian', "Jacobian", "Residuals", "candidate_params", "params0");
else
    ResJac = load('ResidualsAndJacobian.mat');
    Residuals
end

if isempty(Jacobian)
    error('Jacobian calculation failed.');
end

%% 1. Setup and Matrix Calculation
% Assume: 
% - params0: current parameter vector (e.g., from parameter-reference.md)
% - ExpFeatures: 1x12 vector of experimental observations
% - EvaluateModel: function returning 1x12 vector of simulated features
% params = params0;
% p = getParams(params0, ones(length(candidate_params), 1), false, true);
%% 2. Reconstruct arrays: params.arr_2 = 3 -> params.arr(2) = 3
paramVals = zeros(size(candidate_params));
for i_canPan = 1:length(candidate_params)
    if ~contains(candidate_params(i_canPan), '__')
        paramVals(i_canPan) = params0.(candidate_params{i_canPan});
        continue;
    end
    arr = split(candidate_params(i_canPan), '__');
    if isfield(params, arr(1)) ...
        && length(params.(arr{1})) > 1 ... % the array exists
        && ~isempty(str2num(arr{2})) % the index is numeric
            % apply value to an array
            paramVals(i_canPan) = params0.(arr{1})(str2num(arr{2}));
    else 
        disp(['Err...' num2str(i_canPan)]);
    end
end

% candidate_params = {'kah', 'kamh', 'ka', 'kd', 'k1', 'ksr0', 'L_thin', 'gamma_titin'}; % Add all
nParams = length(paramVals);
nFeatures = length(features_model);

% Jacobian = zeros(nFeatures, nParams);
% delta_rel = 1e-4; % 0.01% perturbation

%% 2. Normalization (Elasticity Matrix)
% Transform into % change in feature / % change in parameter
% [Residuals] = evalFeatureCost(features_data, features_model, params0.fn, 1);
featureNames = cellfun(@(s) strtok(s, '|'), params0.fn, 'UniformOutput', false);
features_model0Vals = cellfun(@(n) sum(features_model.(n)), featureNames);

W_feat = diag(1 ./ abs(features_model0Vals + eps));
W_param = diag(paramVals);
S_norm = W_feat * Jacobian * W_param; 

%% 3. Identifiability Analysis (SVD)
[U, Sigma, V] = svd(S_norm, 'econ');
singular_values = diag(Sigma);

figure;
subplot(2,1,1);
bar(singular_values); set(gca, 'YScale', 'log');
title('Information Content (Singular Values)');
ylabel('Strength of Influence');

% Identify the "Null Space" (Parameters that do nothing)
% Parameters with low values in the first few Right Singular Vectors (V) 
% are unidentifiable.

%% 3b. Analyze Parameter Identifiability and Collinearity (V Matrix)
% Determine how many "effective" degrees of freedom the model has
% Let's use a threshold of 1% of the maximum singular value
sv_threshold = 0.01 * max(singular_values);
k_sig = sum(singular_values > sv_threshold);

fprintf('\nThe model has %d parameters, but only %d effective degrees of freedom.\n', nParams, k_sig);

% Calculate Identifiability Score for each parameter
% A parameter is identifiable if it has strong weights in the top k_sig vectors of V
Identifiability_Score = sum(V(:, 1:k_sig).^2, 2);

% Sort by Identifiability
[sorted_ident, id_idx] = sort(Identifiability_Score, 'descend');
sorted_ident_params = candidate_params(id_idx);

figure('Name', 'SVD Insights', 'Position', [100, 100, 1200, 800]);

% Plot Identifiability
subplot(2,2,1);
bar(sorted_ident);
set(gca, 'XTick', 1:nParams, 'XTickLabel', sorted_ident_params, 'TickLabelInterpreter', 'none');
xtickangle(45);
title('Parameter Identifiability Score (from V)');
ylabel('Influence in Active Subspace');
% Parameters with scores near 0 are "dead" or entirely redundant.

% Plot Parameter Collinearity (First 2 Principal Modes)
subplot(2,2,2);
scatter(V(:,1), V(:,2), 'filled');
text(V(:,1)+0.02, V(:,2), candidate_params, 'Interpreter', 'none', 'FontSize', 8);
title('Parameter Collinearity (Modes 1 vs 2)');
xlabel('Principal Mode 1 Weight');
ylabel('Principal Mode 2 Weight');
grid on;
% Parameters clustered closely together are "twins" - the model cannot distinguish them.

%% 3c. Analyze Feature Archetypes (U Matrix)
% The columns of U define the fundamental "shapes" of behavior the model can make.
% We will plot the top 3 most dominant behaviors.

num_modes_to_plot = min(3, k_sig);
for m = 1:num_modes_to_plot
    subplot(2, 2, 2+m);
    barh(U(:, m));
    set(gca, 'YTick', 1:nFeatures, 'YTickLabel', featureNames, 'TickLabelInterpreter', 'none');
    title(sprintf('Feature Archetype %d (from U)', m));
    xlabel('Relative Change Direction');
end

%% 4. Missing Physics Detection (Orthogonal Projection)
% Residual: The gap between model and experiment
ExpFeatures = cellfun(@(n) sum(features_data.(n)), featureNames)';
features_model0Vals = features_model0Vals'; % ensure column vector

% 1. Calculate raw residual
Residual_raw = ExpFeatures - features_model0Vals; 

% 2. NORMALIZE the residual to match S_norm's percentage space
Residual_norm = W_feat * Residual_raw;

% 3. Project normalized residual onto Reachable Space
Reachable_norm = S_norm * (pinv(S_norm) * Residual_norm);

% 4. Calculate Missing Physics in normalized space
Missing_Physics_norm = Residual_norm - Reachable_norm;

% 5. (Optional) Convert back to physical units for interpretability
Missing_Physics_Signal = diag(abs(features_model0Vals + eps)) * Missing_Physics_norm;

%% 5. Visualization
figure;
subplot(1,2,1);
imagesc(S_norm); colorbar;
set(gca, 'XTick', 1:nParams, 'XTickLabel', candidate_params, 'TickLabelInterpreter', 'none');
title('Feature-Parameter Sensitivity Map');

subplot(1,2,2);
barh(Missing_Physics_Signal);
set(gca, 'YTick', 1:nFeatures, 'YTickLabel', featureNames);
title('Structural Error (Cannot be fixed by tuning)');
xlabel('Magnitude of Unreachable Discrepancy');

 
% 5. Analyze Sensitivities
% Sum of absolute Jacobian values across all features for each parameter log-multiplier
sensitivities = sum(abs(Jacobian), 1);

% Sort by sensitivity
[sorted_sens, sort_idx] = sort(sensitivities, 'descend');
sorted_params = candidate_params(sort_idx);

disp('--- Sensitivity Analysis Results ---');
for i_canPan = 1:length(sorted_params)
    fprintf('%2d. %-35s : %f\n', i_canPan, sorted_params{i_canPan}, sorted_sens(i_canPan));
end

% Print the top 12 for the optimization script
top_K = 12;
fprintf('\nTop %d parameters to use in optimization:\n', top_K);
disp(sorted_params(1:top_K)');

% Save results
save('SensitivityAnalysisResults.mat', 'candidate_params', 'Jacobian', 'Residuals', 'sensitivities', 'sorted_params', 'sorted_sens');
disp('Results saved to SensitivityAnalysisResults.mat');
