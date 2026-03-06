% ResidualAndJacobian.m
% Returns the vector of feature costs (Residuals) and the Jacobian (Sensitivity S).
function [Residuals, Jacobian] = ResidualAndJacobian(g, params0, ignoreJac)
    % Use 'persistent' variables to store values between function calls
    persistent S_Cache call_count history
    
    if nargin < 3 
        ignoreJac = false;
    end

    params0.g = g;
    params0.EvalFeatures = true;
    params0.PlotEachSeparately = false;
    params0.BreakOnODEUnstable = true;
    params0.MaxRunTime = 60;
    params0.FV_velocities = -[0.5, 1, 3, 4];
    
    assert(length(params0.g) == length(params0.mods), "The params and g, woe");

    if isempty(call_count)
        call_count = 0;
        S_Cache = []; % Store the last calculated Jacobian
    end

    call_count = call_count + 1;
    
    % --- Configuration ---
    UPDATE_FREQUENCY = 5; % Only recalculate S_Cache every 5 iterations
    % IMPORTANT: delta must be >> sqrt(ODE_RelTol) to be above numerical noise.
    % ODE RelTol = 1e-2, so noise floor is ~1e-2. delta=1e-3 was BELOW noise (bug!).
    % Using 5e-2 (5% parameter change) safely above the noise floor.
    delta = 5e-2;

    % --- Part 1: COSTLY MODEL RUN (Always required for residuals) ---
    try
        LoadData;
        RunBakersExp; % This is your costly model call!
        [Residuals, weights, cost] = evalFeatureCost(features_data, features_model, params0.fn, 1);
    catch e
        disp(e)
        Residuals = 1e3*ones(size(params0.fn));
        Jacobian = S_Cache;
    end

    if isempty(history)
        history = Residuals;
    else
        history = [history; Residuals];
    end
    assignin('base', 'optim_history', history); % Send to workspace in real-time    

    if ignoreJac
        Jacobian = [];
        return;
    end

    
    % --- Part 2: CONDITIONAL JACOBIAN CALCULATION ---

    % Check if it's time to perform the expensive Jacobian calculation
    if isempty(S_Cache) || mod(call_count, UPDATE_FREQUENCY) == 1
        fprintf('  [OPTIMIZER] Recalculating costly Jacobian (S) at call %d...\n', call_count);
        M = length(g);             % Number of parameters
        N = length(Residuals);             % Number of features

        S = zeros(N, M);
        G_base = g;

        featureNames = cellfun(@(s) strtok(s, '|'), params0.fn, 'UniformOutput', false);
        % Use mean() not sum() — features are vectors (one per slack velocity)
        % sum() inflates magnitude by nVelocities; mean() gives consistent scale
        features_model0Vals = cellfun(@(n) mean(features_model.(n)), featureNames);

        % Perform Finite Differences for M parameters (M model calls)
        for i_jac = 1:M
            G_perturbed = G_base;
            G_perturbed(i_jac) = G_base(i_jac) + delta;

            params0.g = G_perturbed;

            % IMPORTANT: Utilize your knowledge of known zero-dependencies here
            try
                disp(i_jac)
                RunBakersExp; % This is your costly model call!
                % [Residuals_pert, weights, cost] = evalFeatureCost(features_data, features_model, params0.fn, 1);
                features_modelVals = cellfun(@(n) mean(features_model.(n)), featureNames);
                S(:, i_jac) = (features_modelVals - features_model0Vals) / delta;
            catch e
                disp(e)
                % Residuals = [1000];
                % Jacobian = S_Cache;
            end
        end
        
        S_Cache = S;
    else
        % Use the cached Jacobian (NO model calls beyond the one for residuals)
        fprintf('  [OPTIMIZER] Reusing cached Jacobian (S) at call %d.\n', call_count);
    end
    
    % Return the cached/updated Jacobian to lsqnonlin
    Jacobian = S_Cache;
end