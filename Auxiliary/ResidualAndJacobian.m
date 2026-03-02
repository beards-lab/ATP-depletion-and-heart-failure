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
    
    assert(length(params0.g) == length(params0.mods), "The params and g, woe");

    if isempty(call_count)
        call_count = 0;
        S_Cache = []; % Store the last calculated Jacobian
    end

    call_count = call_count + 1;
    
    % --- Configuration ---
    UPDATE_FREQUENCY = 5; % Only recalculate S_Cache every 5 iterations
    delta = 1e-3;       % Finite difference step size

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
        
        try
            % Perform Finite Differences for M parameters (M model calls)
            for j = 1:M
                G_perturbed = G_base;
                G_perturbed(j) = G_base(j) + delta;
                
                params0.g = G_perturbed;
                
                % IMPORTANT: Utilize your knowledge of known zero-dependencies here
                RunBakersExp; % This is your costly model call!
                [Residuals_pert, weights, cost] = evalFeatureCost(features_data, features_model, params0.fn, 1);
    
                S(:, j) = (Residuals - Residuals_pert) / delta;
            end
        catch e
            disp(e)
            % Residuals = [1000];
            Jacobian = S_Cache;
        end
        
        S_Cache = S;
    else
        % Use the cached Jacobian (NO model calls beyond the one for residuals)
        fprintf('  [OPTIMIZER] Reusing cached Jacobian (S) at call %d.\n', call_count);
    end
    
    % Return the cached/updated Jacobian to lsqnonlin
    Jacobian = S_Cache;
end