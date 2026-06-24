% ResidualAndJacobian.m
% Returns the vector of feature costs (Residuals) and the Jacobian (Sensitivity S).
function [Residuals] = ResidualOpt(g, params0)
    % Use 'persistent' variables to store values between function calls

    params0.g = g;
    params0.EvalFeatures = true;
    
    assert(length(params0.g) == length(params0.mods), "The params and g, woe");

    % --- Part 1: COSTLY MODEL RUN (Always required for residuals) ---
    try
        LoadData;
        RunBakersExp; % This is your costly model call!
        [Residuals, weights, cost] = evalFeatureCost(features_data, features_model, params0.fn, 1);
    catch e
        disp(e)
        Residuals = 1e3*ones(size(params0.fn));
    end

end