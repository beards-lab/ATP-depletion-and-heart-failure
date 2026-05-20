% RUNBAKERSEXP  Evaluate all Baker lab experimental protocols.
%
%   This script is a coordinator: it calls each protocol function, collects
%   individual cost terms into E, and merges feature structs for feature-based
%   optimization. It expects the following variables in the caller workspace:
%
%     params0      - Parameter struct (built by getParams, with all Run* flags set)
%     ATP_c        - ATP concentration vector (mM)
%     FV_velocities - Shortening velocities for the FV protocol (ML/s)
%     Data_ATP     - Experimental FV data matrix
%     Ktr_mean     - Target ktr value(s) from data (s^-1)
%
%   After running, the following workspace variables are set:
%     E              - Row vector of cost terms (one per active protocol)
%     features_model - Merged struct of model features (when EvalFeatures is true)
%     features_data  - Merged struct of data features  (when EvalFeatures is true)
%     out, outs      - Last output struct(s) from evaluateModel

E = [];
features_model = struct();
features_data  = struct();

% Resolve the ODE function handle.
% Protocol functions also resolve this internally from params0.modelFcn,
% but we keep a local handle for the inline MinStretch/ForceLengthEstim sections.
modelFcn = @dPUdT_CombinedTransitions;
if isfield(params0, 'modelFcn')
    modelFcn = str2func(params0.modelFcn);
end

%% FORCE VELOCITY
if params0.RunForceVelocity
    LoadData;
    [E_fv, outs_fv, fm_fv, fd_fv] = runFVExperiment(params0, ATP_c, Data_ATP);
    E(end+1) = E_fv;
    features_model = mergeStructs(features_model, fm_fv);
    features_data  = mergeStructs(features_data,  fd_fv);
    % expose outs and F_active for downstream sections (RunForceVelocityTime)
    outs = outs_fv;
    F_active = arrayfun(@(o) o.Force(end), outs_fv);
end

%% KTR EXPERIMENT
if params0.RunKtr
    [E_ktr, out_ktr] = runKtrExperiment(params0, Ktr_mean);
    E(end+1) = E_ktr;
    out = out_ktr;
end

%% RAMP UP
if params0.RunStairs
    [E_stairs, out_stairs] = runStairsExperiment(params0);
    E(end+1) = E_stairs;
    out = out_stairs;
end

%% SLACK
if params0.RunSlack
    [E_slack_vec, out_slack, fm_slack, fd_slack] = runSlackExperiment(params0);
    E = [E, E_slack_vec];
    features_model = mergeStructs(features_model, fm_slack);
    features_data  = mergeStructs(features_data,  fd_slack);
    out = out_slack;
end

%% MINIMAL STRETCH
if isfield(params0, 'RunMinStretch') && params0.RunMinStretch

    fname = ['../data/PassiveCaSrc2/' '20241121' '/02_03_Fmax_stiff_ktr.txt'];
    fmaxdata = readtable(fname, 'ReadVariableNames', false);

    datatable = table2array([1 2.0 1].*fmaxdata(:, 1:3));

    h = [2.5, 5, 7.5]*1e-3;
    dt = 1e-3;
    s = 5e-4;

    velocitytable = [
    -2, 0, 0, 2.0;
    s+ 0.1, h(1)/dt, 2*h(1)*dt, 2.0;
    s+ 0.1 + dt, 0, 0, 2.0 + h(1)*2;
    s+ 0.1+0.2, -h(1)/dt, -2*h(1)*dt, 2.0 + h(1)*2;
    s+ 0.1+0.2+dt, 0, 0, 2.0;

    s+ 0.5, h(2)/dt, 2*h(2)*dt, 2.0;
    s+ 0.5 + dt, 0, 0, 2.0 + h(2)*2;
    s+ 0.5+0.2, -h(2)/dt, -2*h(2)*dt, 2.0 + h(2)*2;
    s+ 0.5+0.2+dt, 0, 0, 2.0;

    s+ 0.9, h(3)/dt, 2*h(3)*dt, 3.0;
    s+ 0.9 + dt, 0, 0, 2.0 + h(3)*2;
    s+ 0.9+0.2, -h(3)/dt, -2*h(3)*dt, 2.0 + h(3)*2;
    s+ 0.9+0.2+dt, 0, 0, 2.0;

    1.5, 0, 0, 2.0;
    ];

    v = 250;
    vt_ktr = [
   -0.0146,  -v
   -0.0146 + 0.2/v,   0
   -0.0046,   v
   -0.0046 + 0.25/v,  0
   -0.0001,  -v
   -0.0001 + 0.05/v,  0
    1.0000,   0    ];

    velocitytable = [velocitytable; [vt_ktr + [1.6, 0], vt_ktr(:, 2)/2, vt_ktr(:, 1)*0]];

    params = params0;
    params.Velocity = velocitytable(:, 2);
    params.SL0 = 2.0;
    params.Slim_l = 1.9;
    params.Slim_r = 2.0;
    params = getParams(params, params.g, true);
    params.UseSuperRelaxed = false;
    params.UseSuperRelaxedADP = false;
    params.UseA2MechanicalRecocking = false;

    [~, out] = evaluateModel(modelFcn, velocitytable(:, 1), params);

    mask = out.t >= 0;
    fmaxdataCurrent = {};
    fmaxdataCurrent{1} = datatable;
    fmaxdataCurrent{end+1} = [out.t(mask)', out.SL(mask)', out.Force(mask)'];
    ProcessFmaxData;
end

%% FORCE VELOCITY RESIMULATION
if params0.RunForceVelocityTime
    if exist('F_active', 'var')
        runFVTimecourseExperiment(params0, Data_ATP, FV_velocities, F_active);
    else
        runFVTimecourseExperiment(params0, Data_ATP, FV_velocities, []);
    end
end

%% FORCE LENGTH ESTIMATED relation
if params0.RunForceLengthEstim
    params = params0;
    params.PlotEachSeparately = true;
    gp = ones(5, 1);
    Es = simulateForceLengthEstim(ones(5, 1), params, modelFcn, params.PlotEachSeparately);
    E(end+1) = Es*1e0;
end

%% PARAMETER BOUNDS AUTO-FN (value registration only — fn is unchanged)
% Populate features_model with the current value of every resolvable Tier A–C param.
% This guarantees that any boundary-check entry the user places in params0.fn
% (e.g. 'ka|50-500|10') can find features_model.ka without triggering a 100-penalty.
% params0.fn is NOT modified here; 'AssertParams|XX' expansion is handled below.
if params0.EvalFeatures
    bounds = parameterBounds();
    bn     = fieldnames(bounds);
    params = getParams(params0, params0.g, false, true);
    for i = 1:numel(bn)
        b   = bounds.(bn{i});
        if b.weight == 0; continue; end          % skip Tier D
        val = resolveModelParam(params, bn{i});
        if isnan(val);   continue; end           % absent / inapplicable param
        features_model.(bn{i}) = val;
        features_data.(bn{i})  = nan;
    end
end

%% ASSERT ACTIVE PARAMS EXPANSION
% Replace any 'AssertParams|XX' placeholders in params0.fn with a single
% grouped assertion string built from the current parametrisation switches
% (see assertActiveParams.m).  XX is a weight scale multiplier.
%
% Example: params0.fn = {'ktr', 'A', 'AssertParams|1'};
% After expansion, 'AssertParams|1' becomes e.g.
%   'ka|50-500|10,kd|5-200|10,kstiff1|2000-30000|10,...'
% which evalFeatureCost renders as one joint-cost tile.
%
% Using 'AssertParams|0.1' lowers global physiology pressure by 10×.
if params0.EvalFeatures && isfield(params0, 'fn')
    keep = true(1, numel(params0.fn));
    for i = 1:numel(params0.fn)
        tok = regexp(params0.fn{i}, '^AssertParams\|(\S+)$', 'tokens', 'once');
        if isempty(tok); continue; end
        sc = str2double(tok{1});
        [fn_grouped, feat_assert] = assertActiveParams(params0, sc);
        if isempty(fn_grouped)
            keep(i) = false;            % nothing active — drop the placeholder
        else
            params0.fn{i} = fn_grouped; % replace placeholder with grouped string
        end
        % fn_assert = fieldnames(feat_assert);
        % for j = 1:numel(fn_assert)
        %     features_model.(fn_assert{j}) = feat_assert.(fn_assert{j});
        %     features_data.(fn_assert{j})  = nan;
        % end
    end
    params0.fn = params0.fn(keep);
end

%% FEATURE-BASED COST
if params0.EvalFeatures && isfield(params0, 'fn') && ~isempty(params0.fn)
    E_feat = evalFeatureCost(features_data, features_model, params0.fn, 1);
    E = [E, E_feat];
end

%% PHYSIOLOGY-BOUND COST (fallback — only when EvalFeatures is off)
% When EvalFeatures=true, bounds cost flows via fn entries + evalFeatureCost above.
% This block exists purely for backward compatibility when EvalFeatures=false.
if ~params0.EvalFeatures && isfield(params0, 'w_phys') && params0.w_phys > 0
    %%
    if ~isfield(params0, 'physiologyBounds') || isempty(params0.physiologyBounds)
        params0.physiologyBounds = parameterBounds();
    end
    E_phys = evalPhysiologyCost(params0, params0.physiologyBounds);
    E = [E, params0.w_phys * E_phys];
end

%% READJUST ERROR FUNCTION
if isfield(params0, 'ErrorMultiplier')
    erl = 1:min(length(E), length(params0.ErrorMultiplier));
    E(erl) = E(erl) .* params0.ErrorMultiplier(erl);
end

%% -----------------------------------------------------------------------
function s = mergeStructs(s, s2)
%MERGESTRUCTS Copy all fields from s2 into s, overwriting on collision.
    fns = fieldnames(s2);
    for i = 1:numel(fns)
        s.(fns{i}) = s2.(fns{i});
    end
end
