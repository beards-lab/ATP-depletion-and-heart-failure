function [Et, E, featCost, featNames] = evaluateBakersExp(g, params0, limitNegative, plotDebug)
% EVALUATEBAKERSEXP  Scalar objective for the Baker problem.
%   Et = EVALUATEBAKERSEXP(g, params0) returns the total cost (sum of E in
%   'time' mode, sum of the weighted feature cost in 'Feats' mode).
%
%   [Et, E, FEATCOST, FEATNAMES] = EVALUATEBAKERSEXP(...) additionally returns
%   the per-experiment error E and the PER-FEATURE cost vector FEATCOST (one
%   weighted cost per params0.fn entry -- grouped entries collapsed to their
%   single slot, i.e. per-feature not per-inner-item) with matching FEATNAMES
%   labels. FEATCOST is computed ONLY when the caller requests >= 3 outputs, so
%   the hot fminsearch/surrogateopt path (which asks for just Et) pays nothing.
%   FEATCOST lets callers record which features got better/worse across a run
%   regardless of OptimizeOn.

if nargin < 3
    limitNegative = true;
end
if nargin < 4
    plotDebug = false;
end

% Initialise all outputs so every return path (early skip, error) is safe even
% when the caller requests the extra outputs.
Et = NaN; E = []; featCost = []; featNames = {};
if nargout >= 3 && isfield(params0, 'fn') && ~isempty(params0.fn)
    featNames = iTrimFn(params0.fn);
    featCost  = nan(1, numel(params0.fn));
end
wantFeat = nargout >= 3;

% Evaluate Bakers' problem
if limitNegative && any(g<0)% || ...
        %params0.kstiff1 < params0.kstiff1_n || params0.kstiff2 < params0.kstiff2_n
    Et = NaN;
    warning(' g < 0! Skipping ...');
    return;
end

% --- Joint (dual-condition) ABSOLUTE ATP objective -------------------------
% Et = w8*E(high ATP, as configured) + w2*E(low ATP, AtpJoint_* overrides).
% Each leg is evaluated by a recursive call with the flag cleared, so all the
% existing machinery (RunBakersExp, evalFeatureCost, featCost decomposition) is
% reused verbatim. Guarded => inert when UseAtpJointFit is false.
if isfield(params0, 'UseAtpJointFit') && params0.UseAtpJointFit
    p8 = params0; p8.UseAtpJointFit = false;

    p2 = p8;
    p2.MgATP               = params0.AtpJoint_MgATP2;
    p2.Pi                  = params0.AtpJoint_Pi2;
    p2.velocitytableonfile = params0.AtpJoint_vt2;
    % The low-ATP day has slack data only: drop the experiments it cannot score,
    % otherwise their features register as missing (cost 100 each) every eval.
    if ~params0.AtpJoint_RunFV2
        p2.RunForceVelocity = false; p2.RunForceVelocityTime = false;
    end
    if ~params0.AtpJoint_RunPS2
        p2.RunSlackPassive = false;
    end
    if ~isempty(params0.AtpJoint_fn2); p2.fn = params0.AtpJoint_fn2; end

    if wantFeat
        [Et8, E8, fc8, fn8] = evaluateBakersExp(g, p8, limitNegative, plotDebug);
        [Et2, E2, fc2, fn2] = evaluateBakersExp(g, p2, limitNegative, plotDebug);
        featCost  = [params0.AtpJoint_w8 * fc8, params0.AtpJoint_w2 * fc2];
        featNames = [strcat('hi:', fn8), strcat('lo:', fn2)];
    else
        [Et8, E8] = evaluateBakersExp(g, p8, limitNegative, plotDebug);
        [Et2, E2] = evaluateBakersExp(g, p2, limitNegative, plotDebug);
    end
    Et = params0.AtpJoint_w8 * Et8 + params0.AtpJoint_w2 * Et2;
    E  = [E8(:); E2(:)]';
    return;
end

if ~plotDebug
    params0.PlotEachSeparately = false;
    params0.PlotFeatureFitting = false;
else
    params0.PlotEachSeparately = true;
end
% important to start with the g!!
% params0 = getParams(params0, g, true, true);
params0.g = g;

try
%     Et = 1;
    LoadData;
    lastwarn('', '');
    RunBakersExp;
    if strcmp(params0.OptimizeOn, 'Feats')
        % Et = sum(E);
        % use feature costs instead
        [cost, ~, ~] = evalFeatureCost(features_data, features_model, params0.fn);
        Et = sum(cost);
        if wantFeat; featCost = cost; end   % weighted per-feature (sums to Et here)
        if plotDebug
            plotFeatures(features_data, features_model, [], params0.fn);
            disp(Et);
        end
    else
        Et = sum(E);
        if wantFeat
            % Diagnostic decomposition even though the objective is time-based:
            % which features fit well/poorly at this point. Missing experiments
            % (not run) show as the evalFeatureCost missing-sentinel; callers
            % filter those out.
            try
                [featCost, ~, ~] = evalFeatureCost(features_data, features_model, params0.fn);
            catch
                featCost = nan(1, numel(params0.fn));
            end
        end
    end

    if Et == 0
        % sum fin wong
        Et = 1e3;
    end
catch e
    try
        if ~isempty(e.cause)
            % Extract the cause message if present
            causeMessage = e.cause{end}.message;
            warning(causeMessage);

            % Parse the error value from the cause message
            Et = 1e3 + str2double(causeMessage);
        else
            Et = 1e9;
        end
    catch ee
        Et = 1e12;
    end
end

Et = abs(Et);
end

% -----------------------------------------------------------------------
function nm = iTrimFn(fn)
% Compact per-feature label: keep the primary field, drop selectors/weights/
% ranges and any grouped tail (everything from the first '[', '|' or ',').
nm = cellfun(@(s) regexprep(s, '[\[|,].*$', ''), fn, 'UniformOutput', false);
end
