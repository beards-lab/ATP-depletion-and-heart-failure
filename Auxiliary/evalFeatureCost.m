
function [cost_total, weight, cost] = evalFeatureCost(feats_data, feats_sim, fn, costExp)
% EVALFEATURECOST  Compute weighted feature-based cost between data and simulation.
%
%   [COST_TOTAL, WEIGHT, COST] = EVALFEATURECOST(FEATS_DATA, FEATS_SIM, FN, COSTEXP)
%   evaluates the normalized absolute error for each feature group in FN.
%
%   FN is a cell array of feature specifiers. Each entry is a string that
%   may contain multiple feature names separated by '|'. The last token, if
%   numeric, is treated as a weight (0 = disabled, 1 = default).
%   Example: 'FV_f|FV_v|0.5' — uses both FV_f and FV_v features with weight 0.5
%            'ktr|0'          — evaluates ktr but does not count it in the total
%
%   Cost for each feature:
%     cost(i) = sum(|data/mean(data) - sim/mean(data)|^costExp) + NaN_penalty
%   where NaN_penalty = NAN_COST per NaN in the simulation output.
%
%   If a feature field is absent from either struct (experiment not run),
%   cost(i) = MISSING_FEATURE_COST (100) — penalised but no error raised.
%
%   NAN_COST = 10: chosen to be large relative to a well-fit feature
%   (typical normalized errors are O(0.1-1)), ensuring NaN outputs are
%   heavily penalized without dominating well-fitting features.
%   MISSING_FEATURE_COST = 100: larger than NAN_COST so absent experiments
%   are clearly distinguishable from poor fits.
%
%   Inputs:
%     FEATS_DATA - struct array with experimental feature values
%     FEATS_SIM  - struct array with simulated feature values
%     FN         - cell array of feature specifier strings
%     COSTEXP    - exponent for the error norm (default 1 = L1 norm)
%
%   Outputs:
%     COST_TOTAL - weighted cost per feature (length = numel(FN))
%     WEIGHT     - weight vector extracted from FN
%     COST       - unweighted cost per feature
%
%   See also: evaluateProblem, extractSlackAttributes, extractForceVelocityAttributes

% Penalty applied per NaN in simulation output.
% Set large enough to dominate a typical normalized residual (O(0.1)).
NAN_COST = 10;

% Penalty when a requested feature is absent (experiment not run or not extracted).
% Chosen to be large relative to NAN_COST so missing experiments are clearly
% distinguishable from poor fits.
MISSING_FEATURE_COST = 100;

if nargin < 4
    costExp = 1;
end
for i_feat = 1:size(fn, 2)

    featnames = split(fn{i_feat}, '|');
    feat_y{i_feat} = featnames{1};

    if isempty(str2double(featnames{end})) || isnan(str2double(featnames{end}))
        weight(i_feat) = 1;
    else
        weight(i_feat) = str2double(featnames{end});
    end

    % If either struct is missing the field, the experiment was not run.
    % Penalise but do not error.
    if ~isfield(feats_data, feat_y{i_feat}) || ~isfield(feats_sim, feat_y{i_feat})
        cost(i_feat) = MISSING_FEATURE_COST;
        continue;
    end

    fd = [feats_data.(feat_y{i_feat})];
    fs = [feats_sim.(feat_y{i_feat})];

    % normalization factor
    n_F = mean(fd);

    cost(i_feat) = nansum(abs(fd/n_F - fs/n_F).^costExp) + sum(isnan(fs)*NAN_COST);

end

cost_total = cost.*weight;
