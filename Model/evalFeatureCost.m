
function [cost_total, weight, cost] = evalFeatureCost(feats_data, feats_sim, fn, costExp)
% EVALFEATURECOST  Compute weighted feature-based cost between data and simulation.
%
%   [COST_TOTAL, WEIGHT, COST] = EVALFEATURECOST(FEATS_DATA, FEATS_SIM, FN, COSTEXP)
%   evaluates the normalized absolute error for each feature group in FN.
%
%   FN is a cell array of feature specifiers. Each entry is a string that
%   supports the following forms:
%
%     'field'                   plain feature, weight=1
%     'field|weight'            plain feature with explicit weight
%     'field_y|field_x'         X-Y scatter (normalised by mean(data_y))
%     'field_y|field_x|weight'  X-Y scatter with weight
%     'field|lb-ub'             boundary check: cost=0 inside [lb,ub], quadratic outside
%     'field|lb-ub|weight'      boundary check with explicit weight
%     'f1|w1|lb1-ub1,f2|lb2-ub2'  grouped boundary check: one joint cost, one bar tile
%     'field|@(X,Y_data) expr'  custom cost function (escape hatch)
%
%   Boundary check cost (per element):
%     viol = max(0,(lb-x)/lb)^2 + max(0,(x-ub)/ub)^2
%   Cost = 0 inside [lb,ub]; rises to 1 when x reaches 0 (lower) or 2*ub (upper).
%   feats_data is NOT accessed for boundary entries.
%
%   Grouped entries (comma at top level, outside any parentheses):
%     Each sub-entry is 'field[|weight][|lb-ub]'. The returned cost is the
%     weighted sum of individual violations; weight for the group slot = 1.
%     For non-boundary sub-entries, cost is computed vs feats_data as normal.
%
%   NAN_COST = 1 per NaN in simulation output.
%   MISSING_FEATURE_COST = 100 when a field is absent from feats_sim.
%
%   Inputs:
%     FEATS_DATA - struct array with experimental feature values
%     FEATS_SIM  - struct array with simulated feature values
%     FN         - cell array of feature specifier strings
%     COSTEXP    - exponent for the error norm (default 2 = L2 norm)
%
%   Outputs:
%     COST_TOTAL - weighted cost per feature (length = numel(FN))
%     WEIGHT     - weight vector extracted from FN
%     COST       - unweighted cost per feature
%
%   See also: evaluateProblem, extractSlackAttributes, extractForceVelocityAttributes

NAN_COST = 1;
MISSING_FEATURE_COST = 100;

if nargin < 4
    costExp = 2;
end

cost   = zeros(1, length(fn));
weight = ones(1,  length(fn));

for i_feat = 1:length(fn)

    fn_str = fn{i_feat};

    % === GROUPED mode: comma at top level (outside any parentheses) ===
    if hasTopLevelComma(fn_str)
        [cost(i_feat), weight(i_feat)] = evalGroupedCost( ...
            fn_str, feats_data, feats_sim, NAN_COST, MISSING_FEATURE_COST, costExp);
        continue;
    end

    % Split off custom cost function at '@'
    at_pos = strfind(fn_str, '@');
    if ~isempty(at_pos)
        prefix  = fn_str(1:at_pos(1)-1);
        cost_fn = str2func(strtrim(fn_str(at_pos(1):end)));
    else
        prefix  = fn_str;
        cost_fn = [];
    end

    tokens = split(prefix, '|');
    tokens = tokens(~cellfun(@isempty, tokens));
    [feat_y, sel_y] = parseFeatSelector(tokens{1});

    if ~isempty(cost_fn)
        % === CUSTOM COST FUNCTION mode ===
        weight(i_feat) = 1;
        if ~isfield(feats_data, feat_y) || ~isfield(feats_sim, feat_y)
            cost(i_feat) = MISSING_FEATURE_COST; continue;
        end
        fd = applyFeatSelector(double([feats_data.(feat_y)]), sel_y);
        fs = applyFeatSelector(double([feats_sim.(feat_y)]), sel_y);
        cost(i_feat) = cost_fn(fs, fd) + sum(isnan(fs)) * NAN_COST;
        continue;
    end

    % Delegate to shared sub-entry evaluator (handles boundary + normal modes)
    [cost(i_feat), weight(i_feat)] = evalOneSubEntry( ...
        feat_y, sel_y, tokens(2:end), feats_data, feats_sim, ...
        NAN_COST, MISSING_FEATURE_COST, costExp);
end

cost_total = cost .* weight;
end


% -----------------------------------------------------------------------
function [jcost, jweight] = evalGroupedCost(fn_str, feats_data, feats_sim, ...
                                             NAN_COST, MISSING_FEATURE_COST, costExp)
% Split on top-level commas and sum weighted sub-costs into one slot.
parts   = splitTopLevelComma(fn_str);
jcost   = 0;
jweight = 1;   % group slot always has weight=1; sub-weights applied internally

for pi = 1:numel(parts)
    part   = strtrim(parts{pi});
    tokens = split(part, '|');
    tokens = tokens(~cellfun(@isempty, tokens));
    if isempty(tokens); continue; end
    [feat, sel] = parseFeatSelector(tokens{1});

    [sub_cost, sub_weight] = evalOneSubEntry( ...
        feat, sel, tokens(2:end), feats_data, feats_sim, ...
        NAN_COST, MISSING_FEATURE_COST, costExp);
    jcost = jcost + sub_weight * sub_cost;
end
end


% -----------------------------------------------------------------------
function [c, w] = evalOneSubEntry(feat, sel, tail_tokens, feats_data, feats_sim, ...
                                   NAN_COST, MISSING_FEATURE_COST, costExp)
% Evaluate cost for a single 'field[|weight][|lb-ub]' sub-entry.
% FEAT/SEL come from parseFeatSelector on the first token.
% TAIL_TOKENS is the remaining tokens (everything after the first pipe).
%
% Returns:
%   C  - unweighted cost scalar
%   W  - weight (default 1)

w         = 1;
range_tok = '';

for ti = 1:numel(tail_tokens)
    tok = tail_tokens{ti};
    if isRangeToken(tok)
        range_tok = tok;
    else
        wv = str2double(tok);
        if ~isnan(wv)
            w = wv;
        end
        % Non-numeric, non-range tokens are X-axis field names (XY plots) — ignored for cost.
    end
end

if ~isempty(range_tok)
    % --- BOUNDARY CHECK: cost against sim only ---
    if ~isfield(feats_sim, feat)
        c = MISSING_FEATURE_COST; return;
    end
    fs   = applyFeatSelector(double([feats_sim.(feat)]), sel);
    [lb, ub] = parseRange(range_tok);
    viol = max(0, (lb - fs) / lb) .^ 2 + max(0, (fs - ub) / ub) .^ 2;
    c    = mean(viol, 'omitnan') + sum(isnan(fs)) * NAN_COST;
else
    % --- NORMAL: normalised difference vs data ---
    if ~isfield(feats_data, feat) || ~isfield(feats_sim, feat)
        c = MISSING_FEATURE_COST; return;
    end
    fd  = applyFeatSelector(double([feats_data.(feat)]), sel);
    fs  = applyFeatSelector(double([feats_sim.(feat)]),  sel);
    n_F = mean(fd, 'omitnan');
    if n_F == 0; n_F = 1; end
    c = sum(abs(fd/n_F - fs/n_F).^costExp, 'omitnan') + sum(isnan(fs)) * NAN_COST;
end
end
