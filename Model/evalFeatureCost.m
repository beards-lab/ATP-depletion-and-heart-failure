
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
%     'field_y|field_x'         X-Y scatter (normalized by mean(data_y))
%     'field_y|field_x|weight'  X-Y scatter with weight
%     'field|lb-ub'             boundary check: cost=0 inside [lb,ub], quadratic outside
%     'field|weight|lb-ub'      boundary check with explicit weight
%     'f1|w1|lb1-ub1,f2|lb2-ub2'  grouped boundary check: one joint cost, one bar tile
%     'field|@(X,Y_data) expr'  custom cost function (escape hatch)
%
%   Boundary check cost (per element):
%     viol = max(0,(lb-x)/span)^2 + max(0,(x-ub)/span)^2   where span = ub-lb
%   Cost = 0 inside [lb,ub], rises to 1 at one span outside.
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
        prefix         = fn_str(1:at_pos(1)-1);
        cost_fn        = str2func(strtrim(fn_str(at_pos(1):end)));
        weight(i_feat) = 1;
    else
        prefix  = fn_str;
        cost_fn = [];
    end

    tokens = split(prefix, '|');
    tokens = tokens(~cellfun(@isempty, tokens));
    feat_y = tokens{1};

    % === BOUNDARY CHECK mode: a 'lb-ub' range token is present ===
    if isempty(cost_fn)
        range_tok = '';
        for ti = 2:numel(tokens)
            if isRangeToken(tokens{ti})
                range_tok = tokens{ti};
                break;
            end
        end

        if ~isempty(range_tok)
            % Parse optional weight (numeric token that is not the range)
            for ti = 2:numel(tokens)
                if ~isRangeToken(tokens{ti})
                    w = str2double(tokens{ti});
                    if ~isnan(w)
                        weight(i_feat) = w;
                    end
                end
            end

            if ~isfield(feats_sim, feat_y)
                cost(i_feat) = MISSING_FEATURE_COST;
                continue;
            end
            fs = double([feats_sim.(feat_y)]);

            [lb, ub] = parseRange(range_tok);
            span = ub - lb;
            viol = max(0, (lb - fs) / span) .^ 2 + max(0, (fs - ub) / span) .^ 2;
            cost(i_feat) = mean(viol, 'omitnan') + sum(isnan(fs)) * NAN_COST;
            continue;
        end
    end

    % === NORMAL mode (existing logic) ===
    if isempty(cost_fn)
        w = str2double(tokens{end});
        if ~isnan(w)
            weight(i_feat) = w;
        else
            weight(i_feat) = 1;
        end
    end

    if ~isfield(feats_data, feat_y) || ~isfield(feats_sim, feat_y)
        cost(i_feat) = MISSING_FEATURE_COST;
        continue;
    end

    fd = double([feats_data.(feat_y)]);
    fs = double([feats_sim.(feat_y)]);

    % Determine x-axis field for XY plots (second non-numeric token, no @)
    feat_x = '';
    if isempty(cost_fn) && numel(tokens) >= 2 && isnan(str2double(tokens{2}))
        feat_x = tokens{2};
    end
    if ~isempty(feat_x) && isfield(feats_sim, feat_x)
        % XY mode: fd/fs are the y-values; we don't reindex but keep same normalization
    end

    if ~isempty(cost_fn)
        cost(i_feat) = cost_fn(fs, fd) + sum(isnan(fs)) * NAN_COST;
    else
        n_F = mean(fd, 'omitnan');
        if n_F == 0; n_F = 1; end
        cost(i_feat) = sum(abs(fd/n_F - fs/n_F).^costExp, 'omitnan') + sum(isnan(fs)) * NAN_COST;
    end

end

cost_total = cost .* weight;
end


% -----------------------------------------------------------------------
function [jcost, jweight] = evalGroupedCost(fn_str, feats_data, feats_sim, ...
                                             NAN_COST, MISSING_FEATURE_COST, costExp)
% Split on top-level commas and sum weighted sub-costs.
parts   = splitTopLevelComma(fn_str);
jcost   = 0;
jweight = 1;   % group slot always has weight=1; sub-weights applied internally

for pi = 1:numel(parts)
    part   = strtrim(parts{pi});
    tokens = split(part, '|');
    tokens = tokens(~cellfun(@isempty, tokens));
    if isempty(tokens); continue; end
    feat = tokens{1};

    sub_weight = 1;
    range_tok  = '';
    for ti = 2:numel(tokens)
        if isRangeToken(tokens{ti})
            range_tok = tokens{ti};
        else
            w = str2double(tokens{ti});
            if ~isnan(w)
                sub_weight = w;
            end
        end
    end

    if ~isfield(feats_sim, feat)
        jcost = jcost + sub_weight * MISSING_FEATURE_COST;
        continue;
    end
    fs = double([feats_sim.(feat)]);

    if ~isempty(range_tok)
        [lb, ub] = parseRange(range_tok);
        span     = ub - lb;
        viol     = max(0, (lb - fs) / span) .^ 2 + max(0, (fs - ub) / span) .^ 2;
        sub_cost = mean(viol, 'omitnan') + sum(isnan(fs)) * NAN_COST;
    else
        % Plain feature in group: compare to data
        if ~isfield(feats_data, feat)
            jcost = jcost + sub_weight * MISSING_FEATURE_COST;
            continue;
        end
        fd   = double([feats_data.(feat)]);
        n_F  = mean(fd, 'omitnan');
        if n_F == 0; n_F = 1; end
        sub_cost = sum(abs(fd/n_F - fs/n_F).^costExp, 'omitnan') + sum(isnan(fs)) * NAN_COST;
    end

    jcost = jcost + sub_weight * sub_cost;
end
end


% -----------------------------------------------------------------------
function tf = hasTopLevelComma(s)
% True if s contains a comma outside any parentheses.
depth = 0;
for k = 1:numel(s)
    c = s(k);
    if c == '('; depth = depth + 1; end
    if c == ')'; depth = depth - 1; end
    if c == ',' && depth == 0; tf = true; return; end
end
tf = false;
end


function parts = splitTopLevelComma(s)
% Split s on commas that are not inside parentheses.
parts = {};
depth = 0;
start = 1;
for k = 1:numel(s)
    c = s(k);
    if c == '('; depth = depth + 1; end
    if c == ')'; depth = depth - 1; end
    if c == ',' && depth == 0
        parts{end+1} = s(start:k-1); %#ok<AGROW>
        start = k + 1;
    end
end
parts{end+1} = s(start:end);
end


function tf = isRangeToken(tok)
% True if tok matches 'lb-ub' where lb and ub parse as finite numbers.
% Uses lookbehind to find the separator dash (preceded by a digit).
dash = regexp(tok, '(?<=\d)-', 'once');
if isempty(dash); tf = false; return; end
lb = str2double(tok(1:dash-1));
ub = str2double(tok(dash+1:end));
tf = isfinite(lb) && isfinite(ub);
end


function [lb, ub] = parseRange(tok)
dash = regexp(tok, '(?<=\d)-', 'once');
lb   = str2double(tok(1:dash-1));
ub   = str2double(tok(dash+1:end));
end
