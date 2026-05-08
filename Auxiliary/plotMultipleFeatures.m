function plotMultipleFeatures(feats_cell, labels, colors, markers, fn, lineStyles, fillMarkers, costMap)
%% PLOTMULTIPLEFEATURES  Plot one tile per feature, one series per dataset.
%
%   PLOTMULTIPLEFEATURES(FEATS_CELL, LABELS, COLORS, MARKERS, FN)
%   PLOTMULTIPLEFEATURES(FEATS_CELL, LABELS, COLORS, MARKERS, FN, LINESTYLES, FILLMARKERS)
%   PLOTMULTIPLEFEATURES(FEATS_CELL, LABELS, COLORS, MARKERS, FN, LINESTYLES, FILLMARKERS, COSTMAP)
%
%   Inputs:
%     FEATS_CELL   - 1×N cell array of feature structs (one per dataset)
%     LABELS       - 1×N cell array of display names
%     COLORS       - N×3 RGB matrix
%     MARKERS      - 1×N cell array of marker chars (e.g. {'o','s','^'})
%     FN           - cell array of field names to plot; supports:
%                      'y|x'              X-Y scatter
%                      'field|lb-ub'      single boundary check → bar chart
%                      'f1|lb1-ub1,...'   grouped boundary checks → bar chart
%                    If omitted, uses union of all field names across structs.
%     LINESTYLES   - 1×N cell array of line style strings (default: all '-')
%     FILLMARKERS  - 1×N logical array; true = filled marker (default: all true)
%     COSTMAP      - containers.Map from fn entry → scalar cost; values ≥ 100
%                    shown as "[miss]". Pass [] or omit to suppress.
%
%   See also: plotFeatures, evalFeatureCost, extractSlackAttributes

nDS = numel(feats_cell);

if nargin < 8; costMap    = []; end
if nargin < 6 || isempty(lineStyles);  lineStyles  = repmat({'-'}, 1, nDS); end
if nargin < 7 || isempty(fillMarkers); fillMarkers = true(1, nDS); end

if nargin < 5 || isempty(fn)
    fn = {};
    for di = 1:numel(feats_cell)
        if ~isempty(feats_cell{di})
            fn = union(fn, fieldnames(feats_cell{di}), 'stable');
        end
    end
end

tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

for i_feat = 1:numel(fn)
    nexttile(); ax = gca; hold(ax, 'on');

    fn_str    = fn{i_feat};
    cost_str  = buildCostStr(costMap, fn_str);

    % ── GROUPED: comma at top level ──────────────────────────────────────
    if hasTopLevelComma(fn_str)
        plotGroupedBoundaryTile(ax, fn_str, feats_cell, cost_str);
        hold(ax, 'off');
        continue;
    end

    % ── Parse prefix / @ function ─────────────────────────────────────────
    at_pos = strfind(fn_str, '@');
    if ~isempty(at_pos)
        prefix      = fn_str(1:at_pos(1)-1);
        fn_body     = strtrim(fn_str(at_pos(1):end));
        has_cost_fn = true;
    else
        prefix      = fn_str;
        fn_body     = '';
        has_cost_fn = false;
    end
    tokens = split(prefix, '|');
    tokens = tokens(~cellfun(@isempty, tokens));
    feat_y = tokens{1};

    % ── SINGLE BOUNDARY: range token present ─────────────────────────────
    if ~has_cost_fn
        range_tok = '';
        for ti = 2:numel(tokens)
            if isRangeToken(tokens{ti})
                range_tok = tokens{ti};
                break;
            end
        end
        if ~isempty(range_tok)
            [lb, ub] = parseRange(range_tok);
            plotSingleBoundaryTile(ax, feat_y, lb, ub, feats_cell, cost_str);
            hold(ax, 'off');
            continue;
        end
    end

    % ── NORMAL scatter plot (existing logic) ─────────────────────────────
    feat_x = '';
    if ~has_cost_fn && numel(tokens) >= 2 && ...
            ~strcmp(tokens{2}, '_') && isnan(str2double(tokens{2}))
        feat_x = tokens{2};
    end

    plotted_any = false;
    for di = 1:nDS
        fs = feats_cell{di};
        if isempty(fs) || ~isfield(fs, feat_y); continue; end
        vy = fs.(feat_y);
        if ~isnumeric(vy) || isempty(vy); continue; end

        if isempty(feat_x)
            vx = 1:numel(vy);
        else
            if ~isfield(fs, feat_x); continue; end
            vx = fs.(feat_x);
            if ~isnumeric(vx) || isempty(vx); continue; end
        end

        clr = colors(di, :);
        mkr = markers{di};
        lst = lineStyles{di};
        mfc = clr;
        if ~fillMarkers(di); mfc = 'none'; end
        plot(ax, vx, vy, [mkr lst], 'Color', clr, 'LineWidth', 1.5, ...
            'MarkerFaceColor', mfc, 'DisplayName', labels{di});
        plotted_any = true;
    end

    if ~plotted_any
        axis(ax, 'off');
        text(ax, 0.5, 0.5, sprintf('%s\n(no data)', strrep(fn_str,'_','\_')), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', [0.6 0.6 0.6]);
        hold(ax, 'off');
        continue;
    end

    if has_cost_fn
        fn_disp = fn_body;
        if numel(fn_disp) > 42; fn_disp = [fn_disp(1:42) '...']; end
        title(ax, {[strrep(feat_y,'_','\_') cost_str], fn_disp}, 'Interpreter', 'none');
        xlabel(ax, 'Segment');
    elseif isempty(feat_x)
        title(ax, [strrep(fn_str,'_','\_') cost_str], 'Interpreter', 'none');
        xlabel(ax, 'Segment');
    else
        title(ax, sprintf('%s vs %s%s', feat_y, feat_x, cost_str), 'Interpreter', 'none');
        xlabel(ax, strrep(feat_x,'_','\_'), 'Interpreter', 'none');
    end
    ylabel(ax, strrep(feat_y,'_','\_'), 'Interpreter', 'none');
    box(ax, 'on');
    if i_feat == 1
        legend(ax, 'Location', 'best', 'FontSize', 7);
    end
    hold(ax, 'off');
end
end


% ── Single boundary tile ──────────────────────────────────────────────────
function plotSingleBoundaryTile(ax, feat_y, lb, ub, feats_cell, cost_str)
% Bar chart: one bar per sim value, colored by in/out of [lb, ub].
% Uses only the first non-empty dataset that has the field (typically sim).
fs = [];
for di = 1:numel(feats_cell)
    fc = feats_cell{di};
    if ~isempty(fc) && isfield(fc, feat_y)
        v = double([fc.(feat_y)]);
        if ~isempty(v)
            fs = v;
            break;
        end
    end
end

cla(ax); hold(ax, 'on');
if isempty(fs)
    text(ax, 0.5, 0.5, sprintf('%s\n(no data)', strrep(feat_y,'_','\_')), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'Color', [0.6 0.6 0.6]);
    axis(ax, 'off');
    title(ax, [strrep(feat_y,'_','\_') cost_str], 'Interpreter', 'none');
    return;
end

for idx = 1:numel(fs)
    v = fs(idx);
    if isnan(v)
        clr = [0.85 0.55 0.1];
    elseif v >= lb && v <= ub
        clr = [0.15 0.65 0.15];
    else
        clr = [0.80 0.15 0.15];
    end
    bar(ax, idx, v, 'FaceColor', clr, 'EdgeColor', 'none');
end

yline(ax, lb, '--', 'Color', [0.1 0.6 0.1], 'LineWidth', 1.2, ...
    'Label', sprintf('lb=%.3g', lb), 'LabelVerticalAlignment', 'bottom', ...
    'FontSize', 7);
yline(ax, ub, '--', 'Color', [0.1 0.6 0.1], 'LineWidth', 1.2, ...
    'Label', sprintf('ub=%.3g', ub), 'LabelVerticalAlignment', 'top', ...
    'FontSize', 7);

xlabel(ax, 'Segment', 'FontSize', 8);
ylabel(ax, strrep(feat_y,'_','\_'), 'Interpreter', 'none', 'FontSize', 8);
title(ax, sprintf('%s  [%.3g – %.3g]%s', strrep(feat_y,'_','\_'), lb, ub, cost_str), ...
    'Interpreter', 'none', 'FontSize', 9);
box(ax, 'on');
hold(ax, 'off');
end


% ── Grouped boundary tile ─────────────────────────────────────────────────
function plotGroupedBoundaryTile(ax, fn_str, feats_cell, cost_str)
% Normalized bar chart (physiology-dashboard style) for grouped boundary entries.
% Each sub-entry 'field[|weight][|lb-ub]' gets one row.
% Non-boundary sub-entries show sim value on the y-axis directly (no normalization).

parts = splitTopLevelComma(fn_str);

% Collect parsed sub-entries
sub_names = {};
sub_lbs   = [];
sub_ubs   = [];
sub_vals  = {};   % cell of value vectors (one per feats_cell dataset with data)
sub_has_bound = logical([]);

for pi = 1:numel(parts)
    part   = strtrim(parts{pi});
    tokens = split(part, '|');
    tokens = tokens(~cellfun(@isempty, tokens));
    if isempty(tokens); continue; end
    feat  = tokens{1};
    range_tok = '';
    for ti = 2:numel(tokens)
        if isRangeToken(tokens{ti})
            range_tok = tokens{ti};
            break;
        end
    end

    % Gather sim values (use first dataset that has the field)
    val = NaN;
    for di = 1:numel(feats_cell)
        fc = feats_cell{di};
        if ~isempty(fc) && isfield(fc, feat)
            v = double([fc.(feat)]);
            if ~isempty(v); val = v; break; end
        end
    end

    sub_names{end+1} = feat; %#ok<AGROW>
    sub_vals{end+1}  = val;  %#ok<AGROW>
    if ~isempty(range_tok)
        [lb, ub] = parseRange(range_tok);
        sub_lbs(end+1)        = lb;  %#ok<AGROW>
        sub_ubs(end+1)        = ub;  %#ok<AGROW>
        sub_has_bound(end+1)  = true; %#ok<AGROW>
    else
        sub_lbs(end+1)        = NaN; %#ok<AGROW>
        sub_ubs(end+1)        = NaN; %#ok<AGROW>
        sub_has_bound(end+1)  = false; %#ok<AGROW>
    end
end

if isempty(sub_names)
    axis(ax, 'off'); return;
end

% Build flat bar list: one bar per (sub-entry × element)
bar_labels = {};
bar_vals   = [];
bar_lbs    = [];
bar_ubs    = [];
bar_has_b  = logical([]);

for pi = 1:numel(sub_names)
    v  = sub_vals{pi};
    lb = sub_lbs(pi);
    ub = sub_ubs(pi);
    hb = sub_has_bound(pi);
    nm = sub_names{pi};
    for ei = 1:numel(v)
        if numel(v) > 1
            bar_labels{end+1} = sprintf('%s(%d)', strrep(nm,'_','\_'), ei); %#ok<AGROW>
        else
            bar_labels{end+1} = strrep(nm,'_','\_'); %#ok<AGROW>
        end
        bar_vals(end+1) = v(ei); %#ok<AGROW>
        bar_lbs(end+1)  = lb;    %#ok<AGROW>
        bar_ubs(end+1)  = ub;    %#ok<AGROW>
        bar_has_b(end+1)= hb;    %#ok<AGROW>
    end
end

n = numel(bar_vals);
cla(ax); hold(ax, 'on');

for bi = 1:n
    v  = bar_vals(bi);
    lb = bar_lbs(bi);
    ub = bar_ubs(bi);
    hb = bar_has_b(bi);

    if isnan(v)
        clr = [0.85 0.55 0.1];
    elseif hb && v >= lb && v <= ub
        clr = [0.15 0.65 0.15];
    elseif hb
        clr = [0.80 0.15 0.15];
    else
        clr = [0.30 0.55 0.80];
    end

    bar(ax, bi, v, 'FaceColor', clr, 'EdgeColor', 'none');

    % Draw bound lines specific to this bar's x-position
    if hb && isfinite(lb)
        plot(ax, bi + [-0.4 0.4], [lb lb], '-', 'Color', [0.1 0.6 0.1], 'LineWidth', 1.5);
        plot(ax, bi + [-0.4 0.4], [ub ub], '-', 'Color', [0.1 0.6 0.1], 'LineWidth', 1.5);
    end

    % Value label
    text(ax, bi, max(v, 0) + range([bar_vals bar_lbs bar_ubs], 'all') * 0.03, ...
        sprintf('%.3g', v), 'FontSize', 6, 'HorizontalAlignment', 'center');
end

set(ax, 'XTick', 1:n, 'XTickLabel', bar_labels, 'FontSize', 7);
ylabel(ax, 'Value', 'FontSize', 8);

% Build title from field names (unique, ordered)
seen = {};
for pi = 1:numel(sub_names)
    if ~ismember(sub_names{pi}, seen); seen{end+1} = sub_names{pi}; end %#ok<AGROW>
end
title(ax, [strjoin(cellfun(@(s) strrep(s,'_','\_'), seen, 'UniformOutput', false), ', ') cost_str], ...
    'Interpreter', 'none', 'FontSize', 9);
box(ax, 'on');
hold(ax, 'off');
end


% ── Helpers (mirrors of evalFeatureCost helpers) ──────────────────────────
function cs = buildCostStr(costMap, key)
cs = '';
if isempty(costMap) || ~isa(costMap, 'containers.Map') || ~isKey(costMap, key)
    return;
end
c = costMap(key);
if c >= 100
    cs = ' [miss]';
elseif isnan(c)
    cs = ' [nan]';
else
    cs = sprintf(' [%.2f]', c);
end
end


function tf = hasTopLevelComma(s)
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
