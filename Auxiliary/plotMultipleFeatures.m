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
%     FN           - cell array of field names to plot; supports 'y|x' syntax
%                    for X-Y scatter (same convention as plotFeatures).
%                    If omitted, uses union of all field names across structs.
%     LINESTYLES   - 1×N cell array of line style strings (default: all '-')
%     FILLMARKERS  - 1×N logical array; true = filled marker, false = empty
%                    (default: all true)
%     COSTMAP      - containers.Map from fn entry → scalar cost; when provided,
%                    each tile title shows the cost in brackets, e.g. "ktr [0.34]".
%                    Values ≥ 100 are shown as "[miss]" (missing feature penalty).
%                    Pass [] or omit to suppress cost display.
%
%   Each tile shows one feature vs segment index (or vs another feature
%   when using the 'y|x' syntax).  Datasets missing a requested field are
%   silently skipped for that tile.
%
%   See also: plotFeatures, extractSlackAttributes

nDS = numel(feats_cell);

if nargin < 8
    costMap = [];
end
if nargin < 6 || isempty(lineStyles)
    lineStyles = repmat({'-'}, 1, nDS);
end
if nargin < 7 || isempty(fillMarkers)
    fillMarkers = true(1, nDS);
end

if nargin < 5 || isempty(fn)
    % Union of all fields present in any dataset, in stable order
    fn = {};
    for di = 1:numel(feats_cell)
        if ~isempty(feats_cell{di})
            fn = union(fn, fieldnames(feats_cell{di}), 'stable');
        end
    end
end

tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

for i_feat = 1:numel(fn)
    nexttile(); hold on;

    % Parse 'y|x' or plain 'y'
    parts  = split(fn{i_feat}, '|');
    feat_y = parts{1};
    feat_x = '';
    if numel(parts) >= 2 && ~strcmp(parts{2}, '_') && isnan(str2double(parts{2}))
        feat_x = parts{2};
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

        clr  = colors(di, :);
        mkr  = markers{di};
        lst  = lineStyles{di};
        mfc  = clr;
        if ~fillMarkers(di); mfc = 'none'; end
        plot(vx, vy, [mkr lst], 'Color', clr, 'LineWidth', 1.5, ...
            'MarkerFaceColor', mfc, 'DisplayName', labels{di});
        plotted_any = true;
    end

    if ~plotted_any
        axis off;
        text(0.5, 0.5, sprintf('%s\n(no data)', strrep(fn{i_feat},'_','\_')), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', [0.6 0.6 0.6]);
        continue;
    end

    % Build cost suffix when a costMap was provided
    cost_str = '';
    if ~isempty(costMap) && isa(costMap, 'containers.Map') && isKey(costMap, fn{i_feat})
        c = costMap(fn{i_feat});
        if c >= 100
            cost_str = ' [miss]';
        elseif isnan(c)
            cost_str = ' [nan]';
        else
            cost_str = sprintf(' [%.2f]', c);
        end
    end

    if isempty(feat_x)
        title([strrep(fn{i_feat}, '_', '\_') cost_str], 'Interpreter', 'none');
        xlabel('Segment');
    else
        title(sprintf('%s vs %s%s', feat_y, feat_x, cost_str), 'Interpreter', 'none');
        xlabel(strrep(feat_x, '_', '\_'), 'Interpreter', 'none');
    end
    ylabel(strrep(feat_y, '_', '\_'), 'Interpreter', 'none');
    box on;
    if i_feat == 1
        legend('Location', 'best', 'FontSize', 7);
    end
end
