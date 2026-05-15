function cost = plotFeatures(feats_data, feats_sim, feats_ghost, fn, params)
%% PLOTFEATURES  Plot data vs simulation features (backward-compatible wrapper).
%
%   Wrapper around PLOTMULTIPLEFEATURES for the classic 3-dataset case:
%   data (blue o), simulation (red x), ghost (grey +).
%
%   COST = PLOTFEATURES(FEATS_DATA, FEATS_SIM, FEATS_GHOST, FN)
%   COST = PLOTFEATURES(..., PARAMS) additionally renders the physiology
%   dashboard (parameterBounds vs current params) in figure 80086.
%
%   To build the sample fn, go:
%     fn = fieldnames(feats_sim);
%     quoted = cellfun(@(s) ['''' s ''''], fn, 'UniformOutput', false);
%     fprintf('fn = {%s};\n', strjoin(quoted, ', '));
%   A sample fn:
%     fn = {'ktr','A','t0','SLslack','v_restretch','peak1_y','peak1_t', ...
%           'peak1_SL','peak1_dSL','peak2','vall_t','vall_y','vall2_dy', ...
%           'vall2_t','ovrsht_dy','ovrsht_t','steady'};
%   X-Y plots: plots A as a function of SLslack
%     fn = {'A|SLslack'}
%
%   See also: plotMultipleFeatures, evalFeatureCost

if nargin < 3
    feats_ghost = [];
    fn = fieldnames(feats_data);
elseif nargin < 4
    fn = fieldnames(feats_data);
end

% Build the multi-dataset inputs
feats_cell = {feats_data, feats_sim};
labels     = {'data', 'sim'};
colors     = [0.00 0.45 0.70; 0.84 0.10 0.11];
markers    = {'o', 'x'};

if ~isempty(feats_ghost)
    feats_cell{end+1} = feats_ghost;
    labels{end+1}     = 'ghost';
    colors(end+1, :)  = [0.50 0.50 0.50];
    markers{end+1}    = '+';
end

% Compute per-feature costs (data vs sim only) and build a lookup map for
% display in tile titles.  Missing or NaN features get their penalty values
% from evalFeatureCost (≥ 100 → shown as "[miss]").
[cost, ~, cost_raw] = evalFeatureCost(feats_data, feats_sim, fn);
% cost_raw is the unweighted per-feature error; cost is the weighted total
cmap = containers.Map(fn, num2cell(cost));

figure(80085); clf;
plotMultipleFeatures(feats_cell, labels, colors, markers, fn, [], [], cmap);

% Append a "culprits" tile: bar chart of per-feature cost, sorted desc.
% Makes the dominant cost contributors immediately visible.
nexttile();
plotCostCulprits(fn, cost, gca);

% Append a "bound violations" tile if params provided: which physiology
% bounds are out, with their individual penalty contributions and the sum.
if nargin >= 5 && ~isempty(params)
    params = getParams(params, params.g, false, true);
    [~, violations] = evalPhysiologyCost(params, parameterBounds());
    nexttile();
    plotBoundViolators(violations, gca);

    % Full physiology dashboard goes into a separate figure so it doesn't
    % squeeze the feature tile layout.
    figure(80086); clf;
    plotPhysiologyDashboard(params, [], gca);
end

end


function plotCostCulprits(fn, cost, ax)
% Sort features by cost descending, plot horizontal bar chart.
[c_sorted, ord] = sort(cost(:), 'descend');
fn_sorted = fn(ord);

% Trim labels for display: drop everything from the first | onwards
labels = cellfun(@(s) regexprep(s, '\|.*$', ''), fn_sorted, 'UniformOutput', false);

% Color: red = missing (>= 100), orange = NaN penalty, blue otherwise
colors = repmat([0.30 0.55 0.80], numel(c_sorted), 1);
for i = 1:numel(c_sorted)
    if c_sorted(i) >= 100
        colors(i, :) = [0.80 0.20 0.20];
    elseif isnan(c_sorted(i))
        colors(i, :) = [0.95 0.55 0.20];
    end
end

cla(ax); hold(ax, 'on');
y = 1:numel(c_sorted);
for i = 1:numel(c_sorted)
    barh(ax, y(i), c_sorted(i), 'FaceColor', colors(i, :), 'EdgeColor', 'none');
end
set(ax, 'YTick', y, 'YTickLabel', labels, 'YDir', 'reverse', ...
    'FontSize', 7, 'YLim', [0.5, numel(c_sorted)+0.5]);
xlabel(ax, 'cost', 'FontSize', 8);
title(ax, sprintf('Cost culprits (total %.2f)', sum(c_sorted, 'omitnan')), ...
    'FontSize', 9);
grid(ax, 'on');
box(ax, 'on');
hold(ax, 'off');
end


function plotBoundViolators(violations, ax)
% Horizontal bar chart of physiology-bound violators, sorted by penalty desc.
cla(ax);
if isempty(violations)
    text(ax, 0.5, 0.5, sprintf('No bound violations\n(physiology OK)'), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 9, ...
        'Color', [0.10 0.55 0.10]);
    title(ax, 'Bound violators', 'FontSize', 9);
    axis(ax, 'off');
    return;
end

[~, ord] = sort([violations.penalty], 'descend');
violations = violations(ord);

names = {violations.name};
pens  = [violations.penalty];
vals  = [violations.value];
lbs   = [violations.lb];
ubs   = [violations.ub];

% Build labels with the value and which bound was crossed
disp_labels = cell(numel(names), 1);
for i = 1:numel(names)
    if vals(i) < lbs(i)
        side = sprintf('<%.3g', lbs(i));
    else
        side = sprintf('>%.3g', ubs(i));
    end
    disp_labels{i} = sprintf('%s = %.3g  (%s)', names{i}, vals(i), side);
end

hold(ax, 'on');
y = 1:numel(pens);
for i = 1:numel(pens)
    barh(ax, y(i), pens(i), 'FaceColor', [0.85 0.30 0.30], 'EdgeColor', 'none');
end
set(ax, 'YTick', y, 'YTickLabel', disp_labels, 'YDir', 'reverse', ...
    'FontSize', 7, 'YLim', [0.5, numel(pens)+0.5]);
xlabel(ax, 'weighted penalty', 'FontSize', 8);
title(ax, sprintf('Bound violators (sum %.3g)', sum(pens)), 'FontSize', 9);
grid(ax, 'on');
box(ax, 'on');
hold(ax, 'off');
end
