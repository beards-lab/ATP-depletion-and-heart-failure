function plotPhysiologyDashboard(params, bounds, ax)
%PLOTPHYSIOLOGYDASHBOARD  Bar visualisation of param values vs physiology bounds.
%
%   PLOTPHYSIOLOGYDASHBOARD(PARAMS) uses default bounds from parameterBounds().
%   PLOTPHYSIOLOGYDASHBOARD(PARAMS, BOUNDS, AX) lets you specify both.
%
%   Each bounded parameter (weight > 0) is shown as one row:
%       - Grey bar  : bound range [lb, ub], normalised to [0, 1]
%       - Marker    : current value position
%           * Green : inside bounds
%           * Red   : outside bounds; position uses log10 fold violation:
%                       below lb → x = -log10(lb / val)
%                       above ub → x = 1 + log10(val / ub)
%             so a 10× violation past ub lands at x = 2, a 100× at x = 3, etc.
%
%   The log-fold encoding keeps extreme violators on screen without
%   compressing the in-bounds range.
%
%   See also: parameterBounds, evalPhysiologyCost

if nargin < 2 || isempty(bounds)
    bounds = parameterBounds();
end
if nargin < 3 || isempty(ax)
    ax = gca;
end

names = fieldnames(bounds);

% Filter to bounds with weight > 0 and a present, finite param value
keep = false(size(names));
vals = nan(size(names));
for i = 1:numel(names)
    b = bounds.(names{i});
    if b.weight == 0; continue; end
    v = resolveParam(params, names{i});
    if isnan(v); continue; end
    keep(i) = true;
    vals(i) = v;
end
names = names(keep);
vals  = vals(keep);

if isempty(names)
    cla(ax);
    text(ax, 0.5, 0.5, 'No physiology bounds active', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    return;
end

% Sort by penalty weight (highest first), so Tier A appears at the top
weights = arrayfun(@(i) bounds.(names{i}).weight, 1:numel(names));
[~, ord] = sort(weights, 'descend');
names = names(ord);
vals  = vals(ord);

n = numel(names);
y = (1:n)';

% Compute marker positions and out-of-bounds magnitudes
xpos    = nan(n, 1);
out_of  = false(n, 1);
for i = 1:n
    b = bounds.(names{i});
    span = max(b.ub - b.lb, eps);

    if vals(i) < b.lb
        % Below lb. Use log-fold if both positive, else linear fallback.
        if b.lb > 0 && vals(i) > 0
            xpos(i) = -log10(b.lb / vals(i));
        else
            xpos(i) = (vals(i) - b.lb) / span;  % linear (negative)
        end
        out_of(i) = true;
    elseif vals(i) > b.ub
        if b.ub > 0 && vals(i) > 0
            xpos(i) = 1 + log10(vals(i) / b.ub);
        else
            xpos(i) = (vals(i) - b.lb) / span;
        end
        out_of(i) = true;
    else
        % Inside bounds: linear normalised position in [0, 1]
        xpos(i) = (vals(i) - b.lb) / span;
    end
end

% Adapt xlim to include all violators, with reasonable padding
xmin = min([xpos; 0]);
xmax = max([xpos; 1]);
pad  = 0.2;
xlim_ = [min(-pad, xmin - pad), max(1 + pad, xmax + pad)];

cla(ax); hold(ax, 'on');

for i = 1:n
    b = bounds.(names{i});

    % Bound range bar (grey) at normalised [0, 1]
    rectangle(ax, 'Position', [0, y(i)-0.35, 1, 0.7], ...
        'FaceColor', [0.85 0.85 0.85], 'EdgeColor', [0.6 0.6 0.6]);

    if out_of(i)
        marker_color = [0.85 0.10 0.10];
    else
        marker_color = [0.10 0.55 0.10];
    end

    plot(ax, xpos(i), y(i), 'o', ...
        'MarkerFaceColor', marker_color, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 7);

    % Annotate value to the right of the row
    if out_of(i)
        suffix = sprintf('  ! %.2gx', fold_violation(vals(i), b.lb, b.ub));
    else
        suffix = '';
    end
    text(ax, xlim_(2) - 0.02, y(i), ...
        sprintf('%.3g (w=%g)%s', vals(i), b.weight, suffix), ...
        'FontSize', 7, 'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'right');
end

set(ax, 'YTick', y, 'YTickLabel', names, 'YDir', 'reverse', ...
    'FontSize', 8, 'XLim', xlim_, 'YLim', [0.5 n+0.5]);
xlabel(ax, 'Position in [lb, ub] (log-fold past bounds)');
title(ax, sprintf('Physiology bounds (%d params, %d violations)', n, sum(out_of)), ...
    'FontSize', 9);
xline(ax, 0, '--k', 'lb', 'LabelVerticalAlignment', 'top', 'FontSize', 7);
xline(ax, 1, '--k', 'ub', 'LabelVerticalAlignment', 'top', 'FontSize', 7);

% Decade gridlines for the log-fold region
for d = 1:floor(xlim_(2))
    if d > 1
        xline(ax, d, ':', sprintf('%dx', 10^(d-1)), 'Color', [0.7 0.5 0.5], ...
            'LabelVerticalAlignment', 'top', 'FontSize', 6);
    end
end
for d = 1:abs(floor(xlim_(1)))
    if d >= 1
        xline(ax, -d, ':', sprintf('1/%dx', 10^d), 'Color', [0.5 0.5 0.7], ...
            'LabelVerticalAlignment', 'top', 'FontSize', 6);
    end
end
box(ax, 'on');
hold(ax, 'off');
end


function fv = fold_violation(val, lb, ub)
% Magnitude of bound violation, expressed as a fold factor.
if val < lb
    fv = lb / max(val, eps);
elseif val > ub
    fv = val / max(ub, eps);
else
    fv = 1;
end
end


function val = resolveParam(params, name)
val = NaN;
if isfield(params, name) && isnumeric(params.(name)) && isscalar(params.(name))
    val = params.(name);
    return;
end
sep = strfind(name, '__');
if ~isempty(sep)
    base = name(1:sep(1)-1);
    idx  = str2double(name(sep(1)+2:end));
    if isfield(params, base) && ~isnan(idx) && isnumeric(params.(base)) ...
            && idx >= 1 && idx <= numel(params.(base))
        val = params.(base)(idx);
    end
end
end
