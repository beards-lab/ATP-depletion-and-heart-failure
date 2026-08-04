function out = comparePUParams(paramSets, paramList, opts)
%COMPAREPUPARAMS  Radar (spider) comparison of cross-bridge parameters across paramsets.
%
%   OUT = COMPAREPUPARAMS(PARAMSETS) draws a circular "fingerprint" overlay of
%   the active kinetic/mechanical parameters read by dPUdT_CombinedTransitions,
%   one coloured polygon per paramset, so several parameter snapshots can be
%   compared at a glance.
%
%   OUT = COMPAREPUPARAMS(PARAMSETS, PARAMLIST) restricts the plot to the named
%   parameters (cellstr). Pass [] or {} to use the curated default list of the
%   active dxdt parameters (grouped: attach/detach, strain-dep, force offsets,
%   stiffness, SRX, mech/passive, reg-availability). See iDefaultList below.
%
%   OUT = COMPAREPUPARAMS(PARAMSETS, PARAMLIST, OPTS) sets options:
%     .norm         'minmax' (default; per-axis min-max across sets, any sign),
%                   'first'  (ratio to the FIRST set: set 1 is the unit circle,
%                            radius = value/value_set1, rings marked 0.5x/1x/2x...,
%                            and value_set1 is printed next to each axis label),
%                   'log'    (shared radial log10 scale, magnitude-aware, +ve only),
%                   'none'   (raw values / global max|.|).
%     .resolve      true (default): run each set through getParams so defaults are
%                   filled and '='-expression fields are resolved to numbers.
%     .hideConstant true (default): when >1 set, drop axes identical across all
%                   sets (keeps the plot focused on what differs). Ignored for 1 set.
%     .labels       cellstr legend names, one per set (default 'set 1'..'set N').
%     .title        figure title.
%     .showTable    true (default): also print a name x set value table to console.
%
%   Inputs:
%     PARAMSETS  - struct array OR cell array of params0 structs (as from
%                  loadParams / getParams). Fields may differ between sets.
%
%   Output OUT (struct):
%     .names   parameters actually plotted (after drop rules)
%     .values  raw resolved values, nParam x nSet
%     .R       radial coordinates used, nParam x nSet
%     .table   the console table (also raw values)
%     .fig,.ax figure/axes handles
%     .dropped / .omitted   names removed as constant / non-plottable
%
%   Example:
%     A = loadParams('opt2state_v2_opt');
%     B = loadParams('params_restretch_best');
%     comparePUParams({A,B}, [], struct('labels',{{'v2opt','restretch'}}));
%
%   See also: getParams, loadParams, dPUdT_CombinedTransitions, tunableParams

if nargin < 2, paramList = {}; end
if nargin < 3 || isempty(opts), opts = struct(); end

% ---- options / defaults ------------------------------------------------------
def = struct('norm','minmax','resolve',true,'hideConstant',true, ...
             'title','Cross-bridge parameter comparison (dxdt)', ...
             'rInner',0.12,'showTable',true);
opts = iFillDefaults(opts, def);
if ~isfield(opts,'labels'), opts.labels = {}; end

% ---- normalise paramSets to a cell array of structs --------------------------
if isstruct(paramSets)
    paramSets = num2cell(paramSets);
elseif ~iscell(paramSets)
    error('comparePUParams:badInput', ...
        'paramSets must be a struct array or a cell array of structs.');
end
nSet = numel(paramSets);
if nSet == 0, error('comparePUParams:noSets','No paramsets given.'); end

% ---- resolve each set through getParams (defaults + '=' expressions) ---------
if opts.resolve
    for i = 1:nSet
        try
            paramSets{i} = getParams(paramSets{i}, [], true, false);
        catch ME
            warning('comparePUParams:resolveFailed', ...
                'getParams failed for set %d (%s); using raw struct.', i, ME.message);
        end
    end
end

% ---- parameter list + grouping ----------------------------------------------
if isempty(paramList)
    [paramList, groupOf, groupNames] = iDefaultList();
else
    if ischar(paramList), paramList = cellstr(paramList); end
    groupOf = ones(1, numel(paramList));
    groupNames = {''};
end
paramList = paramList(:)';
groupOf   = groupOf(:)';

% ---- gather values: nP x nSet -----------------------------------------------
nP0 = numel(paramList);
V = nan(nP0, nSet);
for p = 1:nP0
    for i = 1:nSet
        if isfield(paramSets{i}, paramList{p})
            v = paramSets{i}.(paramList{p});
            if isnumeric(v) && isscalar(v), V(p,i) = double(v); end
        end
    end
end

% ---- drop axes constant across sets (comparison-focused; only if >1 set) -----
dropped = {};
if opts.hideConstant && nSet > 1
    keep = true(1, nP0);
    for p = 1:nP0
        vf = V(p, ~isnan(V(p,:)));
        if isempty(vf) || (max(vf) - min(vf)) <= 1e-6*max(1, max(abs(vf)))
            keep(p) = false;
        end
    end
    dropped = paramList(~keep);
    paramList = paramList(keep); groupOf = groupOf(keep); V = V(keep,:);
end

% ---- radial normalisation ----------------------------------------------------
[R, gridSpec, axisNote] = iNormalize(V, opts);

% ---- drop axes that cannot be placed (missing; <=0 in log; sign-flip in first)-
badRow = any(isnan(R), 2);
omitted = paramList(badRow);
paramList = paramList(~badRow); groupOf = groupOf(~badRow);
V = V(~badRow,:); R = R(~badRow,:); axisNote = axisNote(~badRow);
nP = numel(paramList);

% ---- legend labels -----------------------------------------------------------
labels = opts.labels;
if isempty(labels)
    labels = arrayfun(@(k) sprintf('set %d',k), 1:nSet, 'UniformOutput', false);
end
labels = cellstr(labels); labels = labels(:)';

% ---- draw --------------------------------------------------------------------
fig = figure('Color','w','Position',[100 100 940 820]);
ax  = axes('Parent',fig); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');

if nP >= 1
    ang = pi/2 - 2*pi*(0:nP-1)/nP;          % start at top, go clockwise
    iDrawGrid(ax, ang, paramList, groupOf, groupNames, gridSpec, axisNote);

    cols = lines(nSet);
    hSet = gobjects(1, nSet);
    for i = 1:nSet
        r  = R(:,i);
        xc = [r.*cos(ang(:)); r(1)*cos(ang(1))];
        yc = [r.*sin(ang(:)); r(1)*sin(ang(1))];
        patch('Parent',ax,'XData',xc,'YData',yc, ...
              'FaceColor',cols(i,:),'FaceAlpha',0.10, ...
              'EdgeColor',cols(i,:),'LineWidth',1.8);
        hSet(i) = plot(ax, xc, yc, '-o', 'Color',cols(i,:), ...
              'MarkerFaceColor',cols(i,:),'MarkerSize',4,'LineWidth',1.8);
    end
    legend(hSet, labels, 'Location','northeastoutside','Box','off','Interpreter','none');
else
    text(ax, 0, 0, 'No differing parameters to plot', ...
        'HorizontalAlignment','center');
end

if strcmpi(opts.norm, 'first') && ~isempty(labels)
    normTag = sprintf('  [ratios vs %s]', labels{1});
else
    normTag = sprintf('  [norm: %s]', lower(opts.norm));
end
title(ax, [opts.title normTag], 'FontWeight','bold','Interpreter','none');

% ---- console table -----------------------------------------------------------
T = table();
if nP >= 1
    T = array2table(V, ...
        'RowNames', matlab.lang.makeValidName(paramList), ...
        'VariableNames', matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(labels)));
    if opts.showTable
        fprintf('\n=== comparePUParams: %d params x %d sets ===\n', nP, nSet);
        disp(T);
        if ~isempty(dropped)
            fprintf('Hidden (identical across sets): %s\n', strjoin(dropped, ', '));
        end
        if ~isempty(omitted)
            fprintf('Omitted (missing, or <=0 under log norm): %s\n', strjoin(omitted, ', '));
        end
    end
end

out = struct('names',{paramList}, 'values',V, 'R',R, ...
             'table',{T}, 'fig',fig, 'ax',ax, ...
             'dropped',{dropped}, 'omitted',{omitted});
end

% =============================================================================
function [R, gridSpec, axisNote] = iNormalize(V, opts)
%INORMALIZE  Map values to a radial coordinate per opts.norm. NaN => not placed.
%   Also returns GRIDSPEC (radial rings/labels for iDrawGrid) and AXISNOTE (a
%   per-axis annotation string appended to the param label; used by 'first').
[nP, nSet] = size(V);
R = nan(nP, nSet);
axisNote = repmat({''}, 1, nP);
gridSpec = struct('rings',[0.25 0.5 0.75 1.0], ...
                  'labels',{repmat({''},1,4)}, 'refIdx',NaN);

switch lower(opts.norm)
    case 'first'   % ratio to the first set: set 1 -> unit circle (all r=1/rMax)
        fv = V(:,1);
        ratio = nan(nP, nSet);
        for p = 1:nP
            if fv(p) == 0 || ~isfinite(fv(p)), continue; end  % no baseline -> omit
            rr = V(p,:) / fv(p);
            if any(rr(~isnan(rr)) < 0), continue; end         % sign flip -> omit
            ratio(p,:) = rr;
            axisNote{p} = iFmt(fv(p));                         % the scaling value
        end
        good = ratio(~isnan(ratio));
        maxRatio = 1; if ~isempty(good), maxRatio = max(good); end
        rMax = 1.08 * max(1, maxRatio);                        % headroom past outer pt
        R = ratio / rMax;
        cand = [0.25 0.5 1 2 3 4 5 6 8 10 15 20];
        cand = cand(cand <= maxRatio*1.0001);
        if ~any(cand == 1), cand = sort([cand 1]); end
        gridSpec.rings  = cand / rMax;
        gridSpec.labels = arrayfun(@(x) [num2str(x) char(215)], cand, ...
                                   'UniformOutput', false);    % e.g. '2x'
        gridSpec.refIdx = find(cand == 1, 1);
    case 'log'     % shared radial log10 scale across all positive values
        L = V; L(L <= 0) = NaN; L = log10(L);
        lo = min(L(:)); hi = max(L(:));
        if isempty(lo) || ~isfinite(lo) || hi <= lo
            R(~isnan(L)) = 0.6;
        else
            R = opts.rInner + (1-opts.rInner)*(L - lo)/(hi - lo);
        end
    case 'none'    % raw, scaled only by the global largest magnitude
        mx = max(abs(V(:))); if isempty(mx) || mx == 0, mx = 1; end
        R = V/mx;
    otherwise      % 'minmax' per axis across sets (robust to any sign)
        for p = 1:nP
            v = V(p,:); vf = v(~isnan(v));
            if isempty(vf), continue; end
            lo = min(vf); hi = max(vf);
            if hi <= lo
                R(p, ~isnan(v)) = 0.6;                 % constant axis -> mid ring
            else
                R(p,:) = opts.rInner + (1-opts.rInner)*(v - lo)/(hi - lo);
            end
        end
end
end

% =============================================================================
function iDrawGrid(ax, ang, names, groupOf, groupNames, gridSpec, axisNote)
%IDRAWGRID  Concentric rings (per gridSpec), spokes, rotated param labels
%   (with optional axisNote appended), and group labels.
th = linspace(0, 2*pi, 200);
nP = numel(ang);
% radial-tick label angle: the gap just clockwise of the first spoke (no overlap)
if nP >= 2, angGap = ang(1) - pi/nP; else, angGap = ang(1) - 0.15; end
for k = 1:numel(gridSpec.rings)
    rr = gridSpec.rings(k);
    isRef = ~isnan(gridSpec.refIdx) && k == gridSpec.refIdx;
    if isRef
        plot(ax, rr*cos(th), rr*sin(th), '--', 'Color',[0.45 0.45 0.55], 'LineWidth',1.1);
    else
        plot(ax, rr*cos(th), rr*sin(th), '-', 'Color',[0.85 0.85 0.85], 'LineWidth',0.5);
    end
    if ~isempty(gridSpec.labels) && ~isempty(gridSpec.labels{k})
        text(ax, rr*cos(angGap), rr*sin(angGap), gridSpec.labels{k}, ...
            'FontSize',7, 'Color',[0.35 0.35 0.40], ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'BackgroundColor','w', 'Margin',0.5);
    end
end
for p = 1:nP
    a = ang(p);
    plot(ax, [0 cos(a)], [0 sin(a)], '-', 'Color',[0.90 0.90 0.90], 'LineWidth',0.5);
    [deg, align] = iLabelRot(a);
    lab = names{p};
    if ~isempty(axisNote) && numel(axisNote) >= p && ~isempty(axisNote{p})
        lab = sprintf('%s (%s)', names{p}, axisNote{p});
    end
    text(ax, 1.04*cos(a), 1.04*sin(a), lab, 'Rotation',deg, ...
        'HorizontalAlignment',align, 'VerticalAlignment','middle', ...
        'FontSize',8, 'Interpreter','none');
end
% group labels (only when grouping info is meaningful)
if ~isempty(groupNames) && ~isempty(groupNames{1}) && any(groupOf > 0)
    for g = 1:max(groupOf)
        idx = find(groupOf == g);
        if isempty(idx), continue; end
        am = mean(ang(idx));
        [deg, align] = iLabelRot(am);
        text(ax, 1.52*cos(am), 1.52*sin(am), groupNames{g}, 'Rotation',deg, ...
            'HorizontalAlignment',align, 'VerticalAlignment','middle', ...
            'FontSize',9, 'FontWeight','bold', 'Color',[0.20 0.20 0.55], ...
            'Interpreter','none');
    end
end
set(ax, 'XLim',[-1.9 1.9], 'YLim',[-1.9 1.9]);
end

% =============================================================================
function s = iFmt(v)
%IFMT  Compact numeric label for a parameter value.
if ~isfinite(v), s = 'NaN'; return; end
if v == 0, s = '0'; return; end
a = abs(v);
if a >= 1e4 || a < 1e-3
    s = sprintf('%.3g', v);   % e.g. 2.94e+04, 1.13e-03
else
    s = sprintf('%.4g', v);   % e.g. 215.7, 95.03, 0.01127
end
end

% =============================================================================
function [deg, align] = iLabelRot(a)
%ILABELROT  Radial text rotation that stays upright on the left half.
deg = a*180/pi; align = 'left';
if cos(a) < -1e-9
    deg = deg + 180; align = 'right';
end
end

% =============================================================================
function [names, groupOf, groupNames] = iDefaultList()
%IDEFAULTLIST  Curated active dxdt params, grouped by role in the XB cycle.
G = {
 'Attach/Detach', {'ka','kd','k1','k_1','k2','k3','kah','kamh'}
 'Strain-dep',    {'alpha1','alpha2','alpha2_L','alpha2_R','alpha3','s3','estiff'}
 'Force offset',  {'dr','dr1','dr2','dr3'}
 'Stiffness',     {'kstiff1','kstiff2'}
 'SRX',           {'ksr0','kmsr','kmsrd','ksrd','ksr2srd','ksrd2sr','sigma1','sigma2'}
 'Mech/Passive',  {'kSE','mu','k_pas','gamma','Lsc0'}
 'RegAvail',      {'A0','v_ref_reg','tau_reg'}
};
names = {}; groupOf = []; groupNames = G(:,1)';
for g = 1:size(G,1)
    ps = G{g,2};
    names   = [names, ps];                     %#ok<AGROW>
    groupOf = [groupOf, g*ones(1, numel(ps))]; %#ok<AGROW>
end
end

% =============================================================================
function s = iFillDefaults(s, d)
f = fieldnames(d);
for k = 1:numel(f)
    if ~isfield(s, f{k}) || isempty(s.(f{k}))
        s.(f{k}) = d.(f{k});
    end
end
end
