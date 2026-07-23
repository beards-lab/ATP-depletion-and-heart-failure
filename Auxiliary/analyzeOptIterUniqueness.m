function out = analyzeOptIterUniqueness(src, opts)
%ANALYZEOPTITERUNIQUENESS  Are the optimiser's near-best solutions unique?
%
%   OUT = ANALYZEOPTITERUNIQUENESS(SRC) reads the param_opt_iter cloud produced
%   by optimizeFeatures (state.optIter -- one converged (params,cost) record per
%   round, plus a numbered writeParamsToMFile snapshot per improvement) and asks
%   the question the user cares about: among the solutions clustered around the
%   BEST cost, do the parameters CONVERGE BACK to (essentially) one point
%   [UNIQUE], or do DIFFERENT parameter vectors reach the SAME cost
%   [DEGENERATE / sloppy manifold]?
%
%   It answers this with three views + a printed verdict:
%     (1) RADAR  -- comparePUParams overlay of the near-optimal param sets
%                   (ratios vs the single best set). Tight overlap => unique;
%                   fanned-out polygons => degenerate.
%     (2) IDENTIFIABILITY -- per-parameter fractional spread (relative range and
%                   coefficient of variation) across the near-optimal sets,
%                   sorted. Pinned params (spread ~ 0) are IDENTIFIED; params
%                   that roam at ~constant cost are the SLOPPY directions.
%     (3) FEATURE INFLUENCE -- per-feature cost for each near-optimal solution,
%                   read from the vectors optimizeFeatures now STORES per record
%                   (rec.featCost) so NO re-simulation is needed. Three panels:
%                   cost heatmap (feature x solution), BETTER/WORSE vs the best
%                   (diverging: which features improved/degraded), and the mean
%                   per-feature cost (the MOST INFLUENCING features). Reveals
%                   whether solutions reach the same total cost via DIFFERENT
%                   feature trade-offs (a fingerprint of non-uniqueness). Legacy
%                   states without stored featCost fall back to re-simulation
%                   (opts.runFeatures).
%
%   SRC may be:
%     - a tag string        -> loads params/<tag>_state.mat
%     - a path to *_state.mat
%     - a state struct       (with field .optIter)
%
%   OPTS (all optional):
%     .relTol       near-optimal band = cost <= bestcost*(1+relTol). Default 0.15.
%     .absTol       if set, band = cost <= bestcost + absTol (overrides relTol).
%     .minSets      widen the band until at least this many sets qualify. Default 3.
%     .maxSets      cap plotted sets to the this-many lowest-cost. Default 12.
%     .improvedOnly use only IMPROVED records (converged incumbents). Default false
%                   (false keeps the near-miss stalls -- the strongest evidence of
%                   "same cost, different params").
%     .dedup        drop near-identical param vectors (rel. dist < .dedupTol).
%                   Default true. .dedupTol default 0.01 (1%).
%     .paramList    restrict radar/identifiability to these params. Default {} =
%                   the union of the params actually optimised (records.drawn),
%                   falling back to comparePUParams' curated list.
%     .runFeatures  allow RE-SIMULATION as a fallback for the feature panel when
%                   the state has NO stored per-feature costs (legacy runs).
%                   Default true. New runs store featCost, so the panel is drawn
%                   from stored data with no simulation regardless of this flag;
%                   it only matters for old states.
%     .maxFeatureSims cap re-simulations (lowest-cost first). Default 6.
%     .sloppyCV     rel-spread above which a param is called SLOPPY. Default 0.25.
%     .idCV         rel-spread below which a param is called IDENTIFIED. Default 0.05.
%
%   OUT struct: .records .cost .sets .names .values .relSpread .cv
%               .medPairSpread .sloppy .identified .verdict .featCost .featNames
%               .figs
%
%   Example (after a RunOptimFull run tagged 'optfull_FourthSlackTimebase'):
%     analyzeOptIterUniqueness('optfull_FourthSlackTimebase');
%     % fast, no re-simulation:
%     analyzeOptIterUniqueness('optfull_FourthSlackTimebase', struct('runFeatures',false));
%
%   See also: optimizeFeatures, comparePUParams, evalFeatureCost, plotFeatures

if nargin < 2 || isempty(opts); opts = struct(); end
def = struct('relTol',0.15, 'absTol',[], 'minSets',3, 'maxSets',12, ...
             'improvedOnly',false, 'dedup',true, 'dedupTol',0.01, ...
             'paramList',{{}}, 'runFeatures',true, 'maxFeatureSims',6, ...
             'sloppyCV',0.25, 'idCV',0.05);
opts = iFillDefaults(opts, def);

%% ---- resolve source -> optIter records --------------------------------------
state = iResolveState(src);
if ~isfield(state, 'optIter') || isempty(state.optIter)
    error('analyzeOptIterUniqueness:noCloud', ...
        'No state.optIter records found. Run optimizeFeatures first (it populates the param_opt_iter cloud).');
end
rec = state.optIter;
if opts.improvedOnly
    rec = rec([rec.isImproved]);
    if isempty(rec); error('analyzeOptIterUniqueness:noImproved','No improved records.'); end
end

costs = [rec.cost];
bestCost = min(costs);

%% ---- select the near-optimal band -------------------------------------------
if ~isempty(opts.absTol)
    thr = bestCost + opts.absTol;
else
    thr = bestCost * (1 + opts.relTol);
end
sel = find(costs <= thr);
% widen to at least minSets by taking the lowest-cost points
if numel(sel) < opts.minSets
    [~, ord] = sort(costs, 'ascend');
    sel = ord(1:min(opts.minSets, numel(costs)));
end
% sort selected by cost ascending (best first -> radar baseline = the best set)
[~, so] = sort(costs(sel), 'ascend');
sel = sel(so);
rec = rec(sel);
costs = costs(sel);

%% ---- resolve each set through getParams (defaults + '=' expressions) --------
nSet0 = numel(rec);
setsR = cell(1, nSet0);
for i = 1:nSet0
    try
        setsR{i} = getParams(rec(i).params, [], true, false);
    catch
        setsR{i} = rec(i).params;   % fall back to the raw stored struct
    end
end

%% ---- parameter list (what was actually optimised, else curated default) -----
paramList = opts.paramList;
if isempty(paramList)
    drawnAll = {};
    for i = 1:nSet0
        if isfield(rec, 'drawn') && ~isempty(rec(i).drawn); drawnAll = [drawnAll, rec(i).drawn]; end %#ok<AGROW>
    end
    paramList = unique(drawnAll, 'stable');
end
if isempty(paramList); paramList = {}; end   % empty -> comparePUParams default

%% ---- optional de-duplication of near-identical solutions --------------------
Vfull = iGatherValues(setsR, iEffectiveList(paramList, setsR));
if opts.dedup && nSet0 > 1
    keep = true(1, nSet0);
    med  = median(Vfull, 2, 'omitnan');
    for i = 2:nSet0
        for j = find(keep(1:i-1))
            d = iRelDist(Vfull(:,i), Vfull(:,j), med);
            if d < opts.dedupTol; keep(i) = false; break; end
        end
    end
    setsR = setsR(keep); rec = rec(keep); costs = costs(keep);
end
% honour maxSets (already cost-sorted, best first)
if numel(setsR) > opts.maxSets
    setsR = setsR(1:opts.maxSets); rec = rec(1:opts.maxSets); costs = costs(1:opts.maxSets);
end
nSet = numel(setsR);

%% ---- labels -----------------------------------------------------------------
labels = cell(1, nSet);
for i = 1:nSet
    km = ''; if isfield(rec,'isKick') && rec(i).isKick; km = ' K'; end
    rr = NaN; if isfield(rec,'round'); rr = rec(i).round; end
    labels{i} = sprintf('#%d r%g c=%.3f%s', i, rr, costs(i), km);
end
labels{1} = [labels{1} ' [best]'];

fprintf('\n=== analyzeOptIterUniqueness ===\n');
fprintf('near-optimal band: cost <= %.4f  (best %.4f)  -> %d solution(s)\n', thr, bestCost, nSet);

%% ---- VIEW 1: radar overlay (ratios vs best) ---------------------------------
figs = struct();
cp = comparePUParams(setsR, paramList, struct('norm','first', 'resolve',false, ...
        'labels',{labels}, 'showTable',false, ...
        'title',sprintf('param_opt_iter: %d near-best solutions', nSet)));
figs.radar = cp.fig;

%% ---- VIEW 2: identifiability (per-param spread across near-best sets) --------
names = iEffectiveList(paramList, setsR);
V     = iGatherValues(setsR, names);
[relSpread, cv] = deal(nan(numel(names),1));
for p = 1:numel(names)
    v = V(p,:); v = v(~isnan(v));
    if isempty(v); continue; end
    med = median(v); mu = mean(v);
    relSpread(p) = (max(v) - min(v)) / max(eps, abs(med));
    cv(p)        = std(v)          / max(eps, abs(mu));
end
[rs_sorted, ord] = sort(relSpread, 'descend', 'MissingPlacement','last');
names_sorted = names(ord);

figs.ident = figure('Color','w','Name','Parameter identifiability','Position',[120 90 720 760]);
ax = axes('Parent',figs.ident); hold(ax,'on');
y = 1:numel(names_sorted);
for i = 1:numel(y)
    v = rs_sorted(i);
    if isnan(v); continue; end
    if v >= opts.sloppyCV
        col = [0.82 0.22 0.22];        % SLOPPY (roams at ~const cost)
    elseif v <= opts.idCV
        col = [0.20 0.60 0.25];        % IDENTIFIED (pinned)
    else
        col = [0.40 0.50 0.75];        % partial
    end
    barh(ax, y(i), v, 'FaceColor', col, 'EdgeColor','none');
end
xline(ax, opts.idCV,     '--', 'identified', 'Color',[0.20 0.60 0.25], 'LabelVerticalAlignment','bottom','FontSize',7);
xline(ax, opts.sloppyCV, '--', 'sloppy',     'Color',[0.82 0.22 0.22], 'LabelVerticalAlignment','bottom','FontSize',7);
set(ax, 'YTick', y, 'YTickLabel', names_sorted, 'YDir','reverse', ...
    'FontSize',8, 'YLim',[0.5 numel(y)+0.5], 'TickLabelInterpreter','none');
xlabel(ax, 'relative spread across near-best solutions  (max-min)/|median|', 'FontSize',9);
title(ax, sprintf('Parameter identifiability  (%d solutions within %.1f%% of best cost)', ...
    nSet, 100*(thr/bestCost - 1)), 'FontSize',10);
grid(ax,'on'); box(ax,'on');

%% ---- verdict ----------------------------------------------------------------
medPairSpread = iMedianPairwiseSpread(V);
sloppy     = names(relSpread >= opts.sloppyCV & ~isnan(relSpread));
identified = names(relSpread <= opts.idCV);
if nSet < 2
    verdict = 'INCONCLUSIVE: only one near-optimal solution found (widen relTol / run longer).';
elseif medPairSpread < opts.idCV
    verdict = sprintf('UNIQUE: near-best solutions reconverge (median pairwise spread %.1f%%).', 100*medPairSpread);
elseif medPairSpread > opts.sloppyCV
    verdict = sprintf('DEGENERATE / NON-UNIQUE: different params reach ~the same cost (median pairwise spread %.1f%%).', 100*medPairSpread);
else
    verdict = sprintf('PARTIALLY IDENTIFIED: some directions pinned, some sloppy (median pairwise spread %.1f%%).', 100*medPairSpread);
end
fprintf('%s\n', verdict);
if ~isempty(sloppy);     fprintf('  sloppy (roams >%.0f%%): %s\n', 100*opts.sloppyCV, strjoin(sloppy, ', ')); end
if ~isempty(identified); fprintf('  identified (<%.0f%%): %d params\n', 100*opts.idCV, numel(identified)); end
% annotate the identifiability figure
annotation(figs.ident, 'textbox', [0.02 0.005 0.96 0.055], 'String', verdict, ...
    'EdgeColor','none', 'FontSize',9, 'FontWeight','bold', 'Interpreter','none', ...
    'Color', iVerdictColor(verdict), 'HorizontalAlignment','center');

%% ---- VIEW 3: feature influence (which features got better/worse) ------------
% Prefer the PER-FEATURE cost vectors stored per record by optimizeFeatures
% (rec.featCost + state.featNames) -- no re-simulation needed. Fall back to
% re-simulating only for legacy states that predate stored featCost.
[featCost, featNames] = iStoredFeatCost(rec, state);
featSrc = 'stored';
if isempty(featCost) && opts.runFeatures && nSet >= 1
    featSrc = 're-simulated';
    try
        [featCost, featNames] = iFeatureCosts(setsR, costs, opts);
    catch ME
        fprintf('  (feature panel skipped: %s)\n', ME.message);
    end
end
if ~isempty(featCost)
    fprintf('  feature panel: %d features x %d solutions (%s)\n', size(featCost,1), size(featCost,2), featSrc);
    figs.features = iPlotFeatureInfluence(featCost, featNames, labels);
elseif ~opts.runFeatures
    fprintf('  feature panel: no stored featCost in this state; set runFeatures=true to re-simulate.\n');
end

%% ---- output -----------------------------------------------------------------
out = struct('records',rec, 'cost',costs, 'sets',{setsR}, 'names',{names}, ...
             'values',V, 'relSpread',relSpread, 'cv',cv, ...
             'medPairSpread',medPairSpread, 'sloppy',{sloppy}, ...
             'identified',{identified}, 'verdict',verdict, ...
             'featCost',featCost, 'featNames',{featNames}, 'figs',figs);
end

% =============================================================================
function state = iResolveState(src)
if isstruct(src)
    state = src; return;
end
src = char(src);
if endsWith(src, '.mat') && exist(src, 'file')
    f = src;
else
    root = fullfile(fileparts(mfilename('fullpath')), '..');
    f = fullfile(root, 'params', [src '_state.mat']);
end
if ~exist(f, 'file')
    error('analyzeOptIterUniqueness:noState', 'State file not found: %s', f);
end
S = load(f);
if ~isfield(S, 'state'); error('analyzeOptIterUniqueness:badState','%s has no ''state'' variable.', f); end
state = S.state;
end

% =============================================================================
function names = iEffectiveList(paramList, setsR)
% If paramList is empty, borrow comparePUParams' curated default so the
% identifiability view has a sensible axis set even without records.drawn.
if ~isempty(paramList); names = paramList(:)'; return; end
try
    cp = comparePUParams(setsR, {}, struct('resolve',false,'showTable',false,'hideConstant',false));
    names = cp.names; close(cp.fig);
catch
    names = {'ka','kd','k1','k_1','k2','kstiff1','kstiff2','ksr0','kmsrd','A0','v_ref_reg','tau_reg'};
end
end

% =============================================================================
function V = iGatherValues(setsR, names)
% nParam x nSet matrix of resolved scalar values (handles arr__idx notation).
nP = numel(names); nSet = numel(setsR);
V = nan(nP, nSet);
for p = 1:nP
    for i = 1:nSet
        V(p,i) = iCurrentVal(names{p}, setsR{i});
    end
end
end

% =============================================================================
function v = iCurrentVal(name, p)
v = NaN;
if isfield(p, name)
    val = p.(name);
    if isnumeric(val) && isscalar(val); v = double(val); end
elseif contains(name, '__')
    t = split(name, '__');
    if isfield(p, t{1})
        arr = p.(t{1}); idx = str2double(t{2});
        if isnumeric(arr) && numel(arr) >= idx; v = double(arr(idx)); end
    end
end
end

% =============================================================================
function d = iRelDist(a, b, med)
% mean fractional distance between two param vectors, normalised per-param by med.
den = max(eps, abs(med));
m = ~isnan(a) & ~isnan(b);
if ~any(m); d = 0; return; end
d = mean(abs(a(m) - b(m)) ./ den(m));
end

% =============================================================================
function s = iMedianPairwiseSpread(V)
% median over all pairs of the mean per-param fractional distance.
nSet = size(V,2);
if nSet < 2; s = 0; return; end
med = median(V, 2, 'omitnan');
d = [];
for i = 1:nSet
    for j = i+1:nSet
        d(end+1) = iRelDist(V(:,i), V(:,j), med); %#ok<AGROW>
    end
end
s = median(d, 'omitnan');
end

% =============================================================================
function [C, names] = iStoredFeatCost(rec, state)
% Assemble the per-feature cost matrix (nFeat x nSet) from the values stored in
% the optIter records (rec.featCost, labelled by state.featNames). Rows that are
% structurally missing (>= the evalFeatureCost missing-sentinel across ALL
% solutions -- e.g. experiments not run) or all-NaN are dropped. Returns empty
% if the state predates stored featCost (legacy).
C = []; names = {};
if ~isfield(state, 'featNames') || isempty(state.featNames); return; end
if ~isfield(rec, 'featCost'); return; end
names0 = state.featNames(:).';
nF = numel(names0);
cols = cell(1, numel(rec));
for i = 1:numel(rec)
    fc = rec(i).featCost;
    if numel(fc) == nF && ~all(isnan(fc)); cols{i} = fc(:); else; cols{i} = nan(nF,1); end
end
C = cell2mat(cols);          % nF x nSet
MISS = 100;                  % evalFeatureCost MISSING_FEATURE_COST sentinel
keep = true(nF,1);
for r = 1:nF
    rv = C(r, ~isnan(C(r,:)));
    if isempty(rv) || min(rv) >= MISS; keep(r) = false; end
end
C = C(keep,:); names = names0(keep);
if isempty(C) || all(isnan(C(:))); C = []; names = {}; end
end

% =============================================================================
function [C, keptNames] = iFeatureCosts(setsR, costs, opts)
% Re-simulate each near-optimal set and return per-feature cost (nFeat x nSet).
% Only features structurally present in the simulation output are kept.
nRun = min(opts.maxFeatureSims, numel(setsR));
fprintf('  feature panel: simulating %d solution(s) (set opts.runFeatures=false to skip)...\n', nRun);

% LoadData + RunBakersExp are scripts (as in evaluateBakersExp/captureOptimIter):
% they run in this function's workspace and yield features_data/features_model.
LoadData;

fnRef = {};
if isfield(setsR{1}, 'fn') && ~isempty(setsR{1}.fn); fnRef = setsR{1}.fn; end

Cc = {};
fdRef = []; fnKept = {}; keepMask = [];
for i = 1:nRun
    params0 = setsR{i};
    params0.PlotEachSeparately = false;
    params0.PlotFeatureFitting = false;
    params0.mods = {}; params0.g = [];
    fn_i = fnRef; if isfield(params0,'fn') && ~isempty(params0.fn); fn_i = params0.fn; end
    try
        lastwarn('','');
        RunBakersExp;                   % script: yields features_data / features_model
        if isempty(fdRef); fdRef = features_data; end
        [~, ~, craw] = evalFeatureCost(fdRef, features_model, fn_i);
        if isempty(fnKept)
            present = arrayfun(@(k) iFeatPresent(fn_i{k}, features_model), 1:numel(fn_i));
            fnKept  = fn_i(present);
            keepMask = present;
        end
        Cc{end+1} = craw(keepMask(:).'); %#ok<AGROW>
        fprintf('    set #%d (c=%.3f): ok\n', i, costs(i));
    catch ME
        fprintf('    set #%d (c=%.3f): sim failed (%s)\n', i, costs(i), ME.message);
    end
end
if isempty(Cc); C = []; keptNames = {}; return; end
C = cell2mat(Cc(:).');   % nFeat x nOK
keptNames = cellfun(@(s) regexprep(s, '[|,].*$',''), fnKept, 'UniformOutput', false);
end

% =============================================================================
function tf = iFeatPresent(fnStr, feats)
% Does the primary field of an fn selector exist in the sim feature struct?
tok = regexp(fnStr, '^[^|,]*', 'match', 'once');   % up to first | or ,
name = regexprep(tok, '\[.*$', '');                % strip [idx] selector
name = strtrim(name);
tf = isfield(feats, name);
end

% =============================================================================
function fh = iPlotFeatureInfluence(C, featNames, labels)
% Three views of the stored per-feature cost (features sorted by mean cost):
%   (1) cost heatmap (feature x solution)
%   (2) BETTER/WORSE vs the best solution -- diverging: blue = lower cost
%       (better than best), red = higher cost (worse than best). Directly
%       answers "which features got better and which worse".
%   (3) most influencing features -- mean per-feature cost, whisker = spread.
nOK  = size(C,2);
meanC   = mean(C, 2, 'omitnan');
spreadC = max(C,[],2) - min(C,[],2);
[~, ord] = sort(meanC, 'descend', 'MissingPlacement','last');
C = C(ord,:); featNames = featNames(ord); meanC = meanC(ord); spreadC = spreadC(ord);
y = 1:numel(featNames);

hasDelta = nOK >= 2;
D = C - C(:,1);                       % vs the best (column 1 = lowest cost)

fh = figure('Color','w','Name','Feature influence / trade-offs','Position',[130 70 1240 760]);
nTiles = 2 + hasDelta;
tl = tiledlayout(fh, 1, nTiles, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf('Per-feature cost across %d near-best solutions', nOK), 'FontWeight','bold');

% (1) cost heatmap
ax1 = nexttile(tl);
imagesc(ax1, C, 'AlphaData', ~isnan(C)); colormap(ax1, flipud(hot));
cb = colorbar(ax1); cb.Label.String = 'feature cost';
set(ax1, 'YTick',y, 'YTickLabel',featNames, 'TickLabelInterpreter','none', 'FontSize',8);
set(ax1, 'XTick',1:nOK, 'XTickLabel',labels(1:nOK), 'XTickLabelRotation',35);
title(ax1, 'cost heatmap (feature \times solution)'); xlabel(ax1, 'solution');

% (2) better/worse vs best
if hasDelta
    ax2 = nexttile(tl);
    imagesc(ax2, D, 'AlphaData', ~isnan(D));
    colormap(ax2, iDivergingMap()); cmax = max(abs(D(~isnan(D)))); if isempty(cmax)||cmax==0; cmax=1; end
    set(ax2, 'CLim', [-cmax cmax]);
    cb2 = colorbar(ax2); cb2.Label.String = '\Delta cost vs best';
    set(ax2, 'YTick',y, 'YTickLabel',featNames, 'TickLabelInterpreter','none', 'FontSize',8);
    set(ax2, 'XTick',1:nOK, 'XTickLabel',labels(1:nOK), 'XTickLabelRotation',35);
    title(ax2, 'better(blue) / worse(red) vs best'); xlabel(ax2, 'solution');
end

% (3) most influencing features
ax3 = nexttile(tl); hold(ax3,'on');
barh(ax3, y, meanC, 'FaceColor',[0.30 0.45 0.75], 'EdgeColor','none');
for i = 1:numel(y)
    plot(ax3, [meanC(i)-spreadC(i)/2, meanC(i)+spreadC(i)/2], [y(i) y(i)], '-', ...
        'Color',[0.85 0.35 0.10], 'LineWidth',1.3);
end
set(ax3, 'YTick',y, 'YTickLabel',featNames, 'YDir','reverse', 'TickLabelInterpreter','none', ...
    'FontSize',8, 'YLim',[0.5 numel(y)+0.5]);
xlabel(ax3, 'mean feature cost  (whisker = spread)');
title(ax3, 'most influencing features'); grid(ax3,'on'); box(ax3,'on');
end

% =============================================================================
function m = iDivergingMap()
% Blue (negative / better) - white (0) - red (positive / worse).
n = 128; t = linspace(0,1,n)';
lo = [0.13 0.40 0.67]; mid = [1 1 1]; hi = [0.78 0.20 0.20];
top = (1-t).*lo + t.*mid;   bot = (1-t).*mid + t.*hi;
m = [top; bot];
end

% =============================================================================
function c = iVerdictColor(v)
if     startsWith(v, 'UNIQUE');     c = [0.20 0.55 0.25];
elseif startsWith(v, 'DEGENERATE'); c = [0.80 0.20 0.20];
else;                               c = [0.55 0.45 0.10];
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
