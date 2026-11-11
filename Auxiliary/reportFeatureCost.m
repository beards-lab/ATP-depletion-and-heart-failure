function reportFeatureCost(features_data, features_model, fn, costExp)
%REPORTFEATURECOST  Detailed cost + physiology breakdown for manual fitting.
%
%   reportFeatureCost(FEATURES_DATA, FEATURES_MODEL, FN) prints three blocks:
%
%     (1) cost_vec breakdown — every fn entry with its unweighted cost,
%         weight, and weighted cost, sorted so the dominant cost drivers are
%         at the top. Mirrors what plotFeatures computes (costExp=2).
%
%     (2) OUTPUT-bound group expansion — the grouped 'field|lb-ub|w,...' entry
%         broken into individual sub-bounds (value, [lb,ub], weight, weighted
%         violation, OUT flag). Bounds are read straight from FN, so this can
%         never drift from boundedOutputFn.
%
%     (3) raw features: data vs model — a table of every feature in either
%         struct, so the experimental physiology (features_data) can be read
%         directly alongside the simulated values.
%
%   COSTEXP (default 2) matches plotFeatures / the driver GRAND TOTAL. Pass 1
%   to match the E_feat term assembled inside RunBakersExp.
%
%   See also: evalFeatureCost, boundedOutputFn, plotFeatures, DriverBoundedFit.

if nargin < 4 || isempty(costExp); costExp = 2; end

[cv, wv, cvu] = evalFeatureCost(features_data, features_model, fn, costExp);

%% (1) per-fn-entry breakdown, sorted by weighted cost
fprintf('\n===== cost_vec breakdown (%d entries, costExp=%d, sorted by w*cost) =====\n', ...
        numel(fn), costExp);
fprintf('  %-46s %10s %5s %10s\n', 'fn entry', 'cost', 'w', 'w*cost');
[~, ord] = sort(cv, 'descend', 'MissingPlacement', 'last');
for k = ord(:)'
    lbl = fn{k};
    if numel(lbl) > 46; lbl = [lbl(1:43) '...']; end
    fprintf('  %-46s %10.4g %5g %10.4g\n', lbl, cvu(k), wv(k), cv(k));
end
fprintf('  %-46s %10s %5s %10.4g\n', '(total)', '', '', sum(cv, 'omitnan'));

%% (2) expand grouped OUTPUT-bound entry into sub-bounds (bounds read from fn)
gi = find(cellfun(@hasTopLevelComma, fn), 1);
if ~isempty(gi)
    fprintf('\n----- OUTPUT-bound group fn{%d}: per-output cost (sums to %.4g) -----\n', ...
            gi, cv(gi));
    fprintf('  %-13s %10s   %-16s %4s %10s\n', 'output', 'value', '[lb, ub]', 'w', 'w*viol');
    parts = splitTopLevelComma(fn{gi});
    for p = 1:numel(parts)
        toks = split(strtrim(parts{p}), '|');
        toks = toks(~cellfun(@isempty, toks));
        if isempty(toks); continue; end
        [feat, sel] = parseFeatSelector(toks{1});
        w = 1; rng = '';
        for t = 2:numel(toks)
            if isRangeToken(toks{t})
                rng = toks{t};
            else
                wtmp = str2double(toks{t});
                if ~isnan(wtmp); w = wtmp; end
            end
        end
        if isfield(features_model, feat)
            val = applyFeatSelector(double([features_model.(feat)]), sel);
        else
            val = NaN;
        end
        if isempty(rng); continue; end
        [lb, ub] = parseRange(rng);
        viol = max(0, (lb - val) / lb) .^ 2 + max(0, (val - ub) / ub) .^ 2;
        flag = ''; if val < lb || val > ub; flag = '  <-- OUT'; end
        fprintf('  %-13s %10.4g   [%6.3g, %6.3g] %4g %10.4g%s\n', ...
                feat, val, lb, ub, w, w * viol, flag);
    end
end

%% (3) raw data vs model dump for every feature
fprintf('\n----- raw features: data vs model -----\n');
allf = unique([fieldnames(features_data); fieldnames(features_model)], 'stable');
fprintf('  %-18s %-28s %-28s\n', 'feature', 'data', 'model');
for i = 1:numel(allf)
    f = allf{i};
    fprintf('  %-18s %-28s %-28s\n', f, ...
            fmtvec(getFieldSafe(features_data, f)), ...
            fmtvec(getFieldSafe(features_model, f)));
end
fprintf('\n');
end

% -----------------------------------------------------------------------
function c = getFieldSafe(s, f)
% Return per-segment values as a cell (robust to non-numeric fields, e.g. cfit).
if isfield(s, f); c = {s.(f)}; else; c = {}; end
end

function str = fmtvec(c)
if isempty(c); str = '(none)'; return; end
if ~(isnumeric(c{1}) || islogical(c{1}))
    str = sprintf('<%s>', class(c{1})); return;   % e.g. cfit object — skip
end
try
    v = double([c{:}]);
catch
    str = sprintf('<%s>', class(c{1})); return;
end
if isempty(v); str = '(none)'; return; end
if isscalar(v); str = sprintf('%.4g', v); return; end
if numel(v) > 6; v = v(1:6); suf = ' ..'; else; suf = ''; end
str = ['[' strjoin(arrayfun(@(x) sprintf('%.3g', x), v, 'UniformOutput', 0), ' ') ']' suf];
end
