function [val, dat, sem, lab, skipped] = batteryVector(fm, fd, fn)
% batteryVector  The slack-battery features as one vector, with the standard
% error of each element taken from the objective's own noise floors.
%
%   val, dat   model and data value of every element, each divided by the
%              entry's data mean (the normalisation evalFeatureCost uses)
%   sem        relative standard error of the mean of every element
%   lab        label per element
%   skipped    entries left out, with the reason
%
% Only entries that compare the model with a measurement AND carry a noise floor
% are used. In this objective a floor is fl = nEl * s^2 / n (s = relative
% between-preparation spread, n = preparations; RunRefit9Optim, objective
% deadband design), so the squared standard error of one element is fl / nEl,
% divided by the element's inverse-variance weight where the entry has one
% ('w:FIELD'). Boundary checks, grouped physiology bounds and penalties are
% priors, not data, and are left out.

val = []; dat = []; sem = []; lab = {}; skipped = {};
VAR_FLOOR_FRAC = 0.1;                                   % as in evalFeatureCost
for i = 1:numel(fn)
    s = fn{i};
    if hasTopLevelComma(s), skipped{end+1} = sprintf('%s  (grouped bounds)', s); continue; end %#ok<AGROW>
    if contains(s, '@'), skipped{end+1} = sprintf('%s  (custom cost)', s); continue; end %#ok<AGROW>
    tok = split(s, '|'); tok = tok(~cellfun(@isempty, tok));
    [feat, sel] = parseFeatSelector(tok{1});
    fl = 0; wf = ''; isRange = false;
    for t = 2:numel(tok)
        if startsWith(tok{t}, 'w:'), wf = tok{t}(3:end);
        elseif startsWith(tok{t}, 'fl:'), fl = str2double(tok{t}(4:end));
        elseif isRangeToken(tok{t}), isRange = true;
        end
    end
    if isRange, skipped{end+1} = sprintf('%s  (boundary check)', s); continue; end %#ok<AGROW>
    if fl <= 0, skipped{end+1} = sprintf('%s  (no noise floor)', s); continue; end %#ok<AGROW>
    if ~isfield(fd, feat) || ~isfield(fm, feat), skipped{end+1} = sprintf('%s  (missing)', s); continue; end %#ok<AGROW>
    d = applyFeatSelector(double([fd.(feat)]), sel); m = applyFeatSelector(double([fm.(feat)]), sel);
    if numel(m) ~= numel(d), skipped{end+1} = sprintf('%s  (length mismatch)', s); continue; end %#ok<AGROW>
    nF = mean(d, 'omitnan'); if nF == 0, nF = 1; end
    w = ones(size(d));
    if ~isempty(wf)
        vv = applyFeatSelector(double([fd.(wf)]), sel); vbar = mean(vv, 'omitnan');
        if isfinite(vbar) && vbar > 0
            vv = max(vv, VAR_FLOOR_FRAC * vbar); w = 1 ./ vv; w = reshape(w / mean(w, 'omitnan'), size(d));
        end
    end
    se = sqrt(fl / numel(d) ./ w);
    val = [val, m(:)' / nF]; dat = [dat, d(:)' / nF]; sem = [sem, se(:)']; %#ok<AGROW>
    lab = [lab, arrayfun(@(j) sprintf('%s(%d)', tok{1}, j), 1:numel(d), 'uni', 0)]; %#ok<AGROW>
end
end
