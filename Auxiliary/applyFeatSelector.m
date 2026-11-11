function v_red = applyFeatSelector(v, selector)
% APPLYFEATSELECTOR  Reduce a per-segment vector to a scalar per selector.
%
%   V_RED = APPLYFEATSELECTOR(V, SELECTOR) reduces V according to SELECTOR.
%   If SELECTOR is empty, V is returned unchanged.
%
%   Accepted selectors:
%     '1','2',...,'N' - pick segment N (1-based). Out-of-range -> NaN.
%     'a:b'           - contiguous range of segments, clipped to numel(v).
%                       Returns a VECTOR, so the cost stays a mean over the
%                       kept segments. Used to drop a slack cycle that is a
%                       different manoeuvre (see Analyses/RecoveryOscillation).
%     'i,j,k'         - explicit list of segments, clipped the same way.
%     'end'           - last element.
%     'first'         - first non-NaN element (NaN if none).
%     'last'          - last non-NaN element (NaN if none).
%     'mean'          - mean(v, 'omitnan').
%     'median'        - median(v, 'omitnan').
%     'min'           - min(v, [], 'omitnan').
%     'max'           - max(v, [], 'omitnan').
%
%   See also: parseFeatSelector.

if isempty(selector)
    v_red = v;
    return;
end

% Index forms first: 'N', 'a:b', 'i,j,k' and any mixture. This must precede
% str2double, which reads '1,3,5' as the number 135 (comma = thousands
% separator) and would silently turn a segment LIST into a bad index.
if ~isempty(regexp(selector, '^\s*\d+\s*(:\s*\d+\s*)?(,\s*\d+\s*(:\s*\d+\s*)?)*$', 'once'))
    idx = [];
    for part = split(string(selector), ',')'
        ab = str2double(split(part, ':'));
        if numel(ab) == 1
            idx = [idx, ab];                       %#ok<AGROW>
        else
            idx = [idx, ab(1):ab(end)];            %#ok<AGROW>
        end
    end
    idx = unique(idx(idx >= 1 & idx <= numel(v)), 'stable');
    if isempty(idx); v_red = NaN; else; v_red = v(idx); end
    return;
end

switch selector
    case 'end'
        v_red = v(end);
    case 'first'
        idx = find(~isnan(v), 1, 'first');
        if isempty(idx); v_red = NaN; else; v_red = v(idx); end
    case 'last'
        idx = find(~isnan(v), 1, 'last');
        if isempty(idx); v_red = NaN; else; v_red = v(idx); end
    case 'mean'
        v_red = mean(v, 'omitnan');
    case 'median'
        v_red = median(v, 'omitnan');
    case 'min'
        v_red = min(v, [], 'omitnan');
    case 'max'
        v_red = max(v, [], 'omitnan');
    otherwise
        error('applyFeatSelector:bad', 'Unknown selector "[%s]"', selector);
end
end
