function v_red = applyFeatSelector(v, selector)
% APPLYFEATSELECTOR  Reduce a per-segment vector to a scalar per selector.
%
%   V_RED = APPLYFEATSELECTOR(V, SELECTOR) reduces V according to SELECTOR.
%   If SELECTOR is empty, V is returned unchanged.
%
%   Accepted selectors:
%     '1','2',...,'N' - pick segment N (1-based). Out-of-range -> NaN.
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

n = str2double(selector);
if ~isnan(n) && n == round(n) && n >= 1
    if n <= numel(v)
        v_red = v(n);
    else
        v_red = NaN;
    end
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
