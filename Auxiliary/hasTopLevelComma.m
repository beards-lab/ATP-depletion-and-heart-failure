function tf = hasTopLevelComma(s)
%HASTOPLEVELCOMMA  True if S contains a comma that is not inside parentheses.
%
%   Used by evalFeatureCost and plotMultipleFeatures to detect grouped
%   feature specifier strings such as 'ka|50-500|10,kd|5-200|10'.
%
%   See also: splitTopLevelComma, evalFeatureCost

depth = 0;
for k = 1:numel(s)
    c = s(k);
    if c == '('; depth = depth + 1; end
    if c == ')'; depth = depth - 1; end
    if c == ',' && depth == 0; tf = true; return; end
end
tf = false;
end
