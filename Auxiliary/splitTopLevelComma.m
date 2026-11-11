function parts = splitTopLevelComma(s)
%SPLITTOPLEVELCOMMA  Split S on commas at the top nesting level.
%
%   PARTS = SPLITTOPLEVELCOMMA(S) returns a cell array of substrings.
%   Commas inside parentheses (a custom cost function) or square brackets (a
%   segment list such as 'rsK80[1,3,5]') are not treated as delimiters.
%
%   See also: hasTopLevelComma, evalFeatureCost

parts = {};
depth = 0;
start = 1;
for k = 1:numel(s)
    c = s(k);
    if c == '(' || c == '['; depth = depth + 1; end
    if c == ')' || c == ']'; depth = depth - 1; end
    if c == ',' && depth == 0
        parts{end+1} = s(start:k-1); %#ok<AGROW>
        start = k + 1;
    end
end
parts{end+1} = s(start:end);
end
