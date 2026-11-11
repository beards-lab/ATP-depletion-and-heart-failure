function [matchedNames, matchedVals] = matchStructFields(S, pattern, valueFilter, listValues)
%MATCHSTRUCTFIELDS Return names & values of struct fields matching a wildcard pattern
%
% Usage:
%   [names, vals] = matchStructFields(S, 'Use*')
%   [names, vals] = matchStructFields(S, 'Use*', 1)         % only values == 1
%   [names, vals] = matchStructFields(S, '*_L', [], true)   % print values
%   [names, vals] = matchStructFields(S, 'dr*', 42, true)   % name & value filter
%
% Inputs:
%   S            - scalar struct
%   pattern      - string/char with '*' wildcard(s), e.g. 'Use*', '*_L'
%   valueFilter  - (optional) value to match exactly (isequal), [] to skip filter
%   listValues   - (optional) logical, if true prints names & values
%
% Outputs:
%   matchedNames - cell array of matching field names
%   matchedVals  - cell array of corresponding field values

    if nargin < 4
        listValues = false;
    end
    if nargin < 3
        valueFilter = [];
    end

    if ~isstruct(S) || numel(S) ~= 1
        error('S must be a scalar struct.');
    end

    % Ensure pattern is char
    if isstring(pattern)
        pattern = char(pattern);
    end

    % Convert wildcard to regex
    regexPattern = ['^' strrep(pattern, '*', '.*') '$'];

    % Get field names and values
    f = fieldnames(S);
    vals = struct2cell(S);

    % Match field names
    nameMask = ~cellfun('isempty', regexp(f, regexPattern, 'once'));

    % Optional value match
    if ~isempty(valueFilter)
        valueMask = cellfun(@(v) isequal(v, valueFilter), vals);
    else
        valueMask = true(size(vals));
    end

    % Combined mask
    mask = nameMask & valueMask;

    matchedNames = f(mask);
    matchedVals  = vals(mask);

    % Optionally print results
    if listValues
        if isempty(matchedNames)
            fprintf('No matching fields found.\n');
        else
            fprintf('Matching fields for "%s":\n', pattern);
            for i = 1:numel(matchedNames)
                fprintf('  %s: %s\n', matchedNames{i}, mat2str(matchedVals{i}));
            end
        end
    end
end
