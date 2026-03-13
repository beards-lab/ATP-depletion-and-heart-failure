function merged = mergeOutStructs(outs, dim)
%MERGEOUTSTRUCTS Concatenate arrays in a struct array field-by-field.
%
%   merged = mergeOutStructs(outs)
%   merged = mergeOutStructs(outs, dim)
%
%   Inputs:
%     outs  - 1xN struct array; all elements must share the same field names.
%             Typically these are 'out' structs returned by evaluateModel for
%             consecutive simulation chunks that need to be stitched together.
%     dim   - Concatenation dimension (default = 1).
%
%   Output:
%     merged - Single struct with each field concatenated across all elements
%              of outs along dimension dim.
%
%   Field types handled:
%     numeric / logical  - cat(dim, ...)  with vertcat/horzcat fallback
%     cell arrays        - vertcat or horzcat depending on dim
%     string / char      - vertcat (strings) or cat(dim, ...) (char)
%     struct             - first element is kept (no recursive merge)
%     other              - warning issued, field skipped
%
%   Example:
%     out = mergeOutStructs([outs{:}]);

    if nargin < 2 || isempty(dim)
        dim = 1; % default: stack along rows
    end
    if ~isstruct(outs)
        error('Input must be a struct array.');
    end
    if isempty(outs)
        merged = struct();
        return;
    end

    % Ensure all structs share the same fields
    baseFields = fieldnames(outs(1));
    for i = 2:numel(outs)
        if ~isequal(fieldnames(outs(i)), baseFields)
            error('All structs must have identical field names.');
        end
    end

    merged = struct();
    for f = 1:numel(baseFields)
        name = baseFields{f};
        vals = {outs.(name)};  % collect values from each struct for this field

        % Decide concatenation based on class of the first element
        firstVal = vals{1};

        if iscell(firstVal)
            % For cell arrays: use bracket concatenation in the chosen dimension
            if dim == 1
                merged.(name) = vertcat(vals{:});
            else
                merged.(name) = horzcat(vals{:});
            end

        elseif isnumeric(firstVal) || islogical(firstVal)
            % Numeric/logical arrays: use cat(dim, ...)
            try
                merged.(name) = cat(dim, vals{:});
            catch ME %#ok<NASGU>
                % Fallback: try vertical then horizontal
                try
                    merged.(name) = vertcat(vals{:});
                catch
                    merged.(name) = horzcat(vals{:});
                end
            end

        elseif isstring(firstVal) || ischar(firstVal)
            % Strings/char arrays
            if isstring(firstVal)
                merged.(name) = vertcat(vals{:});   % strings usually stack by rows
            else
                % for char, concatenate along dim if sizes permit
                merged.(name) = cat(dim, vals{:});
            end

        elseif isstruct(firstVal)
            % Generic fallback: keep first element (no recursive merge)
            merged.(name) = firstVal;

        else
            warning('Field %s unknown while merging parallel chunks \n', name)
        end
    end
end
