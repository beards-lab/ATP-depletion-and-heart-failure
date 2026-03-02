function params = resolveParams(params0)
% RESOLVEPARAMS Resolves linked fields marked with '=' prefix.
%   Fields like "=a*10" are evaluated with other fields as local variables.
%   Regular strings (no '=' prefix) and numeric values are passed through.

params = params0;
fields = fieldnames(params0);

% Load all numeric fields into a local struct for expression evaluation
localVars = struct();
for i = 1:numel(fields)
    val = params0.(fields{i});
    if isnumeric(val)
        localVars.(fields{i}) = val;
    end
end

% Resolve linked fields
for i = 1:numel(fields)
    val = params0.(fields{i});
    if ischar(val) || isstring(val)
        val = char(val);
        if startsWith(val, '=')
            expr = val(2:end); % strip '='
            try
                % Assign local vars into this scope for eval
                fnames = fieldnames(localVars);
                for j = 1:numel(fnames)
                    name = fnames{j};                    
                    eval(sprintf('%s = localVars.(name);', name));
                end
                params.(fields{i}) = eval(expr);
            catch e
                error('Failed to resolve field "%s" with expression "%s": %s', ...
                      fields{i}, expr, e.message);
            end
        end
    end
end
end