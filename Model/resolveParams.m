function params = resolveParams(params0)
% RESOLVEPARAMS  Evaluate linked parameter expressions marked with '=' prefix.
%
%   PARAMS = RESOLVEPARAMS(PARAMS0) scans all fields of PARAMS0. Fields
%   whose value is a string starting with '=' are evaluated as MATLAB
%   expressions, with all numeric fields of PARAMS0 in scope.
%   All other fields are passed through unchanged.
%
%   Example:
%     p.kah  = 80;
%     p.kamh = '=0.1*kah';   % will become 8.0
%     p = resolveParams(p);   % p.kamh == 8.0
%
%   Inputs:
%     PARAMS0 - parameter struct (may contain '=...' string fields)
%
%   Outputs:
%     PARAMS  - same struct with '=' fields replaced by evaluated values
%
%   See also: getParams

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