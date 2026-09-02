function writeParamsToMFile(filename, params, modNames, comment, quiet)
% Writes params to m file, so we can easily load. 

    if nargin < 5 || isempty(quiet); quiet = false; end
    params = getParams(params, params.g, false, true);
    
    if nargin < 3 || isempty(modNames)
        modNames = fieldnames(params);
    end

    if nargin < 4
        comment = '';
    end

    if ~contains(filename, '.')
        filename = [filename '.m'];
    end

    
    fid = fopen(filename, 'w');
    fprintf(fid, "%% Generated file to manipulate the simulation parameters directly\r\n"); 
    if ~isempty(comment)
        fprintf(fid, "%% %s\r\n", comment); 
    end

    for i_row = 1:length(modNames)
        par = params.(modNames{i_row});
        nam = modNames{i_row};
        if isa(par, 'function_handle') || isa(par, 'struct')
            if ~quiet; fprintf("Skipping %s\r\n", nam); end
            continue;     
        elseif length(par) == 1        
            % fprintf("Writing %s\n", nam);
            fprintf(fid, 'params0.%s = %g;\r\n', nam, par);            
        elseif ischar(par) && ~isempty(par)
             % fprintf("Writing %s\n", nam);            
             fprintf(fid, "params0.%s = '%s';\r\n", nam, par);
        elseif iscell(par)
            % Serialise THIS field, not params.fn. The old code wrote
            % strjoin(params.fn,...) for every cell-valued field, so any
            % snapshot with more than one cell field silently got fn's
            % contents in all of them. That was invisible while fn was the
            % only meaningful cell field (params0.mods was also written wrong
            % but is reset to {} at the end of this file anyway), and it
            % corrupted AtpJoint_fn2 -- the low-ATP leg's feature list -- the
            % moment a second one existed.
            if isempty(par)
                fprintf(fid, "params0.%s = {};\r\n", nam);
            else
                fprintf(fid, "params0.%s = {'%s'};\r\n", nam, strjoin(par, "', '"));
            end
        elseif any(size(par) > 1) ...
            && nam ~= "s" ...
            && nam ~= "g" ...
            && nam ~= "mods" ...
            && nam ~= "datatable" ...
            && nam ~= "PU0" ...
            && all(size(par(:)) < 100) 
            % matrix
            % fprintf("Writing %s\n", nam);
            par_s = mat2str(par);
            fprintf(fid, "params0.%s = %s;\r\n", nam, par_s); 
        else
            if ~quiet; fprintf("Skipping %s\r\n", nam); end
            continue;
        end
    end
    % reset the mods and g
    fprintf(fid, "%% Reset the modifiers, as they are already applied\r\n"); 
    fprintf(fid, "params0.mods = {};\r\n"); 
    fprintf(fid, "params0.g = [];\r\n"); 
    fclose(fid);
end