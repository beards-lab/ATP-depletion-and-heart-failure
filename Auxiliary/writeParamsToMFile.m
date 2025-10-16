function writeParamsToMFile(filename, params, modNames, comment)
% Writes params to m file, so we can easily load. 

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
        if length(par) == 1        
            fprintf(fid, 'params0.%s = %g;\r\n', nam, par);
        elseif ischar(par) && ~isempty(par)
             fprintf(fid, "params0.%s = '%s';\r\n", nam, par);
        elseif any(size(par) > 1) ...
            && nam ~= "s" ...
            && nam ~= "g" ...
            && nam ~= "mods" ...
            && nam ~= "datatable" ...
            && nam ~= "PU0" ...
            && all(size(par(:)) < 100) 
            % matrix
            par_s = mat2str(par);
            fprintf(fid, "params0.%s = %s;\r\n", nam, par_s); 
        else
            fprintf("Skipping %s\r\n", nam);
            continue;
        end
    end
    % reset the mods and g
    fprintf(fid, "%% Reset the modifiers, as they are already applied\r\n"); 
    fprintf(fid, "params0.mods = {};\r\n"); 
    fprintf(fid, "params0.g = [];\r\n"); 
    fclose(fid);
end