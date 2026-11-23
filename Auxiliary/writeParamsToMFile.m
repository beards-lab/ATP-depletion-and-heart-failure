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
            % %.17g, not %g. %g keeps 6 significant digits, and on 2026-09-23
            % that alone moved refit7's accepted best from 0.3355 to 0.3707
            % (+10 %) on reload -- the parameters differed by 1.1e-6, but the
            % model's re-stretch recovery is numerically rough at that level
            % (Analyses/KtrSegmentSlope/RunNumericalNoise.m). A snapshot must
            % reproduce the point the optimiser accepted, bit for bit.
            fprintf(fid, 'params0.%s = %.17g;\r\n', nam, par);
        elseif ischar(par) && ~isempty(par)
             % fprintf("Writing %s\n", nam);            
             fprintf(fid, "params0.%s = '%s';\r\n", nam, par);
        elseif iscell(par)
            % Write THIS field, not params.fn. It used to hardcode params.fn
            % here, so every cell-valued field came back from disk holding the
            % feature list: that is what silently corrupted params0.leanPool in
            % every optimiser snapshot (lean_opt, refit3_opt, ...), leaving the
            % next campaign with an empty parameter pool and nothing to tune.
            fprintf(fid, "params0.%s = {'%s'};\r\n", nam, strjoin(par, "', '"));
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