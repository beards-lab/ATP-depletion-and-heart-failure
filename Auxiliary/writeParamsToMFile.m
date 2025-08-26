function writeParamsToMFile(filename, params, modNames, comment)
% Writes params to m file, so we can easily load. 

    if nargin < 3
        modNames = fieldnames(params);
    end

    if nargin < 4
        comment = '';
    end

    
    fid = fopen(filename, 'w');
    fprintf(fid, "%% Generated file to manipulate the simulation parameters directly\r\n"); 
    if ~isempty(comment)
        fprintf(fid, "%% %s\r\n", comment); 
    end

    params = getParams(params, params.g, false, true);
    for i_row = 1:length(modNames)
        if length(params.(modNames{i_row})) == 1        
            fprintf(fid, 'params0.%s = %g;\r\n', modNames{i_row}, params.(modNames{i_row}));
        elseif ischar(params.(modNames{i_row})) && ~isempty(params.(modNames{i_row}))
             fprintf(fid, "params0.%s = '%s';\r\n", modNames{i_row}, params.(modNames{i_row}));
        else
            fprintf("Skipping %s\r\n", modNames{i_row});
            continue;
        end
    end
    % reset the mods and g
    fprintf(fid, "%% Reset the modifiers, as they are already applied\r\n"); 
    fprintf(fid, "params0.mods = {};\r\n"); 
    fprintf(fid, "params0.g = [];\r\n"); 
    fclose(fid);
end