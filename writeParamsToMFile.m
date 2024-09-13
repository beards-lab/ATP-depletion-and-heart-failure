function writeParamsToMFile(filename, params, modNames)
% Writes params to m file, so we can easily load. 

    if nargin < 3
        modNames = fieldnames(params);
    end
    
    fid = fopen(filename, 'w');
    fprintf(fid, "%% Generated file to manipulate the simulation parameters directly\n"); 
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

    fclose(fid);
end