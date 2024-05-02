function [UpdatedModNames] = getAllDifferent(params)
% lists params different from the base settings - only the modified

    if isfield(params, 'g')
        % update modNames and gs, if appropriate
        params = getParams(params, params.g, false, true);
        params = rmfield(params, 'g');
    else
        % we get all the fields that are in the params0
        params = getParams(params, [], false, false);
    end
    modNames = fieldnames(params);
    % compare to default
    params0 = getParams();
    UpdatedModNames = {};
    for i_mod = 1:length(modNames)
        if ~isfield(params0, modNames{i_mod}) && ...
            ~isempty(params.(modNames{i_mod})) || ... Does not exist in the base set and is not empty in the second
            length(params.(modNames{i_mod})) == 1 && ... OR is it a different single element
            length(params0.(modNames{i_mod})) == 1 && ...            
            params.(modNames{i_mod}) ~= params0.(modNames{i_mod})            
                disp(['Got '  modNames{i_mod}])
                UpdatedModNames{length(UpdatedModNames)+1} = modNames{i_mod};
        else
            % disp(['Drop '  modNames{i_mod}])
        end
    end
end
