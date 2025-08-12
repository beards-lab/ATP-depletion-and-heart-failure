%% Driver_all
% Evaluates all model params in the folder

% clear;

filenames = list_model_param_files();


LoadData; 
done = [];
errs = [];

for i_pfn = 45:length(filenames)
% for i_pfn = done
    cf = figure(1420);clf;
    params0 = getParams();

    params0.modelFcn = 'dPUdTCaSimpleAlternative2State';

    run(filenames{i_pfn})
% initialize parameters
    
    tic
    try
        fprintf('Ruunning #%g: %s...', i_pfn, filenames{i_pfn});
        RunBakersExp;
        % done = [done, i_pfn];
        fprintf('Ok.\n');
        saveas(cf, [filenames{i_pfn} '.png'])
    catch e
        errs = [errs, i_pfn];
        fprintf('Err %s (%s:%g) \n', e.message, e.stack(1).file, e.stack(1).line);
    end

end


function matches = list_model_param_files()
    allFiles = dir('*.m');
    matches = {};
    for k = 1:length(allFiles)
        [~, name, ext] = fileparts(allFiles(k).name);
        if startsWith(name, 'ModelParams') || startsWith(name, 'ModelOptParams')
            matches{end+1} = name; %#ok<AGROW>
        end
    end
end
