%% Driver_all
% Evaluates all model params in the folder

% clear;
% updateFromGhosts = true;
updateFromGhosts = false;

%% get all ghosts and get their params, saving as modelOptParams
if updateFromGhosts
    % List all MAT files starting with 'Ghost_'
    files = dir('Ghost_*.mat');
    
    for i = 1:numel(files)
        fprintf('Checking file: %s... ', files(i).name);
    
        % Load the file
        S = load(files(i).name);
    
    
        % Check for 'params0'
        if isfield(S, 'params0')
            fprintf('Found params0\n');
            [~, base, ~] = fileparts(files(i).name);
            parts = split(base, '_');
            newname = strjoin(parts(2:end-1), '_');
    
            writeParamsToMFile(['ModelParamsGhost_' newname '.m'], params0)
            % === Your stuff goes here ===
            % Example: disp(S.params0)
            % Example: S.params0.myField
        else
            fprintf('No params0\n');
        end
    end
end
%% get all files
filenames = list_model_param_files();
% filenames = load('filenames.mat').errs;
% filenames = load('filenames.mat').done;

LoadData; 
done = [];
errs = [];

%%
for i_pfn = 9:length(filenames)
% for i_pfn = done
    cf = figure(1420);clf;
    params0 = getParams();

    % params0.modelFcn = 'dPUdTCaSimpleAlternative2State';
    run(filenames{i_pfn})

    params0.RunForceVelocity = true;
    params0.RunSlack = true;
    params0.EvalFitSlackOnset = true;
    params0.RunForceLengthEstim = false;
    params0.RunSlackSegments = 'All';
    params0.drawForceOnset = true;

% initialize parameters
    
    tic
    try
        fprintf('Ruunning #%g: %s...', i_pfn, filenames{i_pfn});
        RunBakersExp;
        done = [done, i_pfn];
        fprintf('Ok.\n');
        saveas(cf, [filenames{i_pfn} '.png'])
    catch e
        errs = [errs, i_pfn];
        fprintf('Err %s (%s:%g) \n', e.message, e.stack(1).file, e.stack(1).line);
    end

end

save filenames done errs


function matches = list_model_param_files()
    allFiles = dir('Model*.m');
    matches = {};
    for k = 1:length(allFiles)
        [~, name, ext] = fileparts(allFiles(k).name);
        if startsWith(name, 'ModelParams') || startsWith(name, 'ModelOptParams')
            matches{end+1} = name; %#ok<AGROW>
        end
    end
end
