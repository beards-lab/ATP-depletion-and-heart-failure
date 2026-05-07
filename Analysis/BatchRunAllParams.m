% BatchRunAllParams.m
% Batch-evaluate all 2026 param files from params/ and ModelOptParams/ folders.
% For each file saves:
%   <name>_Fig1.png     — simulation traces (figure 1)
%   <name>_Fig80085.png — feature cost breakdown + total cost
%   <name>_FigRates.png — strain-dependent rate curves (justPlotStateTransitionsFlag)

clear; clc;
% close all; 

%% Discover param files from 2026 across both folders
root = fullfile(fileparts(mfilename('fullpath')), '..');
d1 = dir(fullfile(root, 'params',         'ModelParams*.m'));
d2 = dir(fullfile(root, 'params',         'ModelOptParams*.m'));
d3 = dir(fullfile(root, 'ModelOptParams', 'Model*.m'));
d  = [d1; d2; d3];

mask = [d.datenum] >= datenum(2026, 1, 1);
d = d(mask);
[~, ord] = sort([d.datenum]);
d = d(ord);

%% For iter-series files, keep only the highest iter per series; keep all non-iter files
% e.g. ModelOptParams_iter_1..27 → keep only _iter_27
names = {d.name};
keep  = true(size(names));
iter_pat = regexp(names, '^(.+_iter_)(\d+)\.m$', 'tokens', 'once');
for i = 1:numel(names)
    if isempty(iter_pat{i}); continue; end   % not an iter file
    prefix_i = iter_pat{i}{1};
    num_i    = str2double(iter_pat{i}{2});
    % suppress if any later file has same prefix and higher number
    for j = i+1:numel(names)
        if isempty(iter_pat{j}); continue; end
        if strcmp(iter_pat{j}{1}, prefix_i) && str2double(iter_pat{j}{2}) > num_i
            keep(i) = false;
            break;
        end
    end
end
d = d(keep);

fprintf('Found %d param files from 2026 (iter-series reduced to last only):\n', numel(d));
for i = 1:numel(d)
    fprintf('  %s  (%s)\n', d(i).name, datestr(d(i).datenum, 'yyyy-mm-dd HH:MM'));
end

%% Output folder
out_dir = fullfile(fileparts(mfilename('fullpath')), 'batch_output_passive');
mkdir(out_dir);
fprintf('\nOutput → %s\n\n', out_dir);

%% Loop
for i = 1:numel(d)
    name = d(i).name(1:end-2);  % strip .m
    param_file = fullfile(d(i).folder, d(i).name);
    fprintf('[%d/%d] %s\n', i, numel(d), name);

    %% Step A: normal run
    params0 = getParams();
    try
        run(param_file);
    catch e
        warning('Failed to load %s: %s', name, e.message);
        continue;
    end


    params0.PlotEachSeparately           = true;
    params0.RunForceVelocity             = true;
    params0.justPlotStateTransitionsFlag = false;
    params0.BreakOnODEUnstable           = false;
    params0.ghostSave                    = '';
    params0.ghostLoad                    = '';
    params0.FV_velocities                = [0 -0.5 -1 -2 -4];
    params0.MaxRunTime                   = 120;
    params0.fn = {'FV_f|FV_v', 'ktr|2', 'A|50', 't0|5', 'ktr_rmse|0', 'restretchSlopeStart|0.1', 'restretchSlopeEnd|0.1', 'peak1_y|10', 'peak1_dSL|0.2', 'vall_y|10', 'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|0.1', 'ovrsht_dy|0.1', 'XTOR|@(X, Y_data) 0.01*sum((max(0, 2 - X) + max(0, X - 8)).^2)'};
params0.OptimizeOn = 'Feats';    
params0.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
params0.RunSlackSegments = 'All';
params0.ShowStatePlots = false;

    features_data  = struct();
    features_model = struct();
    ModelParams_Passive;

    figure(1); clf;
    try
        RunBakersExp;
    catch e
        warning('RunBakersExp failed for %s: %s', name, e.message);
    end
    exportgraphics(figure(1), fullfile(out_dir, [name '_Fig1.png']), 'Resolution', 150);

    %% Step B: feature cost figure with total cost
    cost_vec = [];
    if isfield(params0, 'fn') && ~isempty(params0.fn) && ...
            ~isempty(fieldnames(features_data)) && ~isempty(fieldnames(features_model))
        cost_vec = plotFeatures(features_data, features_model, [], params0.fn);
        figure(80085);
        annotation('textbox', [0 0.96 1 0.04], ...
            'String', sprintf('Total cost: %.3f', sum(cost_vec)), ...
            'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontSize', 11);
        exportgraphics(figure(80085), fullfile(out_dir, [name '_Fig80085.png']), 'Resolution', 150);
        fprintf('  total_cost = %.3f\n', sum(cost_vec));
    else
        fprintf('  (no features to plot)\n');
    end

    %% Step C: rates figure via justPlotStateTransitionsFlag
    figure(24); clf;
    params0.justPlotStateTransitionsFlag = true;
    try
        RunBakersExp;
    catch
        % expected — plotStateTransitions draws fig 24 then throws
    end
    exportgraphics(figure(24), fullfile(out_dir, [name '_FigRates.png']), 'Resolution', 150);
    params0.justPlotStateTransitionsFlag = false;

    fprintf('  saved PNGs for %s\n', name);
end

fprintf('\nDone. Output in:\n  %s\n', out_dir);
