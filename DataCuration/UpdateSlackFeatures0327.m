% UpdateSlackFeatures0327.m
% Re-extract slack features for the 03/27/2026 protocol data files using
% extractSlackAttributes with trendline smoothing (useSmoothing=true), and
% write the result back into features_data in each .mat file.
%
% This is the noisier "new" data referenced by params0.velocitytableonfile
% (e.g. params0.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat').
% runSlackExperiment loads features_data straight from the .mat file when
% present, so re-running the smoothing only needs to happen here, once.
%
% Also exports a feature-extraction figure per file to Docs/figures/ for
% visual QC of the fit.

addpath(genpath('..'));

DATA_FILES = { ...
    '../data/protocol_03_27_2026_8mM_slack.mat', ...
    '../data/protocol_03_27_2026_2mM_slack.mat'};

figDir = '../Docs/figures/';
if ~exist(figDir, 'dir'), mkdir(figDir); end

for i_file = 1:numel(DATA_FILES)
    DATA_FILE = DATA_FILES{i_file};
    [~, baseName] = fileparts(DATA_FILE);

    fprintf('\n=== %s ===\n', DATA_FILE);
    datastruct    = load(DATA_FILE);
    datatable     = datastruct.datatable;
    velocitytable = datastruct.velocitytable;

    data_t  = datatable(:, 1);
    data_y  = datatable(:, 3);
    data_SL = datatable(:, 2);

    % Extract with trendline smoothing, plot for QC
    fig = figure(406); clf;
    set(fig, 'Visible', 'on', 'Position', [100 100 1000 650]);
    features_data = extractSlackAttributes(data_t, data_y, data_SL, ...
        velocitytable, [], [], true, true);
    title(sprintf('Slack feature extraction — %s', baseName), 'Interpreter', 'none');
    drawnow;

    % Carry over t0_crossing alias expected by runSlackExperiment
    if isfield(features_data, 't0') && ~isfield(features_data, 't0_crossing')
        features_data.t0_crossing = features_data.t0;
    end

    % Print extracted feature values
    fn = fieldnames(features_data);
    for k = 1:numel(fn)
        val = features_data.(fn{k});
        if ~isnumeric(val), continue; end
        fprintf('  %-20s %s\n', fn{k}, mat2str(round(val, 4, 'significant')));
    end

    % Save features_data back into the .mat file (preserve other fields)
    tmp = load(DATA_FILE);
    tmp.features_data = features_data;
    save(DATA_FILE, '-struct', 'tmp');
    fprintf('Saved features_data to %s\n', DATA_FILE);

    % Export QC figure
    figPath = fullfile(figDir, sprintf('slack_features_%s.png', baseName));
    exportgraphics(fig, figPath, 'Resolution', 150);
    fprintf('Exported figure to %s\n', figPath);
end

fprintf('\nDone.\n');
