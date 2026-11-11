% FocusSlackAttributes.m
% Focus driver: test and tune extractSlackAttributes on raw experimental data.
%
% Run from the repo root in MATLAB.  Does three things:
%   1. Extracts features twice — without and with trendline smoothing — and
%      prints a side-by-side comparison so you can see the effect.
%   2. Plots Figure 406 (produced inside extractSlackAttributes) showing the
%      smoothed trendline overlaid on the raw signal for each slack segment.
%   3. Saves the smoothed features_data back into the .mat file so that
%      runSlackExperiment loads it from cache (lines 48-52) and never
%      re-runs the extraction during optimisation.
%
% The data file defaults to data/bakers_slack8mM_all.mat (the file used by
% all recent param sets).  Change DATA_FILE to target a different condition.

DATA_FILE = 'data/protocol_03_27_2026_8mM_slack.mat';

addpath(genpath('.'));

%% Load data
fprintf('Loading %s ...\n', DATA_FILE);
datastruct = load(DATA_FILE);
datatable    = datastruct.datatable;
velocitytable = datastruct.velocitytable;
velocitytable(1, 1) = -20;   % same pre-processing as runSlackExperiment

data_t  = datatable(:, 1);
data_y  = datatable(:, 3);
data_SL = datatable(:, 2);

%% Pass A — no smoothing (current behaviour)
fprintf('\nPass A: extractSlackAttributes without smoothing...\n');
feats_raw = extractSlackAttributes(data_t, data_y, data_SL, velocitytable, ...
    [], [], false, false);

%% Pass B — trendline smoothing (new behaviour), with plots
fprintf('Pass B: extractSlackAttributes with trendline smoothing...\n');
feats_smooth = extractSlackAttributes(data_t, data_y, data_SL, velocitytable, ...
    [], [], true, true);

%% Side-by-side comparison
fprintf('\n%-22s  %s\n', 'Feature', 'Raw   -->  Smoothed  (per segment)');
fprintf('%s\n', repmat('-', 1, 70));

fields_of_interest = {'peak1_y', 'peak1_t', 'vall_y', 'vall_t', ...
                      'vall2_dy', 'ovrsht_dy', 'ktr', 'A', 'steady'};

for k = 1:numel(fields_of_interest)
    fn = fields_of_interest{k};
    if ~isfield(feats_raw, fn) || ~isfield(feats_smooth, fn)
        continue;
    end
    vr = feats_raw.(fn);
    vs = feats_smooth.(fn);
    raw_str    = strjoin(arrayfun(@(x) sprintf('%7.3f', x), vr, 'UniformOutput', false), '  ');
    smooth_str = strjoin(arrayfun(@(x) sprintf('%7.3f', x), vs, 'UniformOutput', false), '  ');
    fprintf('%-22s  raw:    %s\n', fn, raw_str);
    fprintf('%-22s  smooth: %s\n', '',  smooth_str);
    fprintf('\n');
end

%% Save smoothed features_data into the .mat file for caching
features_data = feats_smooth;
% Carry over t0_crossing alias expected by runSlackExperiment
if isfield(features_data, 't0') && ~isfield(features_data, 't0_crossing')
    features_data.t0_crossing = features_data.t0;
end

fprintf('Saving features_data to %s ...\n', DATA_FILE);
% Load everything that is already in the file, add/replace features_data,
% then write the whole struct back — this preserves datatable, velocitytable,
% and any other fields that are already there.
tmp = load(DATA_FILE);
tmp.features_data = features_data;
save(DATA_FILE, '-struct', 'tmp');
fprintf('Done.  runSlackExperiment will now load features_data from cache.\n');
