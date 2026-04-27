function cost = plotFeatures(feats_data, feats_sim, feats_ghost, fn)
%% PLOTFEATURES  Plot data vs simulation features (backward-compatible wrapper).
%
%   Wrapper around PLOTMULTIPLEFEATURES for the classic 3-dataset case:
%   data (blue o), simulation (red x), ghost (grey +).
%
%   COST = PLOTFEATURES(FEATS_DATA, FEATS_SIM, FEATS_GHOST, FN)
%
%   To build the sample fn, go:
%     fn = fieldnames(feats_sim);
%     quoted = cellfun(@(s) ['''' s ''''], fn, 'UniformOutput', false);
%     fprintf('fn = {%s};\n', strjoin(quoted, ', '));
%   A sample fn:
%     fn = {'ktr','A','t0','SLslack','v_restretch','peak1_y','peak1_t', ...
%           'peak1_SL','peak1_dSL','peak2','vall_t','vall_y','vall2_dy', ...
%           'vall2_t','ovrsht_dy','ovrsht_t','steady'};
%   X-Y plots: plots A as a function of SLslack
%     fn = {'A|SLslack'}
%
%   See also: plotMultipleFeatures, evalFeatureCost

if nargin < 3
    feats_ghost = [];
    fn = fieldnames(feats_data);
elseif nargin < 4
    fn = fieldnames(feats_data);
end

% Build the multi-dataset inputs
feats_cell = {feats_data, feats_sim};
labels     = {'data', 'sim'};
colors     = [0.00 0.45 0.70; 0.84 0.10 0.11];
markers    = {'o', 'x'};

if ~isempty(feats_ghost)
    feats_cell{end+1} = feats_ghost;
    labels{end+1}     = 'ghost';
    colors(end+1, :)  = [0.50 0.50 0.50];
    markers{end+1}    = '+';
end

% Compute per-feature costs (data vs sim only) and build a lookup map for
% display in tile titles.  Missing or NaN features get their penalty values
% from evalFeatureCost (≥ 100 → shown as "[miss]").
[cost, ~, cost_raw] = evalFeatureCost(feats_data, feats_sim, fn);
% cost_raw is the unweighted per-feature error; cost is the weighted total
cmap = containers.Map(fn, num2cell(cost));

figure(80085); clf;
plotMultipleFeatures(feats_cell, labels, colors, markers, fn, [], [], cmap);

end
