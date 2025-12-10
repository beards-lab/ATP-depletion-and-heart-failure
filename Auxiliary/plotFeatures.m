function cost = plotFeatures(feats_data, feats_sim, feats_ghost, fn)
%% plots the features
% To build the sample fn, go:
%   fn = fieldnames(feats_sim);
%   quoted = cellfun(@(s) ['''' s ''''], fn, 'UniformOutput', false);
%   fprintf('fn = {%s};\n', strjoin(quoted, ', '));
% a sample fn
% fn = {'ktr', 'A', 't0', 'SLslack', 'v_restretch', 'peak1_y', 'peak1_t', 'peak1_SL', 'peak1_dSL', 'peak2', 'vall_t', 'vall_y', 'vall2_dy', 'vall2_t', 'ovrsht_dy', 'ovrsht_t', 'steady'};
% X-Y plots: this plots A as a function of SLslack
% fn = {'A|SLslack'}

if nargin < 3
    feats_ghost = [];
    fn = fieldnames(feats_data);
elseif nargin < 4
    fn = fieldnames(feats_data);
end

figure(80085);clf; 
tiledlayout("flow");
% N_feats = size(feats_data, 2);

cost = evalFeatureCost(feats_data, feats_sim, fn);

for i_feat = 1:size(fn, 2)
    nexttile();

    featnames = split(fn{i_feat}, '|');
    feat_y = featnames{1};
    N_feats = length(feats_data.(feat_y));

    if length(featnames) == 1 || (length(featnames) >= 2 && featnames{2} == "_" || isnumeric(featnames{2}))
        fd_x = 1:N_feats;fs_x = 1:N_feats;fg_x = 1:N_feats;
        feat_x = '';
    else
        feat_x = featnames{2};
        fd_x = [feats_data.(feat_x)];
        fs_x = [feats_sim.(feat_x)];
        if exist("feats_ghost", "var") && ~isempty(feats_ghost)
            fg_x = [feats_ghost.(feat_x)];
        else
            fg_x = NaN(size(1:N_feats));
        end
    end

    fd = [feats_data.(feat_y)];
    fs = [feats_sim.(feat_y)];
    if exist("feats_ghost", "var") && ~isempty(feats_ghost)
        fg = [feats_ghost.(feat_y)];
    else
        fg = NaN(size(1:N_feats));
    end

    plot(fd_x, fd, 'o--', fs_x, fs, 'x--', fg_x, fg, '+--');
    if isempty(feat_x)
        title(sprintf('%s (%.3g)', feat_y, cost(i_feat)), Interpreter="none");
    else
        title(sprintf('%s vs %s (%.3g)', feat_y, feat_x, cost(i_feat)), Interpreter="none");
    end
end

% dt vs dSL with extrapolation
% nexttile; plot([feats_data.t0], [feats_data.SLslack]-2.2);ylim([-inf, 0]);title("dt vs dSL");hold on;
% t = (-1:0.1:3)*1e-3;plot(t, -16*t - 0.125);
% nexttile; plot(diff([feats_data.SLslack]-2.2)./diff([feats_data.t0]));ylim([-inf, 0]);title("dt vs dSL (um/s)")

end