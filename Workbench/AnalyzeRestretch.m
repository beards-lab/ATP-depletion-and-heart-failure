% AnalyzeRestretch.m
% Run the current best-fit and focus on the RESTRETCH features, especially the
% post-restretch undershoot (vall2) + overshoot (ovrsht) discrepancy.
% Prints a data-vs-model table for the restretch subset and exports a per-
% segment zoom of the post-restretch hold (data vs model, relative to steady).

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root));
LoadData;

params0 = getParams();
run(fullfile(root, 'params', 'params_reseeded_regavail_opt2.m'));   % current best
params0 = getParams(params0, [], true, false);
params0.RunForceVelocity = false; params0.RunKtr = false; params0.RunStairs = false;
params0.RunForceVelocityTime = false; params0.RunSlack = true;
params0.RunSlackSegments = 'All'; params0.EvalFeatures = true;
params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0; params0.BreakOnODEUnstable = false;

RunBakersExp;   % -> features_model, features_data, out_slack

%% restretch-subset table
restr = {'peak1_y','peak1_t','peak1_dSL','vall_y','vall_t','peak2', ...
         'restretchSlopeStart','restretchSlopeEnd','steady', ...
         'vall2_dy','vall2_t','ovrsht_dy','ovrsht_t'};
fprintf('\n%-20s %-34s %-34s\n','feature','data','model');
fprintf('%s\n', repmat('-',1,90));
for k = 1:numel(restr)
    f = restr{k};
    d = ''; m = '';
    if isfield(features_data,f);  d = num2str(round(double(features_data.(f)),4,'significant')); end
    if isfield(features_model,f); m = num2str(round(double(features_model.(f)),4,'significant')); end
    fprintf('%-20s %-34s %-34s\n', f, d, m);
end

%% per-segment zoom of the post-restretch hold (vel4 -> vel5)
ds = load(fullfile(root,'data',params0.velocitytableonfile));
vt = ds.velocitytable; vt(1,1) = -20; dt = ds.datatable;
segs = find(vt(:,2) < -1)';
t = out_slack.t(:); F = out_slack.Force(:);

fig = figure(808); clf; set(fig,'Position',[60 60 1500 380]);
tl = tiledlayout(1, numel(segs), 'TileSpacing','compact','Padding','compact');
title(tl,'Post-restretch hold: data (black) vs model (blue), relative to each segment steady');
for i = 1:numel(segs)
    s = segs(i);
    if s+4 > size(vt,1); break; end
    t4 = vt(s+3,1); t5 = vt(s+4,1);
    wd = dt(:,1) >= t4-0.01 & dt(:,1) <= t5;
    wm = t      >= t4-0.01 & t      <= t5;
    sd = features_data.steady(i);  sm = features_model.steady(i);
    ax = nexttile(tl); hold(ax,'on');
    plot(ax, dt(wd,1)-t4, dt(wd,3)-sd, 'k-', 'LineWidth',1);
    plot(ax, t(wm)-t4,   F(wm)-sm,    'b-', 'LineWidth',1.5);
    yline(ax,0,':','steady');
    % mark data undershoot/overshoot
    plot(ax, features_data.vall2_t(i),  features_data.vall2_dy(i),  'kv','MarkerSize',8,'LineWidth',1.5);
    plot(ax, features_data.ovrsht_t(i), features_data.ovrsht_dy(i), 'k^','MarkerSize',8,'LineWidth',1.5);
    title(ax, sprintf('seg %d (v=%.0f)', i, features_data.v_restretch(i)));
    xlabel(ax,'t - t_{restretch} (s)'); if i==1; ylabel(ax,'F - steady (kPa)'); end
    ylim(ax,[-16 6]); grid(ax,'on');
end
figDir = fullfile(root,'Docs','figures'); if ~exist(figDir,'dir'); mkdir(figDir); end
exportgraphics(fig, fullfile(figDir,'restretch_overshoot_bestfit.png'), 'Resolution',150);
fprintf('\nSaved Docs/figures/restretch_overshoot_bestfit.png\n');
