% CompareFeatures0410.m
% How far off is the new 04/10/2026 (Male 2-8) prep?
%
% Two comparison figures, one per ATP level, each overlaying three series on the
% standard slack feature panels:
%     * 03/27/2026 M  data   - the prep every current parameter set was fit to
%     * 04/10/2026    data   - the new prep
%     * model                - ThursdayNightFever driven by the 04/10 protocol
%
% Data features come straight from the .mat files (written by
% CreateProtocolVelocityTable + UpdateSlackFeatures); model features come from
% extractSlackAttributes on the simulated trace, so the two are extracted by the
% same code path and are directly comparable.
%
% CAVEAT on the low-ATP model series: ThursdayNightFever is a HIGH-ATP
% parameter set. Only MgATP is switched to 2 mM here; the calibrated low-ATP
% parameterisation (which additionally re-scales k2 -> K_T2 and kmsrd -> K_srx)
% lives in the lowATP_TNF_* snapshots. The green line in the 2 mM figure is
% therefore "what the 8 mM fit predicts at 2 mM", not a low-ATP fit.

clear; close all; clc;
cd(fileparts(mfilename('fullpath')));
addpath(genpath('..'));

figDir = '../Docs/figures/';
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% ── Feature panels (same selection as CompareProtocols) ────────────────────
fn_slack = {'ktr|SLslack', 'ktr2_k|SLslack', 'ktr2_A|SLslack', 'A|SLslack', ...
            'ktr2_omega', 't0', 'restretchSlopeStart', 'restretchSlopeEnd', ...
            'peak1_y', 'peak1_t', 'peak1_dSL', 'vall_y', 'vall_t', 'peak2', ...
            'steady', 'vall2_dy', 'vall2_t', 'ovrsht_dy', 'ovrsht_t'};

%% ── Baseline parameters ────────────────────────────────────────────────────
params0 = getParams(loadParams('ThursdayNightFever'), [], true, false);
params0.PlotEachSeparately = 0;
params0.PlotFeatureFitting = 0;
params0.ghostSave = '';  params0.ghostLoad = '';
params0.EvalFeatures = 1;
params0.RunSlackSegments = 'All';

%% ── The two conditions ─────────────────────────────────────────────────────
cond = struct( ...
    'name',   {'High ATP (8 mM)',                        'Low ATP (2 mM)'}, ...
    'ref',    {'protocol_03_27_2026_8mM_slack.mat',      'protocol_03_27_2026_2mM_slack.mat'}, ...
    'new',    {'protocol_04_10_2026_8mM_slack.mat',      'protocol_04_10_2026_2mM_slack.mat'}, ...
    'MgATP',  {8,                                        2}, ...
    'fig',    {410,                                      411});

for c = 1:numel(cond)
    fprintf('\n================ %s ================\n', cond(c).name);

    Sref = load(['data/' cond(c).ref]);
    Snew = load(['data/' cond(c).new]);
    fd_ref = Sref.features_data;
    fd_new = Snew.features_data;

    % model on the NEW protocol
    p = params0;
    p.velocitytableonfile = cond(c).new;
    p.MgATP = cond(c).MgATP;
    p = getParams(p, p.g, true, false);
    [~, ~, fm_new] = runSlackExperiment(p);

    feats  = {fd_ref, fd_new, fm_new};
    labels = {'03/27 data', '04/10 data', sprintf('model (TNF, %g mM)', cond(c).MgATP)};
    colors = [0.35 0.35 0.35; 0.00 0.45 0.74; 0.20 0.65 0.20];
    markers   = {'s', 'o', '^'};
    lineStyle = {'-', '-', '--'};
    fills     = [true, true, false];

    figure(cond(c).fig); clf;
    set(gcf, 'Name', ['Slack features — ' cond(c).name], ...
             'Units', 'centimeters', 'Position', [2 2 34 24]);
    sgtitle(sprintf('Slack features — %s   (03/27 vs 04/10 vs model)', cond(c).name), ...
            'Interpreter', 'none');
    plotMultipleFeatures(feats, labels, colors, markers, fn_slack, lineStyle, fills);

    % exportgraphics(gcf, fullfile(figDir, sprintf('features_0410_vs_0327_%dmM.png', ...
    %     cond(c).MgATP)), 'Resolution', 150);

    % ── numeric summary: how far off, per feature ──────────────────────────
    fprintf('%-22s %11s %11s %8s %11s %8s\n', 'feature', '03/27', '04/10', ...
            'new/ref', 'model', 'mdl/new');
    keys = {'ktr','A','Am','t0','peak1_y','vall_y','peak2','steady', ...
            'vall2_dy','ovrsht_dy','restretchSlopeStart','SLslack'};
    for k = 1:numel(keys)
        f = keys{k};
        vr = iMean(fd_ref, f);  vn = iMean(fd_new, f);  vm = iMean(fm_new, f);
        fprintf('%-22s %11.4g %11.4g %8.3f %11.4g %8.3f\n', f, vr, vn, vn/vr, vm, vm/vn);
    end
end

%% ── Within-dataset 2 mM / 8 mM ratio: the ATP effect, per prep ─────────────
% This is the comparison that matters, because it cancels prep-to-prep force
% scale. It also exposes the ATP-ORDER confound: 03/27 ran 8 mM first then
% 2 mM, 04/10 ran 2 mM FIRST then 8 mM. Rundown always favours whichever
% condition was recorded first, so any ratio that flips with the order is
% rundown, not ATP; any ratio that replicates across the two orders is real.
fprintf('\n\n============ ATP effect within each prep (2 mM / 8 mM) ============\n');
fprintf('03/27 order: 8 mM -> 2 mM      04/10 order: 2 mM -> 8 mM (REVERSED)\n\n');
D = struct( ...
    'tag',  {'03/27', '04/10'}, ...
    'hi',   {'protocol_03_27_2026_8mM_slack.mat', 'protocol_04_10_2026_8mM_slack.mat'}, ...
    'lo',   {'protocol_03_27_2026_2mM_slack.mat', 'protocol_04_10_2026_2mM_slack.mat'});
keys = {'A', 'Am', 'steady', 'peak1_y', 'peak2', 'vall_y', 'ktr', 't0', 'restretchSlopeStart'};

r = nan(numel(keys), numel(D));
for d = 1:numel(D)
    H = load(['data/' D(d).hi]);  Lo = load(['data/' D(d).lo]);
    for k = 1:numel(keys)
        r(k, d) = iMean(Lo.features_data, keys{k}) / iMean(H.features_data, keys{k});
    end
end
fprintf('%-22s %10s %10s %12s\n', 'feature (2mM/8mM)', D(1).tag, D(2).tag, 'replicates?');
for k = 1:numel(keys)
    if abs(r(k,1) - r(k,2)) <= 0.15 * abs(r(k,1))
        verdict = 'yes';
    else
        verdict = 'NO - order?';
    end
    fprintf('%-22s %10.3f %10.3f %12s\n', keys{k}, r(k,1), r(k,2), verdict);
end
fprintf('\nA ratio that replicates across the reversed ATP order is an ATP effect;\n');
fprintf('one that does not is contaminated by rundown (first-recorded wins).\n');

fprintf('\nFigures exported to %s\n', figDir);

% ---------------------------------------------------------------------------
function m = iMean(s, f)
%IMEAN Mean of a feature across slacks, or NaN when absent / all-NaN.
    if ~isfield(s, f); m = NaN; return; end
    v = s.(f);
    if ~isnumeric(v) || isempty(v); m = NaN; return; end
    m = mean(v, 'omitnan');
end
