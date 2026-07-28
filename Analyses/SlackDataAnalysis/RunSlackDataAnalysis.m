% RunSlackDataAnalysis.m
% Survey of EVERY slack dataset in data/, purely on the data side (no model).
%
% Three products, all written to results/:
%   1. fit_<tag>.png        - the feature-extraction QC plot for each dataset,
%                             so the fits behind every number can be inspected.
%   2. overview_traces.png  - all force traces, one tile per protocol day,
%                             time-aligned on the first slack release.
%   3. features_all.png     - one comparison across ALL datasets on the standard
%                             slack feature panels. Colour = protocol day,
%                             SOLID = high ATP (8 mM), DASHED = low ATP (2 mM),
%                             DASH-DOT = very low ATP (0.2 mM, Baker only).
%   4. consistency table    - printed, plus features_all.mat
%
% WHY EXTRACT FRESH RATHER THAN READ features_data FROM THE .mat FILES
% The stored features_data were written at different times by different script
% revisions (some smoothed, some not; the 03/27 PNB+Mava file predates a
% transition-refinement fix). A cross-dataset comparison is only meaningful if
% every number comes from ONE code path, so this script re-runs
% extractSlackAttributes with identical settings on every dataset and never
% writes back. The .mat files are left untouched.

clear; close all; clc;
cd(fileparts(which('RunSlackDataAnalysis')));
addpath(genpath('../..'));

% isfolder, NOT exist(...,'dir'): addpath(genpath('../..')) puts every other
% analysis's results/ on the path, so exist('results','dir') finds one of those
% and the local folder never gets created.
resDir = 'results';
if ~isfolder(resDir), mkdir(resDir); end
dataDir = '../../data/';

%% ── Dataset registry ───────────────────────────────────────────────────────
% day     : protocol day (sets colour)
% ATP     : mM (sets line style)
% kind    : 'active' (enters the comparison) | 'passive' (PNB+Mava control)
DS = struct( ...
 'tag',  {'baker_8',                'baker_2',           'baker_02', ...
          'd0327_8',                'd0327_2', ...
          'd0403_8',                'd0403_2', ...
          'd0410_8',                'd0410_2', ...
          'd0327_pnb',              'd0403_pnb',         'd0410_pnb'}, ...
 'day',  {'Baker',                  'Baker',             'Baker', ...
          '03/27 M',                '03/27 M', ...
          '04/03 F',                '04/03 F', ...
          '04/10 M2-8',             '04/10 M2-8', ...
          '03/27 M',                '04/03 F',           '04/10 M2-8'}, ...
 'ATP',  {8, 2, 0.2,   8, 2,   8, 2,   8, 2,   8, 8, 8}, ...
 'kind', {'active','active','active', 'active','active', 'active','active', ...
          'active','active', 'passive','passive','passive'}, ...
 'file', {'bakers_slack8mM_all.mat', 'bakers_slack2mM.mat', 'bakers_slack02mM.mat', ...
          'protocol_03_27_2026_8mM_slack.mat', 'protocol_03_27_2026_2mM_slack.mat', ...
          'protocol_04_03_2026_8mM_slack.mat', 'protocol_04_03_2026_2mM_slack.mat', ...
          'protocol_04_10_2026_8mM_slack.mat', 'protocol_04_10_2026_2mM_slack.mat', ...
          'protocol_03_27_2026_ActivePNBMava_slack.mat', ...
          'protocol_04_03_2026_ActivePNBMava_slack.mat', ...
          'protocol_04_10_2026_ActivePNBMava_slack.mat'});

for i = 1:numel(DS)
    if DS(i).ATP == 8
        DS(i).label = sprintf('%s  %g mM', DS(i).day, DS(i).ATP);
    else
        DS(i).label = sprintf('%s  %g mM', DS(i).day, DS(i).ATP);
    end
    if strcmp(DS(i).kind, 'passive'); DS(i).label = [DS(i).day '  PNB+Mava']; end
end

days      = unique({DS.day}, 'stable');
dayColors = lines(numel(days));

%% ── 1. Load + extract features (one code path for everything) ──────────────
fprintf('=== Feature extraction ===\n');
for i = 1:numel(DS)
    S = load([dataDir DS(i).file]);
    DS(i).datatable     = S.datatable;
    DS(i).velocitytable = S.velocitytable;

    % first slack release: the common time origin for the trace overview
    iS = find(S.velocitytable(:,2) < -1, 1, 'first');
    DS(i).t0slack = S.velocitytable(iS, 1);
    DS(i).nSlack  = numel(find(S.velocitytable(:,2) < -1));

    fig = figure(900); clf; set(fig, 'Position', [80 80 1000 640]);
    DS(i).feats = extractSlackAttributes(S.datatable(:,1), S.datatable(:,3), ...
        S.datatable(:,2), S.velocitytable, [], [], true, true);
    title(sprintf('Slack feature extraction — %s   [%s]', DS(i).label, DS(i).file), ...
          'Interpreter', 'none');
    xlabel('Time (s)'); ylabel('Force (kPa)');
    exportgraphics(fig, fullfile(resDir, sprintf('fit_%s.png', DS(i).tag)), 'Resolution', 150);

    fprintf('  %-16s %-46s nSlack=%d  A=%s\n', DS(i).label, DS(i).file, DS(i).nSlack, ...
        mat2str(round(iGet(DS(i).feats,'A'), 1)));
end

%% ── 2. Trace overview, one tile per day, aligned on the first slack ────────
fig = figure(901); clf;
set(fig, 'Position', [30 30 1500 760], 'Name', 'Slack traces — all datasets');
tiledlayout(2, numel(days), 'TileSpacing', 'compact', 'Padding', 'compact');

for d = 1:numel(days)
    sel = find(strcmp({DS.day}, days{d}) & strcmp({DS.kind}, 'active'));
    nexttile(d); hold on; box on;
    for k = sel
        tt = DS(k).datatable(:,1) - DS(k).t0slack;
        plot(tt, DS(k).datatable(:,3), 'LineWidth', 0.8, ...
             'LineStyle', iStyle(DS(k).ATP), 'Color', dayColors(d,:), ...
             'DisplayName', sprintf('%g mM', DS(k).ATP));
    end
    xlim([-0.35, 3.1]); title(days{d}); ylabel('Force (kPa)');
    legend('Location', 'southwest', 'FontSize', 7); ylim([-10 130]);

    % passive control below
    selp = find(strcmp({DS.day}, days{d}) & strcmp({DS.kind}, 'passive'));
    nexttile(numel(days) + d); hold on; box on;
    for k = selp
        tt = DS(k).datatable(:,1) - DS(k).t0slack;
        plot(tt, DS(k).datatable(:,3), 'LineWidth', 0.8, 'Color', [0.4 0.4 0.4], ...
             'DisplayName', 'PNB+Mava');
    end
    if isempty(selp)
        text(0.5, 0.5, 'no passive control', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
    else
        legend('Location', 'northwest', 'FontSize', 7);
    end
    xlim([-0.35, 3.1]); ylim([-5 30]);
    ylabel('Passive force (kPa)'); xlabel('t - t_{first slack} (s)');
end
sgtitle('Slack force traces — every dataset, aligned on the first slack release');
exportgraphics(fig, fullfile(resDir, 'overview_traces.png'), 'Resolution', 150);

%% ── 3. Combined feature comparison ─────────────────────────────────────────
fn_slack = {'ktr|SLslack', 'ktr2_k|SLslack', 'A|SLslack', 'Am|SLslack', ...
            'ktr2_omega', 't0', 'restretchSlopeStart', 'restretchSlopeEnd', ...
            'peak1_y', 'peak1_t', 'peak1_dSL', 'vall_y', 'vall_t', 'peak2', ...
            'steady', 'vall2_dy', 'vall2_t', 'ovrsht_dy', 'ovrsht_t'};

act = find(strcmp({DS.kind}, 'active'));
feats  = {DS(act).feats};
labels = {DS(act).label};
colors = zeros(numel(act), 3);
styles = cell(1, numel(act));
marks  = cell(1, numel(act));
fills  = false(1, numel(act));
for j = 1:numel(act)
    d = find(strcmp(days, DS(act(j)).day));
    colors(j,:) = dayColors(d,:);
    styles{j}   = iStyle(DS(act(j)).ATP);
    marks{j}    = iMarker(DS(act(j)).ATP);
    fills(j)    = DS(act(j)).ATP == 8;      % filled = high ATP
end

fig = figure(902); clf;
set(fig, 'Units', 'centimeters', 'Position', [1 1 40 26], 'Name', 'Slack features — all datasets');
sgtitle('Slack features — all datasets.  Colour = protocol day.  SOLID/filled = 8 mM, DASHED = 2 mM, DASH-DOT = 0.2 mM');
plotMultipleFeatures(feats, labels, colors, marks, fn_slack, styles, fills);
exportgraphics(fig, fullfile(resDir, 'features_all.png'), 'Resolution', 150);

%% ── 3b. Same comparison, force features normalised per prep ────────────────
% Absolute force varies ~20% CV between preps (fibre cross-section, damage,
% rundown), which swamps the panels. Dividing every force-dimension feature by
% that dataset's own mean `steady` asks a different, more useful question: is
% the SHAPE of the slack response reproducible once prep strength is divided
% out? Whatever collapses here is transferable between preps and is what a
% cost function should target; whatever stays scattered is prep-specific.
forceFeats = {'A', 'Am', 'steady', 'peak1_y', 'peak2', 'vall_y', 'vall2_dy', ...
              'restretchSlopeStart', 'ovrsht_dy'};
featsN = feats;
for j = 1:numel(act)
    sc = mean(iGet(DS(act(j)).feats, 'steady'), 'omitnan');
    for k = 1:numel(forceFeats)
        if isfield(featsN{j}, forceFeats{k})
            featsN{j}.(forceFeats{k}) = featsN{j}.(forceFeats{k}) / sc;
        end
    end
end

fn_norm = {'A|SLslack', 'Am|SLslack', 'steady', 'peak1_y', 'peak2', 'vall_y', ...
           'vall2_dy', 'restretchSlopeStart', 'ktr|SLslack', 't0'};
fig = figure(903); clf;
set(fig, 'Units', 'centimeters', 'Position', [1 1 36 20], 'Name', 'Slack features — normalised');
sgtitle('Slack features with force features divided by each prep''s own mean steady force (ktr, t0 left absolute)');
plotMultipleFeatures(featsN, labels, colors, marks, fn_norm, styles, fills);
exportgraphics(fig, fullfile(resDir, 'features_normalised.png'), 'Resolution', 150);

%% ── 4. Consistency tables ──────────────────────────────────────────────────
keys = {'A', 'Am', 'steady', 'peak1_y', 'peak2', 'vall_y', 'ktr', 'ktr2_k', ...
        't0', 'restretchSlopeStart', 'vall2_dy', 'ovrsht_dy', 'SLslack'};

fprintf('\n\n============ ABSOLUTE LEVELS (mean over the 5 slacks) ============\n');
fprintf('%-22s', 'feature');
for j = 1:numel(act); fprintf('%14s', DS(act(j)).label); end
fprintf('\n');
Vals = nan(numel(keys), numel(act));
for k = 1:numel(keys)
    fprintf('%-22s', keys{k});
    for j = 1:numel(act)
        Vals(k,j) = mean(iGet(DS(act(j)).feats, keys{k}), 'omitnan');
        fprintf('%14.4g', Vals(k,j));
    end
    fprintf('\n');
end

% within-day 2 mM / 8 mM ratio -- cancels prep-to-prep force scale
pairDays = {}; ratios = [];
for d = 1:numel(days)
    ih = find(strcmp({DS.day}, days{d}) & [DS.ATP] == 8 & strcmp({DS.kind}, 'active'));
    il = find(strcmp({DS.day}, days{d}) & [DS.ATP] == 2 & strcmp({DS.kind}, 'active'));
    if isempty(ih) || isempty(il); continue; end
    pairDays{end+1} = days{d}; %#ok<SAGROW>
    for k = 1:numel(keys)
        ratios(k, numel(pairDays)) = mean(iGet(DS(il).feats, keys{k}), 'omitnan') ...
                                   / mean(iGet(DS(ih).feats, keys{k}), 'omitnan'); %#ok<SAGROW>
    end
end

fprintf('\n\n============ ATP EFFECT WITHIN EACH PREP (2 mM / 8 mM) ============\n');
fprintf('ATP order:  Baker ?    03/27 = 8->2    04/03 = 8->2    04/10 = 2->8 (REVERSED)\n\n');
fprintf('%-22s', 'feature');
for d = 1:numel(pairDays); fprintf('%12s', pairDays{d}); end
fprintf('%12s%12s\n', 'mean', 'CV');
for k = 1:numel(keys)
    fprintf('%-22s', keys{k});
    for d = 1:numel(pairDays); fprintf('%12.3f', ratios(k,d)); end
    m = mean(ratios(k,:), 'omitnan');
    fprintf('%12.3f%11.0f%%\n', m, 100*std(ratios(k,:), 'omitnan')/abs(m));
end

% between-prep spread at fixed ATP -- how reproducible is a prep?
fprintf('\n\n============ BETWEEN-PREP SPREAD AT FIXED ATP (CV over days) ============\n');
fprintf('%-22s %10s %10s %10s %10s\n', 'feature', 'CV 8mM', 'min 8mM', 'max 8mM', 'CV 2mM');
i8 = find([DS(act).ATP] == 8);  i2 = find([DS(act).ATP] == 2);
for k = 1:numel(keys)
    v8 = Vals(k, i8);  v2 = Vals(k, i2);
    fprintf('%-22s %9.0f%% %10.4g %10.4g %9.0f%%\n', keys{k}, ...
        100*std(v8,'omitnan')/abs(mean(v8,'omitnan')), min(v8), max(v8), ...
        100*std(v2,'omitnan')/abs(mean(v2,'omitnan')));
end

%% ── 5. Protocol geometry: is it even the same experiment? ──────────────────
% Amplitude features are only comparable between datasets if the protocol that
% produced them is the same. Slack DEPTH (SLslack) is shared by construction,
% but the recovery time the fibre is given before the restretch is not, and
% that is what A / Am measure.
fprintf('\n\n============ PROTOCOL GEOMETRY ============\n');
fprintf('%-16s %-26s %-22s %-24s %s\n', 'dataset', 'inter-slack dt (s)', ...
        'slack hold (s)', 'release vel (ML/s)', 'restretch vel (ML/s)');
for j = 1:numel(act)
    vt = DS(act(j)).velocitytable;
    is = find(vt(:,2) < -1);  ir = find(vt(:,2) > 1);
    fprintf('%-16s %-26s %-22s %-24s %s\n', DS(act(j)).label, ...
        mat2str(round(diff(vt(is,1))', 2)), ...
        mat2str(round(vt(is+2,1)' - vt(is+1,1)', 3)), ...
        mat2str(round(vt(is,2)', 0)), ...
        mat2str(round(vt(ir,2)', 1)));
end

%% ── 6. Recovery-truncation diagnostic ──────────────────────────────────────
% A  = fitted asymptote of F(t) = A(1-exp(-ktr(t-t0))) over the slack hold.
% Am = the LAST MEASURED force in that same window.
% So Am/A is the fraction of the asymptote the fibre actually reached before
% the restretch cut the recovery short. A hold that is long relative to 1/ktr
% gives Am/A -> 1 and an unbiased amplitude; a short one biases the MEASURED
% amplitude down, and biases it down HARDER wherever ktr is slower -- i.e.
% preferentially at low ATP. Any low-ATP "force loss" seen in A/Am but not in
% steady/peak2 (which are read after the restretch, long after) is this artifact.
fprintf('\n\n============ RECOVERY TRUNCATION (Am/A, and ktr x hold) ============\n');
fprintf('%-16s %-34s %-34s %s\n', 'dataset', 'Am/A per slack', ...
        'ktr*(hold-t0)  [>3 = complete]', 'hold (s)');
for j = 1:numel(act)
    f  = DS(act(j)).feats;
    vt = DS(act(j)).velocitytable;
    is = find(vt(:,2) < -1);
    hold_s = vt(is+2,1)' - vt(is+1,1)';
    AmA    = iGet(f,'Am') ./ iGet(f,'A');
    kT     = iGet(f,'ktr') .* max(hold_s - iGet(f,'t0'), 0);
    fprintf('%-16s %-34s %-34s %s\n', DS(act(j)).label, ...
        mat2str(round(AmA, 2)), mat2str(round(kT, 1)), mat2str(round(hold_s, 3)));
end

feats_all = rmfield(DS, {'datatable', 'velocitytable'});
save(fullfile(resDir, 'features_all.mat'), 'feats_all', 'keys', 'Vals', 'ratios', 'pairDays');
fprintf('\nSaved %s\n', fullfile(resDir, 'features_all.mat'));

% ---------------------------------------------------------------------------
function s = iStyle(atp)
%ISTYLE Line style encodes ATP: solid = high, dashed = low, dash-dot = very low.
    if     atp >= 8;  s = '-';
    elseif atp >= 2;  s = '--';
    else;             s = '-.';
    end
end

function m = iMarker(atp)
    if     atp >= 8;  m = 'o';
    elseif atp >= 2;  m = 's';
    else;             m = 'd';
    end
end

function v = iGet(s, f)
%IGET Feature vector, or NaN when the field is absent / non-numeric.
    if ~isfield(s, f); v = NaN; return; end
    v = s.(f);
    if ~isnumeric(v) || isempty(v); v = NaN; end
end
