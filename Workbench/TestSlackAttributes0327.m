%% TestSlackAttributes0327.m
% Compare extractSlackAttributes on:
%   A) Baker lab reference data  (bakers_slack8mM_all.mat)
%   B) 03/27/2026 merged 8mM    (02_Merged_8mM_Active.txt + protocol_03_27_2026_velocitytable.mat)
%
% Velocity table columns: [t(s), vel(Lo/s or um/s), vel2, L]
% extractSlackAttributes uses col2 < 0 for segment anchors and linear
% indexing into col1 for phase boundaries — Nx2 and Nx4 both work.
%
% SL units differ: Baker = um, new data = Lo.
% Converted for comparison using Lo_ref = 2.2 um (Baker isometric SL).

close all; clc;
addpath(genpath('..'));

LO_REF_UM = 2.2;  % reference SL (um) for Lo -> um conversion

%% ── A: Baker reference ──────────────────────────────────────────────────
bs    = load('../data/bakers_slack8mM_all.mat');
vt_A  = bs.velocitytable;               % 23x4, SL in um, vel in um/s
dt_A  = bs.datatable;                   % 4000x3: [t, SL(um), F(kPa)]

figure(401); clf; set(gcf, 'Name', 'Baker – fitting', 'Units', 'centimeters', 'Position', [1 1 20 12]);
feats_A = extractSlackAttributes(dt_A(:,1), dt_A(:,3), dt_A(:,2), vt_A, struct(), [], true);
title('Baker 8mM – slack fits');

%% ── B: 03/27/2026 merged 8mM ────────────────────────────────────────────
vts   = load('../data/protocol_03_27_2026_velocitytable.mat');
vt_full = vts.velocitytable;            % 69x4, L in Lo, vel in Lo/s

% Drop near-zero noise negatives (|vel| < 0.1 Lo/s) between real events
noise_rows = vt_full(:,2) < 0 & vt_full(:,2) > -0.1;
vt_clean   = vt_full(~noise_rows, [1 2]);
vt_B       = vt_clean(vt_clean(:,1) >= 74.4, :);   % slack window only

merged  = readmatrix('../data/03 27 2026 M/02_Merged_8mM_Active.txt');
mask_B  = merged(:,1) >= 74.0 & merged(:,1) <= 77.5;
t_B     = merged(mask_B, 1);
L_B     = merged(mask_B, 2);   % Lo
F_B     = merged(mask_B, 3);   % kPa

figure(402); clf; set(gcf, 'Name', '03/27 8mM – fitting', 'Units', 'centimeters', 'Position', [1 1 20 12]);
feats_B = extractSlackAttributes(t_B, F_B, L_B, vt_B, struct(), [], true);
title('03/27/2026 8mM – slack fits');

%% ── Trim to matched valid segments (first 4) ────────────────────────────
% Baker has 5 segments, new data has 5 but seg 5 is corrupt (spurious anchor).
% Use first 4 of each for clean comparison.
N = 4;
trim = @(x) x(1:min(N, numel(x)));

fn = fieldnames(feats_A);
for k = 1:numel(fn)
    feats_A.(fn{k}) = trim(feats_A.(fn{k}));
end
fn = fieldnames(feats_B);
for k = 1:numel(fn)
    feats_B.(fn{k}) = trim(feats_B.(fn{k}));
end

% Convert new-data SL from Lo -> um for apples-to-apples comparison
SLslack_B_um  = feats_B.SLslack  * LO_REF_UM;
SLdiff_B_um   = feats_B.SLdiff   * LO_REF_UM;
peak1SL_B_um  = feats_B.peak1_SL * LO_REF_UM;
peak1dSL_B_um = feats_B.peak1_dSL* LO_REF_UM;

%% ── Comparison figure ───────────────────────────────────────────────────
clr_A = [0.00, 0.45, 0.70];   % blue  – Baker
clr_B = [0.84, 0.10, 0.11];   % red   – 03/27

seg_A = 1:numel(feats_A.ktr);
seg_B = 1:numel(feats_B.ktr);

figure(403); clf;
set(gcf, 'Name', 'Feature comparison: Baker vs 03/27', 'Units', 'centimeters', 'Position', [2 2 24 20]);
tl = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Slack feature comparison — Baker 8mM vs 03/27/2026 8mM');

specs = { ...
    feats_A.ktr,        feats_B.ktr,        '1/s',  'ktr'; ...
    feats_A.A,          feats_B.A,          'kPa',  'A (fit amplitude)'; ...
    feats_A.t0*1000,    feats_B.t0*1000,    'ms',   't_0 (onset delay)'; ...
    feats_A.SLslack,    SLslack_B_um,       'µm',   'SL_{slack}'; ...
    feats_A.SLdiff,     SLdiff_B_um,        'µm',   'SL_{diff}'; ...
    feats_A.peak1_y,    feats_B.peak1_y,    'kPa',  'Peak 1 force'; ...
    feats_A.peak1_t*1000, feats_B.peak1_t*1000, 'ms', 'Peak 1 time'; ...
    feats_A.vall2_dy,   feats_B.vall2_dy,   'kPa',  'Valley 2 \Deltaforce'; ...
    feats_A.v_restretch, feats_B.v_restretch, 'units', 'v_{restretch}'; ...
};

for ii = 1:size(specs,1)
    ax = nexttile(ii); hold on;
    plot(ax, seg_A, specs{ii,1}, 'o-', 'Color', clr_A, 'LineWidth', 1.5, 'MarkerFaceColor', clr_A, 'DisplayName', 'Baker');
    plot(ax, seg_B, specs{ii,2}, 's-', 'Color', clr_B, 'LineWidth', 1.5, 'MarkerFaceColor', clr_B, 'DisplayName', '03/27');
    ylabel(ax, specs{ii,3}); title(ax, specs{ii,4}); box(ax, 'on');
    xticks(ax, 1:max(numel(seg_A), numel(seg_B))); xlabel(ax, 'Segment');
    if ii == 1; legend(ax, 'Location', 'best', 'FontSize', 7); end
end

%% ── Print table ─────────────────────────────────────────────────────────
fprintf('\n%s\n', repmat('-', 1, 72));
fprintf('%-22s  %-28s  %s\n', 'Feature', 'Baker (segs 1-4)', '03/27 (segs 1-4)');
fprintf('%s\n', repmat('-', 1, 72));

compare_fields = {'ktr','A','t0','Am','SLslack','SLdiff','peak1_y','peak1_t','vall2_dy','v_restretch'};
for k = 1:numel(compare_fields)
    fname = compare_fields{k};
    vA = feats_A.(fname);
    vB = feats_B.(fname);
    suffix = '';
    if strcmp(fname,'SLslack'); vB = SLslack_B_um; suffix = ' [Lo->um]'; end
    if strcmp(fname,'SLdiff');  vB = SLdiff_B_um;  suffix = ' [Lo->um]'; end
    fprintf('%-22s  %-28s  %s%s\n', fname, mat2str(round(vA,4,'significant')), mat2str(round(vB,4,'significant')), suffix);
end
fprintf('%s\n', repmat('-', 1, 72));
