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

LO_REF_UM = 2.0;  % reference SL (um) for Lo -> um conversion

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

% The first near-zero-neg row stops the movement (keep).
% A second consecutive near-zero-neg is a duplicate → remove.
% Exact-zero sentinels are always preserved.
is_nzn_B   = vt_full(:,2) < 0 & vt_full(:,2) > -0.1;
noise_rows = is_nzn_B & [false; is_nzn_B(1:end-1)];   % second-of-pair only
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

%% ── C: 04/03/2026 merged 8mM ────────────────────────────────────────────
vts_C      = load('../data/protocol_04_03_2026_velocitytable.mat');
vt_full_C  = vts_C.velocitytable;            % 69x4, L in Lo, vel in Lo/s

is_nzn_C     = vt_full_C(:,2) < 0 & vt_full_C(:,2) > -0.1;
noise_rows_C = is_nzn_C & [false; is_nzn_C(1:end-1)];   % second-of-pair only
vt_clean_C   = vt_full_C(~noise_rows_C, [1 2]);

% Skip the initial slack-constant vel-slack preamble (3 shortening entries at
% 74.43-74.47 s that use a different protocol structure; first real slack at 74.88)
vt_C         = vt_clean_C(vt_clean_C(:,1) >= 74.8, :);   % slack window only

merged_C = readmatrix('../data/04 03 2026 F/02_Merged_8mM_Active.txt');
mask_C   = merged_C(:,1) >= 74.0 & merged_C(:,1) <= 78;
t_C      = merged_C(mask_C, 1);
L_C      = merged_C(mask_C, 2);   % Lo
F_C      = merged_C(mask_C, 3);   % kPa

figure(404); clf; set(gcf, 'Name', '04/03 8mM – fitting', 'Units', 'centimeters', 'Position', [1 1 20 12]);
feats_C = extractSlackAttributes(t_C, F_C, L_C, vt_C, struct(), [], true);
title('04/03/2026 8mM – slack fits');

%% ── Convert Lo -> um for SL features (Baker uses um, new data uses Lo) ──
SLslack_B_um = feats_B.SLslack * LO_REF_UM;
SLdiff_B_um  = feats_B.SLdiff  * LO_REF_UM;
SLslack_C_um = feats_C.SLslack * LO_REF_UM;
SLdiff_C_um  = feats_C.SLdiff  * LO_REF_UM;

%% ── Universal comparison figure & table ─────────────────────────────────
clr_A = [0.00, 0.45, 0.70];   % blue  – Baker
clr_B = [0.84, 0.10, 0.11];   % red   – 03/27
clr_C = [0.47, 0.67, 0.19];   % green – 04/03
datasets = { feats_A, feats_B, feats_C };
labels   = { 'Baker', '03/27', '04/03' };
colors   = { clr_A,   clr_B,   clr_C  };
markers  = { 'o',     's',     '^'    };

% Lo->um fields: B and C are in Lo, A is already in um
lo_fields = {'SLslack', 'SLdiff'};

% Collect all field names present in any dataset
all_fields = fieldnames(feats_A);
for di = 2:numel(datasets)
    fn = fieldnames(datasets{di});
    all_fields = union(all_fields, fn, 'stable');
end

% Separate plottable (numeric, non-empty) from non-plottable fields
plot_fields = {};
for k = 1:numel(all_fields)
    fname = all_fields{k};
    for di = 1:numel(datasets)
        if isfield(datasets{di}, fname)
            v = datasets{di}.(fname);
            if isnumeric(v) && ~isempty(v)
                plot_fields{end+1} = fname; %#ok<SAGROW>
                break;
            end
        end
    end
end

% Auto-size figure grid
ncols = ceil(sqrt(numel(plot_fields)));
nrows = ceil(numel(plot_fields) / ncols);

figure(403); clf;
set(gcf, 'Name', 'Feature comparison: Baker vs 03/27 vs 04/03', 'Units', 'centimeters', ...
    'Position', [2 2 ncols*8 nrows*6]);
tl = tiledlayout(nrows, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Slack feature comparison — Baker 8mM vs 03/27/2026 vs 04/03/2026');

for ii = 1:numel(plot_fields)
    fname = plot_fields{ii};
    ax = nexttile(ii); hold on;
    for di = 1:numel(datasets)
        if ~isfield(datasets{di}, fname); continue; end
        v = datasets{di}.(fname);
        if ~isnumeric(v) || isempty(v); continue; end
        % Apply Lo->um scaling for B and C on SL fields
        if di > 1 && ismember(fname, lo_fields)
            v = v * LO_REF_UM;
        end
        segs = 1:numel(v);
        plot(ax, segs, v, [markers{di} '-'], 'Color', colors{di}, ...
            'LineWidth', 1.5, 'MarkerFaceColor', colors{di}, 'DisplayName', labels{di});
    end
    title(ax, strrep(fname, '_', '\_')); box(ax, 'on');
    xlabel(ax, 'Segment');
    if ismember(fname, lo_fields); ylabel(ax, 'µm'); end
    if ii == 1; legend(ax, 'Location', 'best', 'FontSize', 7); end
end

%% ── Print table (all fields) ────────────────────────────────────────────
COL = 28;
sep = repmat('-', 1, 22 + (COL+2)*numel(datasets) + 2);
fprintf('\n%s\n', sep);
fprintf('%-22s', 'Feature');
for di = 1:numel(datasets); fprintf('  %-*s', COL, labels{di}); end
fprintf('\n%s\n', sep);

for k = 1:numel(all_fields)
    fname = all_fields{k};
    fprintf('%-22s', fname);
    for di = 1:numel(datasets)
        if ~isfield(datasets{di}, fname)
            fprintf('  %-*s', COL, '—');
            continue;
        end
        v = datasets{di}.(fname);
        if ischar(v) || isstring(v)
            fprintf('  %-*s', COL, v);
        elseif isnumeric(v) && ~isempty(v)
            if di > 1 && ismember(fname, lo_fields)
                v = v * LO_REF_UM;
            end
            s = mat2str(round(v, 4, 'significant'));
            fprintf('  %-*s', COL, s);
        else
            fprintf('  %-*s', COL, '[]');
        end
    end
    if ismember(fname, lo_fields); fprintf('  [Lo->um]'); end
    fprintf('\n');
end
fprintf('%s\n', sep);
