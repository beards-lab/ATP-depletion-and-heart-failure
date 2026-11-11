%% Load, merge and plot all Log_*.txt files from a data directory
%
% The merge/alignment itself lives in mergeLogsAndBursts.m so it can be run for
% any protocol day without editing this script; what remains here is the QC
% plotting and the exploratory rundown analysis in the graveyard below.
close all; clear; clc;

% ── Which protocol day to process ────────────────────────────────────────────
% dataDir = '../data/03 27 2026 M/';
% dataDir = '../data/04 03 2026 F/';
dataDir = '../data/04 10 2026 Male 2-8/';

mergeOpts = struct( ...
    'rewriteData',  false, ...   % true => overwrite NN_Merged_*.txt in dataDir
    'saveAllData',  false, ...   % true => write AllDataMerged.mat
    'plotAll',      false, ...   % per-burst alignment diagnostics
    'plotFineAlign', true);      % fine-alignment overlays + MSE curves

% The 04-series slack file opens with an extra FV2 isovelocity ramp at ~74.4 s,
% before the first slack; the default xcorr window would clip it.
if contains(dataDir, '04 03') || contains(dataDir, '04 10')
    mergeOpts.winSlack = [74.2, 78.5];
end

[data, hi_res_ranges, S] = mergeLogsAndBursts(dataDir, mergeOpts);

n      = length(S);
colors = lines(n);
labels = cellfun(@(nm) strrep(nm(1:end-4), '_', ' '), {S.name}', 'UniformOutput', false);

%% Plot Merged Traces vs Original Logs
figure(101);clf;
for i = 1:n
    figure(101);
    ax1_1 = nexttile(1);hold on;
    plot(data(i).t_fineShifted, data(i).F_fineShifted, 'LineWidth', 1);

    ax1_2 = nexttile(2);hold on;
    plot(data(i).t_fineShifted, data(i).L_fineShifted, 'LineWidth', 1);
    linkaxes([ax1_1 ax1_2], 'x');


    figure(i+10); clf;

    ax1 = subplot(2,1,1); hold on;
    % Plot merged trace (thin gray)
    plot(ax1, data(i).t_fineShifted, data(i).F_fineShifted, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    % Plot FULL original log on top (thick colored) — always visible
    plot(ax1, data(i).t_shifted_orig, data(i).F_shifted_orig, 'Color', colors(i,:), 'LineWidth', 2);
    ylabel(ax1, 'Force (kPa)');

    outNameRaw = strrep(S(i).name, 'Log_', 'Merged_');
    title_str = ['Dataset: ', strrep(strrep(outNameRaw, '_', ' '), '.txt', '')];
    title(ax1, title_str);

    ax2 = subplot(2,1,2); hold on;
    % Plot merged trace (thin gray)
    plot(ax2, data(i).t_fineShifted, data(i).L_fineShifted, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    % Plot FULL original log on top (thick colored)
    plot(ax2, data(i).t_shifted_orig, data(i).L_shifted_orig, 'Color', colors(i,:), 'LineWidth', 2);
    ylabel(ax2, 'Length (Lo)');
    xlabel(ax2, 'Time (s)');

    % Mark hi-res injection boundaries with vertical lines + label
    for jj = 1:length(hi_res_ranges{i})
        hr = hi_res_ranges{i}{jj};
        xline(ax1, hr.t_start, '--g', 'LineWidth', 1.5);
        xline(ax1, hr.t_end,   '--g', 'LineWidth', 1.5);
        xline(ax2, hr.t_start, '--g', 'LineWidth', 1.5);
        xline(ax2, hr.t_end,   '--g', 'LineWidth', 1.5);
        % Annotate the burst name at the top of F axis
        yl = ylim(ax1);
        text(ax1, (hr.t_start + hr.t_end)/2, yl(2), ...
            strrep(hr.name, '_', ' '), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 7, 'Color', [0 0.5 0], 'Interpreter', 'none');
    end

    legend(ax2, {'Merged (with Hi-Res)', 'Original Log'}, 'Location', 'best');
    linkaxes([ax1 ax2], 'x');
end

return
%% graveyeard below

%% compensate for run-down


refData = data(2);
t_ref = refData.t_fineShifted;
F_ref = refData.F_fineShifted;
L_ref = refData.L_fineShifted;
t_ref = refData.t_shifted_orig;
F_ref = refData.F_shifted_orig;
L_ref = refData.L_shifted_orig;


% run-down data
rundownData = data(4);
F_rd = interp1(rundownData.t_fineShifted, rundownData.F_fineShifted, t_ref);
L_rd = interp1(rundownData.t_fineShifted, rundownData.L_fineShifted, t_ref);

% passive reference
F_pas = 0; %interp1(data(5).t_fineShifted, data(5).F_fineShifted, t_ref);
F_rdact = F_rd - F_pas;
F_refact = F_ref - F_pas;

L0 = 1.0;
k  = -0.6;   % SL slope: increase until high-SL portion of corrected trace aligns with F_refact
r0 = (68)/(56);
r0=1.214286;
% r828 = 80.09/71;
F_rdact_comp = scaleRundown(F_rdact, L_ref, r0, L0, k);
t_rd_comp = rundownData.t_fineShifted;
F_rd_comp = interp1(t_ref, F_rdact_comp + F_pas, t_rd_comp);
%%
figure(297);clf;
% plot(t_ref, F_refact, t_ref, F_rdact*r0, t_ref, F_rdact_comp,'LineWidth', 1);
ax1 = nexttile(1);plot(t_ref, F_refact + F_pas, t_rd_comp, F_rd_comp,'LineWidth', 1);
% ax2 = nexttile(2);plot(t_ref, L_ref, t_ref, L_rd,'LineWidth', 1);
% linkaxes([ax1, ax2], 'x')
xlim([66.6300   78.5624])
fprintf('With k = %g;r0=%f; cost = %g\n', k, r0, nansum(F_rdact_comp - F_refact));
legend('First', 'Third rescaled with f(L)');
xy_F = linspace(1, 100, 10);%[1:100];
xy_L = linspace(0.6, 1.2, 10);%:0.1:1.2];
sr = scaleRundown(xy_F, xy_L', r0, L0, k)';
% nexttile(2);cla;surf(xy_F, xy_L, sr);
% nexttile(2); plot(xy_F', scaleRundown(xy_F, L_ref, r0, L0, k)');
% plot(t_ref, F_refact, t_ref, F_rdact)
% F_rdact
%% rundown corrig
t = refData.t_shifted;
F = refData.F_shifted;
L = refData.L_shifted;
zones = [71.5, 72.4;77.4, 78.4];
% at 2.0
mask = any(t >= zones(:,1)' & t <= zones(:,2)', 2);
p1 = polyfit(t(mask), F(mask), 1);

% at 2.2
L_smooth = movmean(L, 20, 'omitnan');
falling = find(diff(L_smooth >= 1.1) == -1) + 1;  % first sample below threshold
t_drops = t(falling);
mask = false(size(t));
for i = 1:numel(t_drops)
    mask = mask | (t >= t_drops(i) - 0.15 & t < t_drops(i));
end
p2 = polyfit(t(mask), F(mask), 1);
clf;
plot(t, F, t(mask), F(mask), t, polyval(p1, t), t, polyval(p2, t), linewidth = 1)

%%

%%
% plot(rundownData.t_shifted, rundownData.F_shifted*r828, 'LineWidth', 1);

ax2 = nexttile(2);hold on;
plot(refData.t_shifted, refData.L_shifted, 'LineWidth', 1);
plot(data(3).t_shifted, data(3).L_shifted, 'LineWidth', 1);

%% Report Figure A — Comparison 2, upscaled-4, 3, 5
clr2 = [0.00, 0.45, 0.70];  % blue:   8mM Active
clr3 = [0.84, 0.10, 0.11];  % red:    2mM Active
clr4 = [0.49, 0.18, 0.56];  % purple: 8mM repeat upscaled
clr5 = [0.20, 0.63, 0.17];  % green:  PNB+Mava reference

% Upscale i=4 to match i=2 using scaleRundown
t_ref2 = refData.t_fineShifted;
F_ref2 = refData.F_fineShifted;
L_ref2 = refData.L_fineShifted;
F_rd4  = interp1(rundownData.t_fineShifted, rundownData.F_fineShifted, t_ref2, 'linear', NaN);
F_pas5 = interp1(data(5).t_fineShifted, data(5).F_fineShifted, t_ref2, 'linear', NaN);
F_rdact4      = F_rd4 - F_pas5;
r0_up = 68/56;  L0_up = 1.0;  k_up = -0.8;
F_rdact4_comp = scaleRundown(F_rdact4, L_ref2, r0_up, L0_up, k_up);
F_4_upscaled  = F_rdact4_comp + F_pas5;  % total force, rundown-compensated

zoom_xlims  = {[66, 79], [69.5, 74], [74, 77.5], [66, 79]};
zoom_titles = {'Overview (66–79 s)', 'ktr (69.5–74 s)', 'Slack (74–77.5 s)', 'Length protocol'};

fA = figure('Name', 'Report_Comparison_235', 'Units', 'centimeters', 'Position', [2 2 20 24]);
tlA = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
axA = gobjects(4, 1);

for ii = 1:4
    axA(ii) = nexttile(ii); hold on;
    if ii < 4
        plot(refData.t_fineShifted, refData.F_fineShifted, 'Color', clr2, 'LineWidth', 1.2, 'DisplayName', '8mM Active (i=2)');
        plot(t_ref2, F_4_upscaled,                         '--', 'Color', clr4, 'LineWidth', 1.0, 'DisplayName', '8mM repeat upscaled (i=4)');
        plot(data(3).t_fineShifted, data(3).F_fineShifted, 'Color', clr3, 'LineWidth', 1.2, 'DisplayName', '2mM Active (i=3)');
        plot(data(5).t_fineShifted, data(5).F_fineShifted, 'Color', clr5, 'LineWidth', 1.2, 'DisplayName', 'PNB+Mava (i=5)');
        ylabel('Force (kPa)');
        if ii == 1; legend('Location', 'northwest', 'FontSize', 7); end
    else
        plot(data(5).t_fineShifted, data(5).L_fineShifted, 'Color', clr5, 'LineWidth', 1.2);
        ylabel('Length (Lo)'); xlabel('Time (s)');
    end
    xlim(zoom_xlims{ii});
    title(zoom_titles{ii});
    xline(72.0, '--k', 'Alpha', 0.3, 'HandleVisibility', 'off');
    xline(74.5, '--k', 'Alpha', 0.3, 'HandleVisibility', 'off');
    box on;
end

title(tlA, '03/27/2026 — Force comparison: 8mM vs 2mM vs PNB+Mava');
exportgraphics(fA, fullfile(dataDir, 'Report_Comparison_235.png'), 'Resolution', 300);
fprintf('Saved Report_Comparison_235.png\n');

%% Report Figure B — Rundown linear fits i=2 and i=3
ds_order = [3, 2, 4];
ds_labels = {'2mM Active (i=%g)', '8mM Active (i=%g)', '8mM Active repeat (i=%g)', 'passive(i=%g)'};

fB = figure('Name', 'Report_Rundown_23', 'Units', 'centimeters', 'Position', [2 2 18 16]);

for ii_loop = 1:length(ds_order)
    ds_idx = ds_order(ii_loop);
    t = data(ds_idx).t_fineShifted;
    F = data(ds_idx).F_fineShifted;
    L = data(ds_idx).L_fineShifted;

    % Mask 1: SL=2.0 plateau zones (stable force before each L-drop)
    zones = [71.5, 72.4; 77.4, 78.4];
    mask1 = any(t >= zones(:,1)' & t <= zones(:,2)', 2);
    p1 = polyfit(t(mask1), F(mask1), 1);

    % Mask 2: pre-drop points at SL=2.2 (L >= 1.1, 150ms window before each drop)
    L_smooth = movmean(L, 20, 'omitnan');
    falling  = find(diff(L_smooth >= 1.1) == -1) + 1;
    t_drops  = t(falling);
    mask2 = false(size(t));
    for jj = 1:numel(t_drops)
        mask2 = mask2 | (t >= t_drops(jj) - 0.15 & t < t_drops(jj));
    end
    p2 = polyfit(t(mask2), F(mask2), 1);

    t_line = linspace(71, 79, 500);

    nexttile;hold on;
    win = t >= 70 & t <= 79.5;
    plot(t(win), F(win), 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'DisplayName', 'Raw trace');
    plot(t(mask1), F(mask1), 'o', 'Color', [0.2 0.2 0.8], 'MarkerSize', 3, 'LineStyle', 'none', 'DisplayName', 'SL=2.0 pts');
    plot(t_line, polyval(p1, t_line), '--', 'Color', [0.2 0.2 0.8], 'LineWidth', 1.5, 'DisplayName', sprintf('SL=2.0: %.3f kPa/s', p1(1)));
    plot(t(mask2), F(mask2), 's', 'Color', [0.8 0.2 0.2], 'MarkerSize', 3, 'LineStyle', 'none', 'DisplayName', 'SL=2.2 pts');
    plot(t_line, polyval(p2, t_line), ':', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5, 'DisplayName', sprintf('SL=2.2: %.3f kPa/s', p2(1)));
    xlim([70 79.5]); ylabel('Force (kPa)'); xlabel('Time (s)');
    title(sprintf('%s — Rundown', ds_labels{ii_loop}));
    legend('Location', 'northeast', 'FontSize', 7); box on;
    yl = ylim;
    text(70.2, yl(2) - 0.05*(yl(2)-yl(1)), sprintf('slope_{2.0} = %.3f kPa/s', p1(1)), ...
        'Color', [0.2 0.2 0.8], 'FontSize', 8, 'VerticalAlignment', 'top');
    text(70.2, yl(2) - 0.15*(yl(2)-yl(1)), sprintf('slope_{2.2} = %.3f kPa/s', p2(1)), ...
        'Color', [0.8 0.2 0.2], 'FontSize', 8, 'VerticalAlignment', 'top');
end

sgtitle('Rundown fits: 2mM Active vs 8mM Active');
exportgraphics(fB, fullfile(dataDir, 'Report_Rundown_23.png'), 'Resolution', 300);
fprintf('Saved Report_Rundown_23.png\n');

% -------------------------------------------------------------------------
function F_corrected = scaleRundown(F_rdact, L_ref, r0, L0, k, f0)
% scaleRundown  Scale run-down active force with SL-dependent compensation.
%
%   F_corrected = scaleRundown(F_rdact, L_ref, r0, L0, k)
%
%   Inputs:
%     F_rdact - run-down active force vector (same length as L_ref)
%     L_ref   - sarcomere length vector (same time axis as F_rdact)
%     r0      - base scale factor at pivot length L0 (e.g. 75/66)
%     L0      - pivot sarcomere length (e.g. mean(L_ref))
%     k       - SL slope [1/length_unit]: positive = more compensation at higher SL
%
%   The scale factor varies linearly with SL:
%     scale(L) = r0 + k*(L - L0)
%
%   Start with k=0 to recover flat-ratio behaviour, then increase k until
%   the corrected trace aligns with the reference active force.

    scale = r0 + k .* (L_ref - L0);
    F_corrected = (F_rdact) .* scale;
end
