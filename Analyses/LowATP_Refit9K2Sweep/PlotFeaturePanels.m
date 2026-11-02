% PlotFeaturePanels.m -- every feature, 8 mM vs 2 mM, data vs the model variants.
%   Figure 1: slack features over the five slack-restretch cycles (pooled data
%             band = mean +- the SD the objective uses).
%   Figure 2: ktr-protocol features, force-velocity features and the model
%             read-outs (XTOR, attached / SRX fractions) with the objective's bounds.
% Model runs are the saved ones (no simulation): 8 mM incumbent; 2 mM sets
%   A  R2 x0.54 | B  + ksr0 x0.4 | C  + kmsrd x0.4 (gates at 2 mM only).
rootF = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(rootF); addpath(genpath(rootF));
rs = fullfile(rootF,'Analyses','LowATP_Refit9K2Sweep','results');
Lr = load(fullfile(rs,'levers_ref.mat'));  fm8 = Lr.fm8L;
K2 = load(fullfile(rs,'k2sweep.mat'), 'RK2'); fmA = K2.RK2(abs([K2.RK2.g]-0.54) < 1e-9).fm;
G  = load(fullfile(rs,'ksr0grid.mat'), 'GR'); fmB = G.GR(abs([G.GR.g]-0.54) < 1e-9 & abs([G.GR.f]-0.4) < 1e-9).fm;
E2 = load(fullfile(rs,'extB2.mat'), 'RL');   fmC = E2.RL(abs([E2.RL.f]-0.4) < 1e-9).fm;
BF = {load(fullfile(rs,'bestfit_k2only.mat')), load(fullfile(rs,'bestfit_k2ksr0.mat')), load(fullfile(rs,'bestfit_k2ksr0kmsrd.mat'))};
PT8 = pooledTarget8mM(); PT2 = pooledTarget2mM(); KT8 = pooledKtrTarget(8); KT2 = pooledKtrTarget(2); FV2 = fv2Target();
LoadData; D = Data_ATP;

navy = [0 39 76]/255; red = [0.75 0.2 0.12];
MOD  = {fm8, fmA, fmB, fmC};
MCOL = {navy, [0.85 0.33 0.10], [0.93 0.69 0.13], [0.49 0.18 0.56]};
MLAB = {'8 mM model', '2 mM A: R_2\times0.54', '2 mM B: + k_{sr0}\times0.4', '2 mM C: + k_{msrd}\times0.4'};
MMK  = {'o', 's', 'd', '^'};

%% ---- figure 1: slack features over cycles ---------------------------------------
SL = { 'steady','Steady force (kPa)',1; 'A','Redeveloped force (kPa)',1; 'peak1_y','Re-stretch peak (kPa)',1
       'vall_y','Valley after peak (kPa)',1; 'peak2','Second peak (kPa)',1; 'vall2_dy','Post-peak drop (kPa)',1
       'ktrfit','Slack k_{tr} (s^{-1})',1; 'rsK80','Post-restretch rate rsK80 (s^{-1})',1; 't0','Redevelopment delay t_0 (ms)',1000
       'rsTfitE','Recovery time rsTfitE (ms)',1000; 'vall_tq','Valley time (ms)',1000; 'peak1_tq','Peak time (ms)',1000 };
f1 = figure('Color','w','Position',[20 20 2000 1150],'Visible','off');
tl = tiledlayout(f1, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:size(SL,1)
    nexttile; hold on; box on; grid on; nm = SL{i,1}; sc = SL{i,3}; x = 1:5;
    for c = 1:2
        if c == 1; P = PT8; col = navy; else; P = PT2; col = red; end
        j = find(strcmp(P.name, nm), 1); if isempty(j); continue; end
        v = P.val(j,:)*sc; sd = abs(v)*P.s(j);
        fill([x fliplr(x)], [v+sd fliplr(v-sd)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(x, v, '-', 'Color', col, 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end
    for k = 1:4
        m = iV(MOD{k}, nm)*sc; if isempty(m); continue; end
        plot(x(1:numel(m)), m, ['--' MMK{k}], 'Color', MCOL{k}, 'MarkerFaceColor', MCOL{k}, 'LineWidth', 1.6, 'MarkerSize', 6, 'HandleVisibility', 'off');
    end
    set(gca, 'XTick', 1:5, 'XLim', [0.7 5.3], 'FontSize', 12); title(SL{i,2}, 'FontSize', 13);
    if i > 8; xlabel('slack–restretch cycle'); end
end
iLegend(tl, navy, red, MCOL, MLAB, MMK);
title(tl, 'Slack–restretch features over the five cycles: pooled data (band = \pm objective SD, n = 3) vs model', 'FontSize', 16);
exportgraphics(f1, fullfile(rs,'features_slack.png'), 'Resolution', 170); close(f1);

%% ---- figure 2: ktr protocol, F-V, read-outs -----------------------------------------
f2 = figure('Color','w','Position',[20 20 2000 1150],'Visible','off');
tl = tiledlayout(f2, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
% ktr protocol features: data (8, 2 mM) with SD, models as markers
KP = {'ktrP_k80','k_{tr} protocol k80 (s^{-1})'; 'ktrP_k63','k_{tr} protocol k63 (s^{-1})'; 'ktrP_F20','k_{tr} protocol force at 20 ms (F/F_{iso})'};
for i = 1:3
    nexttile; j = find(strcmp(KT8.name, KP{i,1}), 1);
    iCond(KT8.val(j), KT8.s(j)*KT8.val(j), KT2.val(j), KT2.s(j)*KT2.val(j), cellfun(@(f) iS(f, KP{i,1}), MOD), MCOL, MMK, navy, red);
    title(KP{i,2}, 'FontSize', 13);
end
% FV2 in-prep point
nexttile; iCond(FV2.m8, FV2.s8, FV2.m2, FV2.s2, [BF{1}.fvp8, cellfun(@(b) b.fvp2, BF)], MCOL, MMK, navy, red);
title('FV2 in-prep P/P_0 at 2 ML/s (n = 2)', 'FontSize', 13);
% F-V curve and power
nexttile([1 2]); hold on; box on; grid on;
plot(D(:,1), D(:,2)/D(1,2), 'o', 'Color', navy, 'MarkerFaceColor', navy, 'MarkerSize', 7, 'HandleVisibility', 'off');
plot(D(:,1), D(:,4)/D(1,4), 'o', 'Color', red, 'MarkerFaceColor', red, 'MarkerSize', 7, 'HandleVisibility', 'off');
fvs = {BF{1}.fv8, BF{1}.fv2, BF{2}.fv2, BF{3}.fv2};
for k = 1:4; plot(fvs{k}.FV_v, fvs{k}.FV_fnorm, ['--' MMK{k}], 'Color', MCOL{k}, 'MarkerFaceColor', MCOL{k}, 'LineWidth', 1.6, 'HandleVisibility', 'off'); end
xlim([0 6]); ylim([0 1.05]); set(gca, 'FontSize', 12); xlabel('shortening velocity (ML/s)'); ylabel('F / F_0');
title('Force–velocity (Baker: 8 mM scored, 2 mM not scored)', 'FontSize', 13);
nexttile([1 2]); hold on; box on; grid on;
plot(D(:,1), D(:,1).*D(:,2)/D(1,2), 'o', 'Color', navy, 'MarkerFaceColor', navy, 'MarkerSize', 7, 'HandleVisibility', 'off');
plot(D(:,1), D(:,1).*D(:,4)/D(1,4), 'o', 'Color', red, 'MarkerFaceColor', red, 'MarkerSize', 7, 'HandleVisibility', 'off');
for k = 1:4; plot(fvs{k}.FV_v, fvs{k}.FV_v(:).*fvs{k}.FV_fnorm(:), ['--' MMK{k}], 'Color', MCOL{k}, 'MarkerFaceColor', MCOL{k}, 'LineWidth', 1.6, 'HandleVisibility', 'off'); end
xlim([0 6]); set(gca, 'FontSize', 12); xlabel('shortening velocity (ML/s)'); ylabel('F/F_0 \cdot v (ML/s)');
title('Power', 'FontSize', 13);
% read-outs with the objective's bounds
RO = {'XTOR','Isometric ATP turnover XTOR (s^{-1})',[3 15]; 'XTOR_vmax','Turnover at v_{max} (s^{-1})',[5 25]
      'attached_ss','Attached fraction (isometric)',[0.10 0.45]; 'SRX_ss','Super-relaxed fraction',[0.01 0.40]};
for i = 1:4
    nexttile; hold on; box on; grid on; b = RO{i,3};
    patch([0.5 4.5 4.5 0.5], [b(1) b(1) b(2) b(2)], [0.5 0.5 0.5], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    if strcmp(RO{i,1}, 'attached_ss')
        patch([0.5 4.5 4.5 0.5], [0.12 0.12 0.15 0.15], [0.2 0.6 0.2], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
        text(0.6, 0.135, 'literature', 'FontSize', 10, 'Color', [0.1 0.4 0.1]);
    end
    for k = 1:4
        v = iV(MOD{k}, RO{i,1}); if isempty(v); continue; end
        bar(k, mean(v), 0.6, 'FaceColor', MCOL{k}, 'EdgeColor', 'none');
        plot(k*ones(size(v)), v, '.', 'Color', [0.2 0.2 0.2], 'MarkerSize', 10);
    end
    set(gca, 'XTick', 1:4, 'XTickLabel', {'8 mM','A','B','C'}, 'FontSize', 12, 'XLim', [0.5 4.5]);
    yl = ylim; ylim([0 max(yl(2), b(2))*1.1]);
    title(RO{i,2}, 'FontSize', 13);
end
iLegend(tl, navy, red, MCOL, MLAB, MMK);
title(tl, 'k_{tr} protocol, force–velocity and model read-outs (grey = objective bound; no data for read-outs)', 'FontSize', 16);
exportgraphics(f2, fullfile(rs,'features_ktr_fv.png'), 'Resolution', 170); close(f2);
fprintf('wrote features_slack.png, features_ktr_fv.png\n');

function v = iV(st, nm)
    v = []; if isfield(st, nm); v = double(st.(nm)); v = v(:)'; end
end
function v = iS(st, nm)
    v = NaN; if isfield(st, nm); x = double(st.(nm)); x = x(isfinite(x)); if ~isempty(x); v = mean(x); end; end
end
function iCond(d8, s8, d2, s2, m, MCOL, MMK, navy, red)
    % x = 1: 8 mM (data + 8 mM model), x = 2: 2 mM (data + sets A, B, C)
    hold on; box on; grid on;
    errorbar(1, d8, s8, 'o', 'Color', navy, 'MarkerFaceColor', navy, 'MarkerSize', 9, 'LineWidth', 2, 'CapSize', 12);
    errorbar(2, d2, s2, 'o', 'Color', red, 'MarkerFaceColor', red, 'MarkerSize', 9, 'LineWidth', 2, 'CapSize', 12);
    plot(1.25, m(1), MMK{1}, 'Color', MCOL{1}, 'MarkerFaceColor', MCOL{1}, 'MarkerSize', 9);
    for k = 2:4; plot(2 + 0.15*(k-1), m(k), MMK{k}, 'Color', MCOL{k}, 'MarkerFaceColor', MCOL{k}, 'MarkerSize', 9); end
    set(gca, 'XTick', [1 2], 'XTickLabel', {'8 mM','2 mM'}, 'XLim', [0.6 2.8], 'FontSize', 12);
    yl = ylim; ylim(yl + [-1 1]*0.08*diff(yl));      % keep edge markers visible
end
function iLegend(tl, navy, red, MCOL, MLAB, MMK)
    ax = nexttile(tl, 1); hold(ax, 'on'); h = gobjects(0);
    h(end+1) = patch(ax, nan, nan, navy, 'FaceAlpha', 0.3, 'EdgeColor', navy, 'LineWidth', 2);
    h(end+1) = patch(ax, nan, nan, red, 'FaceAlpha', 0.3, 'EdgeColor', red, 'LineWidth', 2);
    for k = 1:4; h(end+1) = plot(ax, nan, nan, ['--' MMK{k}], 'Color', MCOL{k}, 'MarkerFaceColor', MCOL{k}, 'LineWidth', 1.6); end %#ok<AGROW>
    lg = legend(h, [{'8 mM pooled data \pm SD', '2 mM pooled data \pm SD'}, MLAB], 'Orientation', 'horizontal', 'FontSize', 13);
    lg.Layout.Tile = 'south';
end
