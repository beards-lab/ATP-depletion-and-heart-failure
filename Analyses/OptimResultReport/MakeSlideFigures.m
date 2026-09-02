% MakeSlideFigures.m — the four presentation figures, sized and styled for a
% PowerPoint slide rather than for a report page.
%
% LAYOUT. Four figures in a 2 x 2 grid, each 3:2 landscape, so the block is
% itself 3:2 and drops into the INNER content container of a 16:9 slide (below
% the title, inside the margins) rather than bleeding to the edges. That makes
% each figure ~4.3 x 2.9 in on the slide, which is small -- hence everything
% below is deliberately oversized: 14 pt axes, 15 pt labels, heavy lines.
% MATLAB font sizes are POINTS relative to the figure's paper size, so setting
% the figure in inches makes the on-slide point size exactly what is asked for.
%
% ONE MODEL, ONE PREPARATION. Only the current fit is drawn, labelled simply
% "model" -- the previous parametrisation belongs in the report's comparison
% panels, not on a slide. Data is 03/27 8 mM, the trace the model is fitted to.
%
% LINE WEIGHT IS PER-FIGURE, ON PURPOSE. On the two figures whose data is a
% long noisy trace (whole slack protocol, staircase) a heavy model line covers
% the data completely and the comparison disappears; there the model is drawn
% THIN over thick data. On the two sparse figures (F-V points, one re-stretch
% transient) the model is the heavy line.
%
% Reads results/fit_report_<tag>.mat -- run RunFitReport.m first. No model is
% evaluated here.
%
% Figures -> results/slide_*.png   (300 dpi)
%   slide_fv            force-velocity, normalised POWER (the scored quantity)
%   slide_restretch_all whole slack protocol, all five release/re-stretch cycles
%   slide_restretch_det the PRE-LAST cycle's re-stretch transient (7 ML/s)
%   slide_stairs        staircase prediction (never scored during the fit)

clearvars -except TAG_NEW TAG_BASE; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results');

if ~exist('TAG_NEW','var');  TAG_NEW  = 'spreadw_opt';  end

L = load(fullfile(resDir, sprintf('fit_report_%s.mat', TAG_NEW)));
R = L.R(2);                       % the current fit only

%% ---- slide style --------------------------------------------------------
FIG_W = 4.6; FIG_H = FIG_W*2/3;      % inches, 3:2 -- one cell of the 2x2 block
DPI   = 300;
FS    = 14; FSL = 15; FST = 15; FSLG = 13;
cDat  = [0.10 0.12 0.16];            % data
cMod  = [0.78 0.10 0.14];            % model
LW_SPARSE = 2.8;                     % model line where data is sparse
LW_DENSE  = 1.2;                     % model line over a dense noisy trace
newFig = @() figure('Units','inches','Position',[1 1 FIG_W FIG_H],'Color','w');

%% ---- 1. force-velocity, normalised power -------------------------------
fd = R.fdata;
vp = fd.FV_powerv(:)'; mu = fd.FV_normpowerAvg(:)'; sd = sqrt(fd.FV_normpowerVar(:)');
% v = 0 drawn as the anchor the curve leaves from; it is NOT scored (power is
% identically 0 there by construction in every dataset).
vp0 = [0 vp]; mu0 = [0 mu]; sd0 = [0 sd];

fig = newFig(); ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
fill([vp0 fliplr(vp0)], [mu0+sd0 fliplr(mu0-sd0)], cDat, 'FaceAlpha',0.16, ...
     'EdgeColor','none', 'DisplayName',['data ' char(177) ' SD']);
plot(vp0, mu0, '-', 'Color',cDat, 'LineWidth',2.6, 'DisplayName','data');
plot(vp,  mu,  'o', 'Color',cDat, 'MarkerFaceColor',cDat, 'MarkerSize',6, 'HandleVisibility','off');
plot(vp0, [0 R.fmodel.FV_normpowerAvg(:)'], '-', 'Color',cMod, ...
     'LineWidth',LW_SPARSE, 'DisplayName','model');
plot(vp,  R.fmodel.FV_normpowerAvg(:)', 's', 'Color',cMod, 'MarkerFaceColor',cMod, ...
     'MarkerSize',6, 'HandleVisibility','off');
plot(0, 0, 'o', 'Color',cDat, 'MarkerFaceColor','w', 'LineWidth',1.6, 'MarkerSize',7, ...
     'HandleVisibility','off');
xlim([-0.15 max(vp)+0.25]); ylim([-0.03 0.72]);
xlabel('shortening velocity (ML/s)'); ylabel('normalised power');
title('Force–velocity (power)');
legend('Location','southeast');
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(resDir,'slide_fv.png'), 'Resolution',DPI);

%% ---- 2. whole slack protocol -------------------------------------------
% Dense noisy data: model drawn THIN on top so the trace stays visible.
dtb = R.out_slack.datatable;
Sv  = load(['data/' R.p.velocitytableonfile]);
rel = Sv.velocitytable(Sv.velocitytable(:,2) < -1, 1);  nC = numel(rel);
XL  = [rel(1)-0.20, rel(end)+0.60];

fig = newFig(); ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
m = dtb(:,1) >= XL(1) & dtb(:,1) <= XL(2);
plot(dtb(m,1)-rel(1), dtb(m,3), '-', 'Color',[cDat 0.55], 'LineWidth',2.4, ...
     'DisplayName','data 03/27');
o = R.out_slack; mm = o.t >= XL(1) & o.t <= XL(2);
plot(o.t(mm)-rel(1), o.Force(mm), '-', 'Color',cMod, 'LineWidth',LW_DENSE, 'DisplayName','model');
xlim(XL - rel(1)); ylim([-16 112]);
xlabel('time from first release (s)'); ylabel('force (kPa)');
title('Slack protocol — all five cycles');
legend('Location','southwest');
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(resDir,'slide_restretch_all.png'), 'Resolution',DPI);

%% ---- 3. re-stretch detail, the PRE-LAST cycle --------------------------
cD  = nC - 1;                         % pre-last = cycle 4 of 5 (7 ML/s)
WIN = [0.292 0.40];
fig = newFig(); ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
m = dtb(:,1)-rel(cD) >= WIN(1) & dtb(:,1)-rel(cD) <= WIN(2);
plot(1000*(dtb(m,1)-rel(cD)-WIN(1)), dtb(m,3), '-', 'Color',cDat, 'LineWidth',1.6, ...
     'DisplayName','data 03/27');
o = R.out_slack; mm = o.t-rel(cD) >= WIN(1) & o.t-rel(cD) <= WIN(2);
plot(1000*(o.t(mm)-rel(cD)-WIN(1)), o.Force(mm), '-', 'Color',cMod, ...
     'LineWidth',LW_SPARSE, 'DisplayName','model');
xlim([0 1000*(WIN(2)-WIN(1))]);
xlabel('time from re-stretch (ms)'); ylabel('force (kPa)');
title(sprintf('Re-stretch transient — cycle %d (%.0f ML/s)', cD, fd.v_restretch(cD)));
legend('Location','southeast');
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(resDir,'slide_restretch_det.png'), 'Resolution',DPI);

%% ---- 4. staircase prediction -------------------------------------------
% Dense noisy data again: thin model line.
Sd  = load(['data/' R.p.stairs_velocitytableonfile]);
dt2 = Sd.datatable; vt2 = Sd.velocitytable;
t0s = vt2(find(vt2(:,2)~=0,1),1);

fig = newFig(); ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(dt2(:,1)-t0s, dt2(:,3), '-', 'Color',[cDat 0.55], 'LineWidth',2.0, ...
     'DisplayName','data 03/27');
o = R.out_stairs;
F0 = interp1(o.t, o.Force, t0s, 'linear','extrap'); if F0 == 0; F0 = 1; end
plot(o.t-t0s, o.Force/F0, '-', 'Color',cMod, 'LineWidth',LW_DENSE, 'DisplayName','model');
xlim([-0.10 dt2(end,1)-t0s]); ylim([0.90 1.58]);
xlabel('time from first step (s)'); ylabel('F / F_{iso}');
title('Staircase — PREDICTION (not fitted)');
legend('Location','southeast');
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(resDir,'slide_stairs.png'), 'Resolution',DPI);

fprintf('wrote 4 slide figures (%.2f x %.2f in @ %d dpi) to %s\n', FIG_W, FIG_H, DPI, resDir);
fprintf('block of 2x2 = %.2f x %.2f in — fits the inner container of a 16:9 slide\n', 2*FIG_W, 2*FIG_H);

%% ======================= helpers ========================================
function iStyleAx(ax, FS, FSL, FST, FSLG)
% Slide styling: oversized type, heavier axes, tight margins. Applied last so
% it overrides anything the plotting calls set.
    set(ax, 'FontSize',FS, 'LineWidth',1.1, 'GridAlpha',0.16, 'Box','on', ...
            'TickDir','out', 'TickLength',[0.015 0.015]);
    ax.XLabel.FontSize = FSL;  ax.XLabel.FontWeight = 'bold';
    ax.YLabel.FontSize = FSL;  ax.YLabel.FontWeight = 'bold';
    ax.Title.FontSize  = FST;  ax.Title.FontWeight  = 'bold';
    lg = get(ax,'Legend');
    if ~isempty(lg); lg.FontSize = FSLG; lg.Box = 'on'; end
    % squeeze the axes into the figure: on a slide every wasted margin is
    % font size the reader does not get.
    ax.LooseInset = max(ax.TightInset, [0.02 0.02 0.01 0.01]);
end
