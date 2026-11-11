% MakeSlideFigures.m — the four presentation figures, sized and styled for a
% PowerPoint slide rather than for a report page.
%
% LAYOUT. Three figures are 3:2 landscape (one grid cell); the whole slack
% protocol is TWO CELLS WIDE, because five release/re-stretch cycles squeezed
% into a 3:2 box is a picture of nothing. So the block is 2 columns x 3 rows
% with the wide figure spanning a row -- see MakeSlide.m for the placement.
% Each figure is small on the slide, hence everything below is deliberately
% oversized: 14 pt axes, 15 pt labels, heavy lines. MATLAB font sizes are POINTS
% relative to the figure's paper size, so setting the figure in inches makes the
% on-slide point size exactly what is asked for.
%
% LENGTH ON A SECOND AXIS. Slack, slack-detail and staircase carry the commanded
% length (ML) on the right axis. It is scaled to sit in the BOTTOM THIRD of the
% panel, under the force trace, in a muted blue-grey: the length is context for
% reading the force, not a co-equal signal, and it must never be mistaken for
% one. iLengthAxis() does the scaling so all three panels place it identically.
%
% ONE MODEL, ONE PREPARATION. Only the current fit is drawn, labelled simply
% "model". Data is 03/27 8 mM, the trace the model is fitted to.
%
% LINE WEIGHT IS PER-FIGURE, ON PURPOSE. On the figures whose data is a long
% noisy trace (whole slack protocol, staircase) a heavy model line covers the
% data completely and the comparison disappears; there the model is drawn THIN
% over thick data. On the sparse figures the model is the heavy line.
%
% Reads results/fit_report_<tag>.mat -- run RunFitReport.m first. The F-V panel
% comes from RunFVFine.m (results/fv_fine.mat) when that exists, so the fine
% velocity grid is not thrown away by re-running this. No model is evaluated.
%
% Figures -> results/slide_*.png   (300 dpi)
%   slide_fv            force-velocity, normalised POWER (the scored quantity)
%   slide_restretch_all whole slack protocol, all five cycles      [2 cells wide]
%   slide_restretch_det the PRE-LAST cycle, slack release -> 100 ms post-stretch
%   slide_stairs        staircase prediction (never scored during the fit)

clearvars -except TAG_NEW; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results');

if ~exist('TAG_NEW','var');  TAG_NEW  = 'spreadw_opt';  end

L = load(fullfile(resDir, sprintf('fit_report_%s.mat', TAG_NEW)));
R = L.R(2);                       % the current fit only

%% ---- slide style --------------------------------------------------------
FIG_W = 4.6; FIG_H = FIG_W*2/3;      % inches, 3:2 -- one cell
GAP   = 0.14;                        % the gap the slide leaves between cells
FIG_W2 = 2*FIG_W + GAP;              % two cells wide, same height
DPI   = 300;
FS    = 14; FSL = 15; FST = 15; FSLG = 13;
cDat  = [0.10 0.12 0.16];            % data
cMod  = [0.78 0.10 0.14];            % model
cLen  = [0.42 0.53 0.62];            % commanded length — muted, deliberately
LW_SPARSE = 2.8;                     % model line where data is sparse
LW_DENSE  = 1.2;                     % model line over a dense noisy trace
newFig  = @(w) figure('Units','inches','Position',[1 1 w FIG_H],'Color','w');

%% ---- 1. force-velocity, normalised power -------------------------------
FVF = fullfile(resDir,'fv_fine.mat');
fig = newFig(FIG_W); ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
if isfile(FVF)
    % fine grid from RunFVFine.m, with its plotting exclusions applied there
    C = load(FVF);
    VPM = 6; VDR = 0.25;
    mM = C.vFine <= VPM+1e-9; for vv = VDR; mM = mM & abs(C.vFine-vv) > 1e-9; end
    mD = C.vD <= VPM+1e-9;
    fill([C.vD(mD) fliplr(C.vD(mD))], [C.pMu(mD)+C.pSd(mD) fliplr(C.pMu(mD)-C.pSd(mD))], ...
         cDat, 'FaceAlpha',0.16,'EdgeColor','none','DisplayName',['data ' char(177) ' SD']);
    plot(C.vD(mD), C.pMu(mD), 'o','Color',cDat,'MarkerFaceColor',cDat,'MarkerSize',7, ...
         'LineStyle','none','DisplayName','data');
    plot(C.vFine(mM), C.Pw(mM), '-','Color',cMod,'LineWidth',LW_SPARSE,'DisplayName','model');
    xlim([-0.2 VPM+0.2]); ylim([-0.03 max([C.Pw(mM) C.pMu(mD)+C.pSd(mD)])*1.12]);
else
    fd = R.fdata;
    vp = fd.FV_powerv(:)'; mu = fd.FV_normpowerAvg(:)'; sd = sqrt(fd.FV_normpowerVar(:)');
    fill([[0 vp] fliplr([0 vp])], [[0 mu+sd] fliplr([0 mu-sd])], cDat, 'FaceAlpha',0.16, ...
         'EdgeColor','none','DisplayName',['data ' char(177) ' SD']);
    plot([0 vp], [0 mu], 'o-','Color',cDat,'MarkerFaceColor',cDat,'LineWidth',2.6,'DisplayName','data');
    plot([0 vp], [0 R.fmodel.FV_normpowerAvg(:)'], '-','Color',cMod, ...
         'LineWidth',LW_SPARSE,'DisplayName','model');
    xlim([-0.15 max(vp)+0.25]);
end
xlabel('shortening velocity (ML/s)'); ylabel('normalised power');
title('Force–velocity (power)');
legend('Location','northeast');
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(resDir,'slide_fv.png'), 'Resolution',DPI);

%% ---- 2. whole slack protocol — TWO CELLS WIDE --------------------------
% Dense noisy data: model drawn THIN on top so the trace stays visible.
dtb = R.out_slack.datatable;                 % col2 = SL in microns, col3 = kPa
Sv  = load(['data/' R.p.velocitytableonfile]);
rel = Sv.velocitytable(Sv.velocitytable(:,2) < -1, 1);  nC = numel(rel);
XL  = [rel(1)-0.20, rel(end)+0.60];

fig = newFig(FIG_W2); ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
m = dtb(:,1) >= XL(1) & dtb(:,1) <= XL(2);
o = R.out_slack; mm = o.t >= XL(1) & o.t <= XL(2);
% ML = 1 is pinned to force = 0, so the zero gridline serves both axes. The
% floor is -24 rather than -16 to leave room for the 0.94-to-1 part of the ML
% range once the pin fixes the scale.
YLF = [-24 112];
iLengthAxis(ax, dtb(m,1)-rel(1), dtb(m,2)/2, cLen, FS, FSL, 0.34, YLF, 0);  % SL um -> ML
plot(dtb(m,1)-rel(1), dtb(m,3), '-', 'Color',[cDat 0.55], 'LineWidth',2.4, ...
     'DisplayName','data');
plot(o.t(mm)-rel(1), o.Force(mm), '-', 'Color',cMod, 'LineWidth',LW_DENSE, 'DisplayName','model');
xlim(XL - rel(1)); ylim(YLF);
xlabel('time from first release (s)'); ylabel('force (kPa)');
title('Slack protocol — all five cycles');
% Two columns, not two rows: the clear band between the ML trace (tops out
% ~27 kPa) and the force plateaus (~55 kPa at the last cycle) is only ~28 kPa
% deep, and a stacked legend needs ~33. Laid out horizontally it fits.
lg = legend('Location','southwest','NumColumns',2);
iStyleAx(ax, FS, FSL, FST, FSLG);
% Park it at the right, inside that band — the only region of this panel no
% trace crosses. iStyleAx settles the axes position first, so this comes after.
iPlaceLegend(ax, lg, 0.70, 0.345);
exportgraphics(fig, fullfile(resDir,'slide_restretch_all.png'), 'Resolution',DPI);

%% ---- 3. the PRE-LAST cycle, whole cycle --------------------------------
% From 20 ms BEFORE the slack release to 100 ms after the re-stretch, so the
% panel carries the release, the redevelopment and the re-stretch transient in
% one frame. t = 0 is the slack release.
cD  = nC - 1;                                 % pre-last = cycle 4 of 5 (7 ML/s)
WIN = [-0.020 0.400];
fig = newFig(FIG_W2); ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
m = dtb(:,1)-rel(cD) >= WIN(1) & dtb(:,1)-rel(cD) <= WIN(2);
o = R.out_slack; mm = o.t-rel(cD) >= WIN(1) & o.t-rel(cD) <= WIN(2);
% Same pin as the whole-protocol panel: ML = 1 on the force-zero gridline.
% topFrac is 0.40 here, not 0.34: this cycle only reaches ML 0.96, and with the
% pin fixed the sub-1 span scales as 0.04/span — at 0.34 the "0.96" and "1"
% labels land 0.17 in apart and collide at 13 pt. A taller band spreads them.
YLD = [-25 100];
iLengthAxis(ax, 1000*(dtb(m,1)-rel(cD)), dtb(m,2)/2, cLen, FS, FSL, 0.40, YLD, 0);
plot(1000*(dtb(m,1)-rel(cD)), dtb(m,3), '-', 'Color',cDat, 'LineWidth',1.4, ...
     'DisplayName','data');
plot(1000*(o.t(mm)-rel(cD)), o.Force(mm), '-', 'Color',cMod, ...
     'LineWidth',2.2, 'DisplayName','model');
xlim(1000*WIN); ylim(YLD);
xlabel('time from slack release (ms)'); ylabel('force (kPa)');
title(sprintf('Slack cycle %d  (re-stretch %.0f ML/s)', cD, R.fdata.v_restretch(cD)));
legend('Location','north');   % the plateau's mid-span is the only clear region
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(resDir,'slide_restretch_det.png'), 'Resolution',DPI);

%% ---- 4. staircase prediction -------------------------------------------
% Dense noisy data again: thin model line.
Sd  = load(['data/' R.p.stairs_velocitytableonfile]);
dt2 = Sd.datatable; vt2 = Sd.velocitytable;            % col2 = SL in microns
t0s = vt2(find(vt2(:,2)~=0,1),1);

fig = newFig(FIG_W2); ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
o = R.out_stairs;
F0 = interp1(o.t, o.Force, t0s, 'linear','extrap'); if F0 == 0; F0 = 1; end
% No pin here: F/F_iso has no zero to align to. The force RISES across this
% panel, so the ML axis is extended further (band 0.14 rather than 0.34) —
% that flattens the length trace against the floor and keeps clear air between
% the two all the way to the right-hand end.
YLS = [0.80 1.60];
iLengthAxis(ax, dt2(:,1)-t0s, dt2(:,2)/2, cLen, FS, FSL, 0.14, YLS);
plot(dt2(:,1)-t0s, dt2(:,3), '-', 'Color',[cDat 0.55], 'LineWidth',2.0, ...
     'DisplayName','data');
plot(o.t-t0s, o.Force/F0, '-', 'Color',cMod, 'LineWidth',LW_DENSE, 'DisplayName','model');
xlim([-0.10 dt2(end,1)-t0s]); ylim(YLS);
xlabel('time from first step (s)'); ylabel('F / F_{iso}');
title('Staircase — PREDICTION (not fitted)');
legend('Location','northeast');
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(resDir,'slide_stairs.png'), 'Resolution',DPI);

fprintf('wrote 4 slide figures @ %d dpi to %s\n', DPI, resDir);
fprintf('  three cells  %.2f x %.2f in\n', FIG_W,  FIG_H);
fprintf('  wide  (2x1)  %.2f x %.2f in  (slack protocol)\n', FIG_W2, FIG_H);

%% ======================= helpers ========================================
function iLengthAxis(ax, t, ml, c, FS, FSL, topFrac, leftYlim, pinLeft)
% Commanded length on the right axis, kept in the lower band of the panel so
% force stays the figure's subject and the length reads as context.
%
%   topFrac   where the HIGHEST ML sits, as a fraction of the panel height
%   leftYlim  the force axis limits (set here, so the mapping below is exact)
%   pinLeft   optional: the LEFT-axis value that ML = 1 must line up with.
%             Passing 0 makes the force-zero gridline double as the ML = 1
%             line, which is what ties the two axes together visually.
%
% With a pin the right limits are fully determined — two points fix a linear
% map — so `topFrac` and `pinLeft` together decide everything. The left axis
% therefore needs enough room BELOW pinLeft to hold the sub-1 part of the ML
% range at the resulting scale; the assert below says so if it does not.
    if nargin < 7 || isempty(topFrac); topFrac = 1/3; end
    yyaxis(ax,'right');
    lo = min(ml); hi = max(ml); r = max(hi-lo, eps);
    plot(ax, t, ml, '-', 'Color', c, 'LineWidth', 1.3, 'HandleVisibility','off');

    if nargin >= 9 && ~isempty(pinLeft)
        fPin = (pinLeft - leftYlim(1)) / (leftYlim(2) - leftYlim(1));   % height of ML=1
        if topFrac <= fPin
            error('iLengthAxis:badPin','topFrac (%.3f) must exceed the pin height (%.3f)', ...
                  topFrac, fPin);
        end
        span = (hi - 1) / (topFrac - fPin);
        r0   = 1 - fPin*span;
        if lo < r0
            warning('iLengthAxis:clipped', ...
                ['ML %.3f falls below the axis floor — give the force axis more room ' ...
                 'under %.4g, or lower topFrac'], lo, pinLeft);
        end
    else
        pad  = 0.12*r;
        span = (r + pad) / topFrac;
        r0   = lo - pad;
    end
    ylim(ax, [r0, r0 + span]);
    ylabel(ax, 'ML (1)');
    ax.YAxis(2).Color = c;
    ax.YAxis(2).FontSize = FS - 1;
    ax.YAxis(2).Label.FontSize = FSL - 1;
    ax.YAxis(2).Label.FontWeight = 'bold';
    % the extremes of the actual trace, plus 1 as the shared reference
    ax.YAxis(2).TickValues = unique(round([lo, 1, hi], 2));

    yyaxis(ax,'left');
    ylim(ax, leftYlim);
end

function iPlaceLegend(ax, lg, fx, fy)
% Put the legend's lower-left corner at (fx, fy) as a fraction of the AXES box.
% legend Position is normalised to the FIGURE, so it has to be mapped through
% the axes rectangle; MATLAB has no axes-relative legend placement.
    lg.Units = 'normalized';
    ap = ax.Position;  lp = lg.Position;
    lg.Position = [ap(1) + fx*ap(3), ap(2) + fy*ap(4), lp(3), lp(4)];
end

function iStyleAx(ax, FS, FSL, FST, FSLG)
% Slide styling: oversized type, heavier axes, tight margins. Applied last so
% it overrides anything the plotting calls set.
    if numel(ax.YAxis) > 1; yyaxis(ax,'left'); end
    set(ax, 'FontSize',FS, 'LineWidth',1.1, 'GridAlpha',0.16, 'Box','on', ...
            'TickDir','out', 'TickLength',[0.015 0.015]);
    ax.YAxis(1).Color = [0.15 0.15 0.15];
    ax.XLabel.FontSize = FSL;  ax.XLabel.FontWeight = 'bold';
    ax.YAxis(1).Label.FontSize = FSL;  ax.YAxis(1).Label.FontWeight = 'bold';
    ax.Title.FontSize  = FST;  ax.Title.FontWeight  = 'bold';
    lg = get(ax,'Legend');
    if ~isempty(lg); lg.FontSize = FSLG; lg.Box = 'on'; end
    ax.LooseInset = max(ax.TightInset, [0.02 0.02 0.01 0.01]);
end
