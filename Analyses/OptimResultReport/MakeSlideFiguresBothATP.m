% MakeSlideFiguresBothATP.m -- the three time-series slide panels with BOTH
% [ATP] on the same axes: 8 mM (the spreadw fit) and 2 mM (one ladder rung, RUNG).
%
% Same geometry, styling and ML-axis treatment as MakeSlideFigures.m -- read
% that file for why the panels are two cells wide, why the model line is thin
% over dense data, and how iLengthAxis pins ML = 1 to the force-zero gridline.
% What is new here:
%
% COLOURS, matched to Analyses/AtpEffectHypotheses figs_fvfine/fig_fv_fine_bothatp:
%   model 8 mM  blue      model 2 mM  red
%   data  8 mM  black     data  2 mM  dark grey
%
% ONE LENGTH TRACE, NOT TWO. The 8 mM and 2 mM protocols on 03/27 are the same
% commanded protocol -- releases at 0 / 0.600 / 1.200 / 1.800 / 2.400 s, release
% velocities -36/-44/-53/-62/-70 ML/s, re-stretch ladder 3.9/4.9/5.9/6.8/2.9
% ML/s, SL 1.879-2.201 um -- checked below, not assumed. Drawing the second one
% would add a line that lies exactly on the first.
%
% ABSOLUTE FORCE ON THE SLACK PANELS, on purpose. The measured 2 mM elevation
% (steady 91.3 vs 77.3 kPa on 03/27) is the ATP effect these figures exist to
% show; normalising each condition by its own steady would delete it. Both
% models are fitted in kPa, so the comparison is like for like.
%
% THE STAIRCASE HAS NO 2 mM DATA. All three protocol days recorded the staircase
% at 8 mM only (DataCuration/CreateStairsProtocolVelocityTable.m drives every
% one from a *_Merged_8mM_Active.txt file). That panel is therefore F/F_iso --
% each trace over its OWN isometric -- so the 2 mM curve makes the only claim it
% can support: whether low ATP changes the staircase SHAPE. It is a prediction
% twice over: never fitted, and at an [ATP] never measured.
%
% Inputs:  results/fit_report_<tag>.mat      (8 mM model; RunFitReport.m)
%          results/bothatp/<RUNG>/lowatp_runs.mat  (2 mM; RunBothAtpSlideRuns.m
%                                             in the ATP-wt-lowatp worktree)
% Figures -> results/bothatp/<RUNG>/slide2_*.png  (300 dpi)

% RUNG selects which rung of the ATP hypothesis ladder supplies the 2 mM model.
% The 8 mM model is the ladder's frozen core and is the SAME for all four --
% every gate is exactly 1 at MgATP_ref = 8 -- so only the red curve moves
% between rungs. Each rung reads and writes its own subfolder.
clearvars -except TAG_NEW RUNG RUNGS iR; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results');

if ~exist('TAG_NEW','var');  TAG_NEW  = 'spreadw_opt';  end
if ~exist('RUNG','var');     RUNG     = 'H4';           end

outDir = fullfile(resDir,'bothatp',RUNG); if ~isfolder(outDir); mkdir(outDir); end

L = load(fullfile(resDir, sprintf('fit_report_%s.mat', TAG_NEW)));
R = L.R(2);                                   % 8 mM model
LOW = fullfile(outDir,'lowatp_runs.mat');
if ~isfile(LOW)
    error('MakeSlideFiguresBothATP:noLowAtp', ...
        ['missing %s -- run Analyses/AtpEffectHypotheses/RunBothAtpSlideRuns.m ' ...
         'in the ATP-wt-lowatp worktree first (the ATP gate code lives only there)'], LOW);
end
S = load(LOW);
gateTxt = strjoin(strrep(fieldnames(S.gates),'pAtp_',''), ', ');
if isempty(gateTxt); gateTxt = '(none)'; end
fprintf('%s | 8 mM: %s | 2 mM: %s (core %s, gates %s) built %s\n', RUNG, TAG_NEW, ...
        S.snap, S.core, gateTxt, S.built);
% NO rung subtitle on the panels: which rung a figure belongs to is carried by
% the folder it sits in. It is still printed to stdout above, which is where a
% provenance check belongs.
%
% NO LEGEND on these panels either. The colour key -- black/grey data, blue/red
% model, 8/2 mM -- is identical across every figure in this folder, so it is
% stated once (on the steady-force feature panel) and the rest of the set reads
% against it. Repeating a four-entry legend on ten panels costs a fifth of each
% panel's height to say the same thing ten times.

%% ---- slide style --------------------------------------------------------
FIG_W = 4.6; FIG_H = FIG_W*2/3;      % inches, 3:2 -- one cell
GAP   = 0.14;
FIG_W2 = 2*FIG_W + GAP;              % two cells wide, same height
DPI   = 300;
FS    = 14; FSL = 15; FST = 15; FSLG = 12;
C8 = [0.00 0.35 0.70];               % model 8 mM
C2 = [0.80 0.20 0.10];               % model 2 mM
K8 = [0.04 0.05 0.07];               % data  8 mM -- black
K2 = [0.44 0.45 0.48];               % data  2 mM -- dark grey
cLen = [0.42 0.53 0.62];             % commanded length -- muted, deliberately
LW_DENSE = 1.2;                      % model over a dense noisy trace
newFig = @(w) figure('Units','inches','Position',[1 1 w FIG_H],'Color','w');

%% ---- the two 03/27 slack traces ----------------------------------------
d8 = R.out_slack.datatable;                     % col2 = SL um, col3 = kPa
D2 = load(['data/' S.velocitytableonfile]); d2 = D2.datatable;
Sv = load(['data/' R.p.velocitytableonfile]);
rel  = Sv.velocitytable(Sv.velocitytable(:,2) < -1, 1);  nC = numel(rel);
rel2 = D2.velocitytable(D2.velocitytable(:,2) < -1, 1);
% The two protocols must BE the same protocol for a shared length trace and a
% shared time origin -- checked here rather than trusted.
dRel = (rel2 - rel2(1)) - (rel - rel(1));
if max(abs(dRel)) > 2e-3
    error('MakeSlideFiguresBothATP:protocolMismatch', ...
          'release times differ by up to %.4f s between 8 and 2 mM', max(abs(dRel)));
end

%% ---- 1. whole slack protocol -- TWO CELLS WIDE -------------------------
% Pre-window shortened to -0.10 s: what the reader needs is a glimpse of the
% ML=1.1 baseline sitting steady into the first slack, and 0.1 s is enough.
XL = [rel(1)-0.10, rel(end)+0.60];
LAB_FS = 11; LAB_COL = [0.24 0.30 0.38];   % darker than cLen, for legibility
fig = newFig(FIG_W2); ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
m8 = d8(:,1) >= XL(1) & d8(:,1) <= XL(2);
m2 = d2(:,1) >= XL(1)-rel(1)+rel2(1) & d2(:,1) <= XL(2)-rel(1)+rel2(1);
o8 = R.out_slack;   q8 = o8.t >= XL(1) & o8.t <= XL(2);
o2 = S.slack;       q2 = o2.t >= XL(1)-rel(1)+rel2(1) & o2.t <= XL(2)-rel(1)+rel2(1);
% Force floor raised from -18 to -5 kPa: below zero the slack force barely
% dips, so the extra headroom was blank. The pin (ML=1 = force=0) is DROPPED
% here -- keeping it at this floor would force the ML band into ~5 % of the
% panel height. Without it, iLengthAxis' natural mode gives the ML band the
% topFrac fraction of the panel (0.30 here, up from 0.22), which is what the
% "stretch the ML line" ask needs. ML=1 no longer aligns with the force-zero
% gridline; the ML axis stands on its own ticks instead.
% Top lowered 165 -> 145 now that the legend has gone: 165 existed to hold a
% one-row legend clear of the 2 mM re-stretch peaks (~131 kPa). 145 leaves
% those peaks 14 kPa of air and hands the rest back to the force trace.
YLF = [-5 145];
% Extended ML trace: forced to include a point at t=-0.10 s with ML=1.1, so the
% pre-slack baseline is unambiguously drawn from that point. The data itself
% goes back to -0.435 s at ML ~1.1004, but the muted line was too faint against
% the axis edge to read as "the muscle is being held here"; a synthesised
% anchor at exactly (-0.10, 1.1) fixes the start of the line the user sees.
t8ML = [-0.10; d8(m8,1) - rel(1)]; ml8ML = [1.1; d8(m8,2)/2];
[t8ML, is] = sort(t8ML); ml8ML = ml8ML(is);
iLengthAxis(ax, t8ML, ml8ML, cLen, FS, FSL, 0.30, YLF);
plot(d8(m8,1)-rel(1),  d8(m8,3), '-','Color',[K8 0.85],'LineWidth',2.4,'DisplayName','data 8 mM');
plot(d2(m2,1)-rel2(1), d2(m2,3), '-','Color',[K2 0.85],'LineWidth',2.4,'DisplayName','data 2 mM');
plot(o8.t(q8)-rel(1),  o8.Force(q8), '-','Color',C8,'LineWidth',LW_DENSE,'DisplayName','model 8 mM');
plot(o2.t(q2)-rel2(1), o2.Force(q2),'-','Color',C2,'LineWidth',LW_DENSE,'DisplayName','model 2 mM');
xlim(XL - rel(1)); ylim(YLF);
xlabel('time from first release (s)'); ylabel('force (kPa)');
title('Slack protocol - all five cycles');
iStyleAx(ax, FS, FSL, FST, FSLG);
% Annotations: the slack ML on each low plateau, and the restretch speed to the
% right of each restretch under the 1.1 hold -- the ladder used to be in the
% panel title, but on the actual trace each label reads as a caption of the
% segment it sits by, and the reader no longer has to count cycles.
iSlackAnnot(ax, rel-rel(1), Sv.velocitytable(Sv.velocitytable(:,2)>1,1)-rel(1), ...
            Sv.velocitytable(Sv.velocitytable(:,2)>1,2), d8(:,1)-rel(1), d8(:,2), ...
            1, LAB_FS, LAB_COL);
exportgraphics(fig, fullfile(outDir,'slide2_restretch_all.png'), 'Resolution',DPI);

%% ---- 2. the PRE-LAST cycle ---------------------------------------------
% Pre-window extended from -0.020 to -0.100 s: same 0.1 s of ML=1.1 baseline
% before the release as the whole-protocol panel above.
cD = nC - 1; WIN = [-0.100 0.400];
trest_all = Sv.velocitytable(Sv.velocitytable(:,2)>1,1);
vres_all  = Sv.velocitytable(Sv.velocitytable(:,2)>1,2);
fig = newFig(FIG_W2); ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
m8 = d8(:,1)-rel(cD)  >= WIN(1) & d8(:,1)-rel(cD)  <= WIN(2);
m2 = d2(:,1)-rel2(cD) >= WIN(1) & d2(:,1)-rel2(cD) <= WIN(2);
q8 = o8.t-rel(cD)  >= WIN(1) & o8.t-rel(cD)  <= WIN(2);
q2 = o2.t-rel2(cD) >= WIN(1) & o2.t-rel2(cD) <= WIN(2);
% Same design change as the whole-protocol panel: floor raised to -5 kPa,
% ML pin dropped, ML band takes topFrac = 0.30 of the panel.
YLD = [-5 145];
tDML = [-0.100; d8(m8,1) - rel(cD)]; mlDML = [1.1; d8(m8,2)/2];
[tDML, is] = sort(tDML); mlDML = mlDML(is);
iLengthAxis(ax, 1000*tDML, mlDML, cLen, FS, FSL, 0.30, YLD);
plot(1000*(d8(m8,1)-rel(cD)),  d8(m8,3), '-','Color',K8,'LineWidth',1.4,'DisplayName','data 8 mM');
plot(1000*(d2(m2,1)-rel2(cD)), d2(m2,3), '-','Color',K2,'LineWidth',1.4,'DisplayName','data 2 mM');
plot(1000*(o8.t(q8)-rel(cD)),  o8.Force(q8), '-','Color',C8,'LineWidth',2.2,'DisplayName','model 8 mM');
plot(1000*(o2.t(q2)-rel2(cD)), o2.Force(q2),'-','Color',C2,'LineWidth',2.2,'DisplayName','model 2 mM');
xlim(1000*WIN); ylim(YLD);
xlabel('time from slack release (ms)'); ylabel('force (kPa)');
% Redundancy note: with the ML label on the plateau and the ML/s label by the
% restretch, the title no longer needs to repeat the ladder value.
title(sprintf('Slack cycle %d', cD));
% This panel keeps a legend where the other two dropped theirs: it is the one
% most likely to be shown on its own, and it has somewhere to put it. Top-LEFT
% rather than top-centre -- the only thing in the upper band is the re-stretch
% peak (~132 kPa at t = 300 ms), which is hard right, so the left quarter above
% the pre-slack hold is genuinely empty.
legend('Location','northwest','NumColumns',2);
iStyleAx(ax, FS, FSL, FST, FSLG);
% Slack annotations dropped from this panel for now -- one cycle across two slide
% cells has enough space that the reader can read the title and the ML axis
% ticks without them, and the plateau/restretch values were already given by the
% whole-protocol panel two rows up.
exportgraphics(fig, fullfile(outDir,'slide2_restretch_det.png'), 'Resolution',DPI);

%% ---- 3. staircase ------------------------------------------------------
% F/F_iso, each trace over its own isometric: with no 2 mM staircase recorded,
% shape is the only comparison this panel can support.
Sd = load(['data/' R.p.stairs_velocitytableonfile]);
dt2 = Sd.datatable; vt2 = Sd.velocitytable;
t0s = vt2(find(vt2(:,2)~=0,1),1);
fig = newFig(FIG_W2); ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
n8 = iIso(R.out_stairs, t0s); n2 = iIso(S.stairs, t0s);
YLS = [0.80 1.60];
iLengthAxis(ax, dt2(:,1)-t0s, dt2(:,2)/2, cLen, FS, FSL, 0.14, YLS);
plot(dt2(:,1)-t0s, dt2(:,3), '-','Color',[K8 0.85],'LineWidth',2.0,'DisplayName','data 8 mM');
plot(R.out_stairs.t-t0s, R.out_stairs.Force/n8, '-','Color',C8,'LineWidth',LW_DENSE, ...
     'DisplayName','model 8 mM');
plot(S.stairs.t-t0s,     S.stairs.Force/n2,     '-','Color',C2,'LineWidth',LW_DENSE, ...
     'DisplayName','model 2 mM');
xlim([-0.10 dt2(end,1)-t0s]); ylim(YLS);
xlabel('time from first step (s)'); ylabel('F / F_{iso}');
title('Staircase - PREDICTION (8 mM data only)');

iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(outDir,'slide2_stairs.png'), 'Resolution',DPI);

fprintf('wrote 3 both-ATP slide panels @ %d dpi -> %s\n', DPI, outDir);
fprintf('  wide (2x1)  %.2f x %.2f in\n', FIG_W2, FIG_H);

%% ======================= helpers ========================================
function F0 = iIso(o, t0s)
% Model isometric force at the onset of the first length change.
    F0 = interp1(o.t, o.Force, t0s, 'linear','extrap');
    if F0 == 0; F0 = 1; end
end

function iSlackAnnot(ax, tRel, tRest, vRes, tData, slData, xScale, fs, col)
% ML label on each slack plateau + restretch speed to the right of each
% restretch, on the right-axis (ML) coordinates.
%
% Placement, learned from a round of iteration:
%   ML value  ABOVE the plateau line -- in the band between plateau and ML=1.
%             The old placement below wasted the empty strip between plateau
%             and axis floor; above it, force keeps that room instead.
%   Velocity  right of the restretch, tucked under the ML=1.1 hold line
%             (ml=1.07), i.e. in the empty band between the hold trace and the
%             force plateau. Same rule for every cycle including the last one
%             -- pushing the last label below the plateau put it off the plot.
%
% Inputs are all in the panel's OWN x units: seconds for the whole protocol
% (xScale = 1), milliseconds for the detail panel (xScale = 1000). The offset
% from restretch to velocity label is scaled the same way.
    n = numel(tRel);
    dtR = 0.030 * xScale;               % gap between restretch and its label
    yyaxis(ax,'right');
    for c = 1:n
        tS = tRel(c); tE = tRest(c);
        % ML at slack plateau: median over the middle-third of the plateau,
        % skipping the release ramp at the front and the pre-restretch settle
        % at the back so a fast late oscillation cannot skew it. Rounded to
        % 2 dp as the caller asked; the rounding itself is the label.
        span = tE - tS;
        wS = tS + 0.15*span; wE = tE - 0.05*span;
        win = tData >= wS & tData <= wE;
        if ~any(win); continue; end
        ml = round(median(slData(win))/2, 2);
        text(ax, (tS + tE)/2, ml + 0.02, sprintf('%.2f', ml), ...
             'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
             'FontSize',fs, 'Color',col, 'FontWeight','bold');
        text(ax, tE + dtR, 1.07, sprintf('%.0f ML/s', vRes(c)), ...
             'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
             'FontSize',fs, 'Color',col);
    end
    yyaxis(ax,'left');
end

function iLengthAxis(ax, t, ml, c, FS, FSL, topFrac, leftYlim, pinLeft)
% Commanded length on the right axis, kept in the lower band of the panel so
% force stays the figure's subject and the length reads as context. See
% MakeSlideFigures.m for the full argument; this is a verbatim copy so the two
% figure families place the second axis identically.
    if nargin < 7 || isempty(topFrac); topFrac = 1/3; end
    yyaxis(ax,'right');
    lo = min(ml); hi = max(ml); r = max(hi-lo, eps);
    plot(ax, t, ml, '-', 'Color', c, 'LineWidth', 1.3, 'HandleVisibility','off');
    if nargin >= 9 && ~isempty(pinLeft)
        fPin = (pinLeft - leftYlim(1)) / (leftYlim(2) - leftYlim(1));
        if topFrac <= fPin
            error('iLengthAxis:badPin','topFrac (%.3f) must exceed the pin height (%.3f)', ...
                  topFrac, fPin);
        end
        span = (hi - 1) / (topFrac - fPin);
        r0   = 1 - fPin*span;
        if lo < r0
            warning('iLengthAxis:clipped', ...
                'ML %.3f falls below the axis floor - give the force axis more room under %.4g', ...
                lo, pinLeft);
        end
    else
        pad = 0.12*r; span = (r + pad) / topFrac; r0 = lo - pad;
    end
    ylim(ax, [r0, r0 + span]);
    ylabel(ax, 'ML (1)');
    ax.YAxis(2).Color = c;
    ax.YAxis(2).FontSize = FS - 1;
    ax.YAxis(2).Label.FontSize = FSL - 1;
    ax.YAxis(2).Label.FontWeight = 'bold';
    ax.YAxis(2).TickValues = unique(round([lo, 1, hi], 2));
    yyaxis(ax,'left'); ylim(ax, leftYlim);
end

function iPlaceLegend(ax, lg, fx, fy)
% legend Position is normalised to the FIGURE, so an axes-relative placement
% has to be mapped through the axes rectangle.
    lg.Units = 'normalized';
    ap = ax.Position; lp = lg.Position;
    lg.Position = [ap(1) + fx*ap(3), ap(2) + fy*ap(4), lp(3), lp(4)];
end

function iStyleAx(ax, FS, FSL, FST, FSLG)
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
