% MakeFVFiguresBothATP.m -- force-velocity and power-velocity at BOTH [ATP], on
% the fine simulated grid, in the same slide language as the other panels here.
%
% Ported from Analyses/AtpEffectHypotheses/RunFVFineBothATP.m in the
% ATP-wt-lowatp worktree, with two deliberate changes so the four figure
% families in these folders read as one set:
%
%   DATA COLOURS follow the slack/staircase panels -- black for 8 mM, dark grey
%   for 2 mM. The original used dark blue / dark red markers, which tied the
%   data to the model hue and made "which of these two lines is the
%   measurement" a question the reader had to work out.
%
%   NO RUNG SUBTITLE, as on every other panel in these folders: which rung a
%   figure belongs to is carried by the folder it sits in, and is printed to
%   stdout, which is where a provenance check belongs.
%
%   LEGEND ON FORCE-VELOCITY ONLY. The rest of the set drops it (the colour key
%   is stated once, on the steady-force feature panel), but F-V keeps one --
%   this panel usually travels on its own rather than as part of the folder,
%   and the falling curves leave its top-right corner empty, so the legend
%   costs nothing here. The power panel beside it needs no second copy.
%
% The BAND is the spread between the two independent source curves -- Baker2022
% and IsovelocityForFilip2021 -- each normalised by ITS OWN isometric force
% before averaging, because absolute force is not comparable between
% preparations (56.4 vs 67.4 kPa at 8 mM). Both sets carry a 2 mM column, so the
% same construction works at both concentrations. Note this is a 2-preparation
% spread, NOT the 3-preparation slack envelope used by the feature panels.
%
% NORMALISED POWER = (F/F_iso) * v, so its unit is ML/s = s^-1, not watts. The
% un-normalised product F[kPa]*v[s^-1] is stress x strain rate = kW/m^3, power
% per unit muscle volume; dividing by F_iso expresses that in units of the
% preparation's own isometric stress. Multiply back by F_iso to recover kW/m^3.
%
% Input:   results/bothatp/<RUNG>/fv_fine.mat   (RunFVFineAllRungs.m, in the
%                                                ATP-wt-lowatp worktree)
% Figures -> results/bothatp/<RUNG>/slide2_fv.png        force + power, 2 cells
%            results/bothatp/<RUNG>/slide2_fv_force.png  force only,   1 cell
%            results/bothatp/<RUNG>/slide2_fv_power.png  power only,   1 cell

clearvars -except RUNG RUNGS iR; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));

if ~exist('RUNG','var'); RUNG = 'H4'; end
outDir = fullfile(here,'results','bothatp',RUNG);
FVF = fullfile(outDir,'fv_fine.mat');
if ~isfile(FVF)
    error('MakeFVFiguresBothATP:noSweep', ...
        ['missing %s -- run Analyses/AtpEffectHypotheses/RunFVFineAllRungs.m ' ...
         'in the ATP-wt-lowatp worktree first (the ATP gate code lives only there)'], FVF);
end
C = load(FVF);
gateTxt = strjoin(strrep(fieldnames(C.gates),'pAtp_',''), ', ');
if isempty(gateTxt); gateTxt = '(none)'; end
fprintf('%s | %s | gates: %s | built %s\n', C.rung, C.snap, gateTxt, C.built);
fprintf('  F_iso: 8 mM %.2f | 2 mM %.2f kPa  (2/8 = %.3f)\n', ...
        C.F8iso, C.F2iso, C.F2iso/C.F8iso);

%% ---- style, shared with the other panels in this folder ----------------
FIG_W = 4.6; FIG_H = FIG_W*2/3; GAP = 0.14; DPI = 300;
FS = 14; FSL = 15; FST = 15; FSLG = 11;
C8 = [0.00 0.35 0.70];  C2 = [0.80 0.20 0.10];   % models
K8 = [0.04 0.05 0.07];  K2 = [0.44 0.45 0.48];   % data, as on the slack panels
B8 = [0.62 0.74 0.88];  B2 = [0.92 0.72 0.66];   % bands
LW = 2.8; MS = 7;

% 0.25 ML/s is a protocol artifact, not a model prediction: runFVExperiment
% takes the isometric reference at SL 2.0 but starts every shortening run at
% SL 2.2, so at the slowest speed the 0.4 s ramp spends most of its time at
% longer, stronger lengths and the point lands ABOVE isometric. It stays in the
% .mat; this only keeps it out of a figure where it would read as a claim about
% force-velocity rather than about the protocol.
VPM = 6; VDR = 0.25;
mM = C.vFine <= VPM + 1e-9;
for vv = VDR(:)'; mM = mM & abs(C.vFine - vv) > 1e-9; end
mD = C.vD <= VPM + 1e-9;

% Vmax: where each model curve reaches zero force. A real prediction -- low ATP
% lowers Vmax -- but past it F/F0 goes very slightly negative (~-1e-4), which
% draws a cliff below zero and a spurious kink in power. Clamp for drawing and
% END each curve AT its Vmax, rather than running it flat along zero, which
% would imply the model still predicts something out there.
Fn8p = max(C.Fn8,0); Pw8p = Fn8p .* C.vFine;
Fn2p = max(C.Fn2,0); Pw2p = Fn2p .* C.vFine;
i8 = find(C.Fn8 <= 0,1); v8max = NaN; if ~isempty(i8); v8max = C.vFine(i8); end
i2 = find(C.Fn2 <= 0,1); v2max = NaN; if ~isempty(i2); v2max = C.vFine(i2); end
m8 = mM; if isfinite(v8max); m8 = m8 & C.vFine <= v8max + 1e-9; end
m2 = mM; if isfinite(v2max); m2 = m2 & C.vFine <= v2max + 1e-9; end
fprintf('  Vmax (F/F0 -> 0): 8 mM %s | 2 mM %s ML/s\n', ...
        iNum(v8max,C.VMAX), iNum(v2max,C.VMAX));
fprintf('  peak power: 8 mM %.4f at %.2f | 2 mM %.4f at %.2f ML/s\n', ...
        max(C.Pw8), C.vFine(find(C.Pw8==max(C.Pw8),1)), ...
        max(C.Pw2), C.vFine(find(C.Pw2==max(C.Pw2),1)));

close all force;

%% ---- combined: force + power, two cells wide ---------------------------
fig = figure('Units','inches','Position',[1 1 2*FIG_W+GAP FIG_H],'Color','w');
tl = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
ax = nexttile; iForcePanel(ax, C, mD, m8, m2, Fn8p, Fn2p, v2max, VPM, ...
                           C8,C2,K8,K2,B8,B2, LW,MS, FS,FSL,FST,FSLG, true);
ax = nexttile; iPowerPanel(ax, C, mD, m8, m2, Pw8p, Pw2p, VPM, ...
                           C8,C2,K8,K2,B8,B2, LW,MS, FS,FSL,FST,FSLG, false);

exportgraphics(fig, fullfile(outDir,'slide2_fv.png'), 'Resolution',DPI);

%% ---- the same two panels, standalone, one cell each --------------------
% The slide layout in MakeSlide.m places a single 3:2 cell for F-V; these are
% the pieces for that, so a deck can take either the pair or one of them.
fig = figure('Units','inches','Position',[1 1 FIG_W FIG_H],'Color','w');
ax = axes(fig); iForcePanel(ax, C, mD, m8, m2, Fn8p, Fn2p, v2max, VPM, ...
                            C8,C2,K8,K2,B8,B2, LW,MS, FS,FSL,FST,FSLG, true);

exportgraphics(fig, fullfile(outDir,'slide2_fv_force.png'), 'Resolution',DPI);

fig = figure('Units','inches','Position',[1 1 FIG_W FIG_H],'Color','w');
ax = axes(fig); iPowerPanel(ax, C, mD, m8, m2, Pw8p, Pw2p, VPM, ...
                            C8,C2,K8,K2,B8,B2, LW,MS, FS,FSL,FST,FSLG, false);

exportgraphics(fig, fullfile(outDir,'slide2_fv_power.png'), 'Resolution',DPI);

fprintf('  wrote slide2_fv.png, slide2_fv_force.png, slide2_fv_power.png -> %s\n', outDir);

%% ======================= helpers ========================================
function iForcePanel(ax, C, mD, m8, m2, Fn8p, Fn2p, v2max, VPM, ...
                     C8,C2,K8,K2,B8,B2, LW,MS, FS,FSL,FST,FSLG, withLegend)
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    iBand(ax, C.vD(mD), C.mu8(mD), C.sd8(mD), B8);
    iBand(ax, C.vD(mD), C.mu2(mD), C.sd2(mD), B2);
    p1 = plot(ax, C.vFine(m8), Fn8p(m8), '-','Color',C8,'LineWidth',LW);
    p2 = plot(ax, C.vFine(m2), Fn2p(m2), '-','Color',C2,'LineWidth',LW);
    d1 = plot(ax, C.vD(mD), C.mu8(mD), 'o','Color',K8,'MarkerFaceColor',K8,'MarkerSize',MS);
    d2 = plot(ax, C.vD(mD), C.mu2(mD), 'o','Color',K2,'MarkerFaceColor',K2,'MarkerSize',MS);
    if isfinite(v2max) && v2max <= VPM
        xline(ax, v2max, ':','Color',C2,'LineWidth',2.0,'HandleVisibility','off');
        % Anchored near the FLOOR, not mid-axis: at slide point sizes the
        % rotated label spans a third of the panel height. Kept low even now the
        % legend has gone: the only thing it crosses down here is the tail of
        % the curves, at ~0.03.
        text(ax, v2max-0.06, 0.08, sprintf('V_{max} %.2f', v2max), 'Color',C2, ...
             'FontSize',FSLG-1,'FontWeight','bold','Rotation',90, ...
             'HorizontalAlignment','left');
    end
    xlim(ax,[0 VPM]); ylim(ax,[0 1.02]);
    xlabel(ax,'shortening velocity (ML/s)'); ylabel(ax,'F / F_{iso}');
    title(ax,'Force-velocity');
    if withLegend
        legend(ax, [d1 p1 d2 p2], ...
            {'data 8 mM','model 8 mM','data 2 mM','model 2 mM'}, ...
            'Location','northeast','NumColumns',2);
    end
    iStyleAx(ax, FS, FSL, FST, FSLG);
end

function iPowerPanel(ax, C, mD, m8, m2, Pw8p, Pw2p, VPM, ...
                     C8,C2,K8,K2,B8,B2, LW,MS, FS,FSL,FST,FSLG, withLegend)
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    iBand(ax, C.vD(mD), C.pm8(mD), C.ps8(mD), B8);
    iBand(ax, C.vD(mD), C.pm2(mD), C.ps2(mD), B2);
    p1 = plot(ax, C.vFine(m8), Pw8p(m8), '-','Color',C8,'LineWidth',LW);
    p2 = plot(ax, C.vFine(m2), Pw2p(m2), '-','Color',C2,'LineWidth',LW);
    d1 = plot(ax, C.vD(mD), C.pm8(mD), 'o','Color',K8,'MarkerFaceColor',K8,'MarkerSize',MS);
    d2 = plot(ax, C.vD(mD), C.pm2(mD), 'o','Color',K2,'MarkerFaceColor',K2,'MarkerSize',MS);
    xlim(ax,[0 VPM]);
    top = max([Pw8p(m8), Pw2p(m2), C.pm8(mD)+C.ps8(mD), C.pm2(mD)+C.ps2(mD)]);
    ylim(ax,[0 top*1.10]);           % no legend here -- ordinary headroom
    xlabel(ax,'shortening velocity (ML/s)');
    % The unit is ML/s = s^-1 (see the header): (F/F_iso) is dimensionless and
    % ML is a length ratio, so this is a strain rate, not watts.
    ylabel(ax,'normalised power (s^{-1})');
    title(ax,'Force-velocity (power)');
    if withLegend
        legend(ax, [d1 p1 d2 p2], ...
            {'data 8 mM','model 8 mM','data 2 mM','model 2 mM'}, ...
            'Location','northeast','NumColumns',2);
    end
    iStyleAx(ax, FS, FSL, FST, FSLG);
end

function iBand(ax, x, mu, sd, col)
    g = ~isnan(mu) & ~isnan(sd);
    if ~any(g); return; end
    fill(ax, [x(g) fliplr(x(g))], [mu(g)-sd(g) fliplr(mu(g)+sd(g))], col, ...
         'FaceAlpha',0.45,'EdgeColor','none','HandleVisibility','off');
end

function s = iNum(v, vmax)
    if isnan(v); s = sprintf('>%g', vmax); else; s = sprintf('%.2f', v); end
end

function iStyleAx(ax, FS, FSL, FST, FSLG)
    set(ax, 'FontSize',FS, 'LineWidth',1.1, 'GridAlpha',0.16, 'Box','on', ...
            'TickDir','out', 'TickLength',[0.015 0.015]);
    ax.XLabel.FontSize = FSL;  ax.XLabel.FontWeight = 'bold';
    ax.YLabel.FontSize = FSL;  ax.YLabel.FontWeight = 'bold';
    ax.Title.FontSize  = FST;  ax.Title.FontWeight  = 'bold';
    lg = get(ax,'Legend');
    if ~isempty(lg); lg.FontSize = FSLG; lg.Box = 'on'; end
    ax.LooseInset = max(ax.TightInset, [0.02 0.02 0.01 0.01]);
end
