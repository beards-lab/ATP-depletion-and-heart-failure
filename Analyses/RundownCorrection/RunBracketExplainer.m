% RunBracketExplainer.m
% What "the bracket" means, drawn from the raw timecourses.
%
% THE POINT OF CONFUSION THIS FIGURE EXISTS TO REMOVE
% "Bracket" refers to TIME, not to force. On 03/27 the SAME 8 mM protocol was
% run twice -- once at the start of the session and once at the end -- so the
% 2 mM run sits BETWEEN them in time. Force itself only ever goes DOWN:
%
%       8 mM  @ t=0      68.4 kPa
%       2 mM  @ t=132 s  83.6 kPa   <- higher, but that is the ATP EFFECT
%       8 mM  @ t=694 s  56.5 kPa   <- the same condition, 17.5% weaker
%
% The two 8 mM measurements are the only two points in the whole session
% recorded under IDENTICAL conditions, so the line between them is pure
% rundown, uncontaminated by ATP. Everything else is measured against it.
%
% The ATP effect is then the VERTICAL GAP between the measured 2 mM force and
% that rundown line interpolated to the moment the 2 mM run happened -- not the
% raw 2 mM / 8 mM ratio, which compares a fresher fibre with a tireder one.
%
% Nothing here "inflates" a force. Correcting the later run upward is a
% statement about what it WOULD have measured on an undamaged fibre, not a
% claim that force recovered.
%
% Outputs -> results/bracket_explained.png

clear; close all; clc;
cd(fileparts(which('RunBracketExplainer')));
addpath(genpath('../..'));
resDir = 'results';  if ~isfolder(resDir), mkdir(resDir); end
D = '../../data/03 27 2026 M/';

PHI = 0.452;                     % permanent fraction (RunRundownCorrection.m)
zA  = [71.50 72.35];             % isometric plateau at ML 1.0, before the staircase
zB  = [77.70 78.30];             % isometric plateau at ML 1.0, after the last slack

RUN = struct( ...
 'lab',  {'8 mM (fresh)','2 mM','PNB+Mava','8 mM (repeat)'}, ...
 'file', {'02_Merged_8mM_Active.txt','03_Merged_2mM_Active.txt', ...
          '06_Merged_8mM_Active_PNB_Mava.txt','04_Merged_8mM_Active_repeat.txt'}, ...
 'sess', {0, 132, 348, 694}, ...                        % session time (s) of each activation
 'col',  {[0 0.45 0.74],[0.85 0.33 0.10],[0.5 0.5 0.5],[0 0.30 0.50]});

for i = 1:numel(RUN)
    M = readmatrix([D RUN(i).file]);
    RUN(i).t = M(:,1);  RUN(i).F = M(:,3);
    RUN(i).Fplat = mean(M(M(:,1)>=zA(1) & M(:,1)<=zA(2), 3), 'omitnan');
    m = (M(:,1)>=zA(1)&M(:,1)<=zA(2)) | (M(:,1)>=zB(1)&M(:,1)<=zB(2));
    p = polyfit(M(m,1), M(m,3), 1);  RUN(i).slope = p(1);
    Fs = movmean(M(:,3),51); thr = 0.2*prctile(Fs,99);
    on = find(Fs>thr,1); off = find(Fs>thr,1,'last');
    RUN(i).ton = M(on,1);  RUN(i).Tact = M(off,1)-M(on,1);
end

i1 = 1; i2 = 2; i4 = 4;                       % fresh 8 mM, 2 mM, repeat 8 mM
fprintf('=== the three numbers ===\n');
for i = [i1 i2 i4]
    fprintf('  %-14s session t = %3d s   plateau force = %5.2f kPa\n', ...
            RUN(i).lab, RUN(i).sess, RUN(i).Fplat);
end
fprintf('\n  Force NEVER rises. The 8 mM condition falls %.2f -> %.2f kPa (-%.1f%%).\n', ...
        RUN(i1).Fplat, RUN(i4).Fplat, 100*(1-RUN(i4).Fplat/RUN(i1).Fplat));
fprintf('  The 2 mM run is higher than both only because low ATP raises force.\n');

%% ── Where does the 8 mM line sit at the moment of the 2 mM run? ────────────
% On EFFECTIVE ACTIVATED time: damage accrues only while a run is generating
% force, at phi * |its own slope| kPa per activated second.
dose = @(i) PHI*abs(RUN(i).slope)*RUN(i).Tact;      % kPa of permanent damage per run
d1 = dose(i1);                       % damage done during the fresh 8 mM run
d2 = dose(i2);                       % ... during the 2 mM run
dTot = d1 + d2;                      % PNB+Mava does none (blebbistatin)

F8_at_2mM = RUN(i1).Fplat - d1;      % 8 mM, interpolated to the 2 mM moment
ratio_raw = RUN(i2).Fplat / RUN(i1).Fplat;
ratio_cor = RUN(i2).Fplat / F8_at_2mM;

fprintf('\n=== the bracket, closed ===\n');
fprintf('  damage during the fresh 8 mM run : %.2f kPa\n', d1);
fprintf('  damage during the 2 mM run       : %.2f kPa\n', d2);
fprintf('  predicted total by the repeat    : %.2f kPa\n', dTot);
fprintf('  MEASURED loss (8 mM fresh->rpt)  : %.2f kPa   <- this is what closes it\n', ...
        RUN(i1).Fplat-RUN(i4).Fplat);
fprintf('\n  8 mM interpolated to the 2 mM moment: %.2f kPa\n', F8_at_2mM);
fprintf('  ATP effect  raw       = %.2f / %.2f = x%.3f\n', RUN(i2).Fplat, RUN(i1).Fplat, ratio_raw);
fprintf('  ATP effect  corrected = %.2f / %.2f = x%.3f\n', RUN(i2).Fplat, F8_at_2mM, ratio_cor);

%% ── Figure ─────────────────────────────────────────────────────────────────
fig = figure(960); clf; set(fig,'Position',[20 20 1560 820]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% (a) the four raw timecourses, on their own protocol clock
nexttile([1 2]); hold on; box on;
for i = 1:numel(RUN)
    w = RUN(i).t >= 63 & RUN(i).t <= 82;
    plot(RUN(i).t(w), RUN(i).F(w), '-', 'Color', RUN(i).col, 'LineWidth', 1.1, ...
         'DisplayName', sprintf('%s   (session t = %d s)', RUN(i).lab, RUN(i).sess));
end
patch([zA(1) zA(2) zA(2) zA(1)],[0 0 105 105],[0 0 0],'FaceAlpha',0.06,'EdgeColor','none','HandleVisibility','off');
text(mean(zA), 101, 'plateau window', 'HorizontalAlignment','center','FontSize',7);
xlabel('protocol time within a run (s)'); ylabel('force (kPa)'); ylim([-5 105]);
title('(a) The four 03/27 runs — identical protocol, four different fibre states');
legend('Location','southwest','FontSize',8);

% (b) placed on SESSION time: the actual experiment, start to finish
nexttile; hold on; box on;
for i = 1:numel(RUN)
    w = RUN(i).t >= 63 & RUN(i).t <= 82;
    plot(RUN(i).sess + (RUN(i).t(w)-RUN(i).ton), RUN(i).F(w), '-', ...
         'Color', RUN(i).col, 'LineWidth', 1.0);
end
plot([RUN(i1).sess RUN(i4).sess], [RUN(i1).Fplat RUN(i4).Fplat], 'k--o', ...
     'LineWidth', 2, 'MarkerFaceColor','k','MarkerSize',8);
plot(RUN(i2).sess, RUN(i2).Fplat, 's', 'Color', RUN(i2).col, 'MarkerFaceColor', RUN(i2).col, 'MarkerSize',10);
text(360, 88, 'the BRACKET', 'FontSize',9,'HorizontalAlignment','center','FontWeight','bold');
text(360, 82, 'same condition, measured twice', 'FontSize',7,'HorizontalAlignment','center','Color',[.3 .3 .3]);
xlabel('session time (s)'); ylabel('force (kPa)'); ylim([-5 105]); xlim([-40 740]);
title('(b) The whole session on one clock');

% (c) the damage STAIRCASE: the fibre ages only while it is activated
nexttile([1 2]); hold on; box on;
% Build the staircase: flat in relaxing solution, dropping during each activation.
tS = [RUN(i1).sess, RUN(i1).sess+RUN(i1).Tact, RUN(i2).sess, RUN(i2).sess+RUN(i2).Tact, RUN(i4).sess];
fS = [RUN(i1).Fplat, RUN(i1).Fplat-d1, RUN(i1).Fplat-d1, RUN(i1).Fplat-d1-d2, RUN(i1).Fplat-d1-d2];
plot(tS, fS, 'k-', 'LineWidth', 2.5);
for a = [1 3]      % shade the two force-generating activations
    patch([tS(a) tS(a+1) tS(a+1) tS(a)], [40 40 95 95], [0 0 0], ...
          'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
end
plot(RUN(i1).sess, RUN(i1).Fplat, 'o','Color',RUN(i1).col,'MarkerFaceColor',RUN(i1).col,'MarkerSize',13);
plot(RUN(i4).sess, RUN(i4).Fplat, 'o','Color',RUN(i4).col,'MarkerFaceColor',RUN(i4).col,'MarkerSize',13);
plot(RUN(i2).sess, RUN(i2).Fplat, 's','Color',RUN(i2).col,'MarkerFaceColor',RUN(i2).col,'MarkerSize',13);
plot(RUN(i2).sess, F8_at_2mM, 'kd','MarkerFaceColor','w','MarkerSize',13,'LineWidth',1.8);

plot([RUN(i2).sess RUN(i2).sess], [F8_at_2mM RUN(i2).Fplat], '-', 'Color',[0.85 0.33 0.10], 'LineWidth', 4);
text(RUN(i2).sess+25, mean([F8_at_2mM RUN(i2).Fplat])+1, ...
     sprintf('ATP EFFECT\n%.1f \\rightarrow %.1f kPa\n= \\times%.2f', F8_at_2mM, RUN(i2).Fplat, ratio_cor), ...
     'FontSize',10,'Color',[0.85 0.33 0.10],'FontWeight','bold');
text(RUN(i1).sess-25, RUN(i1).Fplat+4, sprintf('8 mM fresh  %.1f', RUN(i1).Fplat),'FontSize',9);
text(RUN(i4).sess-40, RUN(i4).Fplat-4.5, sprintf('8 mM repeat  %.1f  (-%.0f%%)', ...
     RUN(i4).Fplat, 100*(1-RUN(i4).Fplat/RUN(i1).Fplat)),'FontSize',9,'HorizontalAlignment','right');
text(RUN(i2).sess-8, F8_at_2mM-3.5, sprintf('8 mM would have been %.1f', F8_at_2mM), ...
     'FontSize',8,'HorizontalAlignment','right');
text(430, 92, 'shaded = fibre ACTIVATED (15.3 s each); the staircase only drops there', ...
     'FontSize',8,'Color',[.3 .3 .3],'HorizontalAlignment','center');
text(430, 88, 'flat elsewhere: no force generated, no damage', ...
     'FontSize',8,'Color',[.3 .3 .3],'HorizontalAlignment','center');
xlabel('session time (s)'); ylabel('plateau force at SL 2.0 (kPa)');
ylim([40 95]); xlim([-60 760]);
title('(c) The damage staircase — \phi\cdot|slope|\cdot T_{act} per activation. Force only ever falls.');

% (d) raw vs corrected ratio
nexttile; hold on; box on;
b = bar([ratio_raw ratio_cor]); b.FaceColor='flat'; b.CData=[.6 .6 .6; .2 .5 .8];
set(gca,'XTick',1:2,'XTickLabel',{'raw','vs rundown line'});
ylabel('ATP force ratio 2 mM / 8 mM'); ylim([1 1.45]);
for k=1:2
    v=[ratio_raw ratio_cor]; text(k, v(k)+0.015, sprintf('%.3f',v(k)),'HorizontalAlignment','center','FontSize',10);
end
title(sprintf('(d) 03/27 alone: %.3f -> %.3f', ratio_raw, ratio_cor));

sgtitle('What "the bracket" is: the same 8 mM protocol run twice, so rundown can be seen without ATP');
exportgraphics(fig, fullfile(resDir,'bracket_explained.png'), 'Resolution', 150);
fprintf('\nSaved %s\n', fullfile(resDir,'bracket_explained.png'));
