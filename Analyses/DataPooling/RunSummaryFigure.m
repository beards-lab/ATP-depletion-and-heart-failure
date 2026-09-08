% RunSummaryFigure.m -- one page carrying the DataPooling argument.
% Run RunNormalizationLadder, RunPoolingComparison and RunFVConsistency first.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results');

L  = load(fullfile(resDir,'normalisation_ladder.mat'));
P  = load(fullfile(resDir,'pooling_comparison.mat'));
FVc= load(fullfile(resDir,'fv_consistency.mat'));
T  = buildPooledTarget();

f = figure('Position',[40 40 1420 880],'Color','w');
tiledlayout(f,2,3,'Padding','compact','TileSpacing','compact');

% ---- 1. the nuisance, and what a scale removes -------------------------
nexttile; hold on; box on
pick = {'none/prep','Fpeak/prep+rd','Fiso20b/prep+rd','MoR/prep+rd','Fiso20a/prep+rd', ...
        'Fsteady/prep+rd','Fss22/prep+rd'};
sb = arrayfun(@(q) 100*L.res(strcmp({L.res.name},q{1})).sBar, pick);
sm = arrayfun(@(q) 100*L.res(strcmp({L.res.name},q{1})).sMed, pick);
b = barh([sb; sm]'); b(1).FaceColor=[.72 .74 .80]; b(2).FaceColor=[.16 .32 .58];
set(gca,'YTick',1:numel(pick),'YTickLabel',strrep(pick,'/prep+rd',''),'YDir','reverse');
xlabel('between-preparation relative SD (%)  [leave-one-feature-out]');
legend({'RMS over features','median'},'Location','southeast');
title('(a) any sensible force anchor removes half the spread');

% ---- 2. where along the trace ------------------------------------------
nexttile; hold on; box on
plot(P.TG, 100*median(P.Sraw,1,'omitnan'), 'LineWidth',1.6);
plot(P.TG, 100*median(P.Snrm,1,'omitnan'), 'LineWidth',1.6);
yline(11,'k:','HandleVisibility','off'); yline(2.2,'k:','HandleVisibility','off');
text(0.13, 13.5, 'slack hold: 11 %','FontSize',8);
text(0.36, 4.5, 'recovery at the reference length: 2.2 %','FontSize',8);
xlabel('time from slack release (s)'); ylabel('between-prep rel. SD (%)');
xlim([-0.02 0.57]); ylim([0 34]);
legend({'raw kPa','/ F_{ss22}'},'Location','northeast');
title('(b) what is left is a LENGTH effect, not a scale');

% ---- 3. one scale per prep is the right model --------------------------
nexttile; hold on; box on
bar(100*L.sv(1:min(3,numel(L.sv))), 'FaceColor',[.16 .32 .58]);
set(gca,'XTick',1:3,'XTickLabel',{'PC1','PC2','PC3'});
ylabel('% of between-prep log-variance'); ylim([0 100]);
text(1, 100*L.sv(1)+5, sprintf('|cos| with the flat\ndirection = %.3f', L.cosFlat), ...
     'HorizontalAlignment','center','FontSize',9);
title('(c) the between-prep difference IS one common scale');

% ---- 4. pooling route ---------------------------------------------------
nexttile; hold on; box on
gp = 100*P.gap; gp(~isfinite(gp)) = 0;
bb = bar(gp); bb.FaceColor = 'flat';
big = abs(gp) > 20; bb.CData(big,:) = repmat([.75 .18 .18], sum(big),1);
bb.CData(~big,:) = repmat([.72 .74 .80], sum(~big),1);
set(gca,'XTick',1:numel(P.NM),'XTickLabel',P.NM,'XTickLabelRotation',45);
ylabel('features-of-mean  vs  mean-of-features (%)'); yline(0,'k-');
title('(d) averaging TRACES biases the re-stretch rate by 24 %');

% ---- 5. the FV point ----------------------------------------------------
nexttile; hold on; box on
for j = 1:2; plot(FVc.BK{j,2}/FVc.BK{j,2}(1), FVc.AllV,'o-','LineWidth',1.4); end
for j = 3:4; plot(FVc.BK{j,2}/FVc.BK{j,2}(1), FVc.AllV,'s--','LineWidth',1.1); end
plot(FVc.f8, [2 2],'k^','MarkerFaceColor','k','MarkerSize',10);
plot(FVc.f2, [2 2],'kv','MarkerFaceColor',[.65 .65 .65],'MarkerSize',10);
set(gca,'YDir','reverse'); xlabel('P / P_0'); ylabel('shortening velocity (ML/s)');
legend([cellfun(@(s) s, FVc.BK(:,1)','uni',0) {'FV2 8 mM','FV2 2 mM'}], ...
       'Location','southeast','FontSize',7);
title('(e) the in-prep F-V point lands on the borrowed 8 mM curve');

% ---- 6. the target ------------------------------------------------------
nexttile; hold on; box on
[~, o] = sort(T.s8);
bar(100*T.s8raw(o), 'FaceColor',[.80 .82 .86]);
plot(1:numel(o), 100*T.s8(o), 'o-','Color',[.16 .32 .58],'LineWidth',1.6,'MarkerFaceColor','w');
zz = T.z(o) >= 2;
plot(find(zz), 100*T.s8(o(zz)), 'o','Color',[.75 .18 .18],'MarkerFaceColor',[.75 .18 .18]);
set(gca,'XTick',1:numel(o),'XTickLabel',T.name(o),'XTickLabelRotation',55,'FontSize',8);
ylabel('relative SD (%)');
legend({'raw (2 dof)','moderated','carries the ATP effect (z\geq2)'},'Location','northwest','FontSize',8);
title('(f) the pooled target and its honest error bars');

exportgraphics(f, fullfile(resDir,'datapooling_summary.png'), 'Resolution',150);
fprintf('saved -> %s\n', fullfile(resDir,'datapooling_summary.png'));
