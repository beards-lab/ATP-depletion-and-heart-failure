% PlotRealign.m — the compliant-realignment result in three panels.
%   (a) the authority/bill trade curve, and where the mechanism saturates
%   (b) the LOAD weighting is the selectivity — the load-blind control breaks ktr,
%       and the shared-pool dwell (M2a) breaks it far worse
%   (c) rsK per cycle: baseline, best compensated M1, data

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
res = fullfile(here, 'results');

Lb = load(fullfile(res,'baseline.mat')); B = Lb.B; rsK_d = Lb.rsK_d;
Lr = load(fullfile(res,'loadrealign.mat'));
Lc = load(fullfile(res,'compensate.mat'));

get1 = @(S, nm) S.C{find(cellfun(@(c) ~isempty(c) && c.ok && strcmp(c.name, nm), S.C), 1)};

f = figure('Position',[80 80 1500 460],'Color','w');
tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');

%% (a) authority vs bill ----------------------------------------------------
nexttile; hold on;
gain = {'M1 k=1  tau=20 F=25','M1 k=3  tau=20 F=25','M1 k=6  tau=20 F=25','M1 k=12 tau=20 F=25'};
gv = [1 3 6 12]; gx = nan(size(gv)); gy = nan(size(gv));
for i=1:numel(gain); c = get1(Lr, gain{i}); gx(i) = c.rsK_x; gy(i) = c.L2; end
plot(gx, gy, 'o-', 'LineWidth', 2, 'MarkerFaceColor','w', 'DisplayName','M1 gain k_{cr}, \tau=20 ms');
for i=1:numel(gv); text(gx(i), gy(i), sprintf('  k=%d', gv(i)), 'FontSize',8); end
cmp = {'M1 k=3 + kA2re=30','M1 k=3 + kA2re30+etaM','M1 k=6 + kA2re=30','M1 k=6 + kA2re30+etaM'};
cx = nan(1,4); cy = nan(1,4);
for i=1:numel(cmp); c = get1(Lc, cmp{i}); cx(i) = c.rsK_x; cy(i) = c.L2; end
plot(cx, cy, 's--', 'LineWidth', 2, 'MarkerFaceColor','w', 'DisplayName','+ kA2re / \eta_M compensation');
plot(B.rsK_x, B.L2, 'kp', 'MarkerSize', 16, 'MarkerFaceColor','k', 'DisplayName','baseline');
xline(1, 'k--', 'data', 'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
set(gca,'YScale','log'); grid on; legend('Location','northeast');
xlabel('rsK  (model / data)'); ylabel('feature cost L2');
title(sprintf('(a) full authority over rsK\nx2.35 -> x1.12, bill is restretch SHAPE'));

%% (b) selectivity ----------------------------------------------------------
nexttile;
sel = {'baseline', B.ktr_x, B.rsK_x
       'M1 load-weighted',  get1(Lr,'M1 k=6  tau=20 F=25').ktr_x,  get1(Lr,'M1 k=6  tau=20 F=25').rsK_x
       'M1 LOAD-BLIND',     get1(Lr,'M1 k=6  tau=20 LOADBLIND').ktr_x, get1(Lr,'M1 k=6  tau=20 LOADBLIND').rsK_x
       'M2a SRD dwell .30', get1(Lr,'M2a SRD dwell K=.30').ktr_x,  get1(Lr,'M2a SRD dwell K=.30').rsK_x
       'M2a SRD dwell .10', get1(Lr,'M2a SRD dwell K=.10').ktr_x,  get1(Lr,'M2a SRD dwell K=.10').rsK_x
       'M2b D0 dwell .10',  get1(Lr,'M2b recoil+D0 dwell .10').ktr_x, get1(Lr,'M2b recoil+D0 dwell .10').rsK_x};
bar(cell2mat(sel(:,2:3))); set(gca,'XTickLabel', sel(:,1), 'XTickLabelRotation', 20);
yline(1,'k--','data'); ylabel('model / data'); grid on;
legend({'ktr  (must stay 1)','rsK  (must fall to 1)'}, 'Location','northwest');
title(sprintf('(b) the LOAD weighting IS the selectivity\nremove it and ktr breaks; a shared-pool dwell breaks it worse'));

%% (c) per-cycle ------------------------------------------------------------
nexttile; hold on;
bestn = 'M1 k=6 + kA2re30+etaM'; best = get1(Lc, bestn);
plot(1:5, B.rsK_m, 'o-', 'LineWidth', 2, 'DisplayName','baseline');
plot(1:5, best.rsK_m, 's-', 'LineWidth', 2, 'DisplayName', bestn);
plot(1:5, rsK_d, 'k^--', 'LineWidth', 2, 'MarkerFaceColor','k', 'DisplayName','data 03/27 8 mM');
xlabel('re-stretch cycle'); ylabel('rsK  (s^{-1})'); grid on; legend('Location','northeast');
xlim([0.7 5.3]);
title(sprintf('(c) the level lands on the data\nbut the model is now too FLAT across cycles'));

exportgraphics(f, fullfile(res,'realign.png'), 'Resolution', 140);
fprintf('best: %s -> rsK x%.2f (cost %.3f), ktr x%.2f, steady %.1f, L2 %.1f\n', ...
        bestn, best.rsK_x, best.ct2(find(strcmp(cellfun(@(x) regexprep(strtok(x,'|'),'\[.*$',''), B.fn,'uni',0),'rsK'),1)), ...
        best.ktr_x, best.steady, best.L2);
fprintf('saved results/realign.png\n');
