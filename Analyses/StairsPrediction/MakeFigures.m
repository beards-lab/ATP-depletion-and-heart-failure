% MakeFigures.m — figures for the staircase-prediction / estiff analysis.
clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resdir = fullfile(here,'results');

Lv = load(fullfile(resdir,'levers.mat'));
G2 = load(fullfile(resdir,'grid2d.mat'));
AD = load(fullfile(resdir,'estiff_alldays.mat'));
DS = load(fullfile(resdir,'data_spread.mat'));
Sc = load(fullfile(resdir,'estiff_scan.mat'));
base = Lv.base;

d  = load('data/protocol_03_27_2026_velocitytable_stairs.mat');
dt = d.datatable; vt = d.velocitytable;
tMov = vt(find(vt(:,2)~=0,1),1);

%% ===== FIG 2: the verdict ================================================
f2 = figure('Position',[40 40 1560 900],'Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

% (a) staircase cost vs estiff, all three preparations
nexttile; hold on; box on;
cols = [0 0.35 0.75; 0.85 0.33 0.1; 0.2 0.6 0.2];
for i = 1:numel(AD.AD)
    plot(AD.EV, [AD.AD(i).S.Cfeat], 'o-','LineWidth',1.7,'Color',cols(i,:), ...
        'MarkerFaceColor',cols(i,:),'DisplayName',AD.days{i});
end
xline(base.estiff,'k--','incumbent','LabelVerticalAlignment','bottom','HandleVisibility','off');
set(gca,'YScale','log'); xlabel('estiff'); ylabel('staircase feature cost');
title('(a) every preparation worsens as estiff rises'); legend('Location','northwest');

% (b) block decomposition vs estiff (03/27)
nexttile; hold on; box on;
S2 = Lv.S2; ok = ~cellfun(@isempty,{S2.parts}); P = cat(1,S2(ok).parts); v = [S2(ok).val];
okc = ~cellfun(@isempty,{Sc.R.parts}); Pc = cat(1,Sc.R(okc).parts); vc = [Sc.R(okc).estiff];
vv = [v vc]; PP = [P; Pc]; [vv,ix] = sort(vv); PP = PP(ix);
bl  = {'st_peak','st_nadir','st_rec','st_plat'};
lbl = {'elastic peak','phase-2 nadir','delayed recovery','envelope'};
for j = 1:numel(bl)
    plot(vv, [PP.(bl{j})], 'o-','LineWidth',1.6,'DisplayName',lbl{j});
end
xline(base.estiff,'k--','HandleVisibility','off'); set(gca,'YScale','log');
xlabel('estiff'); ylabel('block cost'); title('(b) which part of the step response pays');
legend('Location','northwest');

% (c) the trade-off: staircase vs the paid-for objective
nexttile; hold on; box on;
ok2 = ~isnan([S2.Et]);
plot([S2(ok2).Cfeat],[S2(ok2).Et],'o-','LineWidth',1.6,'Color',[0 0.35 0.75], ...
     'DisplayName','estiff (compensated)');
for k = find(ok2)
    text(S2(k).Cfeat, S2(k).Et, sprintf(' %.2f',S2(k).val),'FontSize',8,'Color',[0 0.35 0.75]);
end
S3 = Lv.S3; ok3 = ~isnan([S3.Et]);
plot([S3(ok3).Cfeat],[S3(ok3).Et],'s-','LineWidth',1.6,'Color',[0.85 0.33 0.1], ...
     'DisplayName','kSE (series stiffness)');
for k = find(ok3)
    text(S3(k).Cfeat, S3(k).Et, sprintf(' x%.2g',S3(k).val/base.kSE),'FontSize',8,'Color',[0.85 0.33 0.1]);
end
plot(5.415, 0.9786, 'kp','MarkerSize',16,'MarkerFaceColor','y','DisplayName','incumbent');
set(gca,'XScale','log','YScale','log');
xlabel('staircase feature cost'); ylabel('paid-for objective E_t');
title('(c) what the staircase gain costs'); legend('Location','northeast');

% (d)/(e) 2-D grid heatmaps
Cf = arrayfun(@(g) g.Cfeat, G2.G); Et = arrayfun(@(g) g.Et, G2.G);
nexttile; imagesc(log10(Cf)); colorbar; box on;
set(gca,'XTick',1:numel(G2.KV),'XTickLabel',compose('x%.2g',G2.KV), ...
        'YTick',1:numel(G2.EV),'YTickLabel',compose('%.3g',G2.EV));
xlabel('kSE multiplier'); ylabel('estiff'); title('(d) log_{10} staircase cost');
for a = 1:size(Cf,1)
    for b = 1:size(Cf,2)
        text(b,a,sprintf('%.1f',Cf(a,b)),'HorizontalAlignment','center','Color','w','FontWeight','bold');
    end
end

nexttile; imagesc(log10(Et)); colorbar; box on;
set(gca,'XTick',1:numel(G2.KV),'XTickLabel',compose('x%.2g',G2.KV), ...
        'YTick',1:numel(G2.EV),'YTickLabel',compose('%.3g',G2.EV));
xlabel('kSE multiplier'); ylabel('estiff'); title('(e) log_{10} objective E_t');
for a = 1:size(Et,1)
    for b = 1:size(Et,2)
        text(b,a,sprintf('%.1f',Et(a,b)),'HorizontalAlignment','center','Color','w','FontWeight','bold');
    end
end

% (f) model-free amplitude read
nexttile; hold on; box on;
xd = [1 2];
bar(xd, [mean(DS.one) mean(DS.thr)], 0.6, 'FaceColor',[0.82 0.82 0.82],'EdgeColor','k');
errorbar(xd, [mean(DS.one) mean(DS.thr)], [std(DS.one) std(DS.thr)], 'k','LineStyle','none','LineWidth',1.3);
plot(1+0.06*randn(size(DS.one)), DS.one, 'ko','MarkerFaceColor','k','MarkerSize',4);
plot(2+0.06*randn(size(DS.thr)), DS.thr, 'ko','MarkerFaceColor','k','MarkerSize',4);
plot(xd, [mean(DS.one) mean(DS.one)*3], 'r--o','LineWidth',1.6);
set(gca,'XTick',xd,'XTickLabel',{'1% step','3% step'}); ylabel('elastic peak (F_0)');
title(sprintf('(f) measured ratio %.2f, Hookean would give 3.00', DS.R));
legend({'data mean','between-prep SD','individual steps','','linear-spring expectation'}, ...
       'Location','northwest');

exportgraphics(f2, fullfile(resdir,'estiff_verdict.png'),'Resolution',140);

%% ===== FIG 3: traces =====================================================
variants = { base.estiff, 1.0, 'incumbent  e=0.965'; ...
             1.50,        1.0, 'e=1.50 (the hypothesis)'; ...
             0.75,        1.0, 'e=0.75'; ...
             0.85,        1.8, 'e=0.85, kSE x1.8 (best)' };
ref = Lv.ref;
f3 = figure('Position',[40 40 1560 820],'Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile([1 2]); hold on; box on;
plot(dt(:,1), dt(:,3),'-','Color',[.82 .82 .82],'HandleVisibility','off');
plot(dt(:,1), sgolayfilt(dt(:,3),2,101),'k-','LineWidth',1.8,'DisplayName','data');
cl = lines(4); TR = cell(size(variants,1),1);
for i = 1:size(variants,1)
    p = base; p.estiff = variants{i,1};
    sc = (ref.Fiso-ref.Fp1)/(base.kstiff2*ref.pw(variants{i,1}));
    p.kstiff2 = base.kstiff2*sc; p.kstiff2_n = base.kstiff2_n*sc;
    p.kSE = base.kSE*variants{i,2};
    p.RunStairs = 1; p.PlotEachSeparately = false; p.ghostSave=''; p.ghostLoad='';
    [~, os] = runStairsExperiment(p);
    Fm = interp1(os.t, os.Force, tMov,'linear','extrap');
    TR{i} = [os.t(:), os.Force(:)/Fm];
    plot(TR{i}(:,1), TR{i}(:,2),'-','LineWidth',1.5,'Color',cl(i,:),'DisplayName',variants{i,3});
end
xlim([dt(1,1) dt(end,1)]); ylabel('F / F_0'); xlabel('t (s)');
legend('Location','northwest'); title('staircase, whole record');

nexttile; hold on; box on;
plot(dt(:,1), sgolayfilt(dt(:,3),2,101),'k-','LineWidth',1.8);
for i = 1:numel(TR); plot(TR{i}(:,1),TR{i}(:,2),'-','LineWidth',1.5,'Color',cl(i,:)); end
xlim([72.45 72.66]); xlabel('t (s)'); ylabel('F / F_0');
title('first two 1 % steps: peak, undershoot, recovery');

nexttile; hold on; box on;
plot(dt(:,1), sgolayfilt(dt(:,3),2,101),'k-','LineWidth',1.8);
for i = 1:numel(TR); plot(TR{i}(:,1),TR{i}(:,2),'-','LineWidth',1.5,'Color',cl(i,:)); end
xlim([73.14 73.40]); xlabel('t (s)'); ylabel('F / F_0');
title('the closing 3 % step');

exportgraphics(f3, fullfile(resdir,'traces.png'),'Resolution',140);
save(fullfile(resdir,'traces.mat'),'TR','variants');
fprintf('figures -> %s\n', resdir);
