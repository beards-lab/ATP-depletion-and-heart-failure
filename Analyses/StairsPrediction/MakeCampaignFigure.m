% MakeCampaignFigure.m — the two figures the estiff story actually needs.
%
%   estiff_wall.png     why a 60-parameter optimiser with estiff IN its pool
%                       never moved it: the wall against the valley floor.
%   estiff_campaign.png matched-budget refits at seven frozen exponents,
%                       objective and held-out staircase side by side.
%
% Pure plotting: reads results/estiff_wall.mat (RunEstiffWall.m) and
% results/campaign_table.mat (CollectCampaign.m). Deliberately does NO
% simulation and no `clear functions`, so it can be re-run cheaply.
%
% Run: matlab -batch "run('Analyses/StairsPrediction/MakeCampaignFigure.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
res = fullfile(here,'results');

%% ============================ 1. the wall ================================
W = load(fullfile(res,'estiff_wall.mat'));
f = figure('Position',[80 80 1100 420],'Color','w');

subplot(1,2,1); hold on; grid on; box on;
e = W.e0 + W.DE;
plot(e, W.raw, 'o-','LineWidth',1.8,'Color',[0.80 0.20 0.10],'MarkerFaceColor','w');
plot(e, W.rid, 's-','LineWidth',1.8,'Color',[0.00 0.35 0.70],'MarkerFaceColor','w');
set(gca,'YScale','log');
yline(W.rid(W.DE==0),'k:','incumbent','LabelHorizontalAlignment','left');
xlabel('estiff'); ylabel('objective E_t  (log)');
title({'The wall and the valley floor','\pm0.05 in estiff about the incumbent'});
legend({'estiff alone (kstiff2 fixed)','estiff + iso-force kstiff2'}, ...
       'Location','north','FontSize',8);

subplot(1,2,2); hold on; grid on; box on;
plot(100*(W.KS-1), W.ks, 'o-','LineWidth',1.8,'Color',[0.35 0.10 0.55],'MarkerFaceColor','w');
set(gca,'YScale','log');
xlabel('kstiff2 error (%)'); ylabel('objective E_t  (log)');
title({'kstiff2 alone','the coordinate estiff must track'});

exportgraphics(f, fullfile(res,'estiff_wall.png'), 'Resolution', 150);
fprintf('wrote estiff_wall.png\n');

%% ========================= 2. the campaign ===============================
C = load(fullfile(res,'campaign_table.mat')); r = C.r;
E = [r.estiff]; Et = [r.optEt]; St = [r.optStair];
Et0 = [r.seedEt]; St0 = [r.seedStair];

g = figure('Position',[80 80 1150 430],'Color','w');

subplot(1,2,1); hold on; grid on; box on;
plot(E, Et0,'o--','LineWidth',1.2,'Color',[0.60 0.60 0.60],'MarkerFaceColor','w');
plot(E, Et ,'o-' ,'LineWidth',2.0,'Color',[0.00 0.35 0.70],'MarkerFaceColor','w');
set(gca,'YScale','log');
yline(0.9786,'k:','incumbent 0.979','LabelHorizontalAlignment','left');
xlabel('estiff (frozen)'); ylabel('objective E_t  (log)');
title('Scored objective, matched 0.5 h budget');
legend({'iso-force seed','after refit'},'Location','north','FontSize',8);

subplot(1,2,2); hold on; grid on; box on;
plot(E, St0,'o--','LineWidth',1.2,'Color',[0.60 0.60 0.60],'MarkerFaceColor','w');
plot(E, St ,'o-' ,'LineWidth',2.0,'Color',[0.80 0.20 0.10],'MarkerFaceColor','w');
% The two 1.5 h points. These are the whole reason budget matters here: at 0.5 h
% the red curve is monotone and says "go down", and that reverses when e = 0.90 is
% given the same budget the e = 1.25 point had (3.42 -> 4.00, while e = 1.25
% reaches 1.595).
plot([0.90 1.25], [4.000 1.595], 'p-','MarkerSize',14,'LineWidth',1.5, ...
     'Color',[0.10 0.50 0.10],'MarkerFaceColor',[0.60 0.90 0.60]);
text(0.90, 4.000, '  0.90 @ 1.5 h','FontSize',8,'Color',[0.10 0.40 0.10]);
text(1.25, 1.595, '  1.25 @ 1.5 h','FontSize',8,'Color',[0.10 0.40 0.10]);
set(gca,'YScale','log');
yline(5.415,'k:','incumbent 5.42','LabelHorizontalAlignment','left');
xlabel('estiff (frozen)'); ylabel('staircase feature cost  (log)');
title('Held-out staircase, 03/27');
legend({'iso-force seed','after 0.5 h refit','after 1.5 h refit'}, ...
       'Location','north','FontSize',8);

exportgraphics(g, fullfile(res,'estiff_campaign.png'), 'Resolution', 150);
fprintf('wrote estiff_campaign.png\n');

%% ========================= 3. the table ==================================
fprintf('\n%8s %10s %10s %11s %11s %12s\n', ...
        'estiff','seed Et','opt Et','seed stair','opt stair','kstiff2');
for i = 1:numel(r)
    fprintf('%8.3f %10.4f %10.4f %11.3f %11.3f %12.5g\n', ...
        r(i).estiff, r(i).seedEt, r(i).optEt, r(i).seedStair, r(i).optStair, r(i).kstiff2);
end
