% RunLowAtpFitReport.m -- the low-ATP fit against the 3-preparation 2 mM envelope.
%
% Same shape as Analyses/OptimResultReport in the Further-optim worktree, but at
% 2 mM: the model is judged against "the data" (mean +- SD over the three 2026
% protocol days) rather than against the one preparation it was fitted to. The
% fit target, 03/27, is drawn separately in every panel, because a deviation
% from the 3-prep mean is only a defect if it is ALSO a deviation from 03/27.
%
% Absolute force is not transferable between preparations (mean `steady` differs
% by ~20 %), so every trace and every force-dimension feature is divided by ITS
% OWN mean steady before averaging. The band is then the spread of the SHAPE.
% The model is normalised the same way, so both sides are scale-free; absolute
% level is judged separately, in its own panel.
%
% Baseline: the SAME parametrisation as the H0-H4 hypothesis ladder
% (params/lowATP_H2_chem_opt.m = the ladder's core + the H2 gate set).
%
% Panels: slack (whole protocol), slack detail, force-velocity, power-velocity,
%         staircase, and every scored feature.
%
% Run: matlab -batch "run('Analyses/AtpEffectHypotheses/RunLowAtpFitReport.m')"

clear; clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if isempty(gcp('nocreate')); parpool('local',5); end
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end
figDir = fullfile(resDir,'figs_lowatp'); if ~isfolder(figDir); mkdir(figDir); end

if ~exist('SNAP','var'); SNAP = 'lowATP_H2_chem_opt'; end
ATPMM = 2;

DS = { '03/27 M','protocol_03_27_2026_2mM_slack.mat'
       '04/03 F','protocol_04_03_2026_2mM_slack.mat'
       '04/10 M','protocol_04_10_2026_2mM_slack.mat' };
IREF = 1;                              % the preparation the model was fitted to
nP = size(DS,1);

C2  = [0.80 0.20 0.10];   % model, 2 mM
CB  = [0.55 0.55 0.60];   % envelope band
CR  = [0.10 0.10 0.12];   % 03/27 reference trace

TPRE = 0.05; TPOST = 0.545; DT = 0.00025;
tg = (-TPRE:DT:TPOST)';

%% ================= 1. the 2 mM data envelope ==============================
P = struct();
for i = 1:nP
    S = load(['data/' DS{i,2}]);
    P(i).t  = S.datatable(:,1); P(i).ML = S.datatable(:,2); P(i).F = S.datatable(:,3);
    P(i).vt = S.velocitytable;
    P(i).rel = S.velocitytable(S.velocitytable(:,2) < -1, 1);
    P(i).ft = extractSlackAttributes(P(i).t, P(i).F, P(i).ML, P(i).vt, [], [], false, true);
    P(i).sc = mean(double(P(i).ft.steady),'omitnan');
    fprintf('%-9s  steady = %5.1f kPa | %d releases\n', DS{i,1}, P(i).sc, numel(P(i).rel));
end
nC = min(arrayfun(@(x) numel(x.rel), P));

FEAT = { 'steady','S'; 'A','F'; 'peak1_y','F'; 'vall_y','F'; 'peak2','F'
         'vall2_dy','F'; 'ktr','R'; 'rsK','R'; 'restretchSlopeStart','R'
         'peak1_dSL','T'; 'vall_t','T'; 'rsR2','Q' };
ENV = struct('prep',{DS(:,1)'},'scale',[P.sc],'iRef',IREF);
for j = 1:size(FEAT,1)
    nm = FEAT{j,1}; V = nan(nP,nC); ok = true;
    for i = 1:nP
        if ~isfield(P(i).ft,nm); ok=false; break; end
        v = double(P(i).ft.(nm)(:))';
        if numel(v) < nC; ok=false; break; end
        v = v(1:nC);
        if FEAT{j,2}=='F'; v = v / P(i).sc; end
        V(i,:) = v;
    end
    if ~ok; continue; end
    ENV.feat.(nm) = struct('cls',FEAT{j,2},'raw',V, ...
        'mean',mean(V,1,'omitnan'),'sd',std(V,0,1,'omitnan'));
end

% timecourse, per cycle and whole-protocol, aligned on release 1
Fc = nan(nP,numel(tg),nC);
for c = 1:nC
    for i = 1:nP
        tt = P(i).t - P(i).rel(c); m = tt>=tg(1)-2*DT & tt<=tg(end)+2*DT;
        Fc(i,:,c) = interp1(tt(m), P(i).F(m)/P(i).sc, tg, 'linear', NaN);
    end
end
tA = (-TPRE:DT:(nC-1)*0.600+TPOST)';
FA = nan(nP,numel(tA)); LA = nan(nP,numel(tA));
for i = 1:nP
    tt = P(i).t - P(i).rel(1); m = tt>=tA(1)-2*DT & tt<=tA(end)+2*DT;
    FA(i,:) = interp1(tt(m), P(i).F(m)/P(i).sc, tA, 'linear', NaN);
    LA(i,:) = interp1(tt(m), P(i).ML(m),        tA, 'linear', NaN);
end

%% ================= 2. the model ===========================================
p = iCfg(getParams(loadParams(SNAP),[],true,false), ATPMM);

fprintf('\n== slack ==\n');
ps = p; ps.RunSlack = true; ps.RunSlackSegments = 'All';
[~, o2, fm2, fd2] = runSlackExperiment(ps);

fprintf('== force-velocity ==\n');
pf = p; pf.RunSlack = false; pf.RunForceVelocity = true;
[~,~, fvm, fvd] = runFVExperiment(pf, ATP_c, Data_ATP);

fprintf('== staircase ==\n');
sto = []; stErr = NaN;
try
    pst = p; pst.RunSlack = false; pst.RunStairs = true;
    pst.stairs_velocitytableonfile = 'protocol_03_27_2026_velocitytable_stairs.mat';
    [stErr, sto] = runStairsExperiment(pst);
    fprintf('  stairs cost %.4f\n', stErr);
catch ME
    fprintf(2,'  stairs FAILED: %s\n', ME.message);
end

mSc = mean(double(fm2.steady),'omitnan');   % model's own force scale

%% ================= 3. figures =============================================
close all force;

% ---- (a) whole slack protocol, normalised, model vs envelope -------------
fig = figure('Units','pixels','Position',[20 20 1750 820],'Color','w');
tl = tiledlayout(3,1,'Padding','compact','TileSpacing','tight');
ax1 = nexttile([2 1]); hold on; box on; grid on
mu = mean(FA,1,'omitnan'); sd = std(FA,0,1,'omitnan'); gd = ~isnan(mu);
% model is aligned on ITS OWN first release, exactly as the data are
tRelM = fd2.t_seg(1);
fill([tA(gd);flipud(tA(gd))],[mu(gd)'-sd(gd)';flipud(mu(gd)'+sd(gd)')], CB, ...
     'FaceAlpha',0.30,'EdgeColor','none','DisplayName','3-prep mean \pm SD');
plot(tA, mu, '-','Color',CB*0.65,'LineWidth',1.4,'DisplayName','3-prep mean');
plot(tA, FA(IREF,:), '-','Color',CR,'LineWidth',1.0,'DisplayName',[DS{IREF,1} ' (fit target)']);
plot(o2.t - tRelM, o2.Force(:)/mSc, '-','Color',C2,'LineWidth',1.9,'DisplayName','model');
xlim([tA(1) tA(end)]); ylabel('F / F_{steady}','FontSize',13);
legend('Location','southeast','FontSize',11); set(ax1,'FontSize',11);
title(sprintf('2 mM slack protocol, scale-free (each trace / its own steady). Model steady = %.1f kPa; data %.1f\\pm%.1f kPa', ...
      mSc, mean([P.sc]), std([P.sc])),'FontSize',13,'FontWeight','bold');
ax2 = nexttile; hold on; box on; grid on
plot(tA, LA(IREF,:),'-','Color',CR,'LineWidth',1.2);
plot(o2.t - tRelM, o2.SL(:), '-','Color',C2,'LineWidth',1.4);
ylabel('SL (\mum)','FontSize',12); xlabel('time from first release (s)','FontSize',13);
xlim([tA(1) tA(end)]); set(ax2,'FontSize',11); linkaxes([ax1 ax2],'x');
exportgraphics(fig, fullfile(figDir,'fig_slack_full.png'),'Resolution',150);
fprintf('saved fig_slack_full.png\n');

% ---- (b) slack detail: every cycle, model vs envelope --------------------
fig = figure('Units','pixels','Position',[20 20 1850 900],'Color','w');
tl = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
for c = 1:nC
    nexttile; hold on; box on; grid on
    m = squeeze(mean(Fc(:,:,c),1,'omitnan')); s = squeeze(std(Fc(:,:,c),0,1,'omitnan'));
    g = ~isnan(m);
    fill([tg(g);flipud(tg(g))],[m(g)'-s(g)';flipud(m(g)'+s(g)')],CB, ...
         'FaceAlpha',0.30,'EdgeColor','none');
    plot(tg, m,'-','Color',CB*0.65,'LineWidth',1.3);
    plot(tg, squeeze(Fc(IREF,:,c)),'-','Color',CR,'LineWidth',1.0);
    tRel = fd2.t_seg(min(c,numel(fd2.t_seg)));
    tm = o2.t - tRel; k = tm>=tg(1) & tm<=tg(end);
    plot(tm(k), o2.Force(k)/mSc,'-','Color',C2,'LineWidth',1.9);
    xlim([tg(1) tg(end)]); title(sprintf('cycle %d', c),'FontSize',12,'FontWeight','bold');
    if c>3; xlabel('t from release (s)','FontSize',11); end
    if mod(c,3)==1; ylabel('F / F_{steady}','FontSize',11); end
    set(gca,'FontSize',10);
end
nexttile; axis off
text(0,0.75,{'\bf2 mM slack detail','','grey band = 3-prep mean \pm SD', ...
     'black = 03/27 M (the fit target)','red = model', ...
     '','all traces scaled by their own', 'steady force -- shape only'}, ...
     'FontSize',12,'VerticalAlignment','top','Interpreter','tex');
title(tl,'2 mM: each slack cycle, model against the 3-preparation envelope', ...
      'FontSize',15,'FontWeight','bold');
exportgraphics(fig, fullfile(figDir,'fig_slack_detail.png'),'Resolution',150);
fprintf('saved fig_slack_detail.png\n');

% ---- (c) force-velocity + power-velocity, 2 mM ---------------------------
fig = figure('Units','pixels','Position',[20 20 1600 660],'Color','w');
tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile; hold on; box on; grid on
plot(fvd.FV_v, fvd.FV_fnorm,'o','Color',CR,'MarkerSize',10,'LineWidth',1.8, ...
     'MarkerFaceColor','w','DisplayName','data 2 mM (Baker2022)');
plot(fvm.FV_v, fvm.FV_fnorm,'-s','Color',C2,'LineWidth',2.5,'MarkerSize',8, ...
     'MarkerFaceColor',C2,'DisplayName','model 2 mM');
xlabel('shortening velocity (ML/s)','FontSize',13); ylabel('F / F_{iso}','FontSize',13);
title('Force-velocity, 2 mM','FontSize',15,'FontWeight','bold');
legend('Location','northeast','FontSize',12); set(gca,'FontSize',12);
nexttile; hold on; box on; grid on
errorbar(fvd.FV_powerv, fvd.FV_normpowerAvg, sqrt(max(fvd.FV_normpowerVar,0)), ...
   'o','Color',CR,'MarkerSize',10,'LineWidth',1.8,'MarkerFaceColor','w','CapSize',10, ...
   'DisplayName','data 2 mM (2 sets, bars = spread)');
plot(fvm.FV_powerv, fvm.FV_normpower,'-s','Color',C2,'LineWidth',2.5,'MarkerSize',8, ...
     'MarkerFaceColor',C2,'DisplayName','model 2 mM');
xlabel('shortening velocity (ML/s)','FontSize',13);
ylabel('normalised power (F/F_{iso})\timesv','FontSize',13);
title('Power-velocity, 2 mM (the scored FV term)','FontSize',15,'FontWeight','bold');
legend('Location','northeast','FontSize',12); set(gca,'FontSize',12);
exportgraphics(fig, fullfile(figDir,'fig_fv.png'),'Resolution',160);
fprintf('saved fig_fv.png\n');

% ---- (d) staircase --------------------------------------------------------
if ~isempty(sto)
    S = load('data/protocol_03_27_2026_velocitytable_stairs.mat');
    fig = figure('Units','pixels','Position',[20 20 1500 760],'Color','w');
    tl = tiledlayout(3,1,'Padding','compact','TileSpacing','tight');
    a1 = nexttile([2 1]); hold on; box on; grid on
    dt = S.datatable;
    plot(dt(:,1), dt(:,3)/max(dt(:,3)),'-','Color',CR,'LineWidth',1.3, ...
         'DisplayName','data 8 mM (no 2 mM staircase exists)');
    plot(sto.t, sto.Force/max(sto.Force),'-','Color',C2,'LineWidth',2.0, ...
         'DisplayName','model 2 mM (prediction)');
    ylabel('F / F_{max}','FontSize',13); set(a1,'FontSize',11);
    legend('Location','southeast','FontSize',11);
    xlim(a1,[min(dt(:,1))-0.1 max(dt(:,1))+0.1]);
    title(sprintf('Staircase (ramp-up). Model cost vs the 8 mM trace = %.3f', stErr), ...
          'FontSize',14,'FontWeight','bold');
    a2 = nexttile; hold on; box on; grid on
    plot(dt(:,1), dt(:,2),'-','Color',CR,'LineWidth',1.2);
    plot(sto.t, sto.SL,'-','Color',C2,'LineWidth',1.4);
    ylabel('SL (\mum)','FontSize',12); xlabel('time (s)','FontSize',13); set(a2,'FontSize',11);
    xlim(a2,[min(dt(:,1))-0.1 max(dt(:,1))+0.1]);
    linkaxes([a1 a2],'x');
    exportgraphics(fig, fullfile(figDir,'fig_staircase.png'),'Resolution',150);
    fprintf('saved fig_staircase.png\n');
end

% ---- (e) features vs the envelope ----------------------------------------
FL = { 'steady','steady (kPa, absolute)'; 'A','A / steady'; 'peak1_y','peak1 / steady'
       'vall_y','valley / steady'; 'peak2','peak2 / steady'; 'vall2_dy','vall2\_dy / steady'
       'ktr','k_{tr} (s^{-1})'; 'rsK','rsK (s^{-1})'
       'restretchSlopeStart','restretch slope'; 'peak1_dSL','peak1\_dSL (\mum)'
       'vall_t','valley time (s)'; 'rsR2','rsR2' };
fig = figure('Units','pixels','Position',[20 20 1900 1120],'Color','w');
tl = tiledlayout(3,4,'Padding','compact','TileSpacing','compact');
Zs = nan(1,size(FL,1));
for j = 1:size(FL,1)
    nexttile; hold on; box on; grid on
    nm = FL{j,1};
    if ~isfield(ENV.feat,nm) || ~isfield(fm2,nm); axis off; continue; end
    E = ENV.feat.(nm); x = 1:nC;
    mv = double(fm2.(nm)(:))'; mv = mv(1:min(nC,numel(mv)));
    if E.cls=='F'; mv = mv / mSc; end
    errorbar(x, E.mean, E.sd,'o-','Color',CB*0.7,'MarkerSize',8,'LineWidth',1.6, ...
             'MarkerFaceColor','w','CapSize',8);
    plot(x, E.raw(IREF,:),'s-','Color',CR,'MarkerSize',7,'LineWidth',1.3);
    plot(1:numel(mv), mv,'x--','Color',C2,'MarkerSize',12,'LineWidth',2.0);
    z = mean((mv - E.mean(1:numel(mv)))./max(E.sd(1:numel(mv)),eps),'omitnan');
    Zs(j) = z;
    xlim([0.6 nC+0.4]); xticks(x);
    title({FL{j,2}, sprintf('%+.1f SD from 3-prep mean', z)}, ...
          'FontSize',11,'FontWeight','bold','Interpreter','tex');
    if j>8; xlabel('slack cycle','FontSize',11); end
    set(gca,'FontSize',10);
end
lg = legend({'3-prep mean \pm SD','03/27 M (fit target)','model'}, ...
            'Orientation','horizontal','FontSize',13);
lg.Layout.Tile = 'north';
title(tl,'2 mM features: model against the 3-preparation envelope', ...
      'FontSize',16,'FontWeight','bold');
exportgraphics(fig, fullfile(figDir,'fig_features.png'),'Resolution',150);
fprintf('saved fig_features.png\n');

save(fullfile(resDir,'lowatp_fitreport.mat'),'ENV','fm2','fd2','fvm','fvd','Zs','FL','mSc','stErr');
fprintf('\nLOW-ATP FIT REPORT DONE -> %s\n', figDir);

%% ==========================================================================
function p = iCfg(p, atp)
p.UseAtpJointFit=false; p.EvalFeatures=true; p.BreakOnODEUnstable=false;
p.OptimizeOn='Feats'; p.MaxRunTime=600;
p.PlotEachSeparately=0; p.PlotFeatureFitting=0; p.PlotFullscreen=0;
p.RunKtr=false; p.RunStairs=false; p.RunForceVelocityTime=false;
p.RunSlackPassive=false; p.RunForceVelocity=false; p.RunSlack=true;
p.MgATP=atp; p.MgADP=0;
if atp==8
    p.Pi=0; p.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
else
    p.Pi=0.8; p.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
end
end
