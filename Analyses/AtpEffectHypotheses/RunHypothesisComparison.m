% RunHypothesisComparison.m -- H1 vs H2 vs H3 side by side, against the 2 mM data.
%
% The pathway ladder (H0..H4) was reported as costs. This shows what those costs
% look like as CURVES: three hypotheses, each refitted, each drawn against the
% same 3-preparation 2 mM envelope, so "is the extra freedom buying anything you
% can see?" gets a visual answer rather than only a numeric one.
%
%   H1  ATP chemistry only          k2, kah, kamh                      3 gates
%   H2  + SRX nucleotide states     + kmsrd, ksr0, ksr2srd, ksrd2sr     7 gates
%   H3  + mechanical                + eta_M, kSE, mu_neg               10 gates
%
% H0 (no ATP effect) and H4 (all 21) appear in the cost ladder panel but not as
% curves: H0 is off-scale by design, and H4 never converged (one sweep of its
% budget), so drawing it as a fitted model would overstate it.
%
% Parameters come from the ladder's own saved state, so these ARE the fitted
% rungs, not re-derived approximations.
%
% Run: matlab -batch "run('Analyses/AtpEffectHypotheses/RunHypothesisComparison.m')"

clear; clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if isempty(gcp('nocreate')); parpool('local',5); end
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end
figDir = fullfile(resDir,'figs_hyp'); if ~isfolder(figDir); mkdir(figDir); end

LAD = fullfile(root,'Analyses','LowATP_Realign4Transfer','results','atp_ladder2.mat');
L2  = load(LAD);
L1  = load(fullfile(root,'Analyses','LowATP_Realign4Transfer','results','atp_ladder.mat'));

% cost ladder, for the summary panel (2 mM objective incl. force-velocity)
LADDER = { 'H0  none',        0, L1.res(1).c
           'H1  chemistry',   3, L2.OUT(1).c
           'H2  + SRX',       7, L2.OUT(2).c
           'H3  + mechanical',10, L2.OUT(3).c
           'H4  all 21',      21, L2.OUT(4).c };

HY = { 'H1', 'chemistry only (3)',   L2.OUT(1).p, [0.85 0.55 0.10]
       'H2', '+ SRX states (7)',     L2.OUT(2).p, [0.80 0.20 0.10]
       'H3', '+ mechanical (10)',    L2.OUT(3).p, [0.45 0.15 0.55] };
nH = size(HY,1);

CB = [0.55 0.55 0.60]; CR = [0.10 0.10 0.12];

%% ================= 1. 2 mM data envelope ==================================
DS = { '03/27 M','protocol_03_27_2026_2mM_slack.mat'
       '04/03 F','protocol_04_03_2026_2mM_slack.mat'
       '04/10 M','protocol_04_10_2026_2mM_slack.mat' };
IREF = 1; nP = size(DS,1);
P = struct();
for i = 1:nP
    S = load(['data/' DS{i,2}]);
    P(i).ft = extractSlackAttributes(S.datatable(:,1), S.datatable(:,3), ...
                                     S.datatable(:,2), S.velocitytable, [], [], false, true);
    P(i).sc = mean(double(P(i).ft.steady),'omitnan');
end
nC = min(arrayfun(@(x) numel(x.ft.steady), P));
FEAT = { 'steady','S'; 'A','F'; 'peak1_y','F'; 'vall_y','F'; 'peak2','F'
         'vall2_dy','F'; 'ktr','R'; 'rsK','R'; 'restretchSlopeStart','R'
         'peak1_dSL','T'; 'vall_t','T'; 'rsR2','Q' };
ENV = struct();
for j = 1:size(FEAT,1)
    nm = FEAT{j,1}; V = nan(nP,nC); ok = true;
    for i = 1:nP
        if ~isfield(P(i).ft,nm); ok=false; break; end
        v = double(P(i).ft.(nm)(:))';
        if numel(v)<nC; ok=false; break; end
        v = v(1:nC); if FEAT{j,2}=='F'; v = v/P(i).sc; end
        V(i,:) = v;
    end
    if ~ok; continue; end
    ENV.(nm) = struct('cls',FEAT{j,2},'raw',V,'mean',mean(V,1,'omitnan'), ...
                      'sd',std(V,0,1,'omitnan'));
end

%% ================= 2. run each hypothesis at 2 mM =========================
CACHE = fullfile(resDir,'hyp_runs.mat');
if isfile(CACHE) && ~(exist('FORCE_RERUN','var') && FORCE_RERUN)
    fprintf('loading cached hypothesis runs from %s\n', CACHE);
    Q = load(CACHE); R = Q.R;
else
    R = struct('tag',{},'lbl',{},'col',{},'fm',{},'fvm',{},'fvd',{},'sc',{});
    for k = 1:nH
        fprintf('\n== %s: %s ==\n', HY{k,1}, HY{k,2});
        p = iCfg(getParams(HY{k,3}, [], true, false));
        ps = p; ps.RunSlack = true; ps.RunForceVelocity = false; ps.RunSlackSegments = 'AllPar';
        [~,~,fm,~] = runSlackExperiment(ps);
        pf = p; pf.RunSlack = false; pf.RunForceVelocity = true;
        [~,~,fvm,fvd] = runFVExperiment(pf, ATP_c, Data_ATP);
        R(k) = struct('tag',HY{k,1},'lbl',HY{k,2},'col',HY{k,4},'fm',fm, ...
                      'fvm',fvm,'fvd',fvd,'sc',mean(double(fm.steady),'omitnan'));
    end
    save(CACHE,'R'); fprintf('cached -> %s\n', CACHE);
end
close all force;

%% ================= 3. cost ladder =========================================
fig = figure('Units','pixels','Position',[20 20 1250 620],'Color','w');
tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on; grid on
dof = cell2mat(LADDER(:,2)); cst = cell2mat(LADDER(:,3));
plot(dof, cst,'o-','Color',[0.2 0.2 0.2],'LineWidth',2.0,'MarkerSize',11, ...
     'MarkerFaceColor','w');
% H2/H3 sit 3 DOF apart on a flat part of the curve: stagger so they don't collide
dy = [1.10 1.10 1.14 0.84 1.10];
for k = 1:numel(dof)
    text(dof(k)+0.7, cst(k)*dy(k), strtrim(LADDER{k,1}),'FontSize',11,'FontWeight','bold');
end
set(gca,'YScale','log'); xlim([-1.5 26]); ylim([min(cst)*0.75 max(cst)*1.6]);
xlabel('free ATP-gate parameters (DOF)','FontSize',13);
ylabel('2 mM cost (log scale)','FontSize',13);
title('The whole ladder','FontSize',14,'FontWeight','bold'); set(gca,'FontSize',12);

nexttile; hold on; box on; grid on
d2 = dof(2:end); c2 = cst(2:end);
b = bar(1:numel(d2), c2, 0.6,'FaceColor','flat');
for k = 1:numel(d2)
    if k <= nH; b.CData(k,:) = HY{k,4}; else; b.CData(k,:) = [0.6 0.6 0.6]; end
    text(k, c2(k)+0.35, sprintf('%.2f',c2(k)),'HorizontalAlignment','center', ...
         'FontSize',12,'FontWeight','bold');
end
% single-line labels: a literal newline here splits into EXTRA tick labels
xticks(1:numel(d2));
xticklabels(cellfun(@(a,b) sprintf('%s (%d DOF)', strtrim(a), b), ...
    LADDER(2:end,1), num2cell(d2), 'uni', 0));
xtickangle(12);
ylabel('2 mM cost','FontSize',13); ylim([0 max(c2)*1.18]);
title('Refitted rungs (H4 unconverged)','FontSize',14,'FontWeight','bold');
set(gca,'FontSize',11);
title(tl,'ATP pathway ladder: 3 chemically-defensible gates capture 97.8 % of the effect', ...
      'FontSize',15,'FontWeight','bold');
exportgraphics(fig, fullfile(figDir,'fig_hyp_cost.png'),'Resolution',160);
fprintf('saved fig_hyp_cost.png\n');

%% ================= 4. force-velocity / power, H1 vs H2 vs H3 ==============
fig = figure('Units','pixels','Position',[20 20 1680 700],'Color','w');
tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile; hold on; box on; grid on
plot(R(1).fvd.FV_v, R(1).fvd.FV_fnorm,'o','Color',CR,'MarkerSize',11,'LineWidth',2.0, ...
     'MarkerFaceColor','w','DisplayName','data 2 mM');
for k = 1:nH
    plot(R(k).fvm.FV_v, R(k).fvm.FV_fnorm,'-','Color',R(k).col,'LineWidth',2.6, ...
         'DisplayName',sprintf('%s %s',R(k).tag,R(k).lbl));
end
xlabel('shortening velocity (ML/s)','FontSize',13); ylabel('F / F_{iso}','FontSize',13);
title('Force-velocity, 2 mM','FontSize',15,'FontWeight','bold');
legend('Location','northeast','FontSize',12); set(gca,'FontSize',12);

nexttile; hold on; box on; grid on
errorbar(R(1).fvd.FV_powerv, R(1).fvd.FV_normpowerAvg, ...
         sqrt(max(R(1).fvd.FV_normpowerVar,0)),'o','Color',CR,'MarkerSize',11, ...
         'LineWidth',2.0,'MarkerFaceColor','w','CapSize',9,'DisplayName','data 2 mM');
for k = 1:nH
    plot(R(k).fvm.FV_powerv, R(k).fvm.FV_normpower,'-','Color',R(k).col,'LineWidth',2.6, ...
         'DisplayName',sprintf('%s %s',R(k).tag,R(k).lbl));
end
xlabel('shortening velocity (ML/s)','FontSize',13);
ylabel('normalised power (F/F_{iso})\timesv','FontSize',13);
title('Power-velocity, 2 mM','FontSize',15,'FontWeight','bold');
legend('Location','northeast','FontSize',12); set(gca,'FontSize',12);
title(tl,'H1 vs H2 vs H3 at 2 mM: what the extra parameters buy on force-velocity', ...
      'FontSize',15,'FontWeight','bold');
exportgraphics(fig, fullfile(figDir,'fig_hyp_fv.png'),'Resolution',160);
fprintf('saved fig_hyp_fv.png\n');

%% ================= 5. features, H1 vs H2 vs H3 ============================
FL = { 'steady','steady (kPa)'; 'A','A / steady'; 'peak1_y','peak1 / steady'
       'vall_y','valley / steady'; 'peak2','peak2 / steady'; 'vall2_dy','vall2\_dy / steady'
       'ktr','k_{tr} (s^{-1})'; 'rsK','rsK (s^{-1})'
       'restretchSlopeStart','restretch slope'; 'peak1_dSL','peak1\_dSL (\mum)'
       'vall_t','valley time (s)'; 'rsR2','rsR2' };
fig = figure('Units','pixels','Position',[20 20 1900 1120],'Color','w');
tl = tiledlayout(3,4,'Padding','compact','TileSpacing','compact');
Z = nan(nH, size(FL,1));
for j = 1:size(FL,1)
    nexttile; hold on; box on; grid on
    nm = FL{j,1};
    if ~isfield(ENV,nm); axis off; continue; end
    E = ENV.(nm); x = 1:nC;
    errorbar(x, E.mean, E.sd,'o-','Color',CB*0.7,'MarkerSize',8,'LineWidth',1.7, ...
             'MarkerFaceColor','w','CapSize',7);
    plot(x, E.raw(IREF,:),'s-','Color',CR,'MarkerSize',6,'LineWidth',1.1);
    sub = strings(1,0);
    for k = 1:nH
        if ~isfield(R(k).fm,nm); continue; end
        mv = double(R(k).fm.(nm)(:))'; mv = mv(1:min(nC,numel(mv)));
        if E.cls=='F'; mv = mv / R(k).sc; end
        plot(1:numel(mv), mv,'x--','Color',R(k).col,'MarkerSize',11,'LineWidth',2.0);
        Z(k,j) = mean((mv - E.mean(1:numel(mv)))./max(E.sd(1:numel(mv)),eps),'omitnan');
        sub(end+1) = sprintf('%s %+.1f', R(k).tag, Z(k,j)); %#ok<SAGROW>
    end
    xlim([0.6 nC+0.4]); xticks(x);
    title({FL{j,2}, ['SD: ' char(strjoin(sub,'  '))]}, ...
          'FontSize',10,'FontWeight','bold','Interpreter','tex');
    if j>8; xlabel('slack cycle','FontSize',11); end
    set(gca,'FontSize',10);
end
lg = legend([{'3-prep mean \pm SD','03/27 M (fit target)'}, ...
             arrayfun(@(k) sprintf('%s %s',R(k).tag,R(k).lbl), 1:nH,'uni',0)], ...
            'Orientation','horizontal','FontSize',11);
lg.Layout.Tile = 'north';
title(tl,'H1 vs H2 vs H3 at 2 mM: every scored feature against the 3-preparation envelope', ...
      'FontSize',15,'FontWeight','bold');
exportgraphics(fig, fullfile(figDir,'fig_hyp_features.png'),'Resolution',150);
fprintf('saved fig_hyp_features.png\n');

%% ================= 6. where the rungs differ, in SD =======================
fprintf('\n%-22s %8s %8s %8s\n','feature','H1','H2','H3');
for j = 1:size(FL,1)
    if all(isnan(Z(:,j))); continue; end
    fprintf('%-22s %8.1f %8.1f %8.1f\n', FL{j,1}, Z(1,j), Z(2,j), Z(3,j));
end
fprintf('\nmean |SD| over features:  H1 %.2f   H2 %.2f   H3 %.2f\n', ...
        mean(abs(Z(1,:)),'omitnan'), mean(abs(Z(2,:)),'omitnan'), mean(abs(Z(3,:)),'omitnan'));

save(fullfile(resDir,'hyp_comparison.mat'),'R','ENV','Z','FL','LADDER');
fprintf('\nHYPOTHESIS COMPARISON DONE -> %s\n', figDir);

%% ==========================================================================
function p = iCfg(p)
p.UseAtpJointFit=false; p.EvalFeatures=true; p.BreakOnODEUnstable=false;
p.OptimizeOn='Feats'; p.MaxRunTime=600;
p.PlotEachSeparately=0; p.PlotFeatureFitting=0; p.PlotFullscreen=0;
p.RunKtr=false; p.RunStairs=false; p.RunForceVelocityTime=false;
p.RunSlackPassive=false;
p.MgATP=2; p.MgADP=0; p.Pi=0.8;
p.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
end
