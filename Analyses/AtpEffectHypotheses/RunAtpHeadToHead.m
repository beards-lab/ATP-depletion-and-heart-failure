% RunAtpHeadToHead.m -- 8 mM vs 2 mM, data and model together, on one axis.
%
% Two figures, each carrying FOUR things so the ATP effect can be read directly
% rather than inferred from two separate reports:
%
%   fig_fv_headtohead.png    force-velocity and power-velocity:
%                            8 mM data, 2 mM data, 8 mM model, 2 mM model
%   fig_feat_headtohead.png  every scored slack feature, same four series,
%                            data as 3-preparation mean +- SD
%   fig_ratio_headtohead.png the ATP effect itself: 2 mM / 8 mM ratio,
%                            model against data, with the data's own spread
%
% Baseline: the SAME parametrisation as the H0-H4 hypothesis ladder. The gates
% are exactly 1 at MgATP_ref, so ONE parameter file gives both conditions --
% the 8 mM and 2 mM curves below differ only by [MgATP].
%
% Run: matlab -batch "run('Analyses/AtpEffectHypotheses/RunAtpHeadToHead.m')"

clear; clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if isempty(gcp('nocreate')); parpool('local',5); end
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end
figDir = fullfile(resDir,'figs_h2h'); if ~isfolder(figDir); mkdir(figDir); end

if ~exist('SNAP','var'); SNAP = 'lowATP_H2_chem_opt'; end

C8 = [0.00 0.35 0.70];   C2 = [0.80 0.20 0.10];
D8 = [0.05 0.20 0.42];   D2 = [0.48 0.10 0.05];

DS = { '03/27 M','protocol_03_27_2026_%s_slack.mat'
       '04/03 F','protocol_04_03_2026_%s_slack.mat'
       '04/10 M','protocol_04_10_2026_%s_slack.mat' };
IREF = 1; nP = size(DS,1);

%% ================= 1. data envelopes, both conditions =====================
FEAT = { 'steady','S'; 'A','F'; 'peak1_y','F'; 'vall_y','F'; 'peak2','F'
         'vall2_dy','F'; 'ktr','R'; 'rsK','R'; 'restretchSlopeStart','R'
         'peak1_dSL','T'; 'vall_t','T'; 'rsR2','Q' };
ENV = struct();
for a = [8 2]
    tag = sprintf('%dmM', a); P = struct();
    for i = 1:nP
        S = load(['data/' sprintf(DS{i,2}, tag)]);
        ft = extractSlackAttributes(S.datatable(:,1), S.datatable(:,3), ...
                                    S.datatable(:,2), S.velocitytable, [], [], false, true);
        P(i).ft = ft; P(i).sc = mean(double(ft.steady),'omitnan');
    end
    nC = min(arrayfun(@(x) numel(x.ft.steady), P));
    E = struct('scale',[P.sc],'nC',nC);
    for j = 1:size(FEAT,1)
        nm = FEAT{j,1}; V = nan(nP,nC); ok = true;
        for i = 1:nP
            if ~isfield(P(i).ft,nm); ok=false; break; end
            v = double(P(i).ft.(nm)(:))';
            if numel(v)<nC; ok=false; break; end
            v = v(1:nC);
            if FEAT{j,2}=='F'; v = v/P(i).sc; end
            V(i,:) = v;
        end
        if ~ok; continue; end
        E.feat.(nm) = struct('cls',FEAT{j,2},'raw',V, ...
            'mean',mean(V,1,'omitnan'),'sd',std(V,0,1,'omitnan'));
    end
    ENV.(sprintf('a%d',a)) = E;
    fprintf('%s envelope: %d cycles, steady %.1f +- %.1f kPa\n', tag, nC, mean([P.sc]), std([P.sc]));
end

%% ================= 2. model, both conditions ==============================
M = struct();
for a = [8 2]
    p = iCfg(getParams(loadParams(SNAP),[],true,false), a);
    fprintf('\n== model %d mM: slack ==\n', a);
    ps = p; ps.RunSlack = true; ps.RunSlackSegments = 'AllPar';
    [~,~,fm,~] = runSlackExperiment(ps);
    fprintf('== model %d mM: force-velocity ==\n', a);
    pf = p; pf.RunSlack = false; pf.RunForceVelocity = true;
    [~,~,fvm,fvd] = runFVExperiment(pf, ATP_c, Data_ATP);
    M.(sprintf('a%d',a)) = struct('fm',fm,'fvm',fvm,'fvd',fvd, ...
                                  'sc',mean(double(fm.steady),'omitnan'));
end
close all force;

%% ================= 3. FV + power, head to head ============================
fig = figure('Units','pixels','Position',[20 20 1680 700],'Color','w');
tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on; grid on
plot(M.a8.fvd.FV_v, M.a8.fvd.FV_fnorm,'o','Color',D8,'MarkerSize',11,'LineWidth',2.0, ...
     'MarkerFaceColor','w','DisplayName','data 8 mM');
plot(M.a2.fvd.FV_v, M.a2.fvd.FV_fnorm,'o','Color',D2,'MarkerSize',11,'LineWidth',2.0, ...
     'MarkerFaceColor','w','DisplayName','data 2 mM');
plot(M.a8.fvm.FV_v, M.a8.fvm.FV_fnorm,'-','Color',C8,'LineWidth',2.8,'DisplayName','model 8 mM');
plot(M.a2.fvm.FV_v, M.a2.fvm.FV_fnorm,'-','Color',C2,'LineWidth',2.8,'DisplayName','model 2 mM');
xlabel('shortening velocity (ML/s)','FontSize',13); ylabel('F / F_{iso}','FontSize',13);
title('Force-velocity, normalised','FontSize',15,'FontWeight','bold');
legend('Location','northeast','FontSize',12); set(gca,'FontSize',12);

nexttile; hold on; box on; grid on
for a = [8 2]
    k = sprintf('a%d',a); d = M.(k).fvd;
    col = D8; if a==2; col = D2; end
    errorbar(d.FV_powerv, d.FV_normpowerAvg, sqrt(max(d.FV_normpowerVar,0)),'o', ...
        'Color',col,'MarkerSize',11,'LineWidth',2.0,'MarkerFaceColor','w','CapSize',9, ...
        'DisplayName',sprintf('data %d mM',a));
end
plot(M.a8.fvm.FV_powerv, M.a8.fvm.FV_normpower,'-','Color',C8,'LineWidth',2.8,'DisplayName','model 8 mM');
plot(M.a2.fvm.FV_powerv, M.a2.fvm.FV_normpower,'-','Color',C2,'LineWidth',2.8,'DisplayName','model 2 mM');
xlabel('shortening velocity (ML/s)','FontSize',13);
ylabel('normalised power (F/F_{iso})\timesv','FontSize',13);
title('Power-velocity','FontSize',15,'FontWeight','bold');
legend('Location','northeast','FontSize',12); set(gca,'FontSize',12);

title(tl,'8 mM vs 2 mM: force-velocity and power, data and model together', ...
      'FontSize',16,'FontWeight','bold');
exportgraphics(fig, fullfile(figDir,'fig_fv_headtohead.png'),'Resolution',160);
fprintf('saved fig_fv_headtohead.png\n');

%% ================= 4. features, head to head ==============================
FL = { 'steady','steady (kPa)'; 'A','A / steady'; 'peak1_y','peak1 / steady'
       'vall_y','valley / steady'; 'peak2','peak2 / steady'; 'vall2_dy','vall2\_dy / steady'
       'ktr','k_{tr} (s^{-1})'; 'rsK','rsK (s^{-1})'
       'restretchSlopeStart','restretch slope'; 'peak1_dSL','peak1\_dSL (\mum)'
       'vall_t','valley time (s)'; 'rsR2','rsR2' };
fig = figure('Units','pixels','Position',[20 20 1900 1120],'Color','w');
tl = tiledlayout(3,4,'Padding','compact','TileSpacing','compact');
for j = 1:size(FL,1)
    nexttile; hold on; box on; grid on
    nm = FL{j,1}; drew = false;
    for a = [8 2]
        k = sprintf('a%d',a); E = ENV.(k);
        if ~isfield(E.feat,nm) || ~isfield(M.(k).fm,nm); continue; end
        F = E.feat.(nm); x = 1:E.nC;
        cd_ = D8; cm = C8; if a==2; cd_ = D2; cm = C2; end
        errorbar(x, F.mean, F.sd,'o-','Color',cd_,'MarkerSize',8,'LineWidth',1.7, ...
                 'MarkerFaceColor','w','CapSize',7);
        mv = double(M.(k).fm.(nm)(:))'; mv = mv(1:min(E.nC,numel(mv)));
        if F.cls=='F'; mv = mv / M.(k).sc; end
        plot(1:numel(mv), mv,'x--','Color',cm,'MarkerSize',12,'LineWidth',2.1);
        drew = true;
    end
    if ~drew; axis off; continue; end
    % the ATP effect, model vs data, as a ratio
    t = FL{j,2};
    try
        r_d = mean(ENV.a2.feat.(nm).mean)/mean(ENV.a8.feat.(nm).mean);
        m2 = double(M.a2.fm.(nm)(:))'; m8 = double(M.a8.fm.(nm)(:))';
        if ENV.a2.feat.(nm).cls=='F'; m2 = m2/M.a2.sc; m8 = m8/M.a8.sc; end
        r_m = mean(m2,'omitnan')/mean(m8,'omitnan');
        t = {t, sprintf('2/8: data %.2f | model %.2f', r_d, r_m)};
    catch
    end
    title(t,'FontSize',11,'FontWeight','bold','Interpreter','tex');
    if j>8; xlabel('slack cycle','FontSize',11); end
    set(gca,'FontSize',10);
end
lg = legend({'data 8 mM (3-prep \pm SD)','model 8 mM','data 2 mM (3-prep \pm SD)','model 2 mM'}, ...
            'Orientation','horizontal','FontSize',12);
lg.Layout.Tile = 'north';
title(tl,'8 mM vs 2 mM: every scored slack feature, data and model together', ...
      'FontSize',16,'FontWeight','bold');
exportgraphics(fig, fullfile(figDir,'fig_feat_headtohead.png'),'Resolution',150);
fprintf('saved fig_feat_headtohead.png\n');

%% ================= 5. the ATP effect as a ratio ===========================
NMS = FL(:,1)'; rd = nan(1,numel(NMS)); rm = nan(1,numel(NMS)); rsd = nan(1,numel(NMS));
for j = 1:numel(NMS)
    nm = NMS{j};
    if ~isfield(ENV.a2.feat,nm) || ~isfield(ENV.a8.feat,nm); continue; end
    if ~isfield(M.a2.fm,nm) || ~isfield(M.a8.fm,nm); continue; end
    R2 = ENV.a2.feat.(nm).raw; R8 = ENV.a8.feat.(nm).raw;
    n = min(size(R2,2),size(R8,2));
    perPrep = mean(R2(:,1:n),2,'omitnan') ./ mean(R8(:,1:n),2,'omitnan');
    rd(j) = mean(perPrep,'omitnan'); rsd(j) = std(perPrep,0,'omitnan');
    m2 = double(M.a2.fm.(nm)(:))'; m8 = double(M.a8.fm.(nm)(:))';
    if ENV.a2.feat.(nm).cls=='F'; m2 = m2/M.a2.sc; m8 = m8/M.a8.sc; end
    rm(j) = mean(m2,'omitnan')/mean(m8,'omitnan');
end
ok = ~isnan(rd) & ~isnan(rm);
fig = figure('Units','pixels','Position',[20 20 1300 700],'Color','w');
hold on; box on; grid on
x = 1:nnz(ok);
errorbar(x, rd(ok), rsd(ok),'o','Color',[0.15 0.15 0.15],'MarkerSize',11,'LineWidth',2.0, ...
         'MarkerFaceColor','w','CapSize',9,'DisplayName','data (3 preps, \pm SD)');
plot(x, rm(ok),'s','Color',C2,'MarkerSize',12,'LineWidth',2.4,'MarkerFaceColor',C2, ...
     'DisplayName','model');
yline(1,'--','no ATP effect','Color',[0.5 0.5 0.5],'LineWidth',1.2,'HandleVisibility','off');
xticks(x); xticklabels(NMS(ok)); xtickangle(35);
ylabel('2 mM / 8 mM','FontSize',13);
title('The ATP effect itself: 2 mM / 8 mM ratio, model vs data','FontSize',16,'FontWeight','bold');
legend('Location','best','FontSize',12); set(gca,'FontSize',12);
exportgraphics(fig, fullfile(figDir,'fig_ratio_headtohead.png'),'Resolution',160);
fprintf('saved fig_ratio_headtohead.png\n');

fprintf('\n%-22s %10s %10s %10s %9s\n','feature','data 2/8','SD','model 2/8','m/d');
for j = find(ok)
    fprintf('%-22s %10.3f %10.3f %10.3f %9.2f\n', NMS{j}, rd(j), rsd(j), rm(j), rm(j)/rd(j));
end

save(fullfile(resDir,'atp_headtohead.mat'),'ENV','M','rd','rm','rsd','NMS','ok');
fprintf('\nHEAD-TO-HEAD DONE -> %s\n', figDir);

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
