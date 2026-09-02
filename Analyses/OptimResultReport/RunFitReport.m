% RunFitReport.m — standard report for one optimiser result, against the
% previous best and against the three-preparation data envelope.
%
% Runs BOTH parametrisations through every protocol, including the two that
% were NOT scored during the fit (ktr, staircase) -- those are therefore
% genuine PREDICTIONS, and the most informative panels in the set.
%
%   TAG_NEW   the result being reported      (default spreadw_opt)
%   TAG_BASE  what it is being compared to   (default realign4_opt)
%
% Requires results/data_envelope.mat -- run RunDataEnvelope.m first.
%
% Figures -> results/
%   fit_slack_full     whole slack protocol, force + length
%   fit_slack_detail   one tile per slack cycle, release -> recovery -> restretch
%   fit_features       every scored feature vs the 3-prep envelope (mean +- SD)
%   fit_fv             force-velocity, F/F0
%   fit_fvpower        normalised POWER -- the quantity actually scored (70 %)
%   fit_ktr_pred       ktr protocol           PREDICTION (RunKtr was 0 in the fit)
%   fit_stairs_pred    staircase protocol     PREDICTION (RunStairs was 0)
%   fit_costshift      per-feature cost, base vs new, and chi-square

clearvars -except TAG_NEW TAG_BASE REUSE; clc;   % keep caller-set knobs
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir), mkdir(resDir); end

if ~exist('TAG_NEW','var');  TAG_NEW  = 'spreadw_opt';   end
if ~exist('TAG_BASE','var'); TAG_BASE = 'realign4_opt';  end
NAMES = {TAG_BASE, TAG_NEW};

EV = load(fullfile(resDir,'data_envelope.mat'));  ENV = EV.ENV;
SP = load(fullfile(root,'Analyses','FeatureWeighting','results','feature_spread.mat'));

if isempty(gcp('nocreate')); parpool('local', 4); end

%% ---- run both -----------------------------------------------------------
% Both are scored on ONE common fn -- TAG_NEW's. realign4_opt's own list still
% carries the duplicated rsR2 boundary entry, so the two lists differ in length
% and a per-feature comparison against its native list would not line up. The
% duplicate costs 0, so this does not move realign4's headline number.
% REUSE = true replays the last saved run instead of re-simulating — for
% figure/axis tweaks, which otherwise cost 4 minutes of ODE solving each time.
MATF = fullfile(resDir, sprintf('fit_report_%s.mat', TAG_NEW));
if exist('REUSE','var') && REUSE && isfile(MATF)
    fprintf('REUSE: replaying %s (no model evaluated)\n', MATF);
    S = load(MATF, 'R', 'NAMES'); R = S.R; NAMES = S.NAMES;
else
    FN = getParams(loadParams(TAG_NEW), [], true, false).fn;
    for k = 1:2
        fprintf('=== running %s ===\n', NAMES{k});
        R(k) = iRunAll(NAMES{k}, FN);                               %#ok<SAGROW>
        fprintf('    feature cost = %.4f\n', R(k).cost);
    end
end
fd = R(2).fdata;      % data features are identical for both

%% ---- colours ------------------------------------------------------------
cDat = [0.15 0.18 0.22];  cBas = [0.55 0.60 0.66];  cNew = [0.72 0.14 0.18];
cBand= [0.15 0.18 0.22];
LW = {1.6, 2.0};  CC = {cBas, cNew};
% 03/27 is the preparation the model is actually fitted to, so it gets its own
% marker everywhere: a miss against the 3-prep MEAN is only a real defect if it
% is also a miss against THIS trace.
cRef = [0.00 0.45 0.74];
IREF = ENV.iRef;  REFN = ENV.refName;

%% ---- 1. slack, whole protocol ------------------------------------------
% NB every datatable here uses col2 = SL in MICRONS (slack 1.88-2.20, ktr
% 1.60-2.10, stairs 2.00-2.20); only velocitytable col4 is in ML. Checked.
dtb = R(2).out_slack.datatable;
Sv  = load(['data/' R(2).p.velocitytableonfile]);
vt  = Sv.velocitytable;
rel = vt(vt(:,2) < -1, 1);  nC = numel(rel);
XL  = [rel(1)-0.35, rel(end)+0.95];      % the protocol, not the 84 s warm-up

fig = figure('Units','pixels','Position',[20 20 1750 820],'Color','w');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

nexttile([2 1]); hold on; box on; grid on
plot(dtb(:,1), dtb(:,3), '-', 'Color',cDat, 'LineWidth',1.4, 'DisplayName','data 03/27 8 mM');
for k = 1:2
    plot(R(k).out_slack.t, R(k).out_slack.Force, '-', 'Color',CC{k}, 'LineWidth',LW{k}, ...
         'DisplayName',sprintf('%s (%.3f)', NAMES{k}, R(k).cost));
end
ylabel('force (kPa)'); legend('Location','southwest','FontSize',10,'Interpreter','none');
xlim(XL);
title(sprintf('Slack protocol — %s vs %s   (feature cost %.4f -> %.4f, %+.1f %%)', ...
      TAG_BASE, TAG_NEW, R(1).cost, R(2).cost, 100*(R(2).cost-R(1).cost)/R(1).cost), ...
      'FontWeight','bold','Interpreter','none');

nexttile; hold on; box on; grid on
plot(dtb(:,1), dtb(:,2), '-', 'Color',cDat, 'LineWidth',1.2, 'DisplayName','data');
for k = 1:2
    plot(R(k).out_slack.t, R(k).out_slack.SL, '-', 'Color',CC{k}, 'LineWidth',LW{k});
end
xlabel('time (s)'); ylabel('SL (\mum)'); xlim(XL); ylim([1.85 2.25]);
exportgraphics(fig, fullfile(resDir,'fit_slack_full.png'), 'Resolution',150);

%% ---- 2. slack, per-cycle detail ----------------------------------------
fig = figure('Units','pixels','Position',[20 20 1750 880],'Color','w');
tl = tiledlayout(2, nC, 'Padding','compact','TileSpacing','compact');
for c = 1:nC
    % top row: the whole cycle
    nexttile(c); hold on; box on; grid on
    iterPlotCycle(dtb, R, rel(c), [-0.02 0.40], cDat, CC, LW);
    title(sprintf('cycle %d — restretch %.0f ML/s', c, fd.v_restretch(c)));
    if c == 1; ylabel('force (kPa)'); end
    % bottom row: the re-stretch transient, where rsK lives
    nexttile(nC + c); hold on; box on; grid on
    iterPlotCycle(dtb, R, rel(c), [0.29 0.40], cDat, CC, LW);
    xlabel('t - release (s)');
    if c == 1; ylabel('force (kPa)'); end
end
title(tl, {'Slack detail — top: full cycle (release, redevelopment, re-stretch); bottom: the re-stretch transient', ...
           sprintf('%s (grey) vs %s (red) vs data (black)', TAG_BASE, TAG_NEW)}, ...
      'FontWeight','bold','FontSize',12,'Interpreter','none');
exportgraphics(fig, fullfile(resDir,'fit_slack_detail.png'), 'Resolution',150);

%% ---- 3. features against the 3-prep envelope ---------------------------
FL = ENV.featNames;
FL = FL(~ismember(FL, {'SLslack'}));
fig = figure('Units','pixels','Position',[20 20 1750 950],'Color','w');
tl = tiledlayout(4,5,'Padding','compact','TileSpacing','compact');
for j = 1:numel(FL)
    nm = FL{j}; E = ENV.feat.(nm);
    nexttile; hold on; box on; grid on
    x = 1:size(E.mean,2);
    fill([x fliplr(x)], [E.mean+E.sd fliplr(E.mean-E.sd)], cBand, ...
         'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
    plot(x, E.mean, 'o-','Color',cDat,'MarkerFaceColor',cDat,'LineWidth',1.6,'DisplayName','3-prep mean');
    plot(x, E.raw(IREF,:), 'o-','Color',cRef,'MarkerFaceColor',cRef,'LineWidth',1.3, ...
         'MarkerSize',4,'DisplayName',[REFN ' (fitted)']);
    for k = 1:2
        if ~isfield(R(k).fmodel, nm); continue; end
        v = double(R(k).fmodel.(nm)(:))';
        if numel(v) < numel(x); continue; end
        v = v(1:numel(x));
        if E.cls == 'F'; v = v / mean(double(R(k).fmodel.steady),'omitnan'); end
        plot(x, v, 's--','Color',CC{k},'MarkerFaceColor',CC{k},'LineWidth',LW{k}, ...
             'MarkerSize',5,'DisplayName',NAMES{k});
    end
    set(gca,'XTick',x); xlim([0.5 numel(x)+0.5]);
    title(nm,'Interpreter','none','FontSize',10);
    if j == 1; legend('Location','best','FontSize',7,'Interpreter','none'); end
end
title(tl, {'Model features against the THREE-PREPARATION envelope (mean \pm SD, 8 mM)', ...
           sprintf(['force-dimension features are scale-free: each side divided by its own mean steady.  ' ...
                    'Blue = %s, the preparation the model is fitted to'], REFN)}, ...
      'FontWeight','bold','FontSize',13);
exportgraphics(fig, fullfile(resDir,'fit_features.png'), 'Resolution',150);

%% ---- 4/5. force-velocity ------------------------------------------------
AllV = [0, 0.5, 1, 2, 3, 4, 5, 6];
S1 = [56.40, 51.8120, 37.4459, 17.8025, 11.4430, 6.2643, 3.2759, 2.2120];
S2 = [67.3942, 60.9885, 42.7059, 18.8539, 12.5956, 6.4426, 3.5136, 1.8556];

fig = figure('Units','pixels','Position',[20 20 1450 560],'Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile; hold on; box on; grid on
plot(AllV, S1/S1(1), 'o--','Color',cDat,'MarkerFaceColor',cDat,'LineWidth',1.3,'DisplayName','Baker2022');
plot(AllV, S2/S2(1), '^--','Color',[0.45 0.5 0.55],'MarkerFaceColor',[0.45 0.5 0.55],'LineWidth',1.3,'DisplayName','Isovel.2021');
for k = 1:2
    plot(R(k).fmodel.FV_v, R(k).fmodel.FV_fnorm, '-','Color',CC{k},'LineWidth',LW{k}+0.6,'DisplayName',NAMES{k});
end
xlabel('shortening velocity (ML/s)'); ylabel('F / F_0'); legend('Interpreter','none','FontSize',9);
title('Force-velocity, normalised force','FontWeight','bold');

nexttile; hold on; box on; grid on
% v = 0 is EXCLUDED from the scored feature -- power there is identically 0 by
% construction in every dataset, so it can carry no error and would only dilute
% the inverse-variance weights. It is drawn anyway, as the anchor the curve
% leaves from, with an open marker to flag that it is not a fitted point.
vp = fd.FV_powerv(:)'; mu = fd.FV_normpowerAvg(:)'; sd = sqrt(fd.FV_normpowerVar(:)');
vp0 = [0 vp]; mu0 = [0 mu]; sd0 = [0 sd];
fill([vp0 fliplr(vp0)], [mu0+sd0 fliplr(mu0-sd0)], cBand,'FaceAlpha',0.18,'EdgeColor','none', ...
     'DisplayName','\pm SD between the two sources');
plot(vp0, mu0, '-','Color',cDat,'LineWidth',1.8,'DisplayName','data mean');
plot(vp,  mu,  'o','Color',cDat,'MarkerFaceColor',cDat,'HandleVisibility','off');
for k = 1:2
    plot(vp0, [0 R(k).fmodel.FV_normpowerAvg(:)'], '-','Color',CC{k}, ...
         'LineWidth',LW{k}+0.6,'DisplayName',NAMES{k});
    plot(vp,  R(k).fmodel.FV_normpowerAvg(:)', 's','Color',CC{k},'MarkerFaceColor',CC{k}, ...
         'HandleVisibility','off');
end
plot(0, 0, 'o','Color',cDat,'MarkerFaceColor','w','LineWidth',1.4,'MarkerSize',7, ...
     'DisplayName','v = 0 (not scored: power \equiv 0)');
xlim([-0.15 max(vp)+0.3]); ylim([-0.02 inf]);
xlabel('shortening velocity (ML/s)'); ylabel('normalised power  (F/F_0)\cdotv');
legend('Interpreter','none','FontSize',9,'Location','southeast');
iC = find(strcmp(iTrim(R(2).fnNames),'FV_normpowerAvg'),1);
title(sprintf('Normalised POWER — the scored quantity (70 %% of the objective)\ncost %.4f -> %.4f', ...
      R(1).featCost(iC), R(2).featCost(iC)),'FontWeight','bold');
exportgraphics(fig, fullfile(resDir,'fit_fv.png'), 'Resolution',150);

%% ---- 6/7. predictions: ktr and staircase -------------------------------
for which = {'ktr','stairs'}
    w = which{1};
    if strcmp(w,'ktr')
        Sd = load(['data/' R(2).p.ktr_velocitytableonfile]); fld = 'out_ktr';
        ttl = 'k_{tr} protocol'; fnm = 'fit_ktr_pred.png';
    else
        Sd = load(['data/' R(2).p.stairs_velocitytableonfile]); fld = 'out_stairs';
        ttl = 'Staircase (ramp-up) protocol'; fnm = 'fit_stairs_pred.png';
    end
    dt2 = Sd.datatable; vt2 = Sd.velocitytable;
    fig = figure('Units','pixels','Position',[20 20 1450 700],'Color','w');
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
    nexttile([2 1]); hold on; box on; grid on
    plot(dt2(:,1), dt2(:,3), '-','Color',cDat,'LineWidth',1.4,'DisplayName','data (normalised)');
    for k = 1:2
        o = R(k).(fld);
        if strcmp(w,'ktr'); tref = vt2(2,1); else; tref = vt2(find(vt2(:,2)~=0,1),1); end
        F0 = interp1(o.t, o.Force, tref, 'linear','extrap'); if F0 == 0; F0 = 1; end
        plot(o.t, o.Force/F0, '-','Color',CC{k},'LineWidth',LW{k},'DisplayName',NAMES{k});
    end
    ylabel('F / F_{iso}'); legend('Location','best','FontSize',10,'Interpreter','none');
    xlim([dt2(1,1) dt2(end,1)]);
    title({sprintf('%s — PREDICTION (this protocol was NOT scored during the fit)', ttl), ...
           sprintf('%s vs %s', TAG_BASE, TAG_NEW)}, 'FontWeight','bold','Interpreter','none');
    nexttile; hold on; box on; grid on
    plot(dt2(:,1), dt2(:,2), '-','Color',cDat,'LineWidth',1.2);
    for k = 1:2; plot(R(k).(fld).t, R(k).(fld).SL, '-','Color',CC{k},'LineWidth',LW{k}); end
    xlabel('time (s)'); ylabel('SL (\mum)'); xlim([dt2(1,1) dt2(end,1)]);
    exportgraphics(fig, fullfile(resDir,fnm), 'Resolution',150);
end

%% ---- 8. model ON the three-preparation timecourse envelope --------------
% The envelope is the reference the model should be judged against: a thick
% mean line built from all three preps, with a +-SD band, everything scale-free
% (each trace over its own mean steady). The model is normalised the SAME way,
% by its own steady, so a deviation here is shape error and cannot be bought
% by moving the force level.
tg = ENV.tc.t;
fig = figure('Units','pixels','Position',[20 20 1750 900],'Color','w');
tl = tiledlayout(3, nC, 'Padding','compact','TileSpacing','compact');

nexttile([1 nC]); hold on; box on; grid on
tA = ENV.tcAll.t; mu = ENV.tcAll.Fmean(:); sd = ENV.tcAll.Fsd(:);
fill([tA; flipud(tA)], [mu+sd; flipud(mu-sd)], cBand,'FaceAlpha',0.22,'EdgeColor','none', ...
     'DisplayName','3 preps, \pm SD');
plot(tA, mu, '-','Color',cDat,'LineWidth',3.0,'DisplayName','3-prep mean');
plot(tA, ENV.tcAll.F(IREF,:), '--','Color',cRef,'LineWidth',0.9, ...
     'DisplayName',[REFN ' (fitted)']);
for k = 1:2
    [tm, Fm] = iModelNorm(R(k), rel(1));
    plot(tm, Fm, '-','Color',CC{k},'LineWidth',1.7,'DisplayName',NAMES{k});
end
xlim([tA(1) tA(end)]); ylabel('F / own steady'); xlabel('time from first release (s)');
legend('Location','southeast','FontSize',9,'Interpreter','none');
title('Model against the three-preparation envelope — scale-free on both sides','FontWeight','bold');

for c = 1:nC
    for row = 1:2
        nexttile(nC*row + c); hold on; box on; grid on
        m = ENV.tc.Fmean(:,c); s = ENV.tc.Fsd(:,c);
        fill([tg; flipud(tg)], [m+s; flipud(m-s)], cBand,'FaceAlpha',0.22,'EdgeColor','none');
        plot(tg, m, '-','Color',cDat,'LineWidth',2.6);
        plot(tg, squeeze(ENV.tc.F(IREF,:,c)), '--','Color',cRef,'LineWidth',0.9);
        for k = 1:2
            [tm, Fm] = iModelNorm(R(k), rel(c));
            plot(tm, Fm, '-','Color',CC{k},'LineWidth',1.6);
        end
        if row == 1
            xlim([-0.01 0.30]); title(sprintf('cycle %d — %.0f ML/s', c, fd.v_restretch(c)));
            if c == 1; ylabel('recovery'); end
        else
            xlim([0.295 0.40]); xlabel('t - release (s)');
            if c == 1; ylabel('re-stretch'); end
        end
    end
end
title(tl, {'Timecourse envelope — thick black = mean of all three preparations, band = \pm SD', ...
           sprintf(['blue dashed = %s, the preparation the model is fitted to.  ' ...
                    'Middle row: release + redevelopment;  bottom row: the re-stretch transient'], REFN)}, ...
      'FontWeight','bold','FontSize',12);
exportgraphics(fig, fullfile(resDir,'fit_envelope_overlay.png'), 'Resolution',150);

%% ---- 9. cost shift ------------------------------------------------------
nm  = iTrim(R(2).fnNames);
c1  = R(1).featCost; c2 = R(2).featCost;
[~,o] = sort(max(c1,c2),'descend');
fig = figure('Units','pixels','Position',[20 20 1650 660],'Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile; hold on; box on; grid on
b = bar([c1(o); c2(o)]','grouped'); b(1).FaceColor = cBas; b(2).FaceColor = cNew;
set(gca,'XTick',1:numel(o),'XTickLabel',nm(o),'XTickLabelRotation',45, ...
        'TickLabelInterpreter','none','YScale','log'); ylim([1e-4 3]);
ylabel('weighted cost'); legend(NAMES,'Interpreter','none','Location','northeast');
title(sprintf('Per-feature cost:  %.4f -> %.4f', sum(c1), sum(c2)),'FontWeight','bold');
nexttile; hold on; box on; grid on
d = c2(o) - c1(o);
b = bar(d); b.FaceColor = 'flat';
b.CData = repmat([0.20 0.45 0.30], numel(d), 1);
b.CData(d>0,:) = repmat([0.72 0.14 0.18], nnz(d>0), 1);
set(gca,'XTick',1:numel(o),'XTickLabel',nm(o),'XTickLabelRotation',45,'TickLabelInterpreter','none');
ylabel('\Delta cost (new - base)'); yline(0,'k-');
title('green = improved, red = paid','FontWeight','bold');
exportgraphics(fig, fullfile(resDir,'fit_costshift.png'), 'Resolution',150);

%% ---- console summary ----------------------------------------------------
fprintf('\n%-22s %10s %10s %10s\n','feature',TAG_BASE,TAG_NEW,'delta');
for i = o
    fprintf('%-22s %10.4f %10.4f %+10.4f\n', nm{i}, c1(i), c2(i), c2(i)-c1(i));
end
fprintf('%-22s %10.4f %10.4f %+10.4f\n','TOTAL',sum(c1),sum(c2),sum(c2)-sum(c1));

save(fullfile(resDir,sprintf('fit_report_%s.mat',TAG_NEW)), 'R','NAMES','nm','c1','c2','-v7.3');
fprintf('\nwrote results/fit_report_%s.mat and 8 figures\n', TAG_NEW);

%% ======================= helpers ========================================
function R = iRunAll(tag, FN)
% Run every protocol for one snapshot and return everything needed to plot.
    params0 = getParams(loadParams(tag), [], true, false);
    params0.mods = {}; params0.g = ones(1,0);
    if nargin > 1 && ~isempty(FN); params0.fn = FN; end
    params0.RunForceVelocity = true;  params0.RunSlack   = true;
    params0.RunSlackPassive  = true;
    params0.RunKtr           = true;  params0.RunStairs  = true;   % PREDICTIONS
    params0.RunForceVelocityTime = false;
    params0.EvalFeatures = true; params0.OptimizeOn = 'Feats';
    params0.RunSlackSegments = 'AllPar'; params0.BreakOnODEUnstable = false;
    params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0;
    params0.MaxRunTime = 600;

    LoadData; RunBakersExp;                                       %#ok<*NASGU>
    [ct,~,~] = evalFeatureCost(features_data, features_model, params0.fn);

    R = struct('tag',tag,'p',params0,'cost',sum(ct),'featCost',ct, ...
               'fnNames',{params0.fn},'fmodel',features_model,'fdata',features_data, ...
               'out_slack',out_slack,'out_ktr',out_ktr,'out_stairs',out_stairs, ...
               'out_pas',out_pas);
end

function iterPlotCycle(dtb, R, t0, win, cDat, CC, LW)
    m = dtb(:,1)-t0 >= win(1) & dtb(:,1)-t0 <= win(2);
    plot(dtb(m,1)-t0, dtb(m,3), '-','Color',cDat,'LineWidth',1.3);
    for k = 1:2
        o = R(k).out_slack;
        mm = o.t-t0 >= win(1) & o.t-t0 <= win(2);
        plot(o.t(mm)-t0, o.Force(mm), '-','Color',CC{k},'LineWidth',LW{k});
    end
    xlim(win);
end

function s = iTrim(fn)
    s = cellfun(@(x) regexprep(strtok(x,'|'), '\[.*$',''), fn, 'uni', 0);
end

function [tm, Fm] = iModelNorm(Rk, t0)
% Model slack trace, time-referenced to a release and force-normalised by the
% model's OWN mean steady -- the same transform the data envelope uses.
    sc = mean(double(Rk.fmodel.steady), 'omitnan');
    if ~isfinite(sc) || sc == 0; sc = 1; end
    tm = Rk.out_slack.t - t0;
    Fm = Rk.out_slack.Force / sc;
end
