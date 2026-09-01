% MakeReportFigures.m -- the four report figures for the re-based + refitted fit.
%
% Writes into results/figs:
%   fig_fv.png            force-velocity + power-velocity, 8 mM model vs both
%                         data sources, with the 2 mM model curve as a PREDICTION
%                         (there is no 2 mM force-velocity data)
%   fig_slack_full.png    the whole slack protocol, both [ATP], model vs data
%   fig_slack_detail.png  first slack segment, both [ATP], scored features marked
%   fig_features.png      every scored slack feature, both [ATP], data vs model
%
% 8 mM  = params/realign4_opt.m               (Further-optim best)
% 2 mM  = params/lowATP_realign4_gates_opt.m  (that baseline + refitted ATP gates)
%
% Run: matlab -batch "run('Analyses/LowATP_Realign4Transfer/MakeReportFigures.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if isempty(gcp('nocreate')); parpool('local',5); end

if ~exist('SNAP8','var'); SNAP8 = 'realign4_opt'; end
if ~exist('SNAP2','var'); SNAP2 = 'lowATP_realign4_gates_opt'; end
if ~exist('SUFFIX','var'); SUFFIX = ''; end
outdir = fullfile(here,'results',['figs' SUFFIX]);
if ~exist(outdir,'dir'); mkdir(outdir); end
CACHE = fullfile(here,'results',['report_cache' SUFFIX '.mat']);

C8 = [0.00 0.35 0.70];   % 8 mM
C2 = [0.80 0.20 0.10];   % 2 mM
CD = [0.15 0.15 0.15];   % data

%% ---------------------------------------------------------------- run ----
if exist(CACHE,'file') && ~(exist('FORCE_RERUN','var') && FORCE_RERUN)
    fprintf('loading cached runs from %s\n', CACHE); S = load(CACHE);
else
    p8 = iBase(getParams(loadParams(SNAP8), [], true, false), 8);
    p2 = iBase(getParams(loadParams(SNAP2), [], true, false), 2);

    fprintf('\n== FV, 8 mM ==\n');
    q = p8; q.RunForceVelocity = true; q.RunSlack = false;
    [~, ~, S.fvm8, S.fvd8] = runFVExperiment(q, ATP_c, Data_ATP);

    fprintf('== FV, 2 mM (model prediction, no data) ==\n');
    q = p2; q.RunForceVelocity = true; q.RunSlack = false;
    [~, ~, S.fvm2, ~] = runFVExperiment(q, ATP_c, Data_ATP);

    fprintf('== slack, 8 mM ==\n');
    [~, S.o8, S.fm8, S.fd8] = runSlackExperiment(p8);
    fprintf('== slack, 2 mM ==\n');
    [~, S.o2, S.fm2, S.fd2] = runSlackExperiment(p2);

    save(CACHE, '-struct', 'S', '-v7.3');
    fprintf('cached -> %s\n', CACHE);
end
close all force;

%% ------------------------------------------------- 1. FV + power ---------
fig = figure('Units','pixels','Position',[20 20 1600 680],'Color','w');
tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on; grid on
plot(S.fvd8.FV_v(:), S.fvd8.FV_fnorm(:), 'o', 'Color',CD, 'MarkerSize',10, ...
     'LineWidth',1.8,'MarkerFaceColor','w','DisplayName','data (Baker2022)');
plot(S.fvm8.FV_v(:), S.fvm8.FV_fnorm(:), '-o', 'Color',C8,'LineWidth',2.6, ...
     'MarkerSize',8,'MarkerFaceColor',C8,'DisplayName','model 8 mM');
plot(S.fvm2.FV_v(:), S.fvm2.FV_fnorm(:), '--s','Color',C2,'LineWidth',2.4, ...
     'MarkerSize',8,'DisplayName','model 2 mM (prediction)');
xlabel('shortening velocity (ML/s)','FontSize',13);
ylabel('F / F_{iso}','FontSize',13);
title('Force-velocity, normalised','FontSize',15,'FontWeight','bold');
legend('Location','northeast','FontSize',12); set(gca,'FontSize',12);

nexttile; hold on; box on; grid on
vp = S.fvd8.FV_powerv(:);
pA = S.fvd8.FV_normpowerAvg(:); pS = sqrt(max(S.fvd8.FV_normpowerVar(:),0));
errorbar(vp, pA, pS, 'o', 'Color',CD,'MarkerSize',10,'LineWidth',1.8, ...
         'MarkerFaceColor','w','CapSize',10, ...
         'DisplayName','data (mean of 2 sets, bars = spread)');
plot(S.fvm8.FV_powerv(:), S.fvm8.FV_normpower(:), '-o','Color',C8,'LineWidth',2.6, ...
     'MarkerSize',8,'MarkerFaceColor',C8,'DisplayName','model 8 mM');
plot(S.fvm2.FV_powerv(:), S.fvm2.FV_normpower(:), '--s','Color',C2,'LineWidth',2.4, ...
     'MarkerSize',8,'DisplayName','model 2 mM (prediction)');
xlabel('shortening velocity (ML/s)','FontSize',13);
ylabel('normalised power (F/F_{iso}) x v','FontSize',13);
title('Power-velocity: the scored FV term','FontSize',15,'FontWeight','bold');
legend('Location','northeast','FontSize',12); set(gca,'FontSize',12);

title(tl, 'Force-velocity and power-velocity: 8 mM fitted, 2 mM predicted', ...
      'FontSize',16,'FontWeight','bold');
exportgraphics(fig, fullfile(outdir,'fig_fv.png'), 'Resolution', 160);
fprintf('saved fig_fv.png\n');

%% ------------------------------------------- 2. full slack protocol ------
fig = figure('Units','pixels','Position',[20 20 1750 950],'Color','w');
tl = tiledlayout(4,1,'Padding','compact','TileSpacing','tight');
iFullSlack(S.o8, S.fd8, '8 mM ATP', C8, true);
iFullSlack(S.o2, S.fd2, '2 mM ATP', C2, false);
title(tl,'The whole slack protocol: model vs data','FontSize',16,'FontWeight','bold');
exportgraphics(fig, fullfile(outdir,'fig_slack_full.png'), 'Resolution', 150);
fprintf('saved fig_slack_full.png\n');

%% ------------------------------------------- 3. slack detail -------------
fig = figure('Units','pixels','Position',[20 20 1750 980],'Color','w');
tl = tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
iDetail(S.o8, S.fm8, S.fd8, '8 mM ATP', C8);
iDetail(S.o2, S.fm2, S.fd2, '2 mM ATP', C2);
title(tl,'First slack segment: the window every scored feature is read from', ...
      'FontSize',16,'FontWeight','bold');
exportgraphics(fig, fullfile(outdir,'fig_slack_detail.png'), 'Resolution', 150);
fprintf('saved fig_slack_detail.png\n');

%% ------------------------------------------- 4. feature grid -------------
F = { 'steady','steady force (kPa)'; 'A','A, amplitude (kPa)';
      'peak1_y','peak1 (kPa)';       'vall_y','valley (kPa)';
      'peak2','peak2 (kPa)';         'vall2_dy','vall2\_dy (kPa)';
      'ktr','k_{tr} (s^{-1})';       'rsK','rsK, post-restretch rate (s^{-1})';
      'restretchSlopeStart','restretch slope (kPa/s)';
      'peak1_dSL','peak1\_dSL (\mum)'; 'vall_t','valley time (s)';
      'coolDownLS','coolDownLS' };

fig = figure('Units','pixels','Position',[20 20 1900 1150],'Color','w');
tl = tiledlayout(3,4,'Padding','compact','TileSpacing','compact');
for i = 1:size(F,1)
    nexttile; hold on; box on; grid on
    nm = F{i,1};
    d8 = iF(S.fd8,nm); m8 = iF(S.fm8,nm); d2 = iF(S.fd2,nm); m2 = iF(S.fm2,nm);
    n = max([numel(d8) numel(m8) numel(d2) numel(m2) 1]);
    if ~isempty(d8); plot(1:numel(d8),d8,'o-','Color',C8,'MarkerSize',9,'LineWidth',1.8,'MarkerFaceColor','w'); end
    if ~isempty(m8); plot(1:numel(m8),m8,'x--','Color',C8,'MarkerSize',12,'LineWidth',2.0); end
    if ~isempty(d2); plot(1:numel(d2),d2,'o-','Color',C2,'MarkerSize',9,'LineWidth',1.8,'MarkerFaceColor','w'); end
    if ~isempty(m2); plot(1:numel(m2),m2,'x--','Color',C2,'MarkerSize',12,'LineWidth',2.0); end
    xlim([0.6 n+0.4]); xticks(1:n);
    t = F{i,2};
    r8 = iRel(d8,m8); r2 = iRel(d2,m2);
    if ~isnan(r8) || ~isnan(r2)
        t = {t, sprintf('8 mM %.0f%%   |   2 mM %.0f%%', 100*r8, 100*r2)};
    end
    title(t,'FontSize',11,'FontWeight','bold','Interpreter','tex');
    if i > 8; xlabel('slack segment','FontSize',11); end
    set(gca,'FontSize',11);
end
lg = legend({'8 mM data','8 mM model','2 mM data','2 mM model'}, ...
            'Orientation','horizontal','FontSize',13);
lg.Layout.Tile = 'north';
title(tl,'Every scored slack feature: data (o, solid) vs model (x, dashed); % = mean relative error', ...
      'FontSize',16,'FontWeight','bold');
exportgraphics(fig, fullfile(outdir,'fig_features.png'), 'Resolution', 150);
fprintf('saved fig_features.png\n');

fprintf('\nREPORT FIGURES DONE -> %s\n', outdir);

%% ========================================================================
function p = iBase(p, atp)
p.EvalFeatures=true; p.BreakOnODEUnstable=false; p.OptimizeOn='Feats';
p.UseAtpJointFit=false;   % figures score one condition at a time, not the joint objective
p.RunSlackSegments='All'; p.MaxRunTime=600;
p.PlotEachSeparately=0; p.PlotFeatureFitting=0; p.PlotFullscreen=0;
p.RunForceVelocity=false; p.RunForceVelocityTime=false;
p.RunStairs=false; p.RunKtr=false; p.RunSlackPassive=false; p.RunSlack=true;
p.UseAtpDetachR2=0; p.UseAtpKmsrd=0; p.UseStrainDepAtpGate=0;
p.UseR2Ceiling=0; p.UseD0State=0;
p.MgATP=atp; p.MgADP=0;
if atp==8
    p.Pi=0; p.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
    p.ktr_velocitytableonfile='protocol_03_27_2026_velocitytable_ktr.mat';
else
    p.Pi=0.8; p.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
    p.ktr_velocitytableonfile='protocol_03_27_2026_velocitytable_ktr_2mM.mat';
end
end

function iFullSlack(out, fd, ttl, col, isTop)
t = out.t(:); dt = out.datatable; td = dt(:,1); Fd = dt(:,end);
ax = nexttile; hold on; box on; grid on
plot(td, Fd, '-', 'Color',[0.45 0.45 0.45],'LineWidth',0.9,'DisplayName','data');
plot(t, out.Force(:), '-','Color',col,'LineWidth',1.7,'DisplayName','model');
if isfield(fd,'t_seg')
    for k = 1:numel(fd.t_seg)
        xline(fd.t_seg(k),':','Color',[0.6 0.6 0.6],'HandleVisibility','off');
    end
end
ylabel('force (kPa)','FontSize',12); set(ax,'FontSize',11);
title(ttl,'FontSize',14,'FontWeight','bold');
if isTop; legend('Location','northeast','FontSize',11); end
xlim([min(td) max(td)]);

ax2 = nexttile; hold on; box on; grid on
plot(t, out.SL(:), '-','Color',col*0.7,'LineWidth',1.5);
ylabel('SL (\mum)','FontSize',12); set(ax2,'FontSize',11);
xlim([min(td) max(td)]);
if ~isTop; xlabel('time (s)','FontSize',13); end
linkaxes([ax ax2],'x');
end

function iDetail(out, fm, fd, ttl, col)
t = out.t(:); dt = out.datatable; td = dt(:,1); Fd = dt(:,end);
t0 = fd.t_seg(1); w = [t0-0.06, t0+0.62];
ax = nexttile; hold on; box on; grid on
plot(td, Fd, '-','Color',[0.4 0.4 0.4],'LineWidth',1.1,'DisplayName','data');
plot(t,  out.Force(:), '-','Color',col,'LineWidth',2.3,'DisplayName','model');
xlim(w); ylabel('force (kPa)','FontSize',13); set(ax,'FontSize',12);
try
    yline(fd.steady(1),'--','steady','Color',[0.20 0.55 0.15],'LineWidth',1.4, ...
          'FontSize',11,'LabelHorizontalAlignment','left','HandleVisibility','off');
    inw = td>=t0 & td<=t0+0.62; tw = td(inw); Fw = Fd(inw); Fs = movmean(Fw,9);
    [~,iR] = max(diff(Fs));
    seg = (1:numel(tw))' > iR & tw < tw(iR)+0.09;
    [y1,i1] = max(Fs(seg)); ix = find(seg); i1 = ix(i1);
    aft = (1:numel(tw))' > i1 & tw < tw(i1)+0.05;
    [yv,iv] = min(Fs(aft)); ix = find(aft); iv = ix(iv);
    af2 = (1:numel(tw))' > iv & tw < tw(iv)+0.030;
    [y2,i2] = max(Fs(af2)); ix = find(af2); i2 = ix(i2);
    cA = [0.85 0.33 0.10];
    mk = {tw(i1),y1,'peak1'; tw(iv),yv,'valley'; tw(i2),y2,'peak2'};
    for q = 1:size(mk,1)
        plot(mk{q,1},mk{q,2},'o','MarkerSize',9,'LineWidth',1.9,'Color',cA, ...
             'MarkerFaceColor','w','HandleVisibility','off');
        text(mk{q,1},mk{q,2},['  ' mk{q,3}],'Color',cA,'FontSize',11,'FontWeight','bold');
    end
    text(t0+0.055, fd.A(1)*0.40, sprintf('k_{tr}: data %.0f | model %.0f s^{-1}', ...
         fd.ktr(1), fm.ktr(1)),'Color',[0 0.45 0.74],'FontSize',12,'FontWeight','bold');
    text(tw(i2)+0.055, fd.steady(1)*0.46, sprintf('rsK: data %.0f | model %.0f s^{-1}', ...
         fd.rsK(1), fm.rsK(1)),'Color',[0 0.45 0.74],'FontSize',12,'FontWeight','bold');
catch ME
    fprintf('  [markers skipped, %s] %s\n', ttl, ME.message);
end
legend('Location','southeast','FontSize',12);
title(ttl,'FontSize',14,'FontWeight','bold');
xlabel('time (s)','FontSize',13);
end

function v = iF(s,nm)
if isfield(s,nm); v = double(s.(nm)(:))'; else; v = []; end
end

function r = iRel(d,m)
if isempty(d)||isempty(m)||numel(d)~=numel(m)||mean(abs(d))<1e-9; r = NaN; return; end
r = mean(abs(m(:)-d(:))./max(abs(d(:)),eps),'omitnan');
end
