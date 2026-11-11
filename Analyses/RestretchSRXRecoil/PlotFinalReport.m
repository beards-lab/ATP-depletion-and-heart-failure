% PlotFinalReport.m — the four report figures on the POWER objective.
%
%   results/fin_slack_full.png    whole slack protocol, data vs model(s)
%   results/fin_slack_cycle1.png  zoom on the first slack cycle + length command
%   results/fin_fv.png            force-velocity
%   results/fin_pv.png            POWER-velocity (normalised), both source datasets
%
% Full battery (FV + slack + passive), scored with the realign4 objective:
% ovrsht_dy and FV_fnorm dropped, FV_normpowerAvg|w:FV_normpowerVar|62.2779,
% passive weights x0.446621.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
res = fullfile(here,'results'); if ~exist(res,'dir'); mkdir(res); end
if isempty(gcp('nocreate')); parpool('local', 5); end

% col, then line style. ThursdayNightFever is carried as a GHOST: the historical
% parametrisation this whole line of work descends from, drawn dashed so it reads
% as a reference rather than a competitor.
if ~exist('CASES','var') || isempty(CASES)
    CASES = { 'TNF (ghost)',     'params/ThursdayNightFever.m', [0.49 0.18 0.56], '--'
              'baseline',        'params/rskR2_w025_opt.m',     [0.00 0.45 0.74], '-'
              'realign3',        'params/realign3_opt.m',       [0.47 0.67 0.19], '-'
              'realign4 (best)', 'params/realign4_opt.m',       [0.85 0.33 0.10], '-' };
end

%% ---- objective: built ONCE from the BASELINE snapshot ---------------------
% Deliberately not from CASES{1}: the ghost is an older parametrisation whose own
% fn differs, and the whole comparison is only meaningful on one fixed objective.
pb0  = getParams(loadParams('params/rskR2_w025_opt.m'), [], true, false);
fn   = pb0.fn;
prim = cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), fn,'uni',0);
fn   = fn(~ismember(prim,{'ovrsht_dy','FV_fnorm'}));
prim = prim(~ismember(prim,{'ovrsht_dy','FV_fnorm'}));
for i = find(startsWith(prim,'PS_'))
    tk = strsplit(fn{i},'|');
    if numel(tk)==1, wOld=1; tk{end+1}=''; else, wOld=str2double(tk{end}); end %#ok<AGROW>
    tk{end} = num2str(wOld*0.446621,'%.6g'); fn{i} = strjoin(tk,'|');
end
fn{end+1} = 'FV_normpowerAvg|w:FV_normpowerVar|62.2779';
FN = fn;
fprintf('objective: %d terms (from the baseline snapshot)\n', numel(FN));

S  = load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat'));
vt = S.velocitytable; dt = S.datatable;
iS = find(vt(:,2) < -1); iR = find(vt(:,2) > 1);
t0 = dt(1,1);

R = struct('name',{},'t',{},'F',{},'feat',{},'ct2',{},'L2',{},'pas',{},'col',{},'ls',{});
for c = 1:size(CASES,1)
    nm = CASES{c,1};
    params0 = getParams(loadParams(CASES{c,2}), [], true, false);
    params0.MaxRunTime = 900;
    params0.RunForceVelocity=1; params0.RunSlack=1; params0.RunSlackPassive=1;
    params0.RunKtr=0; params0.RunStairs=0; params0.RunForceVelocityTime=0;
    params0.EvalFeatures=1; params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0;
    params0.RunSlackSegments='AllPar'; params0.mods={}; params0.g=ones(1,0);
    params0.fn = FN;

    ATP_c = params0.MgATP; FV_velocities = params0.FV_velocities; %#ok<NASGU>
    fprintf('running %-18s ... ', nm); tic; RunBakersExp; fprintf('%.0f s\n', toc);
    ct2 = evalFeatureCost(features_data, features_model, FN, 2);
    if ~exist('FD','var'); FD = features_data; end

    R(end+1) = struct('name', nm, 't', out_slack.t(:), 'F', out_slack.Force(:), ...
        'feat', features_model, 'ct2', ct2(:), 'L2', sum(ct2), ...
        'pas', struct('t',out_pas.t(:),'F',out_pas.Force(:),'dt',out_pas.datatable), ...
        'col', CASES{c,3}, 'ls', CASES{c,4}); %#ok<SAGROW>
    fprintf('   L2 = %.4f\n', sum(ct2));
end
lab = arrayfun(@(r) sprintf('%s  (L2 %.3f)', r.name, r.L2), R, 'uni', 0);
LW = 1.8;

%% FIG 1 — full slack
f1 = figure('Position',[40 40 1800 620],'Color','w'); hold on;
plot(dt(:,1)-t0, dt(:,3), 'Color',[.55 .55 .55],'LineWidth',1.0,'DisplayName','data 03/27 8 mM');
for i=1:numel(R); plot(R(i).t-t0, R(i).F, R(i).ls, 'Color',R(i).col,'LineWidth',LW,'DisplayName',lab{i}); end
xlabel('time (s)'); ylabel('force (kPa)'); grid on; box on; set(gca,'FontSize',11);
legend('Location','southeast','FontSize',11);
xlim([vt(iS(1),1)-t0-0.06, vt(end,1)-t0]);
title('Slack protocol — force, model vs data (03/27, 8 mM ATP)','FontSize',13);
exportgraphics(f1, fullfile(res,'fin_slack_full.png'),'Resolution',150);

%% FIG 2 — first slack cycle
ta = vt(iS(1),1)-0.020; tb = vt(iS(2),1)+0.010;
f2 = figure('Position',[40 40 1500 800],'Color','w');
ax1 = subplot(4,1,1:3); hold on;
m = dt(:,1)>=ta & dt(:,1)<=tb;
plot((dt(m,1)-vt(iS(1),1))*1000, dt(m,3), 'Color',[.55 .55 .55],'LineWidth',1.0,'DisplayName','data');
for i=1:numel(R)
    mm = R(i).t>=ta & R(i).t<=tb;
    plot((R(i).t(mm)-vt(iS(1),1))*1000, R(i).F(mm), R(i).ls, 'Color',R(i).col,'LineWidth',LW,'DisplayName',R(i).name);
end
yl = ylim;
for xv = [(vt(iR(1),1)-vt(iS(1),1))*1000, (vt(iR(1)+1,1)-vt(iS(1),1))*1000]
    plot([xv xv], yl, 'k:','HandleVisibility','off');
end
text((vt(iR(1),1)-vt(iS(1),1))*1000, yl(2),' re-stretch','VerticalAlignment','top','FontSize',10);
text((vt(iR(1)+1,1)-vt(iS(1),1))*1000, yl(1),' ramp end','VerticalAlignment','bottom','FontSize',10);
ylabel('force (kPa)'); grid on; box on; legend('Location','east','FontSize',11);
ttl = sprintf('First slack cycle.   rsK: data %.0f', FD.rsK(1));
for i=1:numel(R); ttl=[ttl sprintf(' | %s %.0f', R(i).name, R(i).feat.rsK(1))]; end %#ok<AGROW>
title(ttl,'FontSize',12); set(gca,'FontSize',11);
ax2 = subplot(4,1,4);
tq = linspace(ta,tb,3000); Lq = interp1(vt(:,1),vt(:,4),tq,'previous','extrap');
plot((tq-vt(iS(1),1))*1000, Lq,'k-','LineWidth',1.6);
xlabel('time after first release (ms)'); ylabel('length (ML)'); grid on; box on; set(gca,'FontSize',11);
linkaxes([ax1 ax2],'x'); xlim([(ta-vt(iS(1),1))*1000,(tb-vt(iS(1),1))*1000]);
exportgraphics(f2, fullfile(res,'fin_slack_cycle1.png'),'Resolution',150);

%% FIG 3 — force-velocity
f3 = figure('Position',[60 60 1100 640],'Color','w'); hold on;
plot(FD.FV_v, FD.FV_fnorm,'ks--','MarkerSize',11,'MarkerFaceColor','k','LineWidth',1.6, ...
     'DisplayName','data (Baker2022, the fitted set)');
for i=1:numel(R)
    plot(R(i).feat.FV_v, R(i).feat.FV_fnorm,['o' R(i).ls],'Color',R(i).col,'MarkerSize',9, ...
        'MarkerFaceColor','w','LineWidth',LW,'DisplayName',lab{i});
end
xlabel('shortening velocity (\mum/s)'); ylabel('force / F_{iso}');
grid on; box on; set(gca,'FontSize',12); legend('Location','northeast','FontSize',11);
title('Force–velocity  (NB: no longer scored directly — folded into power)','FontSize',13);
exportgraphics(f3, fullfile(res,'fin_fv.png'),'Resolution',150);

%% FIG 4 — POWER-velocity, the scored quantity
f4 = figure('Position',[60 60 1150 660],'Color','w'); hold on;
v = FD.FV_powerv;
hb = plot(v, FD.FV_normpower1,'^:','Color',[.35 .35 .35],'MarkerSize',9,'LineWidth',1.3, ...
     'DisplayName','data source 1 (Baker2022)');
hi = plot(v, FD.FV_normpower2,'v:','Color',[.35 .35 .35],'MarkerSize',9,'LineWidth',1.3, ...
     'DisplayName','data source 2 (Isovel2021)');
errorbar(v, FD.FV_normpowerAvg, sqrt(FD.FV_normpowerVar),'ks-','MarkerSize',12, ...
     'MarkerFaceColor','k','LineWidth',2,'CapSize',10,'DisplayName','data mean \pm sd  (SCORED)');
for i=1:numel(R)
    plot(v, R(i).feat.FV_normpowerAvg,['o' R(i).ls],'Color',R(i).col,'MarkerSize',10, ...
        'MarkerFaceColor','w','LineWidth',2.2,'DisplayName',lab{i});
end
for k=1:numel(v)
    text(v(k), FD.FV_normpowerAvg(k)+0.045, sprintf('w=%.2f', ...
        localW(FD.FV_normpowerVar,k)),'HorizontalAlignment','center','FontSize',9,'Color',[.3 .3 .3]);
end
xlabel('shortening velocity (\mum/s)'); ylabel('normalised power  (F/F_{iso} \times v)');
grid on; box on; set(gca,'FontSize',12); legend('Location','southeast','FontSize',10);
title('Power–velocity — the quantity now being fitted (w = inverse-variance weight)','FontSize',13);
exportgraphics(f4, fullfile(res,'fin_pv.png'),'Resolution',150);

%% FIG 5 — feature cost breakdown
nm = cellfun(@(x) regexprep(strtok(x,'|'),'\[.*$',''), FN,'uni',0);
M  = cell2mat(arrayfun(@(r) r.ct2(:), R,'uni',0));      % nFeat x nCase
sel = find(max(M,[],2) > 1e-3);
[~,ord] = sort(max(M(sel,:),[],2),'descend'); sel = sel(ord);   % rank by worst case,
                                                    % so nothing important is buried
f5 = figure('Position',[20 20 1900 780],'Color','w');
b = bar(M(sel,:),'grouped');
for i=1:numel(b)
    b(i).FaceColor = R(i).col;
    if strcmp(R(i).ls,'--')          % mark the ghost: hatched-looking outline, no fill
        b(i).FaceAlpha = 0.35; b(i).EdgeColor = R(i).col; b(i).LineWidth = 1.4;
    end
end
set(gca,'XTick',1:numel(sel),'XTickLabel',nm(sel),'XTickLabelRotation',40, ...
    'YScale','log','FontSize',12);
ylabel('feature cost (L2, log scale)','FontSize',13); grid on; box on;
legend(arrayfun(@(r) sprintf('%s   TOTAL %.3f', r.name, r.L2), R,'uni',0), ...
       'Location','northeast','FontSize',12);
title('Feature cost — full battery (FV + slack + passive), power objective','FontSize',14);
exportgraphics(f5, fullfile(res,'fin_features.png'),'Resolution',150);

%% tables
[~,o] = sort(M(:,end),'descend');
fprintf('\n%-22s', 'feature'); for i=1:numel(R); fprintf(' %16s', R(i).name); end; fprintf('\n');
for j = o(:)'
    if max(M(j,:)) < 1e-3; continue; end
    fprintf('%-22s', nm{j}); fprintf(' %16.4f', M(j,:)); fprintf('\n');
end
fprintf('%-22s', 'TOTAL'); fprintf(' %16.4f', [R.L2]); fprintf('\n');

fprintf('\n%-16s %10s', 'observable','DATA'); for i=1:numel(R); fprintf(' %16s', R(i).name); end; fprintf('\n');
for k = {'rsK','ktr','steady','peak1_y','vall_y','peak2','vall2_dy'}
    if ~isfield(FD,k{1}); continue; end
    fprintf('%-16s %10.3f', k{1}, mean(FD.(k{1}),'omitnan'));
    for i=1:numel(R); fprintf(' %16.3f', mean(R(i).feat.(k{1}),'omitnan')); end; fprintf('\n');
end
fprintf('\n%-16s %10s', 'norm.power','DATA'); for i=1:numel(R); fprintf(' %16s', R(i).name); end; fprintf('\n');
for k = 1:numel(v)
    fprintf('  v=%-12.1f %10.4f', v(k), FD.FV_normpowerAvg(k));
    for i=1:numel(R); fprintf(' %16.4f', R(i).feat.FV_normpowerAvg(k)); end; fprintf('\n');
end
save(fullfile(res,'final_report.mat'),'R','FD','FN','CASES');
fprintf('\nsaved fin_slack_full / fin_slack_cycle1 / fin_fv / fin_pv / fin_features .png\n');

function w = localW(vr, k)
    ww = 1./max(vr, 0.1*mean(vr)); ww = ww/mean(ww); w = ww(k);
end
