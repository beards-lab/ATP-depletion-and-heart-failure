% PlotReportFigures.m — report figures on the FULL battery, scored with the
% OPTIMISER's objective, so every number here is directly comparable to the cost
% the refit is minimising.
%
%   results/fig_slack_full.png     whole slack protocol, data vs model(s)
%   results/fig_slack_cycle1.png   zoom on the first slack cycle + length command
%   results/fig_fv.png             force-velocity, data vs model
%   results/fig_passive.png        passive (PNB/Mava) protocol, data vs model
%   results/fig_features.png       feature-cost breakdown, full width
%
% Battery: FV + slack + passive (RunKtr = 0, RunStairs = 0) — exactly what
% RunOptimRealign.m runs. Objective: params0.fn with DROP_FEATS removed, matching
% the reoptim (ovrsht_dy dropped on request; the same window is still held by
% coolDownLS and rsK).
%
% Set CASES before running (label, snapshot path or '' for the analysis baseline,
% override struct).

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
res = fullfile(here,'results'); if ~exist(res,'dir'); mkdir(res); end
if isempty(gcp('nocreate')); parpool('local', 4); end

SNAP0 = 'params/rskR2_w025_opt.m';
if ~exist('CASES','var') || isempty(CASES)
    CASES = { 'baseline',          '',                      struct()
              'reoptim (realign2)','params/realign2_opt.m', struct() };
end
if ~exist('DROP_FEATS','var'); DROP_FEATS = {'ovrsht_dy'}; end

S  = load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat'));
vt = S.velocitytable; dt = S.datatable;
iS = find(vt(:,2) < -1); iR = find(vt(:,2) > 1);
t0 = dt(1,1);

R = struct('name',{},'t',{},'F',{},'feat',{},'ct2',{},'L2',{},'fv',{},'pas',{});
for c = 1:size(CASES,1)
    nm = CASES{c,1}; snap = CASES{c,2}; over = CASES{c,3};
    if isempty(snap); snap = SNAP0; end
    params0 = getParams(loadParams(snap), [], true, false);
    f = fieldnames(over); for i=1:numel(f); params0.(f{i}) = over.(f{i}); end
    params0 = getParams(params0, [], true, false);

    params0.MaxRunTime = 900;
    params0.RunForceVelocity = 1; params0.RunSlack = 1; params0.RunSlackPassive = 1;
    params0.RunKtr = 0; params0.RunStairs = 0; params0.RunForceVelocityTime = 0;
    params0.EvalFeatures = 1; params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0;
    params0.RunSlackSegments = 'AllPar';
    params0.mods = {}; params0.g = ones(1,0);

    % objective = optimiser's
    prim = cellfun(@(s) regexprep(strtok(s,'|'), '\[.*$',''), params0.fn, 'uni', 0);
    params0.fn = params0.fn(~ismember(prim, DROP_FEATS));

    ATP_c = params0.MgATP;               %#ok<NASGU> read by RunBakersExp
    FV_velocities = params0.FV_velocities; %#ok<NASGU>

    fprintf('running %-20s ... ', nm); tic;
    RunBakersExp;                        % script: leaves out_slack / outs_fv / out_pas
    fprintf('%.0f s\n', toc);

    if c == 1
        keep = cellfun(@(q) isfield(features_model,q) && isfield(features_data,q), ...
                       cellfun(@(s) regexprep(strtok(s,'|'), '\[.*$',''), params0.fn, 'uni', 0));
        FN = params0.fn(keep); FD = features_data;
        if ~all(keep)
            fprintf('  WARNING: %d fn entries unscorable and excluded\n', nnz(~keep));
        end
    end
    ct2 = evalFeatureCost(features_data, features_model, FN, 2);

    R(end+1) = struct('name', nm, 't', out_slack.t(:), 'F', out_slack.Force(:), ...
        'feat', features_model, 'ct2', ct2(:), 'L2', sum(ct2), ...
        'fv', struct('v', features_model.FV_v, 'fnorm', features_model.FV_fnorm), ...
        'pas', struct('t', out_pas.t(:), 'F', out_pas.Force(:), 'dt', out_pas.datatable)); %#ok<SAGROW>
    fprintf('  L2 (optimiser objective) = %.4f\n', sum(ct2));
end

col = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];
LW  = 1.8;
lab = arrayfun(@(r) sprintf('%s  (L2 %.3f)', r.name, r.L2), R, 'uni', 0);

%% ---- FIG 1: full slack ---------------------------------------------------
f1 = figure('Position',[40 40 1800 620],'Color','w'); hold on;
plot(dt(:,1)-t0, dt(:,3), 'Color',[.55 .55 .55], 'LineWidth',1.0, 'DisplayName','data 03/27 8 mM');
for i=1:numel(R); plot(R(i).t-t0, R(i).F, 'Color',col(i,:), 'LineWidth',LW, 'DisplayName',lab{i}); end
xlabel('time (s)'); ylabel('force (kPa)'); grid on; box on; set(gca,'FontSize',11);
legend('Location','southeast','FontSize',11);
xlim([vt(iS(1),1)-t0-0.06, vt(end,1)-t0]);
title('Slack protocol — force, model vs data (03/27, 8 mM ATP)','FontSize',13);
exportgraphics(f1, fullfile(res,'fig_slack_full.png'), 'Resolution', 150);

%% ---- FIG 2: first slack cycle -------------------------------------------
ta = vt(iS(1),1) - 0.020; tb = vt(iS(2),1) + 0.010;
f2 = figure('Position',[40 40 1500 800],'Color','w');
ax1 = subplot(4,1,1:3); hold on;
m = dt(:,1) >= ta & dt(:,1) <= tb;
plot((dt(m,1)-vt(iS(1),1))*1000, dt(m,3), 'Color',[.55 .55 .55], 'LineWidth',1.0, 'DisplayName','data');
for i=1:numel(R)
    mm = R(i).t >= ta & R(i).t <= tb;
    plot((R(i).t(mm)-vt(iS(1),1))*1000, R(i).F(mm), 'Color',col(i,:), 'LineWidth',LW, 'DisplayName',R(i).name);
end
yl = ylim;
for xv = [(vt(iR(1),1)-vt(iS(1),1))*1000, (vt(iR(1)+1,1)-vt(iS(1),1))*1000]
    plot([xv xv], yl, 'k:', 'HandleVisibility','off');
end
text((vt(iR(1),1)-vt(iS(1),1))*1000, yl(2), ' re-stretch', 'VerticalAlignment','top','FontSize',10);
text((vt(iR(1)+1,1)-vt(iS(1),1))*1000, yl(1), ' ramp end', 'VerticalAlignment','bottom','FontSize',10);
ylabel('force (kPa)'); grid on; box on; legend('Location','east','FontSize',11);
ttl = sprintf('First slack cycle.   rsK: data %.0f', FD.rsK(1));
for i=1:numel(R); ttl = [ttl sprintf(' | %s %.0f', R(i).name, R(i).feat.rsK(1))]; end %#ok<AGROW>
title(ttl,'FontSize',12); set(gca,'FontSize',11);
ax2 = subplot(4,1,4);
tq = linspace(ta, tb, 3000); Lq = interp1(vt(:,1), vt(:,4), tq, 'previous','extrap');
plot((tq-vt(iS(1),1))*1000, Lq, 'k-', 'LineWidth',1.6);
xlabel('time after first release (ms)'); ylabel('length (ML)'); grid on; box on; set(gca,'FontSize',11);
linkaxes([ax1 ax2],'x'); xlim([(ta-vt(iS(1),1))*1000, (tb-vt(iS(1),1))*1000]);
exportgraphics(f2, fullfile(res,'fig_slack_cycle1.png'), 'Resolution', 150);

%% ---- FIG 3: force-velocity ----------------------------------------------
f3 = figure('Position',[60 60 1100 640],'Color','w'); hold on;
plot(FD.FV_v, FD.FV_fnorm, 'ks--', 'MarkerSize',11,'MarkerFaceColor','k','LineWidth',1.6, ...
     'DisplayName','data (Baker 2022, 8 mM)');
for i=1:numel(R)
    plot(R(i).fv.v, R(i).fv.fnorm, 'o-', 'Color',col(i,:), 'MarkerSize',9, ...
         'MarkerFaceColor','w','LineWidth',LW, 'DisplayName',lab{i});
end
xlabel('shortening velocity (\mum/s)'); ylabel('force / F_{iso}  (normalised)');
grid on; box on; set(gca,'FontSize',12); legend('Location','northeast','FontSize',11);
title('Force–velocity — the term the reoptim must not break (FV\_fnorm, weight 10)','FontSize',13);
exportgraphics(f3, fullfile(res,'fig_fv.png'), 'Resolution', 150);

%% ---- FIG 4: passive ------------------------------------------------------
f4 = figure('Position',[60 60 1500 620],'Color','w'); hold on;
pd = R(1).pas.dt;
plot(pd(:,1)-pd(1,1), pd(:,3), 'Color',[.55 .55 .55],'LineWidth',1.0,'DisplayName','data (PNB+Mava)');
for i=1:numel(R)
    plot(R(i).pas.t-pd(1,1), R(i).pas.F, 'Color',col(i,:),'LineWidth',LW,'DisplayName',lab{i});
end
xlabel('time (s)'); ylabel('force (kPa)'); grid on; box on; set(gca,'FontSize',11);
legend('Location','northeast','FontSize',11);
title('Passive protocol (active attachment OFF) — PS\_* features','FontSize',13);
exportgraphics(f4, fullfile(res,'fig_passive.png'), 'Resolution', 150);

%% ---- FIG 5: features -----------------------------------------------------
nm = cellfun(@(x) regexprep(strtok(x,'|'), '\[.*$',''), FN, 'uni', 0);
M  = cell2mat(arrayfun(@(r) r.ct2(:), R, 'uni', 0));
sel = find(max(M,[],2) > 1e-3); [~,ord] = sort(M(sel,1),'descend'); sel = sel(ord);
f5 = figure('Position',[20 20 1900 780],'Color','w');
b = bar(M(sel,:),'grouped'); for i=1:numel(b); b(i).FaceColor = col(i,:); end
set(gca,'XTick',1:numel(sel),'XTickLabel',nm(sel),'XTickLabelRotation',40,'YScale','log','FontSize',12);
ylabel('feature cost (L2, log scale)','FontSize',13); grid on; box on;
legend(arrayfun(@(r) sprintf('%s   TOTAL %.3f', r.name, r.L2), R, 'uni',0), ...
       'Location','northeast','FontSize',12);
title('Feature cost — FULL battery (FV + slack + passive), optimiser objective','FontSize',14);
exportgraphics(f5, fullfile(res,'fig_features.png'), 'Resolution', 150);

%% ---- tables --------------------------------------------------------------
fprintf('\n%-20s', 'feature'); for i=1:numel(R); fprintf(' %16s', R(i).name); end; fprintf('\n');
for j = sel(:)'
    fprintf('%-20s', nm{j}); for i=1:numel(R); fprintf(' %16.4f', R(i).ct2(j)); end; fprintf('\n');
end
fprintf('%-20s', 'TOTAL'); for i=1:numel(R); fprintf(' %16.4f', R(i).L2); end; fprintf('\n');

fprintf('\n%-14s %10s', 'observable','DATA'); for i=1:numel(R); fprintf(' %16s', R(i).name); end; fprintf('\n');
for k = {'rsK','ktr','steady','peak1_y','vall_y','peak2','vall2_dy','FV_fnorm'}
    if ~isfield(FD,k{1}); continue; end
    fprintf('%-14s %10.3f', k{1}, mean(FD.(k{1}),'omitnan'));
    for i=1:numel(R); fprintf(' %16.3f', mean(R(i).feat.(k{1}),'omitnan')); end
    fprintf('\n');
end
save(fullfile(res,'report_figs.mat'),'R','FD','FN','CASES');
fprintf('\nsaved fig_slack_full / fig_slack_cycle1 / fig_fv / fig_passive / fig_features .png\n');
