% ReportSnapshot.m — the report: force TIMECOURSE (data vs model) plus the full
% FEATURE table, for any set of parameter snapshots / overrides.
%
% Set CASES before running:
%   CASES = { 'label', <snapshot path or ''>, <override struct> ; ... }
% '' as the path means "the analysis baseline" (params/rskR2_w025_opt.m).
%
% Produces results/report_timecourse.png, results/report_features.png and a
% printed feature table with per-feature L2 cost for every case side by side.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if isempty(gcp('nocreate')); parpool('local', 5); end
res = fullfile(here,'results'); if ~exist(res,'dir'); mkdir(res); end

SNAP0 = 'params/rskR2_w025_opt.m';
if ~exist('CASES','var') || isempty(CASES)
    CASES = { 'baseline',    '', struct()
              'M1 realign',  '', struct('UseLoadRealign',true,'k_cr',6,'tau_cr',0.020, ...
                                        'F_cr',25,'kA2re',30,'eta_M',[]) };  % [] = x0.3, filled below
end

S  = load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat'));
vt = S.velocitytable; dt = S.datatable;
iS = find(vt(:,2) < -1); iR = find(vt(:,2) > 1);
t0 = dt(1,1);                                    % protocol clock offset

R = struct('name',{},'t',{},'F',{},'feat',{},'ct2',{},'L2',{});
for c = 1:size(CASES,1)
    nm = CASES{c,1}; snap = CASES{c,2}; over = CASES{c,3};
    if isempty(snap); snap = SNAP0; end
    p = getParams(loadParams(snap), [], true, false);
    if isfield(over,'eta_M') && isempty(over.eta_M); over.eta_M = 0.3*p.eta_M; end
    f = fieldnames(over); for i=1:numel(f); p.(f{i}) = over.(f{i}); end
    p = getParams(p, [], true, false);
    p.MaxRunTime = 600; p.RunSlackSegments = 'All';   % sequential: real carryover
    p.EvalFeatures = 1; p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0;

    fprintf('running %-14s ... ', nm); tic;
    [~, o, fm, fd] = runSlackExperiment(p);
    fprintf('%.0f s\n', toc);

    % Score ONLY the fn entries this run can actually produce. This is the slack
    % protocol alone, so the FV_* and PS_* entries have no model feature and
    % would each score MISSING_FEATURE_COST = 100 x weight — which looks like an
    % ordinary bad fit in the total and is not one. Same guard as evalRs.
    if c == 1
        prim = cellfun(@(s) regexprep(strtok(s,'|'), '\[.*$',''), p.fn, 'uni', 0);
        scorable = cellfun(@(q) isfield(fm, q) && isfield(fd, q), prim);
        FN = p.fn(scorable); FD = fd;
        fprintf('  scoring %d/%d fn entries (slack protocol only; dropped: %s)\n', ...
                nnz(scorable), numel(p.fn), strjoin(unique(prim(~scorable)), ', '));
    end
    [ct2, ~, ~] = evalFeatureCost(fd, fm, FN, 2);

    R(end+1) = struct('name', nm, 't', o.t(:), 'F', o.Force(:), 'feat', fm, ...
                      'ct2', ct2(:), 'L2', sum(ct2)); %#ok<SAGROW>
end

%% ---- timecourse ----------------------------------------------------------
f1 = figure('Position',[60 60 1500 760],'Color','w');
tiledlayout(f1, 2, 3, 'Padding','compact','TileSpacing','compact');
col = lines(numel(R));

% full protocol
nexttile([1 3]); hold on;
plot(dt(:,1)-t0, dt(:,3), 'Color',[.45 .45 .45], 'LineWidth', 1.2, 'DisplayName','data 03/27 8 mM');
for i=1:numel(R); plot(R(i).t-t0, R(i).F, 'Color', col(i,:), 'LineWidth',1.4, 'DisplayName', R(i).name); end
xlabel('time (s)'); ylabel('force (kPa)'); grid on; legend('Location','southeast');
xlim([vt(iS(1),1)-t0-0.05, vt(end,1)-t0]);
title('slack protocol — full force trace');

% per-cycle zoom on the restretch + recovery, the window rsK scores
for cyc = 1:3
    nexttile; hold on;
    ta = vt(iR(cyc),1) - 0.010; tb = vt(iR(cyc)+1,1) + 0.150;
    m = dt(:,1) >= ta & dt(:,1) <= tb;
    plot((dt(m,1)-vt(iR(cyc),1))*1000, dt(m,3), 'Color',[.45 .45 .45], 'LineWidth',1.4);
    for i=1:numel(R)
        mm = R(i).t >= ta & R(i).t <= tb;
        plot((R(i).t(mm)-vt(iR(cyc),1))*1000, R(i).F(mm), 'Color', col(i,:), 'LineWidth',1.4);
    end
    xlabel('t after re-stretch ramp start (ms)'); ylabel('force (kPa)'); grid on;
    ttl = sprintf('cycle %d  |  rsK data %.0f', cyc, FD.rsK(cyc));
    for i=1:numel(R); ttl = [ttl sprintf('  %s %.0f', R(i).name, R(i).feat.rsK(cyc))]; end %#ok<AGROW>
    title(ttl, 'FontSize', 8);
end
exportgraphics(f1, fullfile(res,'report_timecourse.png'), 'Resolution', 130);

%% ---- feature table -------------------------------------------------------
nm = cellfun(@(x) regexprep(strtok(x,'|'), '\[.*$',''), FN, 'uni', 0);
fprintf('\n%-18s', 'feature');
for i=1:numel(R); fprintf(' %11s', R(i).name); end
fprintf('   (L2 cost per feature)\n');
[~, o] = sort(R(1).ct2, 'descend');
for j = o(:)'
    if max(arrayfun(@(r) r.ct2(j), R)) < 5e-3; continue; end
    fprintf('%-18s', nm{j});
    for i=1:numel(R); fprintf(' %11.3f', R(i).ct2(j)); end
    fprintf('\n');
end
fprintf('%-18s', 'TOTAL L2 (slack)');
for i=1:numel(R); fprintf(' %11.3f', R(i).L2); end
fprintf('\n(slack-protocol features only — FV and passive are not run here)\n');

% headline observables, model vs data
fprintf('\n%-18s %11s', 'observable', 'DATA');
for i=1:numel(R); fprintf(' %11s', R(i).name); end; fprintf('\n');
show = {'rsK','ktr','steady','peak1_y','vall_y','peak2','ovrsht_dy','vall2_dy','A'};
for k = 1:numel(show)
    if ~isfield(FD, show{k}); continue; end
    fprintf('%-18s %11.2f', show{k}, mean(FD.(show{k}),'omitnan'));
    for i=1:numel(R)
        if isfield(R(i).feat, show{k})
            fprintf(' %11.2f', mean(R(i).feat.(show{k}),'omitnan'));
        else; fprintf(' %11s','-'); end
    end
    fprintf('\n');
end

f2 = figure('Position',[60 60 1250 480],'Color','w');
M = cell2mat(arrayfun(@(r) r.ct2(:), R, 'uni', 0));   % nFeat x nCase
sel = find(max(M, [], 2) > 0.02);
[~, ord] = sort(M(sel,1), 'descend'); sel = sel(ord);
bar(M(sel,:));
set(gca,'XTick',1:numel(sel),'XTickLabel',nm(sel),'XTickLabelRotation',45, 'YScale','log');
ylabel('L2 cost (log)'); grid on; legend({R.name}, 'Location','northeast');
title('feature cost breakdown — slack protocol');
exportgraphics(f2, fullfile(res,'report_features.png'), 'Resolution', 130);

save(fullfile(res,'report.mat'), 'R', 'FD', 'FN', 'CASES');
fprintf('\nsaved results/report_timecourse.png, report_features.png, report.mat\n');
