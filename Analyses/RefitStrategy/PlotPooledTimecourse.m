% PlotPooledTimecourse.m -- put the model against the target it is ACTUALLY
% fitted to, instead of against one strong preparation.
%
% The earlier figure plotted 03/27, which is the strongest of the three preps
% (Fss22 79.76 vs a pooled 69.35, range 1.37x), so the model looked ~14 % low
% when it was in fact +4 % on its own target. This builds a pooled TIMECOURSE
% the same way the pooled FEATURES are built:
%
%   1. align each prep on its own slack #1 onset (identical 5-slack ML ladder
%      and 0.25 ms sampling on all three; only the absolute start time differs)
%   2. divide each prep's force by its own pre-slack isometric at SL 2.2, then
%      multiply by the pooled Fref -- the per-PREP anchor of
%      Analyses/DataPooling, which is what makes the level a measurement
%   3. take the MEAN and the MEDIAN across preps
%   4. extract features from each with extractSlackAttributes and keep whichever
%      reproduces pooledTarget8mM more closely -- the pooled target pools
%      FEATURES, not traces, so the two need not agree and the closer one is
%      the honest backdrop
%
%   matlab -batch "run('Analyses/RefitStrategy/PlotPooledTimecourse.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RefitStrategy','results');
c = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_pt_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
c.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(c, 5); end

DAYS = {'03_27','04_03','04_10'};
PT   = pooledTarget8mM();
fprintf('pooled Fref %.2f kPa, per prep %s (CV %.1f %%)\n', ...
        PT.Fref, mat2str(round(PT.FrefPrep,2)), 100*PT.FrefCV);

%% ---- 1-2. align and normalise ------------------------------------------
DT = 2.5e-4; TREL = (-0.40 : DT : 3.55)';        % relative to slack #1 onset
Y = nan(numel(TREL), 3); SLm = nan(numel(TREL), 3); anc = nan(1,3);
for i = 1:3
    S = load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{i})));
    v = S.velocitytable; iS = find(v(:,2)<-1); t0 = v(iS(1),1);
    d = S.datatable;
    % anchor: the prep's own pre-slack isometric at SL 2.2, 300 ms before slack 1
    anc(i) = median(d(d(:,1) > t0-0.30 & d(:,1) < t0-0.01, 3));
    Y(:,i)   = interp1(d(:,1)-t0, d(:,3), TREL, 'linear', NaN) * (PT.Fref/anc(i));
    SLm(:,i) = interp1(d(:,1)-t0, d(:,2), TREL, 'linear', NaN);
    if i == 1; VT1 = v; VT1(:,1) = VT1(:,1) - t0; end
end
fprintf('measured Fss22 anchors %s kPa  (pooledTarget FrefPrep %s)\n', ...
        mat2str(round(anc,2)), mat2str(round(PT.FrefPrep,2)));
POOL_MEAN = mean(Y,2,'omitnan');
POOL_MED  = median(Y,2,'omitnan');
SL_MEAN   = mean(SLm,2,'omitnan');

%% ---- 3. which backdrop reproduces the pooled FEATURES better? ----------
fprintf('\n=== does the pooled TRACE reproduce the pooled FEATURES? ===\n');
FLD = {'A','peak1_y','vall_y','peak2','steady','vall2_dy'};
cand = {'mean', POOL_MEAN; 'median', POOL_MED};
err = zeros(1,2); FT = cell(1,2);
for k = 1:2
    ok = isfinite(cand{k,2});
    f = extractSlackAttributes(TREL(ok), cand{k,2}(ok), SL_MEAN(ok), VT1, struct(), [], false, true);
    FT{k} = f;
    e = [];
    for q = FLD
        j = find(strcmp(PT.name, q{1}),1);
        if isempty(j) || ~isfield(f,q{1}); continue; end
        tgt = PT.val(j,:); mod_ = f.(q{1})(:)';
        n = min(numel(tgt), numel(mod_));
        e = [e, abs(mod_(1:n)-tgt(1:n))./abs(tgt(1:n))]; %#ok
    end
    err(k) = 100*mean(e,'omitnan');
    fprintf('  %-7s trace: mean |error| vs the pooled feature target = %.2f %%\n', cand{k,1}, err(k));
end
[~,kBest] = min(err);
POOL = cand{kBest,2};
fprintf('  -> using the %s trace as the backdrop\n', cand{kBest,1});
for q = FLD
    j = find(strcmp(PT.name,q{1}),1); if isempty(j); continue; end
    fprintf('    %-10s target %s\n', q{1}, mat2str(round(PT.val(j,:),2)));
    fprintf('    %-10s trace  %s\n', '',   mat2str(round(FT{kBest}.(q{1})(:)',2)));
end

%% ---- 4. the models -----------------------------------------------------
TAGS = {'ktrfence_seed','ktrfence_opt'};
M = cell(1,2);
for k = 1:2
    p = getParams(loadParams(TAGS{k}), [], true, false);
    p.PlotEachSeparately=0; p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0;
    p.EvalFeatures=0; p.BreakOnODEUnstable=0; p.MaxRunTime=600;
    p = getParams(p, p.g, true);
    [~, os] = runSlackExperiment(p);
    svt = load(fullfile(root,'data',p.velocitytableonfile)).velocitytable;
    iS = find(svt(:,2)<-1);
    M{k} = struct('t', os.t(:)-svt(iS(1),1), 'F', os.Force(:));
end

%% ---- 5. figure ---------------------------------------------------------
f1=figure('Position',[20 20 1600 880],'Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
CS=[0.85 0.33 0.10]; CO=[0 0.45 0.74]; CP=[0.15 0.15 0.15]; CB=[0.75 0.75 0.75];
iSv = find(VT1(:,2)<-1); iRv = find(VT1(:,2)>1);

nexttile([1 2]); hold on; box on; grid on
for i=1:3; plot(TREL, Y(:,i), '-', 'Color', CB, 'LineWidth', 0.6); end
plot(TREL, POOL, '-', 'Color', CP, 'LineWidth', 2.0);
plot(M{1}.t, M{1}.F, '-', 'Color', CS, 'LineWidth', 1.2);
plot(M{2}.t, M{2}.F, '-', 'Color', CO, 'LineWidth', 1.6);
% the pooled FEATURE targets, drawn where they are measured
for q = {'A','peak1_y','peak2','steady'}
    j = find(strcmp(PT.name,q{1}),1); if isempty(j); continue; end
    for s = 1:5
        switch q{1}
            case 'A';       x = VT1(iSv(s)+1,1) + 0.25;
            case 'peak1_y'; x = VT1(iRv(s),1) + 0.004;
            case 'peak2';   x = VT1(iRv(s),1) + 0.022;
            case 'steady';  x = VT1(iRv(s),1) + 0.24;
        end
        plot(x, PT.val(j,s), 'v', 'Color',[0.85 0 0], 'MarkerFaceColor',[0.85 0 0], 'MarkerSize',5);
    end
end
xlim([-0.3 3.4]); xlabel('t since slack #1 onset (s)'); ylabel('F (kPa)');
title(sprintf(['A   the target the fit ACTUALLY chases: pooled %s of 3 preps, each divided by its own ' ...
   'Fss22 and rescaled to Fref = %.1f kPa'], cand{kBest,1}, PT.Fref));
legend({'03/27 norm.','04/03 norm.','04/10 norm.',sprintf('pooled %s',cand{kBest,1}), ...
        'seed','opt','pooled feature targets'},'Location','southeast','Box','off','NumColumns',2);

% --- cycle 2 zoom ---
nexttile; hold on; box on; grid on
t0=VT1(iSv(2),1); t1=VT1(iRv(2)+1,1)+0.30;
m=TREL>=t0-0.05 & TREL<=t1;
for i=1:3; plot((TREL(m)-t0)*1e3, Y(m,i), '-','Color',CB,'LineWidth',0.6); end
plot((TREL(m)-t0)*1e3, POOL(m), '-','Color',CP,'LineWidth',2.0);
for k=1:2
    mm=M{k}.t>=t0-0.05 & M{k}.t<=t1;
    plot((M{k}.t(mm)-t0)*1e3, M{k}.F(mm),'-','Color',(k==1)*CS+(k==2)*CO,'LineWidth',1.2+0.4*(k==2));
end
xlabel('t since slack #2 release (ms)'); ylabel('F (kPa)');
title('B   slack cycle 2 against the pooled trace');

% --- A/steady ratio, the thing the objective cannot see ---
nexttile; hold on; box on; grid on
jA=find(strcmp(PT.name,'A'),1); jS=find(strcmp(PT.name,'steady'),1);
rt = PT.val(jA,:)./PT.val(jS,:);
rs = [58.826 55.781 52.071 47.975 42.941]./[71.971 72.038 72.136 72.199 55.766];
ro = [56.577 52.364 47.482 42.394 36.390]./[72.952 73.074 73.145 73.223 52.534];
plot(1:5, rt,'o-','Color',CP,'LineWidth',2,'MarkerFaceColor',CP);
plot(1:5, rs,'s-','Color',CS,'LineWidth',1.5);
plot(1:5, ro,'^-','Color',CO,'LineWidth',1.5);
xlabel('slack segment'); ylabel('A / steady'); xlim([0.8 5.2]); xticks(1:5);
title('C   A/steady — a shape the objective never scores (A and steady are separate deadbanded terms)');
legend({'pooled target','seed','opt'},'Location','southwest','Box','off');

exportgraphics(f1, fullfile(resdir,'pooled_timecourse.png'),'Resolution',120);
fprintf('\nwrote %s\n', fullfile(resdir,'pooled_timecourse.png'));
save(fullfile(resdir,'pooled_timecourse.mat'),'TREL','Y','POOL_MEAN','POOL_MED','anc','err','kBest');
fprintf('\nDONE\n');
