% RunPeakAudit.m -- are peak1_dSL and vall_t measurements, or detector jitter?
%
% They are 57 % of the refit3 objective (0.1326 + 0.1112 of 0.4271) and the two
% noisiest features in the set: between-preparation relative SD 52 % and 57 %,
% against 13-25 % for everything else. Neither has ever been audited. The last
% restretch metric that was (rsK -> rsK80) turned out to be measuring the wrong
% thing.
%
% THE SUSPICION, from the definition:
%     peak1_dSL = SL(t_peak) - SL(window end)
%   SL is the COMMANDED motor length x2, not a measurement -- and it is
%   identical across the three preparations by protocol. So every bit of
%   peak1_dSL's between-prep spread has to come from where findpeaks puts
%   t_peak on the raw force trace. During the restretch ramp SL advances
%   ~0.0025 um per 0.25 ms sample, and the feature's whole range is
%   0.007-0.030 um -- so ~3 samples of peak jitter spans the entire feature.
%
% TESTS
%   T1  is SL really protocol-identical across preps?
%   T2  what is d(peak1_dSL)/d(sample) during the ramp?
%   T3  bootstrap: add each trace's OWN measured noise, re-detect, and get the
%       estimator's SD. If that approaches the between-prep SD, the "spread" is
%       the detector, not the preparation.
%   T4  do better-posed alternatives have a tighter spread?
%
% Pure data analysis -- no model runs, so it does not compete with a campaign.
%
%   matlab -batch "run('Analyses/RestretchPeakAudit/RunPeakAudit.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RestretchPeakAudit','results');
if ~exist(resdir,'dir'); mkdir(resdir); end
rng(7);
DAYS = {'03_27','04_03','04_10'};
PEAK_MIN_DISTANCE_S = 0.001;   % as in extractSlackAttributes.m:30
PEAK_MIN_PROMINENCE = 0.1;     % ...:31

S = struct();
for d = 1:3
    L = load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d})));
    S(d).dt = L.datatable; S(d).vt = L.velocitytable; S(d).fd = L.features_data;
end
iS = find(S(1).vt(:,2)<-1); iR = find(S(1).vt(:,2)>1);

%% ---- T1: is SL the same command on every day? --------------------------
fprintf('\n=============== T1: is SL a measurement or the command? ===============\n');
for d = 1:3
    v = S(d).vt; dt = S(d).dt;
    m = dt(:,1) >= v(iR(2),1) & dt(:,1) <= v(iR(2)+1,1);
    t = dt(m,1)-v(iR(2),1); sl = dt(m,2);
    S(d).rampT = t; S(d).rampSL = sl;
    fprintf('  %s  ramp %5.2f ms, SL %.4f -> %.4f um, mean step %.5f um/sample\n', ...
            DAYS{d}, 1e3*(t(end)-t(1)), sl(1), sl(end), mean(diff(sl)));
end
n = min(arrayfun(@(x) numel(x.rampSL), S));
SLm = cell2mat(arrayfun(@(x) x.rampSL(1:n), S, 'UniformOutput', false));
fprintf('  max |SL_day - SL_mean| over the ramp: %.6f um  (sample step %.5f um)\n', ...
        max(max(abs(SLm - mean(SLm,2)))), mean(mean(diff(SLm))));
fprintf('  -> SL is %s\n', ternary(max(max(abs(SLm-mean(SLm,2)))) < 2*mean(mean(diff(SLm))), ...
        'the SAME COMMANDED trajectory on all three days', 'genuinely different per day'));

%% ---- T2: sensitivity of peak1_dSL to one sample of peak jitter ---------
fprintf('\n=============== T2: how much is one sample worth? ===============\n');
dSLds = mean(mean(diff(SLm)));
PT = pooledTarget8mM(); jP = find(strcmp(PT.name,'peak1_dSL'),1);
tgt = PT.val(jP,:);
fprintf('  SL advances %.5f um per 0.25 ms sample during the restretch\n', dSLds);
fprintf('  pooled peak1_dSL target: %s um  (range %.4f)\n', mat2str(round(tgt,4)), max(tgt)-min(tgt));
fprintf('  -> the FULL range of the feature is %.1f samples of peak-detector jitter\n', ...
        (max(tgt)-min(tgt))/dSLds);
fprintf('  -> one sample of jitter = %.0f %% of the mean feature value\n', 100*dSLds/mean(tgt));

%% ---- T3: bootstrap the detector against each trace's own noise ---------
fprintf('\n=============== T3: detector jitter vs between-prep spread ===============\n');
NB = 200;
pk = nan(3,5); vt_ = nan(3,5); pkB = nan(3,5,NB); vtB = nan(3,5,NB);
for d = 1:3
    dt = S(d).dt; v = S(d).vt;
    % noise level: RMS of the high-pass of a quiet isometric stretch
    q = dt(dt(:,1) > v(iS(1),1)-0.30 & dt(:,1) < v(iS(1),1)-0.02, 3);
    sig = std(q - movmean(q,9));
    S(d).sigma = sig;
    for s = 1:5
        m  = dt(:,1) >= v(iR(s),1) & dt(:,1) <= v(iR(s)+1,1);
        t  = dt(m,1); y = dt(m,3); sl = dt(m,2);
        [pk(d,s), vt_(d,s)] = iPeakVall(t, y, sl, PEAK_MIN_DISTANCE_S, PEAK_MIN_PROMINENCE);
        for b = 1:NB
            [pkB(d,s,b), vtB(d,s,b)] = iPeakVall(t, y + sig*randn(size(y)), sl, ...
                                                 PEAK_MIN_DISTANCE_S, PEAK_MIN_PROMINENCE);
        end
    end
    fprintf('  %s  force noise sigma = %.3f kPa\n', DAYS{d}, sig);
end
relBetween = @(X) std(X,0,1)./abs(mean(X,1));
jitP = squeeze(std(pkB,0,3))./abs(pk);     % per prep x segment, relative
jitV = squeeze(std(vtB,0,3))./abs(vt_);
fprintf('\n  %-16s %8s %8s %8s %8s %8s   mean\n','','seg1','seg2','seg3','seg4','seg5');
fprintf('  %-16s %8.3f %8.3f %8.3f %8.3f %8.3f   %.3f\n','peak1_dSL between', relBetween(pk), mean(relBetween(pk)));
fprintf('  %-16s %8.3f %8.3f %8.3f %8.3f %8.3f   %.3f\n','peak1_dSL JITTER', mean(jitP,1), mean(jitP(:)));
fprintf('  %-16s %8.3f %8.3f %8.3f %8.3f %8.3f   %.3f\n','vall_t between',   relBetween(vt_), mean(relBetween(vt_)));
fprintf('  %-16s %8.3f %8.3f %8.3f %8.3f %8.3f   %.3f\n','vall_t JITTER',    mean(jitV,1), mean(jitV(:)));
fprintf('\n  detector jitter explains %.0f %% of peak1_dSL''s spread and %.0f %% of vall_t''s\n', ...
        100*mean(jitP(:))/mean(relBetween(pk)), 100*mean(jitV(:))/mean(relBetween(vt_)));

%% ---- T4: better-posed alternatives -------------------------------------
fprintf('\n=============== T4: alternatives, by between-prep spread ===============\n');
ALT = struct('name',{},'val',{});
p1y = nan(3,5); p1t = nan(3,5); dSLf = nan(3,5); stdy = nan(3,5); t80 = nan(3,5);
for d = 1:3
    dt = S(d).dt; v = S(d).vt;
    for s = 1:5
        m = dt(:,1) >= v(iR(s),1) & dt(:,1) <= v(iR(s)+1,1);
        t = dt(m,1); y = dt(m,3); sl = dt(m,2);
        [~,~,tp,yp] = iPeakVall(t,y,sl,PEAK_MIN_DISTANCE_S,PEAK_MIN_PROMINENCE);
        p1y(d,s) = yp;
        p1t(d,s) = tp - t(1);                                   % peak time from ramp start
        % SL reached when force first crosses 80 % of the peak -- a geometry
        % measure that does not depend on locating a maximum
        i80 = find(y >= 0.8*yp, 1);
        if ~isempty(i80); dSLf(d,s) = sl(i80) - sl(end); end
        t80(d,s) = t(max(i80,1)) - t(1);
        ss = dt(:,1) >= v(iR(s)+1,1)+0.20 & dt(:,1) <= v(iR(s)+1,1)+0.28;
        stdy(d,s) = median(dt(ss,3));
    end
end
tbl = { 'peak1_dSL  (current)',      pk
        'vall_t     (current)',      vt_
        'peak1_t    from ramp start',p1t
        'dSL at 80% of peak',        dSLf
        't at 80% of peak',          t80
        'peak1_y / steady',          p1y./stdy };
fprintf('  %-28s %8s %8s\n','quantity','mean rel SD','vs peak1_dSL');
base = mean(relBetween(pk));
for i = 1:size(tbl,1)
    r = mean(relBetween(tbl{i,2}),'omitnan');
    fprintf('  %-28s %8.3f %11.2fx\n', tbl{i,1}, r, r/base);
end

%% ---- figure -------------------------------------------------------------
f=figure('Position',[30 30 1500 780],'Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
CC=lines(3);
for s = [1 3 5]
    nexttile; hold on; box on; grid on
    for d = 1:3
        dt=S(d).dt; v=S(d).vt;
        m = dt(:,1) >= v(iR(s),1)-0.001 & dt(:,1) <= v(iR(s),1)+0.012;
        plot((dt(m,1)-v(iR(s),1))*1e3, dt(m,3), '-', 'Color', CC(d,:), 'LineWidth', 1.1);
        [~,~,tp,yp] = iPeakVall(dt(dt(:,1)>=v(iR(s),1)&dt(:,1)<=v(iR(s)+1,1),1), ...
                                dt(dt(:,1)>=v(iR(s),1)&dt(:,1)<=v(iR(s)+1,1),3), ...
                                dt(dt(:,1)>=v(iR(s),1)&dt(:,1)<=v(iR(s)+1,1),2), ...
                                PEAK_MIN_DISTANCE_S,PEAK_MIN_PROMINENCE);
        plot((tp-v(iR(s),1))*1e3, yp, 'v', 'Color', CC(d,:), 'MarkerFaceColor', CC(d,:), 'MarkerSize',8);
    end
    xlabel('t since restretch start (ms)'); ylabel('F (kPa)');
    title(sprintf('restretch %d — detected peak1 (triangles)', s));
    if s==1; legend(DAYS,'Box','off','Location','southeast','Interpreter','none'); end
end
nexttile; hold on; box on; grid on
bar([mean(relBetween(pk)) mean(jitP(:)); mean(relBetween(vt_)) mean(jitV(:))]);
set(gca,'XTick',1:2,'XTickLabel',{'peak1\_dSL','vall\_t'}); ylabel('relative SD');
legend({'between preparations','detector jitter alone'},'Box','off'); grid on
title('how much of the "spread" is the detector');
nexttile; hold on; box on; grid on
r = cellfun(@(x) mean(relBetween(x),'omitnan'), tbl(:,2));
barh(r); set(gca,'YTick',1:numel(r),'YTickLabel',tbl(:,1),'YDir','reverse', ...
    'TickLabelInterpreter','none','FontSize',9);
xlabel('mean between-preparation relative SD'); title('alternatives'); grid on
nexttile; hold on; box on; grid on
plot(1:5, tgt, 'ko-','LineWidth',2,'MarkerFaceColor','k');
for d=1:3; plot(1:5, pk(d,:), 's-','Color',CC(d,:),'LineWidth',1.2); end
xlabel('slack segment'); ylabel('peak1\_dSL (\mum)'); xticks(1:5);
yline(dSLds,':k'); text(3.2, dSLds*1.4, 'one sample of jitter', 'FontSize',8);
title('the "arc" against one sample of detector jitter');
legend([{'pooled target'} DAYS],'Box','off','Location','northwest','Interpreter','none');
exportgraphics(f, fullfile(resdir,'peak_audit.png'),'Resolution',120);
fprintf('\nwrote %s\n', fullfile(resdir,'peak_audit.png'));
save(fullfile(resdir,'peak_audit.mat'),'pk','vt_','jitP','jitV','tbl','dSLds');
fprintf('\nDONE\n');

function [dSL, tv, tp, yp] = iPeakVall(t, y, sl, mind, minprom)
    dSL = NaN; tv = NaN; tp = NaN; yp = NaN;
    [py, pt] = findpeaks(y, t, 'MinPeakDistance', mind, 'MinPeakProminence', minprom);
    if isempty(py); return; end
    tp = pt(1); yp = py(1);
    i = find(t >= tp, 1);
    dSL = sl(i) - sl(end);
    after = t > tp;
    if any(after)
        [~, k] = min(y(after)); ta = t(after);
        tv = ta(k) - t(1);
    end
end

function o = ternary(c,a,b); if c; o=a; else; o=b; end; end
