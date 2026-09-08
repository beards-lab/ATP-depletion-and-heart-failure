% RunMetricAudit.m -- is rsK a bad NUMBER or a bad ESTIMATOR?
%
% THE COMPLAINT
% DataPooling found that the post-re-stretch recovery rate `rsK` is the one
% feature that cannot survive trace averaging (-24 %, 1.19 preparation widths)
% and RestretchRecoveryFit found that the model's rsK-vs-amplitude slope has
% the WRONG SIGN. Both could be a property of the muscle -- or an artefact of
% how the rate is put on the trace. This script asks the second question
% properly before any more model structure is blamed.
%
% SEVEN ESTIMATORS (rsEstimators.m), all on the same shape z = (y-B)/(P-B):
%   k_fit   current: 1-exp(-k(t-t0)), k and t0 both free, whole window
%   k_t0z   the same with the clock started AT the valley (t0 == 0)
%   k_band  log-linear on z in [0.15 0.90]: the top of the rise only
%   k_t50   ln2 / t50                       (model-free)
%   k_2080  ln4 / (t80-t20)                 (model-free, toe- and plateau-free)
%   k_mrt   1/MRT, MRT = INT (1-z) dt, finite-window corrected  (LINEAR in y)
%   k_area  1/INT (1-z) dt, uncorrected
%
% SIX TESTS, in increasing order of how much they should decide the choice:
%   T1 IDENTITY   time-dilate a real trace by a known factor c. A rate
%                 estimator MUST return c*k. This is the only test with a
%                 ground truth, and it is run on the actual data.
%   T2 DELAY      inject a pure 5 ms delay. A RATE must not move; only t0 may.
%   T3 SPAN       truncate the window to 60/100/280 ms. A rate must not care.
%   T4 NOISE      add measured-amplitude noise, 200 draws: the estimator's own
%                 variance, i.e. how much of the between-prep SD it manufactures.
%   T5 POOLING    mean-of-features vs features-of-mean-trace, on two
%                 alignments (ramp-end clock, and valley-aligned). This is the
%                 DataPooling failure, re-run per estimator.
%   T6 TARGET     between-preparation SD, within-run repeatability, and the
%                 2 mM/8 mM discriminability z of each estimator. An estimator
%                 that is stable but blind to ATP is worthless.
%
% Output -> results/metric_audit.mat, results/metric_audit.png
%
% Run from anywhere.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
addpath(here);
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

EST  = {'k_fit','k_t0z','k_band','k_t50','k_2080','k_mrt','k_area'};
LBL  = {'current','t0=0','band 15-90','t50','t20-t80','MRT (lin)','1/area'};
nE   = numel(EST);
DAYS = {'03/27','04/03','04/10'};
COND = {'8mM','2mM'};
nP = numel(DAYS); nC = numel(COND); nS = 5;

%% ================= gather every post-re-stretch window =================
% W{p,c,g} = struct('t',..,'y',..) -- the same window extractSlackAttributes
% uses: [velocity_segment(4) velocity_segment(5)], i.e. ramp end -> segment end.
W = cell(nP,nC,nS); SLs = nan(nP,nC,nS); Fref = nan(nP,1);
R = poolLoadRuns();
for p = 1:nP
    i8 = find(strcmp({R.day},DAYS{p}) & strcmp({R.cond},'8mM'));
    Fref(p) = R(i8).A.Fss22;
    for c = 1:nC
        S  = load(fullfile('data', slackFile(DAYS{p},COND{c})));
        t  = S.datatable(:,1); F = S.datatable(:,3); L = S.datatable(:,2);
        vt = S.velocitytable;
        seg = find(vt(:,2) < -1)';
        for g = 1:min(nS,numel(seg))
            vs = vt(seg(g):seg(g)+4, 1);
            w  = t >= vs(4) & t <= vs(5);
            W{p,c,g} = struct('t', t(w)-vs(4), 'y', F(w), 'vs', vs);
            wq = t >= vs(2) & t <= vs(3);
            SLs(p,c,g) = median(L(wq));           % the slack length of this segment
        end
    end
end

%% ---- baseline values, every estimator x every window ------------------
K = nan(nP,nC,nS,nE); AUX = cell(nP,nC,nS);
for p = 1:nP, for c = 1:nC, for g = 1:nS
    if isempty(W{p,c,g}); continue; end
    E = rsEstimators(W{p,c,g}.t, W{p,c,g}.y, true);
    AUX{p,c,g} = E;
    for e = 1:nE; K(p,c,g,e) = E.(EST{e}); end
end, end, end

fprintf('=== baseline: each estimator, 8 mM, mean over 5 segments (1/s) ===\n');
fprintf('%-12s', 'estimator'); fprintf('%9s', DAYS{:}); fprintf('%9s %8s\n','mean','CV%');
for e = 1:nE
    v = squeeze(mean(K(:,1,:,e),3,'omitnan'));
    fprintf('%-12s', LBL{e}); fprintf('%9.1f', v);
    fprintf('%9.1f %7.1f%%\n', mean(v,'omitnan'), 100*std(v,'omitnan')/mean(v,'omitnan'));
end

%% ================= T1  IDENTITY: known time dilation ==================
% t -> t/c compresses the recovery, so every rate must come back c times
% larger. Resampled onto the ORIGINAL grid so the sampling is not the cue.
CF = [0.6 0.8 1.25 1.6667];
Kd = nan(nP,nC,nS,nE,numel(CF));
for p = 1:nP, for c = 1:nC, for g = 1:nS
    if isempty(W{p,c,g}); continue; end
    for j = 1:numel(CF)
        [td, yd] = dilate(W{p,c,g}.t, W{p,c,g}.y, CF(j));
        Ed = rsEstimators(td, yd, true);
        for e = 1:nE; Kd(p,c,g,e,j) = Ed.(EST{e}); end
    end
end, end, end
fprintf('\n=== T1 IDENTITY: time-dilate the real trace by c, recover c ===\n');
fprintf('(median |k(c)/k(1)/c - 1| over all 30 windows; 0 = exact)\n');
fprintf('%-12s', 'estimator'); fprintf('%8s', string(CF)); fprintf('%9s\n','worst');
T1 = nan(nE,numel(CF));
for e = 1:nE
    for j = 1:numel(CF)
        r = Kd(:,:,:,e,j)./K(:,:,:,e)/CF(j) - 1;
        T1(e,j) = median(abs(r(:)),'omitnan');
    end
    fprintf('%-12s', LBL{e}); fprintf('%7.1f%%', 100*T1(e,:));
    fprintf('%8.1f%%\n', 100*max(T1(e,:)));
end

%% ================= T2  DELAY: a pure 5 ms shift =======================
% Insert 5 ms of flat baseline before the rise. The rate is unchanged by
% construction; anything that moves is charging delay to the rate.
Kdl = nan(nP,nC,nS,nE);
for p = 1:nP, for c = 1:nC, for g = 1:nS
    if isempty(W{p,c,g}); continue; end
    [td, yd] = injectDelay(W{p,c,g}.t, W{p,c,g}.y, 0.005);
    Ed = rsEstimators(td, yd, true);
    for e = 1:nE; Kdl(p,c,g,e) = Ed.(EST{e}); end
end, end, end
fprintf('\n=== T2 DELAY: inject a pure +5 ms onset delay ===\n');
fprintf('(median relative change in k; 0 = the rate ignored the delay)\n');
T2 = nan(nE,1);
for e = 1:nE
    r = Kdl(:,:,:,e)./K(:,:,:,e) - 1;
    T2(e) = median(r(:),'omitnan');
    fprintf('%-12s %+7.1f%%\n', LBL{e}, 100*T2(e));
end

%% ================= T3  SPAN ===========================================
SPAN = [0.060 0.100 0.280];
Ks = nan(nP,nC,nS,nE,numel(SPAN));
for p = 1:nP, for c = 1:nC, for g = 1:nS
    if isempty(W{p,c,g}); continue; end
    for j = 1:numel(SPAN)
        m = W{p,c,g}.t <= SPAN(j);
        if nnz(m) < 25; continue; end
        Ed = rsEstimators(W{p,c,g}.t(m), W{p,c,g}.y(m), true);
        for e = 1:nE; Ks(p,c,g,e,j) = Ed.(EST{e}); end
    end
end, end, end
fprintf('\n=== T3 SPAN: truncate the window ===\n');
fprintf('%-12s', 'estimator'); fprintf('%9s', "T="+string(1000*SPAN)+"ms"); fprintf('\n');
T3 = nan(nE,numel(SPAN));
for e = 1:nE
    for j = 1:numel(SPAN)
        r = Ks(:,:,:,e,j)./K(:,:,:,e) - 1;
        T3(e,j) = median(r(:),'omitnan');
    end
    fprintf('%-12s', LBL{e}); fprintf('%+8.1f%%', 100*T3(e,:)); fprintf('\n');
end

%% ================= T4  NOISE ==========================================
% Measured noise amplitude: SD of the residual about a 5 ms movmedian, taken
% on the plateau of every window. 200 draws on 03/27 8 mM segment 1.
nz = [];
for p = 1:nP, for c = 1:nC, for g = 1:nS
    if isempty(W{p,c,g}); continue; end
    yv = W{p,c,g}.y; tv = W{p,c,g}.t;
    m = tv > 0.15;
    if nnz(m) > 20
        nz(end+1) = std(yv(m) - smoothdata(yv(m),'movmedian',25)); %#ok<SAGROW>
    end
end, end, end
sigN = median(nz);
fprintf('\n=== T4 NOISE: sigma = %.3f kPa (measured on the plateaus) ===\n', sigN);
rng(7);
T4 = nan(nE,1); nDraw = 200;
Kn = nan(nDraw,nE);
base = W{1,1,1};
for d = 1:nDraw
    Ed = rsEstimators(base.t, base.y + sigN*randn(size(base.y)), true);
    for e = 1:nE; Kn(d,e) = Ed.(EST{e}); end
end
fprintf('%-12s %10s %10s\n','estimator','CV_noise','bias');
for e = 1:nE
    T4(e) = std(Kn(:,e),'omitnan')/mean(Kn(:,e),'omitnan');
    fprintf('%-12s %9.1f%% %+9.1f%%\n', LBL{e}, 100*T4(e), ...
            100*(mean(Kn(:,e),'omitnan')/K(1,1,1,e)-1));
end

%% ================= T5  POOLING ========================================
% Two alignments. (a) ramp-end clock -- what DataPooling did. (b) valley
% clock -- each trace shifted so its own located minimum sits at t=0, which
% is the obvious repair if the bias is onset jitter.
fprintf('\n=== T5 POOLING: mean-of-features vs features-of-mean-trace ===\n');
TG  = (0 : 2e-4 : 0.280)';
T5  = nan(nE,2);
Kg  = nan(nE,2,nS);
MeanTr = nan(2,nS,numel(TG));
for al = 1:2
    for g = 1:nS
        V = nan(nP,numel(TG));
        for p = 1:nP
            if isempty(W{p,1,g}); continue; end
            tv = W{p,1,g}.t; yv = W{p,1,g}.y / Fref(p);     % normalised
            if al == 2
                yl = smoothdata(W{p,1,g}.y,'movmedian',max(3,round(0.005/median(diff(tv)))));
                nWin = nnz(tv <= tv(1)+0.080);
                [~, iB] = min(yl(1:nWin));
                tv = tv - tv(iB);
            end
            V(p,:) = interp1(tv, yv, TG, 'linear', NaN);
        end
        mv = mean(V,1,'omitnan') * mean(Fref);              % back on a kPa axis
        MeanTr(al,g,:) = mv;
        Eg = rsEstimators(TG, mv(:), true);
        for e = 1:nE; Kg(e,al,g) = Eg.(EST{e}); end
    end
end
fprintf('%-12s %10s %10s %10s %10s %8s\n','estimator','mean-of-f','ramp-clock','valley-clock', ...
        'gap ramp','gap vall');
for e = 1:nE
    mf = mean(K(:,1,:,e),'all','omitnan');
    g1 = mean(Kg(e,1,:),'omitnan'); g2 = mean(Kg(e,2,:),'omitnan');
    T5(e,:) = [g1/mf-1, g2/mf-1];
    fprintf('%-12s %10.1f %10.1f %12.1f %9.1f%% %8.1f%%\n', LBL{e}, mf, g1, g2, ...
            100*T5(e,1), 100*T5(e,2));
end

%% ================= T6  TARGET QUALITY =================================
% s_between : between-preparation relative SD of the per-prep mean (8 mM)
% s_within  : within-run repeatability = residual about a linear fit on
%             SLslack across the five segments, relative
% ratio     : pooled 2 mM / 8 mM, its between-prep SD, and z = |r-1|/sqrt(s^2+sr^2)
fprintf('\n=== T6 TARGET QUALITY ===\n');
fprintf('%-12s %9s %9s %8s %9s %9s %7s\n','estimator','s_betw','s_with','bio%', ...
        'ratio','s_ratio','z');
T6 = nan(nE,6);
for e = 1:nE
    pm = nan(nP,1); pw = nan(nP,1); rt = nan(nP,1);
    for p = 1:nP
        v8 = squeeze(K(p,1,:,e)); v2 = squeeze(K(p,2,:,e));
        pm(p) = mean(v8,'omitnan');
        sl = squeeze(SLs(p,1,:));
        gd = isfinite(v8) & isfinite(sl);
        if nnz(gd) >= 3
            pf = polyfit(sl(gd), v8(gd), 1);
            pw(p) = std(v8(gd) - polyval(pf, sl(gd))) / abs(pm(p));
        end
        rt(p) = mean(v2,'omitnan')/mean(v8,'omitnan');
    end
    sB = std(pm,'omitnan')/mean(pm,'omitnan');
    sW = mean(pw,'omitnan');
    sMeas = sW/sqrt(nS);                       % what reaches a 5-segment mean
    bio = max(0, 1 - (sMeas/sB)^2);
    r = mean(rt,'omitnan'); sr = std(rt,'omitnan')/abs(r);
    zz = abs(r-1)/sqrt(sB^2+sr^2);
    T6(e,:) = [sB sW bio r sr zz];
    fprintf('%-12s %8.1f%% %8.1f%% %7.0f%% %9.3f %8.1f%% %7.2f\n', LBL{e}, ...
            100*sB, 100*sW, 100*bio, r, 100*sr, zz);
end

%% ================= composite ==========================================
fprintf('\n=== VERDICT (lower is better on the first four) ===\n');
fprintf('%-12s %8s %8s %8s %8s %8s %8s %6s\n','estimator','identity','delay', ...
        'span','noise','pool','s_betw','z_ATP');
score = nan(nE,1);
for e = 1:nE
    score(e) = max(T1(e,:)) + abs(T2(e)) + max(abs(T3(e,:))) + T4(e) + abs(T5(e,1));
    fprintf('%-12s %7.1f%% %+7.1f%% %+7.1f%% %7.1f%% %+7.1f%% %7.1f%% %6.2f\n', LBL{e}, ...
        100*max(T1(e,:)), 100*T2(e), 100*max(abs(T3(e,:))), 100*T4(e), ...
        100*T5(e,1), 100*T6(e,1), T6(e,6));
end
[~, best] = min(score);
fprintf('\nlowest total distortion: %s\n', LBL{best});

save(fullfile(resDir,'metric_audit.mat'), 'K','EST','LBL','T1','T2','T3','T4','T5','T6', ...
     'CF','SPAN','sigN','Kn','Kg','MeanTr','TG','DAYS','COND','SLs','Fref','score');

%% ================= figure =============================================
f = figure('Position',[60 60 1280 720],'Color','w');
tiledlayout(f,2,3,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on
cols = lines(nP);
for p = 1:nP
    E = AUX{p,1,1}; wv = W{p,1,1};
    z = (wv.y - E.B)/E.A;
    plot(1000*wv.t, z, 'Color',[cols(p,:) .55], 'LineWidth',1.0, 'DisplayName',DAYS{p});
end
yline(0.15,'k:','HandleVisibility','off'); yline(0.90,'k:','HandleVisibility','off');
yline(0.5,'k--','HandleVisibility','off');
xlabel('time from ramp end (ms)'); ylabel('z = (F-B)/(P-B)');
xlim([0 120]); ylim([-0.15 1.2]); legend('Location','southeast','FontSize',8);
title('(a) the shape every estimator sees (seg 1, 8 mM)');

nexttile; hold on; box on
bar(100*[max(T1,[],2) abs(T2) max(abs(T3),[],2) T4 abs(T5(:,1))]);
set(gca,'XTick',1:nE,'XTickLabel',LBL,'XTickLabelRotation',35);
ylabel('distortion (%)'); legend({'identity','delay','span','noise','pooling'}, ...
    'Location','northeast','FontSize',8);
title('(b) the five distortions, per estimator');

nexttile; hold on; box on
for e = 1:nE
    plot(CF, 100*T1(e,:), 'o-','LineWidth',1.2,'DisplayName',LBL{e});
end
xlabel('imposed time-dilation factor c'); ylabel('|recovered/c - 1| (%)');
legend('Location','northwest','FontSize',7); title('(c) T1: the only test with ground truth');

nexttile; hold on; box on
b = bar(100*T5); set(gca,'XTick',1:nE,'XTickLabel',LBL,'XTickLabelRotation',35);
legend({'ramp-end clock','valley clock'},'Location','southwest','FontSize',8);
ylabel('features-of-mean vs mean-of-features (%)');
title('(d) T5: does trace averaging still bias it?');

nexttile; hold on; box on
scatter(100*T6(:,1), T6(:,6), 90, 'filled');
for e = 1:nE; text(100*T6(e,1)+0.3, T6(e,6), LBL{e}, 'FontSize',8); end
xlabel('between-preparation SD (%)'); ylabel('ATP discriminability z');
title('(e) T6: stable AND still sees ATP'); grid on

nexttile; hold on; box on
for p = 1:nP
    plot(squeeze(SLs(p,1,:)), squeeze(K(p,1,:,1)), 'o--','Color',cols(p,:), ...
         'LineWidth',1.1,'DisplayName',[DAYS{p} ' current']);
    plot(squeeze(SLs(p,1,:)), squeeze(K(p,1,:,best)), 's-','Color',cols(p,:), ...
         'LineWidth',1.4,'DisplayName',[DAYS{p} ' ' LBL{best}]);
end
xlabel('slack length (ML)'); ylabel('recovery rate (1/s)');
legend('Location','best','FontSize',6); title('(f) the amplitude dependence, both metrics');

exportgraphics(f, fullfile(resDir,'metric_audit.png'), 'Resolution',150);
fprintf('\nwrote %s\n', fullfile(resDir,'metric_audit.png'));

%% ---- helpers ----------------------------------------------------------
function [td, yd] = dilate(t, y, c)
%DILATE compress/stretch the time axis by 1/c, resample on the original grid.
    t0 = t(1);
    tsrc = t0 + (t - t0)/c;
    td = t;
    yd = interp1(tsrc, y, td, 'linear');
    keep = isfinite(yd);
    td = td(keep); yd = yd(keep);
end

function [td, yd] = injectDelay(t, y, dly)
%INJECTDELAY prepend DLY seconds of the first sample, keep the same grid.
    dt = median(diff(t));
    n  = round(dly/dt);
    yd = [repmat(y(1), n, 1); y(1:end-n)];
    td = t;
end

function f = slackFile(day, cond)
    tag = struct('x03_27','03_27_2026','x04_03','04_03_2026','x04_10','04_10_2026');
    f = sprintf('protocol_%s_%s_slack.mat', tag.(['x' strrep(day,'/','_')]), cond);
end
