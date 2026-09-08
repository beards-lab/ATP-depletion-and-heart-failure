% RunPoolingComparison.m -- how to combine three preparations into one target,
% and how much of the residual spread is noise rather than biology.
%
% Three questions the normalisation ladder does not answer:
%
% 1. TIMECOURSE. The ladder scores extracted features. Does the same
%    normalisation also collapse the raw traces, and where along the trace does
%    the agreement break down?
%
% 2. POOL WHERE? Two ways to make one target out of three preparations:
%       mean-of-features   extract per preparation, then average the numbers
%       features-of-mean   average the traces, then extract once
%    Feature extraction is NON-linear (exponential fits, peak picking, argmax),
%    so these are not the same number. Averaging traces that differ in t0 and in
%    ktr smears the recovery and biases the fitted rate. This measures the gap.
%
% 3. IS THE SPREAD NOISE OR BIOLOGY? Each run contains five slack segments.
%    Regress a feature on SLslack within a run and the residual is that run's
%    own repeatability -- measurement noise. What is left between preparations
%    after removing it is genuine preparation-to-preparation biology. The split
%    says whether the next experiment should be MORE PREPARATIONS or a BETTER
%    MEASUREMENT, which is not a question the pooled SD alone can answer.
%
% Output -> results/pooling_comparison.mat, results/pooling_traces.png

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

R    = poolLoadRuns();
days = unique({R.day},'stable'); nP = numel(days);
COND = {'8mM','2mM'};
nSeg = 5;

% relative-time grid, one slack cycle: release at 0, restretch at ~0.30
TG = (-0.020 : 2e-4 : 0.570)';   nT = numel(TG);

%% ---- 1. aligned traces -------------------------------------------------
% Tr(p,c,seg,:) = force on the common grid; Fref(p,c) = the Fss22 anchor
Tr   = nan(nP,numel(COND),nSeg,nT);
Fref = nan(nP,numel(COND));
for p = 1:nP
    for c = 1:numel(COND)
        i  = find(strcmp({R.day},days{p}) & strcmp({R.cond},COND{c}));
        S  = load(fullfile('data', slackFile(days{p},COND{c})));
        t  = S.datatable(:,1); F = S.datatable(:,3);
        vt = S.velocitytable;
        rel = find(vt(:,2) < -20);                    % the five slack releases
        for g = 1:min(nSeg,numel(rel))
            Tr(p,c,g,:) = interp1(t - vt(rel(g),1), F, TG, 'linear', NaN);
        end
        Fref(p,c) = R(i).A.Fss22;
    end
end

%% ---- pointwise between-preparation spread, raw vs normalised ----------
fprintf('=== pointwise between-preparation spread of the force trace ===\n');
fprintf('(median over the trace of the across-prep relative SD, 8 mM)\n');
PHASE = {'release+hold  (0 - 0.29 s)', TG >= 0     & TG < 0.290
         'restretch     (0.29-0.32 s)', TG >= 0.290 & TG < 0.320
         'recovery      (0.32-0.57 s)', TG >= 0.320
         'WHOLE TRACE',                 TG >  -0.02};
fprintf('%-28s %10s %10s %10s\n','phase','raw kPa','/Fss22','MoR');
Sraw = nan(nSeg,nT); Snrm = Sraw; Smor = Sraw;
for g = 1:nSeg
    V = squeeze(Tr(:,1,g,:));                          % nP x nT
    Sraw(g,:) = std(V,0,1,'omitnan') ./ abs(mean(V,1,'omitnan'));
    Vn = V ./ Fref(:,1);
    Snrm(g,:) = std(Vn,0,1,'omitnan') ./ abs(mean(Vn,1,'omitnan'));
    % median-of-ratios on the trace itself: one factor per prep, from the whole
    % trace, held out nowhere -- shown for reference only (it is fitted to the
    % very thing it is then scored on, so it is a lower bound, not a result)
    Lg = log(max(V,0.1)); a = exp(median(Lg - mean(Lg,1,'omitnan'), 2, 'omitnan'));
    Vm = V ./ (a/exp(mean(log(a))));
    Smor(g,:) = std(Vm,0,1,'omitnan') ./ abs(mean(Vm,1,'omitnan'));
end
for j = 1:size(PHASE,1)
    m = PHASE{j,2};
    fprintf('%-28s %9.1f%% %9.1f%% %9.1f%%\n', PHASE{j,1}, ...
        100*median(Sraw(:,m),'all','omitnan'), ...
        100*median(Snrm(:,m),'all','omitnan'), ...
        100*median(Smor(:,m),'all','omitnan'));
end

%% ---- 2. mean-of-features vs features-of-mean-trace ---------------------
% Build one grand-average trace per condition (normalised traces, averaged on
% the common grid), splice it back onto a canonical velocity table, and run the
% SAME extraction that produced the per-preparation features.
fprintf('\n=== pooling: mean-of-features vs features-of-mean-trace (8 mM) ===\n');
Sm = load(fullfile('data', slackFile(days{1},'8mM')));
vt0 = Sm.velocitytable; t0v = Sm.datatable(:,1); L0 = Sm.datatable(:,2);
rel0 = find(vt0(:,2) < -20);

% normalised grand-average force, spliced back onto 03/27's own time base so
% the canonical velocity table stays valid.
% NOTE the re-scaling by Fbar0: extractSlackAttributes carries ABSOLUTE
% thresholds in kPa (FORCE_THRESHOLD_KPA = 10 selects the active window), so a
% trace in F/Fss22 units silently returns NaN for every fitted quantity. The
% pooled trace is put back on a kPa axis at the mean anchor level, and the
% force-dimension features are divided by that same number afterwards.
Fbar0 = mean(Fref(:,1));
Fbar  = squeeze(mean(Tr(:,1,:,:) ./ Fref(:,1), 1, 'omitnan'));   % nSeg x nT
Fgrand = Sm.datatable(:,3) * (Fbar0/Fref(1,1));                  % start from 03/27
for g = 1:nSeg
    tt = t0v - vt0(rel0(g),1);
    m  = tt >= TG(1) & tt <= TG(end);
    Fgrand(m) = Fbar0 * interp1(TG, Fbar(g,:), tt(m), 'linear', 'extrap');
end
fPool = extractSlackAttributes(t0v, Fgrand, L0, vt0, [], [], false, true);

NM = {'A','Am','peak1_y','vall_y','peak2','steady','vall2_dy','rsA','ktr','rsK','t0','rsT0','rsR2'};
isFf = ismember(NM, {'A','Am','peak1_y','vall_y','peak2','steady','vall2_dy','rsA'});
fprintf('%-12s %12s %12s %9s %10s\n','feature','mean-of-feat','feat-of-mean','gap','gap/s');
gap = nan(1,numel(NM));
for k = 1:numel(NM)
    per = nan(nP,nSeg);
    for p = 1:nP
        i = find(strcmp({R.day},days{p}) & strcmp({R.cond},'8mM'));
        v = double(R(i).feat.(NM{k})(:))'; v = v(1:min(nSeg,numel(v)));
        if isFf(k); v = v / Fref(p,1); end
        per(p,1:numel(v)) = v;
    end
    mof = mean(per(:),'omitnan');
    fom = mean(double(fPool.(NM{k})),'omitnan');
    if isFf(k); fom = fom / Fbar0; end       % pooled trace lives on a kPa axis
    sK  = mean(std(per,0,1,'omitnan')./abs(mean(per,1,'omitnan')),'omitnan');
    gap(k) = (fom-mof)/abs(mof);
    fprintf('%-12s %12.4g %12.4g %8.1f%% %9.2f\n', NM{k}, mof, fom, 100*gap(k), abs(gap(k))/sK);
end
fprintf('  gap/s = the pooling choice expressed in between-preparation SDs.\n');
fprintf('  |gap/s| > 1 means the two pooling routes disagree by more than the\n');
fprintf('  data can distinguish, i.e. the choice is NOT cosmetic for that feature.\n');

%% ---- 3. noise vs biology ----------------------------------------------
% Within a run the five segments differ only by SLslack (by design), so the
% residual of a linear fit on SLslack is the run's own repeatability. A
% one-way random-effects split then gives
%     s_between^2  =  s_bio^2 + s_meas^2 / nSeg
fprintf('\n=== is the residual spread noise or biology? (8 mM) ===\n');
fprintf('s_meas is the WITHIN-run repeatability; only s_meas/sqrt(%d) of it reaches\n', nSeg);
fprintf('the between-prep SD, because each feature is already a 5-segment mean.\n');
fprintf('%-20s %9s %9s %9s %9s %8s\n','feature','s_total','s_meas','contrib','s_bio','bio share');
VD = struct('name',{},'stot',{},'smeas',{},'sbio',{});
for k = 1:numel(NM)
    per = nan(nP,nSeg); wres = nan(nP,1);
    for p = 1:nP
        i = find(strcmp({R.day},days{p}) & strcmp({R.cond},'8mM'));
        v = double(R(i).feat.(NM{k})(:))'; v = v(1:min(nSeg,numel(v)));
        sl = double(R(i).feat.SLslack(:))'; sl = sl(1:numel(v));
        if isFf(k); v = v / Fref(p,1); end
        per(p,1:numel(v)) = v;
        g = isfinite(v) & isfinite(sl);
        if sum(g) >= 4
            r = v(g) - polyval(polyfit(sl(g),v(g),1), sl(g));
            wres(p) = std(r) * sqrt((sum(g)-1)/(sum(g)-2));   % dof for the fit
        end
    end
    mu    = mean(per(:),'omitnan');
    stot  = mean(std(per,0,1,'omitnan')./abs(mean(per,1,'omitnan')),'omitnan');
    smeas = mean(wres,'omitnan')/abs(mu);
    sbio  = sqrt(max(stot^2 - smeas^2/nSeg, 0));
    fprintf('%-20s %8.1f%% %8.1f%% %8.1f%% %8.1f%% %7.0f%%\n', NM{k}, 100*stot, ...
            100*smeas, 100*smeas/sqrt(nSeg), 100*sbio, 100*sbio^2/max(stot^2,eps));
    VD(end+1) = struct('name',NM{k},'stot',stot,'smeas',smeas,'sbio',sbio); %#ok<SAGROW>
end
fprintf(['  bio share near 100%%  -> the preparations genuinely differ; only MORE\n' ...
         '                          PREPARATIONS shrink the error bar.\n' ...
         '  bio share well below -> the segment-to-segment scatter dominates; a\n' ...
         '                          cleaner extraction or more segments would help.\n']);

%% ---- 4. what is LEFT after the scale? ----------------------------------
% The trace panel says the residual is not uniform: after dividing by Fss22 the
% preparations agree to ~2 % during the post-restretch recovery (at ML 1.10) and
% still differ by ~11 % during the slack HOLD (ML 0.94-1.02). A pure scale
% cannot do that -- a scale is the same at every length. The natural second
% nuisance is therefore a LENGTH OFFSET: the same lesion RundownCorrection
% found (series-elastic creep = a horizontal shift of the length-tension curve,
% d = 0.044-0.060 ML), differing between preparations.
%
% Test: `A` (the recovery asymptote at the slack length) against SLslack, one
% line per preparation. If the three are horizontal translates, one d_p per
% preparation removes the remaining spread; if they differ in slope, it does not.
fprintf('\n=== what is left after the scale: a length offset? (8 mM) ===\n');
LNM  = {'A','Am','vall_y','peak1_y','peak2'};
Lmat = nan(nP,nSeg,numel(LNM)); SLm = nan(nP,nSeg);
for p = 1:nP
    i = find(strcmp({R.day},days{p}) & strcmp({R.cond},'8mM'));
    SLm(p,:) = double(R(i).feat.SLslack(1:nSeg))';
    for k = 1:numel(LNM)
        Lmat(p,:,k) = double(R(i).feat.(LNM{k})(1:nSeg))'/Fref(p,1);
    end
end
% d_p is fitted on `A` ALONE, gauge sum(d) = 0, then applied unchanged to the
% other four features. Those four are a genuine held-out test: if the residual
% really is a shared length offset, their spread must fall too, with no further
% fitting. If it only fell for `A`, the offset would just be absorbing A's noise.
sl0 = linspace(max(min(SLm,[],2)), min(max(SLm,[],2)), 20);
obj = @(d) spreadAt(SLm, Lmat(:,:,1), [0; d(:)] - mean([0; d(:)]), sl0);
dOpt = fminsearch(obj, [0 0], optimset('Display','off','TolX',1e-6));
dOpt = [0; dOpt(:)]; dOpt = dOpt - mean(dOpt);
fprintf('  fitted on `A` only:  d = %s um SL  (gauge: sum = 0)\n', strtrim(sprintf('%+.3f ', dOpt)));
fprintf('  RundownCorrection''s series-elastic creep shift, for scale: 0.088-0.120 um\n\n');
fprintf('%-12s %14s %14s %10s\n','feature','scale only','+ length offset','fitted?');
for k = 1:numel(LNM)
    s0 = spreadAt(SLm, Lmat(:,:,k), zeros(nP,1), sl0);
    s1 = spreadAt(SLm, Lmat(:,:,k), dOpt,        sl0);
    fprintf('%-12s %13.1f%% %13.1f%% %10s\n', LNM{k}, 100*s0, 100*s1, ...
            tern(k==1,'FITTED','held out'));
end
fprintf('  (spread is read off each preparation''s own five points, interpolated\n');
fprintf('   onto a common SL grid after the shift -- no smoothing model.)\n');

%% ---- figure ------------------------------------------------------------
f = figure('Position',[80 60 1180 760],'Color','w'); tl = tiledlayout(f,2,2,'Padding','compact');

nexttile; hold on; box on
cols = lines(nP);
for p = 1:nP
    plot(TG, squeeze(Tr(p,1,1,:))/Fref(p,1), 'Color',cols(p,:), 'LineWidth',1.2);
end
plot(TG, Fbar(1,:), 'k-', 'LineWidth',2);
xlabel('time from slack release (s)'); ylabel('F / F_{ss22}');
title('segment 1, 8 mM: normalised traces + grand average');
legend([days {'pooled'}], 'Location','southeast'); xlim([-0.02 0.57]);

nexttile; hold on; box on
plot(TG, 100*median(Sraw,1,'omitnan'), 'LineWidth',1.5);
plot(TG, 100*median(Snrm,1,'omitnan'), 'LineWidth',1.5);
plot(TG, 100*median(Smor,1,'omitnan'), 'LineWidth',1.5);
xlabel('time from slack release (s)'); ylabel('between-prep rel. SD (%)');
title('where along the trace do the preparations disagree?');
legend({'raw kPa','/ F_{ss22}','MoR (lower bound)'},'Location','northeast');
xlim([-0.02 0.57]); ylim([0 45]);

nexttile; hold on; box on
b = bar(100*[[VD.stot]; [VD.sbio]]');
b(1).FaceColor = [.80 .80 .85]; b(2).FaceColor = [.20 .35 .60];
set(gca,'XTick',1:numel(NM),'XTickLabel',NM,'XTickLabelRotation',45);
ylabel('relative SD (%)'); title('total between-prep spread vs the biological part');
legend({'total','biology (noise removed)'},'Location','northwest');

nexttile; hold on; box on
bar(100*gap); yline(0,'k-');
set(gca,'XTick',1:numel(NM),'XTickLabel',NM,'XTickLabelRotation',45);
ylabel('(feat-of-mean - mean-of-feat) / mean-of-feat  (%)');
title('the pooling choice is not cosmetic for the fitted rates');

exportgraphics(f, fullfile(resDir,'pooling_traces.png'), 'Resolution',150);
save(fullfile(resDir,'pooling_comparison.mat'), 'TG','Tr','Fref','Sraw','Snrm','Smor', ...
     'NM','gap','VD','Fbar','fPool','days','COND');
fprintf('\nsaved -> %s\n', fullfile(resDir,'pooling_comparison.mat'));

%% ======================= helpers ========================================
function s = spreadAt(SLm, V, d, slg)
% Shift each preparation's length axis by d(p), read its own five points onto
% the common grid slg, and return the mean across-prep relative SD there.
    nP = size(SLm,1); Q = nan(nP, numel(slg));
    for p = 1:nP
        [x, j] = sort(SLm(p,:) + d(p));
        y = V(p,j); g = isfinite(x) & isfinite(y);
        if sum(g) < 2; continue; end
        Q(p,:) = interp1(x(g), y(g), slg, 'linear', NaN);
    end
    ok = all(isfinite(Q),1);
    if ~any(ok); s = NaN; return; end
    s = mean(std(Q(:,ok),0,1)./abs(mean(Q(:,ok),1)));
end

function o = tern(c,a,b); if c; o = a; else; o = b; end; end

function f = slackFile(day, cond)
    tag = struct('x03_27','03_27_2026','x04_03','04_03_2026','x04_10','04_10_2026');
    f = sprintf('protocol_%s_%s_slack.mat', tag.(['x' strrep(day,'/','_')]), cond);
end
