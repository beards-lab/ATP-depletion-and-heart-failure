% RunAmplitudeIsTheWeakLink.m -- every estimator failed the same two tests
% (identity at c = 0.6, span at 60 ms) in the same direction. That pattern says
% the failure is NOT in the rate formula but in the one thing they all share:
%
%     A = median of the last 15 % of the window
%
% Both failing tests move the END of the window: dilating by c = 0.6 pulls the
% plateau sample back to where z is still ~0.99, truncating to 60 ms pulls it
% back to z ~ 0.85. A is then too small and every rate comes out too fast.
%
% TEST: repeat T1 and T3 with A supplied from OUTSIDE (opts.Aext), taken from
% the untouched full-length trace. If the distortion collapses, the amplitude
% convention was the whole problem, and the fix is to source A from `steady`
% (which is measured separately, at the reference length, over the last 20 ms
% before the segment ends) rather than from the tail of the fit window.
%
% ALSO: how much does the crossing smoother matter for the model-free
% estimators' noise? t50/t20-t80 looked noisy in T4 (11-13 %), but T4 adds a
% full measurement noise on top of data that already contains one, so it is an
% upper bound. Scanned here.
%
% Output -> results/amplitude_weaklink.mat / .png

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root)); addpath(here);
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

EST = {'k_fit','k_t0z','k_band','k_t50','k_2080'};
LBL = {'current','t0=0','band 15-90','t50','t20-t80'};
nE = numel(EST);
days = {'03/27','04/03','04/10'}; COND = {'8mM','2mM'};
nP = 3; nC = 2; nS = 5;

%% ---- windows + the externally measured amplitude ---------------------
% Aext = steady - B, with `steady` the median of the last 20 ms before the
% segment ends (what extractSlackAttributes calls feats.steady) and B the
% located post-re-stretch minimum. Segment 5 re-stretches to ML 1.00 rather
% than 1.10, so its steady is its own -- this is per segment, not a constant.
W = cell(nP,nC,nS); Aext = nan(nP,nC,nS); SLs = nan(nP,nC,nS);
for p = 1:nP, for c = 1:nC
    S = load(fullfile('data', slackFile(days{p},COND{c})));
    t = S.datatable(:,1); F = S.datatable(:,3); L = S.datatable(:,2);
    v = S.velocitytable; sg = find(v(:,2) < -1)';
    for g = 1:nS
        vs = v(sg(g):sg(g)+4,1);
        w = t >= vs(4) & t <= vs(5);
        W{p,c,g} = struct('t',t(w)-vs(4),'y',F(w));
        steady = median(F(t >= vs(5)-0.02 & t <= vs(5)));
        E0 = rsEstimators(W{p,c,g}.t, W{p,c,g}.y, true);
        Aext(p,c,g) = steady - E0.B;
        SLs(p,c,g) = median(L(t >= vs(2) & t <= vs(3)));
    end
end, end

fprintf('=== the two amplitude conventions ===\n');
fprintf('%-28s %10s %10s %9s\n','','tail-of-window','steady - B','differ');
for g = 1:nS
    a = nan(nP,1); b = nan(nP,1);
    for p = 1:nP
        E0 = rsEstimators(W{p,1,g}.t, W{p,1,g}.y, true);
        a(p) = E0.A; b(p) = Aext(p,1,g);
    end
    fprintf('  segment %d (8 mM, kPa)%9s %10.2f %10.2f %8.1f%%\n', g, '', ...
            mean(a), mean(b), 100*(mean(b)/mean(a)-1));
end

%% ---- T1 identity, both conventions -----------------------------------
CF = [0.6 0.8 1.25 1.6667];
T1 = nan(nE,numel(CF),2);
for mode = 1:2
    for j = 1:numel(CF)
        r = nan(nP,nC,nS,nE);
        for p = 1:nP, for c = 1:nC, for g = 1:nS
            o = struct(); if mode == 2; o.Aext = Aext(p,c,g); end
            E0 = rsEstimators(W{p,c,g}.t, W{p,c,g}.y, true, o);
            [td,yd] = dilate(W{p,c,g}.t, W{p,c,g}.y, CF(j));
            Ed = rsEstimators(td, yd, true, o);
            for e = 1:nE
                r(p,c,g,e) = Ed.(EST{e})/E0.(EST{e})/CF(j) - 1;
            end
        end, end, end
        for e = 1:nE
            v = r(:,:,:,e); T1(e,j,mode) = median(abs(v(:)),'omitnan');
        end
    end
end
fprintf('\n=== T1 IDENTITY: median |recovered/c - 1| ===\n');
fprintf('%-12s %28s | %28s\n','','A = tail of window','A = steady - B (external)');
fprintf('%-12s', ''); fprintf('%7s', string(CF)); fprintf(' |'); fprintf('%7s', string(CF)); fprintf('\n');
for e = 1:nE
    fprintf('%-12s', LBL{e}); fprintf('%6.1f%%', 100*T1(e,:,1));
    fprintf(' |'); fprintf('%6.1f%%', 100*T1(e,:,2)); fprintf('\n');
end

%% ---- T3 span, both conventions ---------------------------------------
SPAN = [0.060 0.100 0.279];
T3 = nan(nE,numel(SPAN),2);
for mode = 1:2
    for j = 1:numel(SPAN)
        r = nan(nP,nC,nS,nE);
        for p = 1:nP, for c = 1:nC, for g = 1:nS
            o = struct(); if mode == 2; o.Aext = Aext(p,c,g); end
            E0 = rsEstimators(W{p,c,g}.t, W{p,c,g}.y, true, o);
            m = W{p,c,g}.t <= SPAN(j);
            if nnz(m) < 25; continue; end
            Ed = rsEstimators(W{p,c,g}.t(m), W{p,c,g}.y(m), true, o);
            for e = 1:nE; r(p,c,g,e) = Ed.(EST{e})/E0.(EST{e}) - 1; end
        end, end, end
        for e = 1:nE
            v = r(:,:,:,e); T3(e,j,mode) = median(v(:),'omitnan');
        end
    end
end
fprintf('\n=== T3 SPAN: median relative change when the window is truncated ===\n');
fprintf('%-12s %24s | %24s\n','','A = tail of window','A = steady - B');
fprintf('%-12s', ''); fprintf('%8s', "T"+string(1000*SPAN)); fprintf(' |');
fprintf('%8s', "T"+string(1000*SPAN)); fprintf('\n');
for e = 1:nE
    fprintf('%-12s', LBL{e}); fprintf('%+7.1f%%', 100*T3(e,:,1));
    fprintf(' |'); fprintf('%+7.1f%%', 100*T3(e,:,2)); fprintf('\n');
end

%% ---- crossing-smoother scan for the model-free estimators ------------
CW = [0.002 0.005 0.010 0.020];
rng(11); nDraw = 120;
% noise level, measured on the plateaus
nz = [];
for p = 1:nP, for c = 1:nC, for g = 1:nS
    yv = W{p,c,g}.y; tv = W{p,c,g}.t; m = tv > 0.15;
    if nnz(m) > 20; nz(end+1) = std(yv(m)-smoothdata(yv(m),'movmedian',25)); end %#ok<SAGROW>
end, end, end
sigN = median(nz);
fprintf('\n=== crossing smoother: noise CV and between-prep SD (sigma=%.2f kPa) ===\n', sigN);
fprintf('%-12s', 'crossWin'); fprintf('%9s', string(1000*CW)+"ms"); fprintf('   |');
fprintf('%9s', string(1000*CW)+"ms"); fprintf('\n');
fprintf('%-12s %36s %38s\n','', 'CV added by re-noising', 'between-preparation SD');
SC = nan(nE,numel(CW),2);
for w = 1:numel(CW)
    o = struct('crossWin',CW(w));
    % (a) noise CV on one window
    Kn = nan(nDraw,nE);
    base = W{1,1,1};
    for d = 1:nDraw
        Ed = rsEstimators(base.t, base.y + sigN*randn(size(base.y)), true, o);
        for e = 1:nE; Kn(d,e) = Ed.(EST{e}); end
    end
    % (b) between-prep SD of the 5-segment mean
    pm = nan(nP,nE);
    for p = 1:nP
        vv = nan(nS,nE);
        for g = 1:nS
            Ed = rsEstimators(W{p,1,g}.t, W{p,1,g}.y, true, o);
            for e = 1:nE; vv(g,e) = Ed.(EST{e}); end
        end
        pm(p,:) = mean(vv,1,'omitnan');
    end
    for e = 1:nE
        SC(e,w,1) = std(Kn(:,e),'omitnan')/mean(Kn(:,e),'omitnan');
        SC(e,w,2) = std(pm(:,e),'omitnan')/mean(pm(:,e),'omitnan');
    end
end
for e = 1:nE
    fprintf('%-12s', LBL{e}); fprintf('%8.1f%%', 100*SC(e,:,1)); fprintf('   |');
    fprintf('%8.1f%%', 100*SC(e,:,2)); fprintf('\n');
end

save(fullfile(resDir,'amplitude_weaklink.mat'),'T1','T3','SC','CW','CF','SPAN', ...
     'Aext','EST','LBL','sigN');

%% ---- figure ----------------------------------------------------------
f = figure('Position',[60 60 1250 430],'Color','w');
tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');
nexttile; hold on; box on
b = bar(100*[max(T1(:,:,1),[],2) max(T1(:,:,2),[],2)]);
set(gca,'XTick',1:nE,'XTickLabel',LBL,'XTickLabelRotation',30);
ylabel('worst |recovered/c - 1| (%)');
legend({'A from window tail','A from steady'},'Location','northeast','FontSize',8);
title('(a) T1 identity: the amplitude was the error');

nexttile; hold on; box on
bar(100*[abs(T3(:,1,1)) abs(T3(:,1,2))]);
set(gca,'XTick',1:nE,'XTickLabel',LBL,'XTickLabelRotation',30);
ylabel('|change| at a 60 ms window (%)');
legend({'A from window tail','A from steady'},'Location','northeast','FontSize',8);
title('(b) T3 span: same story');

nexttile; hold on; box on
for e = 1:nE
    plot(1000*CW, 100*SC(e,:,2), 'o-','LineWidth',1.3,'DisplayName',LBL{e});
end
xlabel('crossing smoother (ms)'); ylabel('between-preparation SD (%)');
legend('Location','northwest','FontSize',7); set(gca,'XScale','log');
title('(c) how much smoothing the crossings want');
exportgraphics(f, fullfile(resDir,'amplitude_weaklink.png'), 'Resolution',150);
fprintf('\nwrote %s\n', fullfile(resDir,'amplitude_weaklink.png'));

%% ---- helpers ---------------------------------------------------------
function [td, yd] = dilate(t, y, c)
    t0 = t(1); tsrc = t0 + (t - t0)/c; td = t;
    yd = interp1(tsrc, y, td, 'linear');
    keep = isfinite(yd); td = td(keep); yd = yd(keep);
end

function f = slackFile(day, cond)
    tag = struct('x03_27','03_27_2026','x04_03','04_03_2026','x04_10','04_10_2026');
    f = sprintf('protocol_%s_%s_slack.mat', tag.(['x' strrep(day,'/','_')]), cond);
end
