% RunSegment5Check.m -- segment 5 is not the same manoeuvre as segments 1-4.
%
% CheckPoolingArtefact localised the whole "-24 % trace averaging" result to a
% SINGLE segment: on the pooled trace, segments 1-4 give rsK = 51/47/52/50
% (against a mean-of-features 53.5) and segment 5 gives 9.9. That is not a
% smearing bias, it is one broken number. Reading the velocity table explains
% why, and the explanation matters far beyond pooling:
%
%   segment   release to      re-stretch      re-stretch    window
%             (ML)            back to (ML)    velocity      vs(4)->vs(5)
%   1-4       1.020 .. 0.960  1.1005          ~3.9-6.8      0.279 s
%   5         0.940           1.0002          ~2.9          0.844 s
%
% The fifth segment re-stretches to a DIFFERENT FINAL LENGTH (ML 1.00, not
% 1.10), at a lower ramp velocity, and its window is three times longer
% because the protocol simply holds until t = 78 s afterwards.
%
% Three consequences, tested here:
%   C1 the pooled-trace failure is a splice hole: the DataPooling grid ends at
%      release+0.570 s but segment 5's window runs to release+1.16 s, so the
%      plateau that sets the fit amplitude is un-averaged data.
%   C2 the 5-segment mean of rsK is not a mean over ONE manoeuvre. The
%      "amplitude dependence" dk/dSLslack is dominated by the one segment that
%      also changes final length and ramp speed.
%   C3 the between-preparation spread of rsK is concentrated in segment 5.
%
% Output -> results/segment5.mat / .png

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root)); addpath(here);
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

days = {'03/27','04/03','04/10'}; COND = {'8mM','2mM'};
nP = 3; nC = 2; nS = 5;
EST = {'k_fit','k_t0z','k_band','k_t50','k_2080'};
LBL = {'current','t0=0','band 15-90','t50','t20-t80'};
nE  = numel(EST);

%% ---- C1: the protocol, segment by segment ----------------------------
fprintf('=== the protocol, segment by segment (03/27; identical on all days) ===\n');
S = load(fullfile('data','protocol_03_27_2026_8mM_slack.mat'));
vt = S.velocitytable; seg = find(vt(:,2) < -1)';
fprintf('%4s %10s %10s %10s %10s %10s\n','seg','slack ML','back to ML','d(ML)','v(ML/s)','window s');
PROT = nan(nS,5);
for g = 1:nS
    r = seg(g):seg(g)+4;
    Lslack = vt(r(2),4); Lback = vt(r(4),4);
    dur = vt(r(4),1) - vt(r(3),1);
    PROT(g,:) = [Lslack Lback Lback-Lslack (Lback-Lslack)/dur vt(r(5),1)-vt(r(4),1)];
    fprintf('%4d %10.4f %10.4f %10.4f %10.2f %10.3f\n', g, PROT(g,:));
end
fprintf('  -> segment 5 re-stretches to a DIFFERENT LENGTH at a DIFFERENT SPEED\n');
fprintf('     and its feature window is %.1fx longer.\n', PROT(5,5)/mean(PROT(1:4,5)));

%% ---- gather every window, both conditions ----------------------------
W = cell(nP,nC,nS); SLs = nan(nP,nC,nS); Fref = nan(nP,1);
R = poolLoadRuns();
for p = 1:nP
    Fref(p) = R(find(strcmp({R.day},days{p}) & strcmp({R.cond},'8mM'),1)).A.Fss22;
    for c = 1:nC
        Sd = load(fullfile('data', slackFile(days{p},COND{c})));
        t = Sd.datatable(:,1); F = Sd.datatable(:,3); L = Sd.datatable(:,2);
        v = Sd.velocitytable; sg = find(v(:,2) < -1)';
        for g = 1:nS
            vs = v(sg(g):sg(g)+4,1);
            w = t >= vs(4) & t <= vs(5);
            W{p,c,g} = struct('t',t(w)-vs(4),'y',F(w));
            SLs(p,c,g) = median(L(t >= vs(2) & t <= vs(3)));
        end
    end
end

% two window conventions: the protocol window, and a COMMON 279 ms window
K  = nan(nP,nC,nS,nE); Kc = K;
for p = 1:nP, for c = 1:nC, for g = 1:nS
    E  = rsEstimators(W{p,c,g}.t, W{p,c,g}.y, true);
    m  = W{p,c,g}.t <= 0.279;
    Ec = rsEstimators(W{p,c,g}.t(m), W{p,c,g}.y(m), true);
    for e = 1:nE; K(p,c,g,e) = E.(EST{e}); Kc(p,c,g,e) = Ec.(EST{e}); end
end, end, end

fprintf('\n=== does capping segment 5 at the common 279 ms window change it? ===\n');
fprintf('%-12s %12s %12s %10s\n','estimator','seg5 native','seg5 capped','change');
for e = 1:nE
    a = mean(K(:,1,5,e),'omitnan'); b = mean(Kc(:,1,5,e),'omitnan');
    fprintf('%-12s %12.1f %12.1f %+9.1f%%\n', LBL{e}, a, b, 100*(b/a-1));
end

%% ---- C2: the amplitude dependence, with and without segment 5 --------
fprintf('\n=== C2 dk/dSLslack (1/s per ML), 8 mM -- the slope the model must match ===\n');
fprintf('%-12s %22s %22s\n','estimator','segments 1-5','segments 1-4 only');
fprintf('%-12s %7s %7s %7s %7s %7s %7s\n','', days{:}, days{:});
SLP = nan(nE,nP,2);
for e = 1:nE
    for p = 1:nP
        y = squeeze(Kc(p,1,:,e)); x = squeeze(SLs(p,1,:));
        gd = isfinite(y) & isfinite(x);
        pf = polyfit(x(gd), y(gd), 1); SLP(e,p,1) = pf(1);
        gd4 = gd; gd4(5) = false;
        pf = polyfit(x(gd4), y(gd4), 1); SLP(e,p,2) = pf(1);
    end
    fprintf('%-12s %7.0f %7.0f %7.0f %7.0f %7.0f %7.0f\n', LBL{e}, ...
            SLP(e,:,1), SLP(e,:,2));
end
fprintf('\n%-12s %14s %14s %10s\n','estimator','mean 1-5','mean 1-4','shrinks by');
for e = 1:nE
    a = mean(SLP(e,:,1)); b = mean(SLP(e,:,2));
    fprintf('%-12s %14.0f %14.0f %9.0f%%\n', LBL{e}, a, b, 100*(1-abs(b)/abs(a)));
end
fprintf(['  A negative slope = the rate is FASTER after a bigger slack (smaller\n' ...
         '  SLslack). RestretchRecoveryFit measured data -126 and model +714 on\n' ...
         '  the 5-segment version of this slope.\n']);

%% ---- C3: where the between-preparation spread lives ------------------
fprintf('\n=== C3 between-preparation relative SD, per segment (8 mM) ===\n');
fprintf('%-12s %8s %8s %8s %8s %8s | %9s %9s\n','estimator','s1','s2','s3','s4','s5', ...
        'mean 1-5','mean 1-4');
SD = nan(nE,nS);
for e = 1:nE
    for g = 1:nS
        v = squeeze(Kc(:,1,g,e));
        SD(e,g) = std(v,'omitnan')/mean(v,'omitnan');
    end
    m5 = std(mean(Kc(:,1,:,e),3,'omitnan'))/mean(mean(Kc(:,1,:,e),3,'omitnan'));
    m4 = std(mean(Kc(:,1,1:4,e),3,'omitnan'))/mean(mean(Kc(:,1,1:4,e),3,'omitnan'));
    fprintf('%-12s %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%% | %8.1f%% %8.1f%%\n', ...
            LBL{e}, 100*SD(e,:), 100*m5, 100*m4);
end

%% ---- and the ATP effect, with and without segment 5 -----------------
fprintf('\n=== the 2 mM / 8 mM ratio, with and without segment 5 ===\n');
fprintf('%-12s %10s %9s %7s | %10s %9s %7s\n','estimator','r (1-5)','sd','z','r (1-4)','sd','z');
ATP = nan(nE,6);
for e = 1:nE
    r5 = nan(nP,1); r4 = nan(nP,1);
    for p = 1:nP
        r5(p) = mean(Kc(p,2,:,e),'omitnan')/mean(Kc(p,1,:,e),'omitnan');
        r4(p) = mean(Kc(p,2,1:4,e),'omitnan')/mean(Kc(p,1,1:4,e),'omitnan');
    end
    s5b = std(mean(Kc(:,1,:,e),3,'omitnan'))/mean(mean(Kc(:,1,:,e),3,'omitnan'));
    s4b = std(mean(Kc(:,1,1:4,e),3,'omitnan'))/mean(mean(Kc(:,1,1:4,e),3,'omitnan'));
    z5 = abs(mean(r5)-1)/sqrt(s5b^2 + (std(r5)/mean(r5))^2);
    z4 = abs(mean(r4)-1)/sqrt(s4b^2 + (std(r4)/mean(r4))^2);
    ATP(e,:) = [mean(r5) std(r5)/mean(r5) z5 mean(r4) std(r4)/mean(r4) z4];
    fprintf('%-12s %10.3f %8.1f%% %7.2f | %10.3f %8.1f%% %7.2f\n', LBL{e}, ...
            ATP(e,1), 100*ATP(e,2), ATP(e,3), ATP(e,4), 100*ATP(e,5), ATP(e,6));
end

save(fullfile(resDir,'segment5.mat'),'PROT','K','Kc','SLs','SLP','SD','ATP','EST','LBL','days');

%% ---- figure ---------------------------------------------------------
f = figure('Position',[60 60 1250 460],'Color','w');
tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');
cols = lines(nP);

nexttile; hold on; box on
for p = 1:nP
    for g = 1:nS
        E = rsEstimators(W{p,1,g}.t, W{p,1,g}.y, true);
        m = W{p,1,g}.t <= 0.15;
        st = '-'; lw = 1.0;
        if g == 5; st = '-'; lw = 2.2; end
        plot(1000*W{p,1,g}.t(m), (W{p,1,g}.y(m)-E.B)/E.A, st, 'LineWidth',lw, ...
             'Color',[cols(p,:) 0.5], 'HandleVisibility','off');
    end
end
xlabel('time from ramp end (ms)'); ylabel('z'); ylim([-0.2 1.3]); xlim([0 150]);
title('(a) all 15 windows, 8 mM (thick = segment 5)');

nexttile; hold on; box on
for p = 1:nP
    plot(squeeze(SLs(p,1,:)), squeeze(Kc(p,1,:,1)), 'o--','Color',cols(p,:), ...
         'LineWidth',1.2,'DisplayName',[days{p} ' current']);
    plot(squeeze(SLs(p,1,:)), squeeze(Kc(p,1,:,4)), 's-','Color',cols(p,:), ...
         'LineWidth',1.6,'DisplayName',[days{p} ' t50']);
end
xline(1.89,'k:','HandleVisibility','off');
text(1.885,80,'seg 5','Rotation',90,'FontSize',8);
xlabel('slack length (ML)'); ylabel('rate (1/s)');
legend('Location','northeast','FontSize',6);
title('(b) the amplitude dependence; segment 5 is the lever');

nexttile; hold on; box on
bar(100*SD'); set(gca,'XTick',1:nS,'XTickLabel',{'s1','s2','s3','s4','s5'});
ylabel('between-preparation SD (%)'); legend(LBL,'Location','northwest','FontSize',7);
title('(c) the spread lives in segment 5');
exportgraphics(f, fullfile(resDir,'segment5.png'), 'Resolution',150);
fprintf('\nwrote %s\n', fullfile(resDir,'segment5.png'));

function f = slackFile(day, cond)
    tag = struct('x03_27','03_27_2026','x04_03','04_03_2026','x04_10','04_10_2026');
    f = sprintf('protocol_%s_%s_slack.mat', tag.(['x' strrep(day,'/','_')]), cond);
end
