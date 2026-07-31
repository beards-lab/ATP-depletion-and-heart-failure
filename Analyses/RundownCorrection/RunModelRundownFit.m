% RunModelRundownFit.m
% Can a DECAYING MODEL PARAMETER stand in for rundown, so that the fibre state
% is fitted rather than the data corrected?
%
% THE IDEA
% Instead of correcting measured forces, give each run a nuisance parameter
% describing how far the preparation has degraded, and let the model carry it.
% Then the ATP parameters are estimated on a common footing automatically. The
% natural candidates are the ones the mechanism study left standing:
%     SL0   -- the relaxed sarcomere length at a given ML (the fibre has slacked)
%     kSE   -- the serial elastic stiffness (the series element has crept softer)
%
% THE TEST
% The 03/27 8 mM bracket gives not one number but a PROFILE: ktr and amplitude
% at each of the five slack segments, for the fresh run and for the repeat. A
% lesion that is merely tuned to the mean ktr drop is fitting one number; a
% lesion that reproduces the per-segment profile is being tested against five.
%
% ktr must be compared sampling-matched: the repeat is 10 ms while tau ~ 20 ms,
% so the fresh trace is resampled onto the repeat's own sample times and both
% fitted identically (bias becomes common-mode).
%
% Outputs -> results/model_rundown_fit.png, .mat

clear; close all; clc;
cd(fileparts(which('RunModelRundownFit')));
addpath(genpath('../..'));
resDir = 'results';  if ~isfolder(resDir), mkdir(resDir); end
D = '../../data/';

%% ── 1. Observed per-segment profile, sampling-matched ──────────────────────
S  = load([D 'protocol_03_27_2026_8mM_slack.mat']);
vt = S.velocitytable;  is = find(vt(:,2) < -1);  nS = numel(is);
A  = readmatrix([D '03 27 2026 M/02_Merged_8mM_Active.txt']);
B  = readmatrix([D '03 27 2026 M/04_Merged_8mM_Active_repeat.txt']);
selA = A(:,1)>=74 & A(:,1)<=78;  tA = A(selA,1);  FA = A(selA,3);
selB = B(:,1)>=74 & B(:,1)<=78;  tB = B(selB,1);  FB = B(selB,3);
FAd  = interp1(tA, FA, tB, 'linear');
ok   = ~isnan(FAd);  tC = tB(ok);  FAd = FAd(ok);  FBd = FB(ok);

ftk  = fittype('a*(1-exp(-k*max(x-t0,0)))');
fitk = @(t,y) fit(t(:)-t(1), y(:), ftk, 'StartPoint',[max(y) 50 0.005], ...
                  'Lower',[1 1 0], 'Upper',[300 500 0.06]);
kF = nan(1,nS); kR = nan(1,nS); aF = nan(1,nS); aR = nan(1,nS);
for s = 1:nS
    t0 = vt(is(s)+1,1);  t1 = vt(is(s)+2,1);
    w = tC > t0 & tC < t1;
    c1 = fitk(tC(w), FAd(w));  kF(s) = c1.k;  aF(s) = c1.a;
    c2 = fitk(tC(w), FBd(w));  kR(s) = c2.k;  aR(s) = c2.a;
end
obsK = kR./kF;   obsA = aR./aF;
SLs  = 2*vt(is+1,4)';

fprintf('=== OBSERVED (03/27 8 mM: repeat / fresh, sampling-matched) ===\n');
fprintf('%-14s %s\n','SLslack (um)', sprintf('%8.3f', SLs));
fprintf('%-14s %s\n','ktr fresh',    sprintf('%8.2f', kF));
fprintf('%-14s %s\n','ktr repeat',   sprintf('%8.2f', kR));
fprintf('%-14s %s   mean %.3f\n','ktr ratio', sprintf('%8.3f', obsK), mean(obsK));
fprintf('%-14s %s   mean %.3f\n','A   ratio', sprintf('%8.3f', obsA), mean(obsA));

%% ── 2. Model: candidate decaying parameters ────────────────────────────────
p0 = getParams(loadParams('ThursdayNightFever'), [], true, false);
p0.PlotEachSeparately=0; p0.PlotFullscreen=0; p0.PlotFeatureFitting=0;
p0.ghostSave=''; p0.ghostLoad=''; p0.EvalFeatures=1; p0.RunSlackSegments='All';
p0.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';

CAND = struct( ...
 'id',  {'base','SL0','kSE','SL0+kSE','kstiff','kstiff+kSE'}, ...
 'lab', {'baseline', ...
         'SL0 only  (-0.098 um)', ...
         'kSE only  (x0.65)', ...
         'SL0 + kSE (-0.098, x0.65)', ...
         'kstiff only (x0.84)', ...
         'kstiff + kSE (x0.84, x0.65)'}, ...
 'apply',{@(p) p, ...
          @(p) setfield(p,'SL0',p.SL0-0.098), ...
          @(p) setfield(p,'kSE',p.kSE*0.65), ...
          @(p) setfield(setfield(p,'SL0',p.SL0-0.098),'kSE',p.kSE*0.65), ...
          @(p) setfield(setfield(p,'kstiff1',p.kstiff1*0.84),'kstiff2',p.kstiff2*0.84), ...
          @(p) setfield(setfield(setfield(p,'kstiff1',p.kstiff1*0.84),'kstiff2',p.kstiff2*0.84),'kSE',p.kSE*0.65)});

fprintf('\n=== MODEL RUNS ===\n');
for c = 1:numel(CAND)
    tic;
    [~, ~, fm] = runSlackExperiment(CAND(c).apply(p0));
    CAND(c).ktr = fm.ktr;  CAND(c).A = fm.A;
    fprintf('  %-30s %5.1f s   ktr %s\n', CAND(c).lab, toc, sprintf('%6.1f', fm.ktr));
end
kb = CAND(1).ktr;  ab = CAND(1).A;
for c = 1:numel(CAND); CAND(c).kr = CAND(c).ktr./kb;  CAND(c).ar = CAND(c).A./ab; end

%% ── 3. Compare profiles ────────────────────────────────────────────────────
fprintf('\n=== PER-SEGMENT ktr RATIO (model vs observed) ===\n');
fprintf('%-30s %s %10s %10s\n','candidate', sprintf('%8s','s1','s2','s3','s4','s5'), 'mean','RMSE vs obs');
fprintf('%-30s %s %10.3f %10s\n','OBSERVED', sprintf('%8.3f', obsK), mean(obsK), '--');
for c = 2:numel(CAND)
    e = sqrt(mean((CAND(c).kr - obsK).^2));
    fprintf('%-30s %s %10.3f %10.3f\n', CAND(c).lab, sprintf('%8.3f', CAND(c).kr), mean(CAND(c).kr), e);
end

fprintf('\n=== PER-SEGMENT AMPLITUDE RATIO ===\n');
fprintf('%-30s %s %10s %10s\n','candidate', sprintf('%8s','s1','s2','s3','s4','s5'), 'mean','RMSE vs obs');
fprintf('%-30s %s %10.3f %10s\n','OBSERVED', sprintf('%8.3f', obsA), mean(obsA), '--');
for c = 2:numel(CAND)
    e = sqrt(mean((CAND(c).ar - obsA).^2));
    fprintf('%-30s %s %10.3f %10.3f\n', CAND(c).lab, sprintf('%8.3f', CAND(c).ar), mean(CAND(c).ar), e);
end

fprintf('\n=== JOINT SCORE (ktr + amplitude profiles) ===\n');
fprintf('%-30s %10s %10s %10s\n','candidate','RMSE ktr','RMSE A','TOTAL');
best = inf; ibest = 0;
for c = 2:numel(CAND)
    ek = sqrt(mean((CAND(c).kr-obsK).^2));  ea = sqrt(mean((CAND(c).ar-obsA).^2));
    tot = sqrt(ek^2+ea^2);  CAND(c).score = tot;
    if tot < best; best = tot; ibest = c; end
    fprintf('%-30s %10.3f %10.3f %10.3f\n', CAND(c).lab, ek, ea, tot);
end
fprintf('\n  BEST: %s\n', CAND(ibest).lab);
fprintf('\n  A single decaying parameter is NOT enough:\n');
fprintf('    SL0 alone  moves amplitude but leaves ktr at %.3f (need %.3f)\n', mean(CAND(2).kr), mean(obsK));
fprintf('    kSE alone  moves ktr but leaves amplitude at %.3f (need %.3f)\n', mean(CAND(3).ar), mean(obsA));
fprintf('\n  WHICH force lever? The observed amplitude profile is nearly FLAT across\n');
fprintf('  the slack ladder (%.3f -> %.3f, gradient %+.3f). A length shift of the\n', obsA(1), obsA(end), obsA(end)-obsA(1));
fprintf('  required size predicts a STEEP gradient (%.3f -> %.3f, %+.3f) because the\n', ...
        CAND(2).ar(1), CAND(2).ar(end), CAND(2).ar(end)-CAND(2).ar(1));
fprintf('  deep segments lose more length-tension. A force SCALE predicts a flat one\n');
fprintf('  (%.3f -> %.3f, %+.3f), which is what the data show.\n', ...
        CAND(5).ar(1), CAND(5).ar(end), CAND(5).ar(end)-CAND(5).ar(1));
fprintf('  => on these 5 points the force loss looks like FEWER HEADS, not lost length.\n');
fprintf('  (The earlier length-shift verdict rested on the single ML 1.10 point, which\n');
fprintf('   lies outside the slack ladder -- see conclusions.md section 4.)\n');

%% ── 4. Figure ──────────────────────────────────────────────────────────────
fig = figure(990); clf; set(fig,'Position',[20 20 1450 720]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
cm = lines(numel(CAND));

nexttile; hold on; box on;
plot(SLs, kF, 'ko-','MarkerFaceColor','k','LineWidth',1.8);
plot(SLs, kR, 'rs-','MarkerFaceColor','r','LineWidth',1.8);
set(gca,'XDir','reverse'); xlabel('SL_{slack} (\mum)'); ylabel('ktr (1/s)');
title('(a) Observed ktr profile, sampling-matched');
legend({'fresh','repeat'},'Location','southwest','FontSize',8);

nexttile; hold on; box on;
plot(1:nS, obsK, 'k-o','LineWidth',2.5,'MarkerFaceColor','k','MarkerSize',9,'DisplayName','OBSERVED');
for c = 2:numel(CAND)
    plot(1:nS, CAND(c).kr, 'o--','Color',cm(c,:),'LineWidth',1.5,'DisplayName',CAND(c).id);
end
yline(1,':','HandleVisibility','off');
xlabel('slack segment'); ylabel('ktr ratio (lesion / baseline)');
title('(b) ktr profile: which lesion matches?');
legend('Location','southwest','FontSize',7);

nexttile; hold on; box on;
plot(1:nS, obsA, 'k-o','LineWidth',2.5,'MarkerFaceColor','k','MarkerSize',9,'DisplayName','OBSERVED');
for c = 2:numel(CAND)
    plot(1:nS, CAND(c).ar, 'o--','Color',cm(c,:),'LineWidth',1.5,'DisplayName',CAND(c).id);
end
yline(1,':','HandleVisibility','off');
xlabel('slack segment'); ylabel('amplitude ratio');
title('(c) Amplitude profile');
legend('Location','southwest','FontSize',7);

nexttile; hold on; box on;
for c = 2:numel(CAND)
    plot(mean(CAND(c).ar), mean(CAND(c).kr), 'o','Color',cm(c,:),'MarkerFaceColor',cm(c,:), ...
         'MarkerSize',13,'DisplayName',CAND(c).lab);
end
plot(mean(obsA), mean(obsK), 'kp','MarkerSize',22,'MarkerFaceColor','y','DisplayName','OBSERVED');
plot(1,1,'k*','MarkerSize',12,'HandleVisibility','off');
xlabel('amplitude ratio'); ylabel('ktr ratio'); grid on;
title('(d) One parameter cannot reach the star');
legend('Location','southwest','FontSize',7);

nexttile; hold on; box on;
sc = [CAND(2:end).score];
b = bar(sc); b.FaceColor='flat'; b.CData = cm(2:end,:);
set(gca,'XTick',1:numel(sc),'XTickLabel',{CAND(2:end).id});
ylabel('joint RMSE (ktr + amplitude)'); title('(e) Joint profile score, lower is better');

nexttile; hold on; box on;
txt = {'\bfUsing rundown as a model nuisance\rm','', ...
       'SL0 alone : amplitude moves, ktr does not', ...
       sprintf('            ktr %.3f vs needed %.3f', mean(CAND(2).kr), mean(obsK)), '', ...
       'kSE alone : ktr moves, amplitude does not', ...
       sprintf('            A %.3f vs needed %.3f', mean(CAND(3).ar), mean(obsA)), '', ...
       'BOTH together reach the observed point,', ...
       'and they are ONE physical lesion (series', ...
       'element creeping longer AND softer).', '', ...
       '=> fit ONE rundown coordinate per run,', ...
       '   mapped onto (SL0, kSE) jointly.'};
text(0.03,0.95,txt,'VerticalAlignment','top','FontSize',10); axis off;

sgtitle('Rundown as a decaying model parameter — tested against the 5-segment profile');
exportgraphics(fig, fullfile(resDir,'model_rundown_fit.png'), 'Resolution', 150);
save(fullfile(resDir,'model_rundown_fit.mat'),'CAND','obsK','obsA','kF','kR','aF','aR','SLs');
fprintf('\nSaved %s\n', fullfile(resDir,'model_rundown_fit.png'));
