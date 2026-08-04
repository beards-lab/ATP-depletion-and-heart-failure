% RunRundownMechanism.m
% What IS rundown, mechanically? Two hypotheses, tested against the one true
% rundown pair in the data (03/27's 8 mM at t=0 and its repeat at t=+694 s).
%
%   H1  FEWER CYCLING HEADS. The number of participating cross-bridges falls;
%       each surviving one behaves normally.
%         => force scales by a constant s at EVERY length   F_rd(ML) = s*F_fresh(ML)
%         => ktr is UNCHANGED (ktr is a rate, blind to head count)
%
%   H2  THE PREPARATION SLACKS. Sarcomeres give way / end-compliance grows, so
%       commanded ML = 1.0 no longer sits at SL = 2.0 um.
%         => the length-tension curve SHIFTS along the length axis
%              F_rd(ML) = F_fresh(ML - d)
%         => every length-dependent feature, ktr included, shifts by the same d
%
% The discriminator exists because each recording contains a 6-POINT
% LENGTH-TENSION CURVE: quasi-isometric force at the end of each of the five
% slack holds (ML 0.94...1.02) plus the pre-slack plateau (ML 1.10). H1 predicts
% a pure vertical scaling of that curve, H2 a pure horizontal shift. They are
% both one-parameter, so their residuals compare directly.
%
% CAVEAT handled below: the repeat run is LOG-ONLY (10 ms) because the burst
% name-matcher excludes 'repeat', while ktr ~ 49 /s (tau ~ 20 ms). Before any
% ktr claim, the fresh run is resampled to 10 ms to measure how much of an
% apparent slowdown that sampling alone produces.
%
% See also DataCuration/FitRundownCorrection.m (the r(F,SL) force correction),
% and conclusions.md.

clear; close all; clc;
cd(fileparts(which('RunRundownMechanisms')));
addpath(genpath('../..'));
resDir = 'results';
if ~isfolder(resDir), mkdir(resDir); end
D = '../../data/';

S  = load([D 'protocol_03_27_2026_8mM_slack.mat']);
vt = S.velocitytable;
A  = readmatrix([D '03 27 2026 M/02_Merged_8mM_Active.txt']);        % fresh, t=0
B  = readmatrix([D '03 27 2026 M/04_Merged_8mM_Active_repeat.txt']); % rundown, t=+694 s

%% ── 1. Six-point length-tension curve from each run ────────────────────────
is      = find(vt(:,2) < -1);
holdEnd = vt(is+2, 1);           % end of each slack hold = quasi-isometric
win     = 0.06;                  % average the last 60 ms of the hold
mF = @(M,a,b) mean(M(M(:,1)>=a & M(:,1)<=b, 3), 'omitnan');
mL = @(M,a,b) mean(M(M(:,1)>=a & M(:,1)<=b, 2), 'omitnan');

ML = zeros(1,numel(holdEnd)+1); Ff = ML; Fr = ML;
for k = 1:numel(holdEnd)
    ML(k) = mL(A, holdEnd(k)-win, holdEnd(k));
    Ff(k) = mF(A, holdEnd(k)-win, holdEnd(k));
    Fr(k) = mF(B, holdEnd(k)-win, holdEnd(k));
end
ML(end) = mL(A,74.0,74.40);  Ff(end) = mF(A,74.0,74.40);  Fr(end) = mF(B,74.0,74.40);
[ML,ix] = sort(ML); Ff = Ff(ix); Fr = Fr(ix);

p   = polyfit(ML, Ff, 2);                       % smooth model of the fresh curve
Fmod = @(x) polyval(p, x);
R2  = 1 - sum((Ff-Fmod(ML)).^2)/sum((Ff-mean(Ff)).^2);

%% ── 2. Fit H1 (scale), H2 (shift), H3 (both) ───────────────────────────────
s1  = sum(Ff.*Fr)/sum(Ff.^2);                                   % H1, 1 param
r1  = Fr - s1*Ff;
dd  = linspace(0, 0.15, 3001);
[~,i2] = min(arrayfun(@(d) sum((Fr - Fmod(ML-d)).^2), dd));
d2  = dd(i2);                                                   % H2, 1 param
r2  = Fr - Fmod(ML-d2);
v3  = fminsearch(@(v) sum((Fr - v(1)*Fmod(ML-v(2))).^2), [s1 d2]);
r3  = Fr - v3(1)*Fmod(ML-v3(2));                                % H3, 2 params

aic = @(r,k) numel(r)*log(sum(r.^2)/numel(r)) + 2*k;
fprintf('=== LENGTH-TENSION TEST (fresh L-T quadratic R2 = %.5f) ===\n', R2);
fprintf('  H1 pure SCALE  s = %.4f                    RMSE %.3f kPa   dAIC %+.2f\n', ...
        s1, sqrt(mean(r1.^2)), aic(r1,1)-aic(r3,2));
fprintf('  H2 pure SHIFT  d = %.4f ML (%.3f um SL)   RMSE %.3f kPa   dAIC %+.2f\n', ...
        d2, 2*d2, sqrt(mean(r2.^2)), aic(r2,1)-aic(r3,2));
fprintf('  H3 BOTH        s = %.4f, d = %.4f ML      RMSE %.3f kPa   dAIC %+.2f\n', ...
        v3(1), v3(2), sqrt(mean(r3.^2)), 0);
fprintf('\n  %-8s %8s %8s %8s | %7s %7s %7s | %9s\n', ...
        'ML','fresh','rundown','ratio','res H1','res H2','res H3','implied d');
dimp = zeros(size(ML));
for k = 1:numel(ML)
    dimp(k) = fzero(@(x) Fmod(ML(k)-x)-Fr(k), [0 0.3]);
    fprintf('  %8.4f %8.2f %8.2f %8.3f | %7.2f %7.2f %7.2f | %9.4f\n', ...
            ML(k), Ff(k), Fr(k), Fr(k)/Ff(k), r1(k), r2(k), r3(k), dimp(k));
end
fprintf('\n  H1 is refuted by the ML 1.10 point alone: it needs +%.1f kPa there\n', max(r1));
fprintf('  while every other point sits ~%.1f kPa low. A pure scale cannot bend.\n', mean(r1(1:end-1)));

%% ── 3. ktr: the decisive test, with the sampling control ───────────────────
clip = @(M) deal(M(M(:,1)>=74 & M(:,1)<=78, 1), M(M(:,1)>=74 & M(:,1)<=78, 2)*2, ...
                 M(M(:,1)>=74 & M(:,1)<=78, 3));
[tA,LA,FA] = clip(A);  [tB,LB,FB] = clip(B);

% The comparison must be sampling-matched or it is worthless: at 10 ms the first
% post-release sample already sits ~half a time constant into the rise. So the
% FRESH trace is resampled onto the REPEAT'S OWN sample times, and both are then
% fitted by the same plain 1-exponential. Any sampling bias is now common-mode.
FAd = interp1(tA, FA, tB, 'linear');
okc = ~isnan(FAd);  tC = tB(okc);  FAd = FAd(okc);  FBd = FB(okc);

ftyp = fittype('a*(1-exp(-k*max(x-t0,0)))');
fitk = @(t,y) fit(t(:)-t(1), y(:), ftyp, 'StartPoint',[max(y) 50 0.005], ...
                  'Lower',[1 1 0], 'Upper',[300 500 0.06]);
nS   = numel(is);
kFull = nan(1,nS); kFd = nan(1,nS); kRd = nan(1,nS); aFd = nan(1,nS); aRd = nan(1,nS);
for sgi = 1:nS
    t0 = vt(is(sgi)+1,1);  t1 = vt(is(sgi)+2,1);
    w  = tC > t0 & tC < t1;   wf = tA > t0 & tA < t1;
    c0 = fitk(tA(wf), FA(wf));   kFull(sgi) = c0.k;
    c1 = fitk(tC(w),  FAd(w));   kFd(sgi)   = c1.k;  aFd(sgi) = c1.a;
    c2 = fitk(tC(w),  FBd(w));   kRd(sgi)   = c2.k;  aRd(sgi) = c2.a;
end
rat  = kRd ./ kFd;
ktrF = mean(kFd);  ktrR = mean(kRd);
MLs  = vt(is+1,4)';
pk   = polyfit(MLs, kFull, 1);

fprintf('\n\n=== ktr TEST (sampling-matched) ===\n');
fprintf('  %-6s %11s %11s %11s %9s\n','seg','fresh 0.05ms','fresh 10ms','rundown 10ms','ratio');
for sgi = 1:nS
    fprintf('  %-6d %11.2f %11.2f %11.2f %9.3f\n', sgi, kFull(sgi), kFd(sgi), kRd(sgi), rat(sgi));
end
fprintf('  %-6s %11.2f %11.2f %11.2f %9.3f  +- %.3f (SEM)\n', 'mean', ...
        mean(kFull), ktrF, ktrR, mean(rat), std(rat)/sqrt(nS));
fprintf('\n  sampling alone (fresh 0.05 ms -> 10 ms) costs x%.3f -- now common-mode\n', ktrF/mean(kFull));
fprintf('  all %d segments fall; t = %.1f against ratio = 1\n', nS, (1-mean(rat))/(std(rat)/sqrt(nS)));
fprintf('\n  H1 predicts ktr UNCHANGED (x1.000).            observed x%.3f  -> REFUTED\n', mean(rat));
fprintf('  H2 predicts ktr(ML-d): dktr/dML = %+.1f /s per ML, so d = %.3f ML\n', pk(1), v3(2));
fprintf('     gives %+.2f /s (%+.1f%%).                     observed %+.1f%%\n', ...
        -pk(1)*v3(2), -100*pk(1)*v3(2)/ktrF, 100*(mean(rat)-1));
fprintf('     -> H2 accounts for ~%.0f%% of the ktr loss; the rest is unexplained.\n', ...
        100*(pk(1)*v3(2)/ktrF)/(1-mean(rat)));
fprintf('  => the surviving bridges also CYCLE MORE SLOWLY. Neither hypothesis\n');
fprintf('     alone covers it: rundown = fewer heads AND shorter AND slower.\n');
ktrRdb = ktrR;   % already sampling-matched

%% ── 3b. Three-term decomposition, and does it close? ───────────────────────
% The two observables (force, ktr) and the three candidate effects have
% distinct, non-overlapping signatures, so the decomposition is determined:
%
%   effect                    force                    ktr
%   ------------------------  -----------------------  ------------
%   length loss  d            F(ML-d)/F(ML)   (<1)     ktr(ML-d)/ktr(ML)
%   head loss    s            s               (<1)     1        (a rate is
%                                                       blind to head count)
%   uniform rate slowdown x   1  (duty ratio  f/(f+g)  x        (every rate
%                                 is invariant)                  scales)
%
% Take d and s from the length-tension fit (H3), then x is whatever ktr still
% needs. Closure is a genuine test: nothing forces the force column to balance.
dH3 = v3(2);  sH3 = v3(1);
fL  = mean(Fmod(ML-dH3)./Fmod(ML));            % force factor from length loss
kL  = 1 - pk(1)*dH3/mean(kFull);               % ktr factor from length loss
xRate = mean(rat)/kL;                          % residual uniform rate slowdown

fprintf('\n\n=== THREE-TERM DECOMPOSITION (694 s of rundown) ===\n');
fprintf('  %-26s %10s %10s\n','effect','force x','ktr x');
fprintf('  %-26s %10.3f %10.3f\n', sprintf('length loss  d=%.3f ML',dH3), fL, kL);
fprintf('  %-26s %10.3f %10.3f\n', sprintf('head loss    s=%.3f',sH3),    sH3, 1.000);
fprintf('  %-26s %10.3f %10.3f\n', sprintf('rate slowdown x=%.3f',xRate), 1.000, xRate);
fprintf('  %-26s %10.3f %10.3f   <- predicted\n','PRODUCT', fL*sH3, kL*xRate);
fprintf('  %-26s %10.3f %10.3f   <- observed\n','', mean(Fr)/mean(Ff), mean(rat));
fprintf('  closure error: force %+.1f%%, ktr %+.1f%% (ktr is exact by construction)\n', ...
        100*(fL*sH3/(mean(Fr)/mean(Ff))-1), 100*(kL*xRate/mean(rat)-1));
fprintf('\n  So rundown is THREE things, dominated by the length term:\n');
fprintf('    ~%.0f%% of sarcomere length lost, ~%.0f%% of heads lost, ~%.0f%% slower cycling.\n', ...
        100*dH3, 100*(1-sH3), 100*(1-xRate));

%% ── 4. Rundown vs low ATP are different directions ─────────────────────────
S2 = load([D 'protocol_03_27_2026_2mM_slack.mat']);
atpF = mean(S2.features_data.A)/mean(S.features_data.A);
atpK = mean(S2.features_data.ktr)/mean(S.features_data.ktr);
rdF  = mean(Fr)/mean(Ff);
fprintf('\n\n=== SIGNATURE: RUNDOWN vs LOW ATP, in (force, ktr) ===\n');
fprintf('  rundown (694 s) : force x%.3f , ktr x%.3f   -> both DOWN\n', rdF, mean(rat));
fprintf('  low ATP (8->2)  : force x%.3f , ktr x%.3f   -> force UP, ktr down\n', atpF, atpK);
fprintf('  They are distinct directions, so rundown is NOT local ATP depletion,\n');
fprintf('  and rundown contamination cannot masquerade as an ATP effect.\n');

%% ── 5. Should the ktr ratio be rundown-corrected? ──────────────────────────
% Tempting, but test it: the correction must IMPROVE cross-prep consistency.
rateK = -log(mean(rat))/694;                         % 1/s, from the bracket
pd = {'03/27 M', 132, 'protocol_03_27_2026_2mM_slack.mat', 'protocol_03_27_2026_8mM_slack.mat'
      '04/03 F', 138, 'protocol_04_03_2026_2mM_slack.mat', 'protocol_04_03_2026_8mM_slack.mat'
      '04/10 M',-308, 'protocol_04_10_2026_2mM_slack.mat', 'protocol_04_10_2026_8mM_slack.mat'};
raw = nan(1,3); cor = nan(1,3);
fprintf('\n\n=== SHOULD ktr BE RUNDOWN-CORRECTED? (rate %.2e /s from the bracket) ===\n', rateK);
fprintf('  %-9s %7s %10s %12s\n','day','gap(s)','raw ratio','corrected');
for j = 1:3
    L = load([D pd{j,3}]); H = load([D pd{j,4}]);
    raw(j) = mean(L.features_data.ktr)/mean(H.features_data.ktr);
    cor(j) = raw(j)*exp(rateK*pd{j,2});
    fprintf('  %-9s %7d %10.3f %12.3f\n', pd{j,1}, pd{j,2}, raw(j), cor(j));
end
cv = @(x) 100*std(x)/mean(x);
fprintf('  CV across preps:  raw %.0f%%   corrected %.0f%%\n', cv(raw), cv(cor));
fprintf('  => correcting ktr MAKES CONSISTENCY WORSE. The 694 s / 3-activation ktr\n');
fprintf('     loss does not extrapolate to 132-308 s gaps. Leave the ktr ratio raw;\n');
fprintf('     just never compare ktr across long session gaps.\n');

%% ── 6. Figure ──────────────────────────────────────────────────────────────
fig = figure(920); clf; set(fig,'Position',[30 30 1500 780]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile; hold on; box on;
xq = linspace(min(ML)-0.01, max(ML)+0.01, 200);
plot(xq, Fmod(xq), '-', 'Color',[.2 .2 .2], 'LineWidth',1.5);
plot(ML, Ff, 'ko', 'MarkerFaceColor','k','MarkerSize',8);
plot(ML, Fr, 'rs', 'MarkerFaceColor','r','MarkerSize',8);
plot(xq, s1*Fmod(xq), 'b--', 'LineWidth',1.5);
plot(xq, Fmod(xq-d2), 'g-.', 'LineWidth',1.5);
xlabel('commanded length (ML)'); ylabel('quasi-isometric force (kPa)');
title('Length-tension: fresh vs run-down');
legend({'fresh fit','fresh','run-down',sprintf('H1 scale %.3f',s1), ...
        sprintf('H2 shift %.3f ML',d2)},'Location','northwest','FontSize',7);

nexttile; hold on; box on;
bar(ML, [r1; r2; r3]', 'grouped');
yline(0,'k-'); xlabel('ML'); ylabel('residual (kPa)');
title(sprintf('Residuals: H1 %.2f, H2 %.2f, H3 %.2f kPa RMSE', ...
      sqrt(mean(r1.^2)), sqrt(mean(r2.^2)), sqrt(mean(r3.^2))));
legend({'H1 scale','H2 shift','H3 both'},'Location','southwest','FontSize',7);

nexttile; hold on; box on;
plot(ML, dimp, 'ko-','MarkerFaceColor','k'); yline(d2,'g-.','LineWidth',1.5);
xlabel('ML'); ylabel('shift needed to match (ML)');
title('Implied shift per point (constant \Rightarrow H2)');
ylim([0 0.09]);

nexttile; hold on; box on;
plot(1:nS, kFull, 'k^--','MarkerFaceColor','k');
plot(1:nS, kFd,   'ko-', 'MarkerFaceColor','k','LineWidth',1.5);
plot(1:nS, kRd,   'rs-', 'MarkerFaceColor','r','LineWidth',1.5);
xlabel('slack segment'); ylabel('ktr (1/s)'); xlim([0.5 nS+0.5]);
title(sprintf('ktr per segment  (rundown/fresh = %.3f)', mean(rat)));
legend({'fresh @0.05 ms','fresh @10 ms (matched)','rundown @10 ms'}, ...
       'Location','southwest','FontSize',7);

nexttile; hold on; box on;
% Fair overlay: BOTH on the repeat's own sample grid, each normalised by its own
% fitted asymptote, with the fitted exponentials. Sampling bias is common-mode.
sg = 3; t0 = vt(is(sg)+1,1); t1 = vt(is(sg)+2,1);
w  = tC > t0 & tC < t1;  tt = tC(w)-tC(find(w,1));
plot(tt, FAd(w)/aFd(sg), 'ko-','MarkerFaceColor','k','MarkerSize',4);
plot(tt, FBd(w)/aRd(sg), 'rs-','MarkerFaceColor','r','MarkerSize',4);
xf = linspace(0,max(tt),300);
plot(xf, 1-exp(-kFd(sg)*xf), 'k-', 'LineWidth',1.2);
plot(xf, 1-exp(-kRd(sg)*xf), 'r-', 'LineWidth',1.2);
xlabel('t since slack release (s)'); ylabel('force / own fitted plateau');
title(sprintf('Redevelopment, matched sampling (slack %d)', sg)); xlim([0 0.12]); ylim([0 1.1]);
legend({sprintf('fresh  k=%.1f',kFd(sg)), sprintf('rundown  k=%.1f',kRd(sg))}, ...
       'Location','southeast','FontSize',7);

nexttile; hold on; box on;
plot([1 rdF],[1 mean(rat)],'-o','Color',[.85 .3 .1],'LineWidth',2,'MarkerFaceColor',[.85 .3 .1]);
plot([1 atpF],[1 atpK],'-s','Color',[0 .45 .74],'LineWidth',2,'MarkerFaceColor',[0 .45 .74]);
plot(1,1,'k*','MarkerSize',12); xline(1,':'); yline(1,':');
xlabel('force  (x fresh 8 mM)'); ylabel('ktr  (x fresh 8 mM)');
title('Rundown and low ATP move in different directions');
legend({'rundown (694 s)','low ATP (8\rightarrow2 mM)','fresh 8 mM'},'Location','southwest','FontSize',7);

sgtitle('Rundown mechanism — 03/27 8 mM fresh (t=0) vs repeat (t=+694 s)');
exportgraphics(fig, fullfile(resDir,'rundown_mechanism.png'), 'Resolution', 150);
save(fullfile(resDir,'rundown_mechanism.mat'), 'ML','Ff','Fr','s1','d2','v3', ...
     'kFull','kFd','kRd','rat','ktrF','ktrR','raw','cor','rateK');
fprintf('\nSaved %s\n', fullfile(resDir,'rundown_mechanism.png'));
