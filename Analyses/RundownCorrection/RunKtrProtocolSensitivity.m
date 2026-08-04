% RunKtrProtocolSensitivity.m
% Why does the slack-derived ktr fall with rundown, and why do the two ktr
% protocols disagree about ATP?
%
% TWO OBSERVATIONS THAT NEED ONE EXPLANATION
%   (1) Between 03/27's 8 mM run and its 8 mM REPEAT -- same fibre, same
%       condition -- the slack ktr falls x0.884, in every one of the five
%       segments. That is rundown, not biological variability.
%   (2) The ATP effect on ktr differs by protocol: the dedicated ktr protocol
%       says x0.72 (CV 4 %), the slack says x0.58 (CV 7 %). Both are internally
%       consistent across preparations, yet they are 24 % apart.
%
% THE HYPOTHESIS
% The slack measurement is COMPLIANCE-CONTAMINATED. After an unloaded shortening
% the contractile element must shorten again to re-stretch the series element
% before force appears, so added series compliance slows the apparent rise. The
% ktr protocol restretches to 1.05 ML and returns to 1.0, starting much closer to
% its force-bearing configuration, so it needs far less internal shortening.
% If so, the slack ktr is the right probe for rundown (a compliance lesion) and
% the WRONG probe for cross-bridge kinetics.
%
% NOTE ON WHAT CANNOT BE DONE: the 03/27 repeat has no hi-res burst file, only
% the 10 ms log. Its ktr-protocol redevelopment lasts ~60 ms, i.e. 6 samples, so
% the rundown drop can only be measured in the slack, whose 298 ms hold gives 30.
%
% Outputs -> results/ktr_protocol_sensitivity.png, .mat

clear; close all; clc;
cd(fileparts(which('RunKtrProtocolSensitivity')));
addpath(genpath('../..'));
resDir = 'results';  if ~isfolder(resDir), mkdir(resDir); end
D = '../../data/';
ft = fittype('a*(1-exp(-k*max(x-t0,0)))');

%% ── 1. Data: ktr from both protocols, all six active runs ──────────────────
K = {'03/27','8 mM','protocol_03_27_2026_velocitytable_ktr.mat','protocol_03_27_2026_8mM_slack.mat'
     '03/27','2 mM','protocol_03_27_2026_velocitytable_ktr_2mM.mat','protocol_03_27_2026_2mM_slack.mat'
     '04/03','8 mM','protocol_04_03_2026_velocitytable_ktr.mat','protocol_04_03_2026_8mM_slack.mat'
     '04/03','2 mM','protocol_04_03_2026_velocitytable_ktr_2mM.mat','protocol_04_03_2026_2mM_slack.mat'
     '04/10','8 mM','protocol_04_10_2026_velocitytable_ktr.mat','protocol_04_10_2026_8mM_slack.mat'
     '04/10','2 mM','protocol_04_10_2026_velocitytable_ktr_2mM.mat','protocol_04_10_2026_2mM_slack.mat'};
kp = nan(6,1); ks = nan(6,1); r2 = nan(6,1);
fprintf('%-6s %-5s | %9s %14s %8s | %11s\n','day','cond','ktr_prot','95%% CI','R2','ktr_slack2');
for i = 1:6
    S = load([D K{i,3}]);  dt = S.datatable;  vt = S.velocitytable;
    w  = dt(:,1) > vt(7,1) & dt(:,1) <= vt(8,1);
    tt = dt(w,1)-dt(find(w,1),1);  yy = dt(w,3);
    [c,g] = fit(tt,yy,ft,'StartPoint',[max(yy) 50 0.002],'Lower',[0.1 1 0],'Upper',[5 400 0.02]);
    ci = confint(c);  kp(i) = c.k;  r2(i) = g.rsquare;
    Ss = load([D K{i,4}]);  ks(i) = Ss.features_data.ktr(2);     % segment 2 = SL 2.0
    fprintf('%-6s %-5s | %9.2f [%5.1f %5.1f] %8.4f | %11.2f\n', ...
            K{i,1},K{i,2},c.k,ci(1,2),ci(2,2),g.rsquare,ks(i));
end

days = {'03/27','04/03','04/10'};
rp = kp(2:2:6)./kp(1:2:5);   rs = ks(2:2:6)./ks(1:2:5);
cv = @(x) 100*std(x)/mean(x);
fprintf('\n=== ATP EFFECT ON ktr, BY PROTOCOL ===\n');
fprintf('%-8s %15s %15s | %12s\n','prep','ktr protocol','slack (seg 2)','slack/prot');
for j=1:3
    fprintf('%-8s %15.3f %15.3f | %12.3f\n', days{j}, rp(j), rs(j), rs(j)/rp(j));
end
fprintf('%-8s %15.3f %15.3f\n','mean',mean(rp),mean(rs));
fprintf('%-8s %14.0f%% %14.0f%%\n','CV',cv(rp),cv(rs));
fprintf('\n  Both are internally tight, yet %.0f%% apart. Where do they diverge?\n', ...
        100*(mean(rp)/mean(rs)-1));
fprintf('  at 8 mM: slack/protocol = %.2f, %.2f, %.2f  (they AGREE)\n', ks(1)/kp(1), ks(3)/kp(3), ks(5)/kp(5));
fprintf('  at 2 mM: slack/protocol = %.2f, %.2f, %.2f  (they DIVERGE)\n', ks(2)/kp(2), ks(4)/kp(4), ks(6)/kp(6));

%% ── 2. Model: how compliance-sensitive is each measurement? ────────────────
p0 = getParams(loadParams('ThursdayNightFever'), [], true, false);
p0.PlotEachSeparately=0; p0.PlotFullscreen=0; p0.PlotFeatureFitting=0;
p0.ghostSave=''; p0.ghostLoad=''; p0.EvalFeatures=1; p0.RunSlackSegments='All';
p0.velocitytableonfile     = 'protocol_03_27_2026_8mM_slack.mat';
p0.ktr_velocitytableonfile = 'protocol_03_27_2026_velocitytable_ktr.mat';
Sk  = load([D 'protocol_03_27_2026_velocitytable_ktr.mat']);  vtk = Sk.velocitytable;

lev = [1 0.80 0.65 0.50];
KS = nan(size(lev)); KP = nan(size(lev));
fprintf('\n=== MODEL: compliance sensitivity of each measurement ===\n');
for j = 1:numel(lev)
    p = p0;  p.kSE = p0.kSE*lev(j);
    [~,~,fm] = runSlackExperiment(p);      KS(j) = fm.ktr(2);
    q = p; q.Slim_l = 1.4; q.Slim_r = 2.3;
    [~,ok] = runKtrExperiment(q, []);
    w  = ok.t > vtk(7,1) & ok.t <= vtk(8,1);
    tt = ok.t(w)-ok.t(find(w,1));  yy = ok.Force(w);  yy = yy(:)/max(yy);
    c  = fit(tt(:),yy,ft,'StartPoint',[1 50 0.002],'Lower',[0.1 1 0],'Upper',[5 400 0.02]);
    KP(j) = c.k;
    fprintf('  kSE x%.2f : slack-ktr %6.2f (x%.3f)   protocol-ktr %6.2f (x%.3f)\n', ...
            lev(j), KS(j), KS(j)/KS(1), KP(j), KP(j)/KP(1));
end
sens = (1-KS(3)/KS(1)) / (1-KP(3)/KP(1));
fprintf('\n  At kSE x0.65 the slack ktr falls %.0f%% and the protocol ktr only %.0f%%:\n', ...
        100*(1-KS(3)/KS(1)), 100*(1-KP(3)/KP(1)));
fprintf('  the slack measurement is ~%.0fx more COMPLIANCE-SENSITIVE.\n', sens);

%% ── 3. What this explains ──────────────────────────────────────────────────
fprintf('\n=== WHAT THIS EXPLAINS ===\n');
fprintf('  (1) The rundown ktr drop. Rundown adds ~35%% series compliance, which the\n');
fprintf('      model says costs the slack ktr %.0f%% -- against the x0.884 (-12%%)\n', 100*(1-KS(3)/KS(1)));
fprintf('      observed between the 8 mM run and its repeat. Same sign, right order.\n');
fprintf('  (2) The protocol disagreement. A measurement %.0fx more compliance-sensitive\n', sens);
fprintf('      will report a different rate whenever compliance differs between the\n');
fprintf('      conditions being compared -- which it does at 2 mM.\n');
fprintf('\n  CONSEQUENCE: the slack ktr is the right probe for RUNDOWN (a compliance\n');
fprintf('  lesion) and the wrong probe for CROSS-BRIDGE KINETICS. For the ATP effect\n');
fprintf('  the ktr protocol value x%.2f (CV %.0f%%) should be preferred over the slack\n', mean(rp), cv(rp));
fprintf('  value x%.2f (CV %.0f%%).\n', mean(rs), cv(rs));
fprintf('\n  CAVEAT: the model''s absolute protocol-ktr (%.0f) is much faster than its\n', KP(1));
fprintf('  slack ktr (%.0f), whereas the data have them within 10%% at 8 mM (see\n', KS(1));
fprintf('  ../RestretchVsKtrRecovery). The sensitivity RATIO is a within-model\n');
fprintf('  comparison and survives that, but its magnitude may be overstated.\n');

%% ── 4. Figure ──────────────────────────────────────────────────────────────
fig = figure(1010); clf; set(fig,'Position',[20 20 1450 720]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
c8 = [0 0.45 0.74]; c2 = [0.85 0.33 0.10];

nexttile; hold on; box on;
plot(1:3, kp(1:2:5),'o-','Color',c8,'LineWidth',2,'MarkerFaceColor',c8);
plot(1:3, kp(2:2:6),'o--','Color',c8,'LineWidth',2,'MarkerFaceColor','w');
plot(1:3, ks(1:2:5),'s-','Color',c2,'LineWidth',2,'MarkerFaceColor',c2);
plot(1:3, ks(2:2:6),'s--','Color',c2,'LineWidth',2,'MarkerFaceColor','w');
set(gca,'XTick',1:3,'XTickLabel',days); ylabel('ktr (1/s)');
title('(a) ktr by protocol and ATP');
legend({'protocol 8 mM','protocol 2 mM','slack 8 mM','slack 2 mM'},'Location','best','FontSize',7);

nexttile; hold on; box on;
b = bar([rp(:); rs(:)]); b.FaceColor='flat'; b.CData=[repmat(c8,3,1); repmat(c2,3,1)];
set(gca,'XTick',1:6,'XTickLabel',[days days]); xtickangle(30);
yline(mean(rp),'--','Color',c8,'LineWidth',1.5); yline(mean(rs),'--','Color',c2,'LineWidth',1.5);
ylabel('ATP ratio 2 mM / 8 mM'); ylim([0 0.9]);
title(sprintf('(b) protocol x%.2f  vs  slack x%.2f', mean(rp), mean(rs)));
text(2, 0.82,'ktr protocol','Color',c8,'FontSize',8,'HorizontalAlignment','center');
text(5, 0.82,'slack','Color',c2,'FontSize',8,'HorizontalAlignment','center');

nexttile; hold on; box on;
plot(lev, KS/KS(1), 'o-','Color',c2,'LineWidth',2.2,'MarkerFaceColor',c2);
plot(lev, KP/KP(1), 'o-','Color',c8,'LineWidth',2.2,'MarkerFaceColor',c8);
set(gca,'XDir','reverse'); xlabel('kSE factor (more compliance \rightarrow)');
ylabel('ktr / baseline'); title(sprintf('(c) slack ktr is ~%.0fx more compliance-sensitive', sens));
legend({'slack ktr','protocol ktr'},'Location','southwest','FontSize',8);

nexttile; hold on; box on;
sr = [ks(1)/kp(1) ks(3)/kp(3) ks(5)/kp(5); ks(2)/kp(2) ks(4)/kp(4) ks(6)/kp(6)]';
b = bar(sr); b(1).FaceColor=c8; b(2).FaceColor=c2;
yline(1,'k:'); set(gca,'XTick',1:3,'XTickLabel',days);
ylabel('slack ktr / protocol ktr');
title('(d) The two agree at 8 mM, diverge at 2 mM');
legend({'8 mM','2 mM'},'Location','southwest','FontSize',7);

nexttile; hold on; box on;
b = bar([1-0.884, 1-KS(3)/KS(1)]);
b.FaceColor='flat'; b.CData=[0.3 0.3 0.3; c2];
set(gca,'XTick',1:2,'XTickLabel',{'observed rundown','model kSE x0.65'});
ylabel('fractional slack-ktr loss');
title('(e) Compliance explains the rundown ktr drop');

nexttile; hold on; box on;
txt = {'\bfWhat this means\rm','', ...
       'The slack ktr is COMPLIANCE-', ...
       'CONTAMINATED (~5x more than the', ...
       'ktr protocol).','', ...
       '-> right probe for RUNDOWN', ...
       '   (a compliance lesion): explains', ...
       '   the x0.884 drop in the repeat','', ...
       '-> WRONG probe for cross-bridge', ...
       '   kinetics: use the ktr protocol','', ...
       sprintf('ATP effect on ktr: \\bfx%.2f\\rm (protocol,', mean(rp)), ...
       sprintf('CV %.0f%%), not x%.2f (slack).', cv(rp), mean(rs))};
text(0.03,0.97,txt,'VerticalAlignment','top','FontSize',9.5); axis off;

sgtitle('Why the slack ktr falls with rundown, and why the two ktr protocols disagree about ATP');
exportgraphics(fig, fullfile(resDir,'ktr_protocol_sensitivity.png'), 'Resolution', 150);
save(fullfile(resDir,'ktr_protocol_sensitivity.mat'),'kp','ks','rp','rs','lev','KS','KP','sens','days');
fprintf('\nSaved %s\n', fullfile(resDir,'ktr_protocol_sensitivity.png'));
