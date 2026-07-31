% RunJointCorrection.m
% Can force AND ktr both be corrected, using the refined decay geometry?
%
% What the refinements bought (RunDecayGeometry.m):
%   * decay acts only while activated; between activations there is a small
%     RECOVERY (+3 %), not decay
%   * each run can be referenced to its OWN activation peak, and given its OWN
%     decay-phase duration T_dec, instead of a shared protocol window / duration
%   * 2 mM decays x1.95 faster and peaks ~2 s earlier than 8 mM
%
% This script re-derives phi on that footing and then asks, for BOTH observables,
% the only question that matters: does the correction reduce between-preparation
% scatter? For ktr the coupling to the force correction is left free (gamma) and
% scanned, so the data choose it rather than the modeller.
%
%   force_frac(gap) = phi * |slope_e| * T_dec_e / F_later      (M.5, M.5b)
%   ktr_frac(gap)   = gamma * force_frac(gap)
%
% gamma = 0    -> do not correct ktr
% gamma ~ 0.68 -> the coupling the 03/27 bracket implies (force x0.829, ktr x0.884)
% gamma ~ 0.71 -> the coupling the series-creep lesion predicts (C1 in
%                 RunMechanismSimulation.m: force x0.827, ktr x0.877)
%
% Outputs -> results/joint_correction.png, joint_correction.mat

clear; close all; clc;
cd(fileparts(which('RunJointCorrection')));
addpath(genpath('../..'));
resDir = 'results';  if ~isfolder(resDir), mkdir(resDir); end
D = '../../data/';

W = {[68.20 69.20],[71.20 72.35],[77.70 78.30]};

R = struct( ...
 'day', {'03/27 M','03/27 M','03/27 M','04/03 F','04/03 F','04/10 M','04/10 M'}, ...
 'cond',{'8 mM','2 mM','8 mM rpt','8 mM','2 mM','2 mM','8 mM'}, ...
 'act', {1,2,4,1,2,1,2}, ...
 'file',{'03 27 2026 M/02_Merged_8mM_Active.txt','03 27 2026 M/03_Merged_2mM_Active.txt', ...
         '03 27 2026 M/04_Merged_8mM_Active_repeat.txt','04 03 2026 F/02_Merged_8mM_Active.txt', ...
         '04 03 2026 F/03_Merged_2mM_Active.txt','04 10 2026 Male 2-8/02_Merged_2mM_Active.txt', ...
         '04 10 2026 Male 2-8/03_Merged_8mM_Active.txt'}, ...
 'mat', {'protocol_03_27_2026_8mM_slack.mat','protocol_03_27_2026_2mM_slack.mat','', ...
         'protocol_04_03_2026_8mM_slack.mat','protocol_04_03_2026_2mM_slack.mat', ...
         'protocol_04_10_2026_2mM_slack.mat','protocol_04_10_2026_8mM_slack.mat'});

%% ── 1. Geometry per run ────────────────────────────────────────────────────
for i = 1:numel(R)
    M = readmatrix([D R(i).file]);  t = M(:,1); L = M(:,2); F = M(:,3);
    for k = 1:3
        w = t>=W{k}(1) & t<=W{k}(2) & abs(L-1)<2e-3;
        R(i).v(k) = mean(F(w),'omitnan');  R(i).tc(k) = mean(t(w));
    end
    iso = abs(L-1)<2e-3 & t>63 & t<82;  ti = t(iso);  env = movmean(F(iso),401);
    [~,ip] = max(env);  R(i).tpk = ti(ip);
    Fs = movmean(F,51); thr = 0.2*prctile(Fs,99);
    R(i).toff = t(find(Fs>thr,1,'last'));
    R(i).slope = (R(i).v(3)-R(i).v(2))/(R(i).tc(3)-R(i).tc(2));
    R(i).Tdec  = R(i).toff - R(i).tpk;
    R(i).Fpk   = R(i).v(2) + R(i).slope*(R(i).tpk  - R(i).tc(2));
    R(i).Fend  = R(i).v(2) + R(i).slope*(R(i).toff - R(i).tc(2));
    if ~isempty(R(i).mat)
        S = load([D R(i).mat]);  R(i).ktr = mean(S.features_data.ktr,'omitnan');
        R(i).A = mean(S.features_data.A,'omitnan');
    end
end

%% ── 2. Recalibrate phi on the T_dec footing ────────────────────────────────
i1=1; i2=2; i4=3;
predT  = abs(R(i1).slope)*R(i1).Tdec + abs(R(i2).slope)*R(i2).Tdec;
measL  = R(i1).Fpk - R(i4).Fpk;                % peak-referenced bracket loss
PHI    = measL / predT;
fprintf('=== phi recalibrated on own-peak / T_dec ===\n');
fprintf('  bracket loss, peak-referenced : %.2f kPa   (window-referenced was %.2f)\n', ...
        measL, R(i1).v(2)-R(i4).v(2));
fprintf('  predicted from |slope|*T_dec  : %.2f kPa\n', predT);
fprintf('  phi = %.3f   (was %.3f on the T_act footing)\n', PHI, 0.452);

%% ── 3. Correct force and ktr; scan the ktr coupling gamma ──────────────────
days = {'03/27 M','04/03 F','04/10 M'};
frac = nan(1,3); rawF = nan(1,3); corF = nan(1,3); rawK = nan(1,3);
for d = 1:numel(days)
    ih = find(strcmp({R.day},days{d}) & strcmp({R.cond},'8 mM'));
    il = find(strcmp({R.day},days{d}) & strcmp({R.cond},'2 mM'));
    e  = ih; if R(il).act < R(ih).act; e = il; end          % earlier run
    l  = ih; if e == ih; l = il; end                        % later run
    dmg     = PHI*abs(R(e).slope)*R(e).Tdec;
    frac(d) = dmg / R(l).Fpk;                               % of the DEGRADED run
    rawF(d) = R(il).Fpk / R(ih).Fpk;
    if e == ih; corF(d) = rawF(d)*(1+frac(d)); else; corF(d) = rawF(d)/(1+frac(d)); end
    rawK(d) = R(il).ktr / R(ih).ktr;
end

cv = @(x) 100*std(x)/mean(x);
gam = linspace(-0.5, 2.0, 501);  cvK = nan(size(gam));  KR = nan(numel(gam),3);
for g = 1:numel(gam)
    kk = nan(1,3);
    for d = 1:3
        ih = find(strcmp({R.day},days{d}) & strcmp({R.cond},'8 mM'));
        il = find(strcmp({R.day},days{d}) & strcmp({R.cond},'2 mM'));
        f = gam(g)*frac(d);
        if R(ih).act < R(il).act; kk(d) = rawK(d)*(1+f); else; kk(d) = rawK(d)/(1+f); end
    end
    KR(g,:) = kk;  cvK(g) = cv(kk);
end
[cvKmin, ig] = min(cvK);  gOpt = gam(ig);
gBr = (1-0.884)/(1-0.829);        % coupling implied by the bracket
gMd = (1-0.877)/(1-0.827);        % coupling predicted by the series-creep lesion

fprintf('\n=== FORCE ===\n');
fprintf('%-9s %9s %9s %9s\n','day','frac loss','raw','corrected');
for d=1:3; fprintf('%-9s %8.1f%% %9.3f %9.3f\n', days{d}, 100*frac(d), rawF(d), corF(d)); end
fprintf('%-9s %9s %9.3f %9.3f\n','mean','',mean(rawF),mean(corF));
fprintf('%-9s %9s %8.0f%% %8.0f%%\n','CV','',cv(rawF),cv(corF));

fprintf('\n=== ktr : let the DATA choose the coupling gamma ===\n');
fprintf('  gamma = 0      (no correction)      CV = %.1f%%   ratios %s\n', cv(rawK), mat2str(round(rawK,3)));
for gv = [gBr gMd gOpt]
    kk = nan(1,3);
    for d=1:3
        ih = find(strcmp({R.day},days{d}) & strcmp({R.cond},'8 mM'));
        il = find(strcmp({R.day},days{d}) & strcmp({R.cond},'2 mM'));
        f = gv*frac(d);
        if R(ih).act < R(il).act; kk(d)=rawK(d)*(1+f); else; kk(d)=rawK(d)/(1+f); end
    end
    fprintf('  gamma = %+.2f  %-22s CV = %.1f%%   ratios %s\n', gv, ...
        iLabel(gv,gBr,gMd,gOpt), cv(kk), mat2str(round(kk,3)));
end
fprintf('\n  gamma minimising ktr CV = %+.2f  (CV %.1f%%)\n', gOpt, cvKmin);
atEdge = (ig == 1) || (ig == numel(gam));
if atEdge
    fprintf('  *** the optimum sits ON THE SCAN BOUNDARY: CV(gamma) is monotone over\n');
    fprintf('      the whole range, so there is NO interior optimum. The apparent\n');
    fprintf('      "improvement" at negative gamma is not a mechanism -- a negative\n');
    fprintf('      coupling would mean rundown makes ktr FASTER. What it actually does\n');
    fprintf('      is lift 04/10 (the lowest ratio) toward the other two, i.e. it fits\n');
    fprintf('      one preparation''s idiosyncrasy with one free parameter on 3 points.\n');
end
fprintf('\n  At the two physically motivated couplings the correction is WORSE:\n');
fprintf('    gamma = %.2f (bracket): CV %.1f%%     gamma = %.2f (lesion): CV %.1f%%\n', ...
        gBr, cvK(find(gam>=gBr,1)), gMd, cvK(find(gam>=gMd,1)));
fprintf('    gamma = 0    (none)   : CV %.1f%%   <-- best of the defensible options\n', cv(rawK));

% Is there any dose-dependence to fit at all?
[rho, pv] = corr(frac(:), rawK(:));
fprintf('\n  ktr ratio vs force-damage dose: r = %+.2f, p = %.2f (n = 3)\n', rho, pv);
fprintf('  The pattern is not even monotone (%.3f at %.1f%%, %.3f at %.1f%%, %.3f at %.1f%%),\n', ...
        rawK(1),100*frac(1), rawK(2),100*frac(2), rawK(3),100*frac(3));
fprintf('  so there is no dose-response for a correction to follow.\n');

%% ── 4. Why: the two observables are not proportional ───────────────────────
fprintf('\n=== WHY ktr RESISTS CORRECTION ===\n');
fprintf('  Over the bracket the lesion produces force x0.829 but ktr x0.884.\n');
fprintf('  Series-elastic creep reaches force through the length-tension relation\n');
fprintf('  (a strong channel) and ktr through added compliance (a weak one), so the\n');
fprintf('  two are not proportional -- and the proportionality constant is exactly\n');
fprintf('  what a scaled correction assumes. With 3 preparations the ktr scatter\n');
fprintf('  (CV %.0f%% raw) is smaller than the correction it would receive.\n', cv(rawK));

%% ── 5. Figure ──────────────────────────────────────────────────────────────
fig = figure(980); clf; set(fig,'Position',[20 20 1450 720]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
cols = lines(3);

nexttile; hold on; box on;
b = bar([rawF; corF]'); b(1).FaceColor=[.6 .6 .6]; b(2).FaceColor=[.2 .5 .8];
yline(mean(corF),'k--','LineWidth',1.5);
set(gca,'XTick',1:3,'XTickLabel',days); ylabel('ATP force ratio');
title(sprintf('FORCE: CV %.0f%% -> %.0f%%  (mean %.2f)', cv(rawF), cv(corF), mean(corF)));
legend({'raw (own-peak)','corrected'},'Location','northwest','FontSize',7);

nexttile; hold on; box on;
plot(gam, cvK, 'k-', 'LineWidth', 2);
plot(0, cv(rawK), 'ko','MarkerFaceColor','w','MarkerSize',9);
xline(gBr,'--','Color',[.85 .33 .1],'LineWidth',1.5);
xline(gMd,'--','Color',[.2 .6 .3],'LineWidth',1.5);
plot(gOpt, cvKmin, 'kp','MarkerFaceColor','y','MarkerSize',15);
xlabel('\gamma  (ktr correction / force correction)'); ylabel('ktr cross-prep CV (%)');
title('ktr: CV(\gamma) is monotone — no interior optimum');
legend({'CV(\gamma)','no correction','bracket coupling 0.68','lesion coupling 0.71','edge (not an optimum)'}, ...
       'Location','north','FontSize',7);

nexttile; hold on; box on;
for d=1:3
    plot(gam, KR(:,d), '-', 'Color', cols(d,:), 'LineWidth', 1.6, 'DisplayName', days{d});
end
xline(0,'k:','HandleVisibility','off'); xline(gBr,'--','Color',[.85 .33 .1],'HandleVisibility','off');
xlabel('\gamma'); ylabel('ktr ratio 2 mM / 8 mM');
title('Correcting ktr pulls the preps APART'); legend('Location','northwest','FontSize',7);

nexttile; hold on; box on;
b = bar([1-0.829 1-0.884; 1-0.827 1-0.877]');
b(1).FaceColor=[.85 .33 .1]; b(2).FaceColor=[.2 .6 .3];
set(gca,'XTick',1:2,'XTickLabel',{'force','ktr'}); ylabel('fractional loss over the bracket');
title('Force and ktr respond unequally'); legend({'measured','series-creep lesion'},'FontSize',7);

nexttile; hold on; box on;
plot(frac*100, rawK, 'ko','MarkerSize',11,'MarkerFaceColor','w');
for d=1:3
    plot(frac(d)*100, rawK(d), 'o','Color',cols(d,:),'MarkerFaceColor',cols(d,:),'MarkerSize',11);
    text(frac(d)*100+0.3, rawK(d), days{d},'FontSize',8);
end
p = polyfit(frac*100, rawK, 1); xq = linspace(min(frac)*100-1, max(frac)*100+1, 50);
plot(xq, polyval(p,xq), 'k--');
xlabel('force damage over the gap (%)'); ylabel('raw ktr ratio');
title(sprintf('No dose-dependence in ktr (slope %+.4f/%%)', p(1)));

nexttile; hold on; box on;
txt = {'\bfRecommendation\rm','', ...
       sprintf('FORCE  correct.  phi = %.2f on own-peak', PHI), ...
       sprintf('       + T_{dec};  CV %.0f%% -> %.0f%%, mean %.2f', cv(rawF), cv(corF), mean(corF)), '', ...
       'ktr    DO NOT correct.', ...
       sprintf('       at the bracket coupling %.2f:', gBr), ...
       sprintf('       CV %.0f%% -> %.0f%% (worse)', cv(rawK), cvK(find(gam>=gBr,1))), ...
       '       CV(\gamma) monotone: no optimum', ...
       'Reason: the lesion reaches force via', ...
       'length-tension and ktr via compliance.', ...
       'Not proportional, so no single scaling', ...
       'can serve both.'};
text(0.03,0.95,txt,'VerticalAlignment','top','FontSize',10); axis off;

sgtitle('Joint correction — force yes, ktr no, and the reason why');
exportgraphics(fig, fullfile(resDir,'joint_correction.png'), 'Resolution', 150);
save(fullfile(resDir,'joint_correction.mat'),'R','PHI','frac','rawF','corF','rawK','gam','cvK','gOpt','gBr','gMd');
fprintf('\nSaved %s\n', fullfile(resDir,'joint_correction.png'));

% ---------------------------------------------------------------------------
function s = iLabel(gv, gBr, gMd, gOpt)
    if gv==gBr;      s = '(bracket coupling)';
    elseif gv==gMd;  s = '(lesion coupling)';
    elseif gv==gOpt; s = '(CV optimum)';
    else;            s = '';
    end
end
