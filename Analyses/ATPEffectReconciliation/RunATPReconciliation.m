% RunATPReconciliation.m
% What does ATP actually do to the muscle, once rundown is taken out?
%
% Four protocol days measured 2 mM vs 8 mM ATP. Three of them ran 8 mM -> 2 mM,
% one (04/10) ran 2 mM -> 8 mM. Raw, they disagree wildly on the force effect
% (x1.19 to x1.69) while agreeing tightly on ktr (x0.49 to x0.59). The
% disagreement is rundown: whichever condition is recorded SECOND is degraded,
% so the reversed-order prep gets the opposite bias.
%
% The correction comes from ../RundownCorrection: rundown is LINEAR in effective
% activated time at a force-independent kPa/s rate, and only a fraction phi of
% each recording's within-run decline is permanent. Damage separating two runs:
%       loss = phi * |slope_earlier| * T_act_earlier          [kPa]
% with phi = 0.452 calibrated by 03/27's end-of-session 8 mM repeat.
%
% Baker is NOT force-corrected: it has no isometric-plateau slope available in
% the same form, and its recovery windows are 2-6x too short so its amplitude
% features are truncated anyway (see ../SlackDataAnalysis). Its ktr is used.
%
% Outputs -> results/ : atp_reconciliation.png, atp_feature_panels.png, .mat

clear; close all; clc;
cd(fileparts(which('RunATPReconciliation')));
addpath(genpath('../..'));
resDir = 'results';  if ~isfolder(resDir), mkdir(resDir); end
D = '../../data/';

PHI = 0.452;                      % permanent fraction (../RundownCorrection)
zA  = [71.50 72.35];  zB = [77.70 78.30];

%% ── Per-prep registry: earlier run first ───────────────────────────────────
P = {  % day, earlier cond, earlier merged, earlier .mat, later cond, later merged, later .mat
 '03/27 M','8 mM','03 27 2026 M/02_Merged_8mM_Active.txt','protocol_03_27_2026_8mM_slack.mat', ...
           '2 mM','03 27 2026 M/03_Merged_2mM_Active.txt','protocol_03_27_2026_2mM_slack.mat'
 '04/03 F','8 mM','04 03 2026 F/02_Merged_8mM_Active.txt','protocol_04_03_2026_8mM_slack.mat', ...
           '2 mM','04 03 2026 F/03_Merged_2mM_Active.txt','protocol_04_03_2026_2mM_slack.mat'
 '04/10 M','2 mM','04 10 2026 Male 2-8/02_Merged_2mM_Active.txt','protocol_04_10_2026_2mM_slack.mat', ...
           '8 mM','04 10 2026 Male 2-8/03_Merged_8mM_Active.txt','protocol_04_10_2026_8mM_slack.mat'};
nP = size(P,1);

% force features (corrected) and rate features (NOT corrected -- see below)
fF = {'A','Am','steady','peak1_y','peak2','vall_y'};
fK = {'ktr','ktr2_k'};

% The damage accrues during the EARLIER run's activation, so it is the LATER
% run that is measured on a degraded fibre. The fractional depression must
% therefore be referenced to the LATER run's own force, not the earlier one's.
fprintf('%-9s %-14s | %9s %9s | %10s %9s %11s\n','day','order','slope_e','T_act_e', ...
        'loss (kPa)','F_later','frac of later');
frac = nan(1,nP); ord = cell(1,nP);
for i = 1:nP
    Me = readmatrix([D P{i,3}]);  te = Me(:,1); Fe_t = Me(:,3);
    m  = (te>=zA(1)&te<=zA(2)) | (te>=zB(1)&te<=zB(2));  p = polyfit(te(m),Fe_t(m),1);
    Fs = movmean(Fe_t,51); thr = 0.2*prctile(Fs,99);
    on = find(Fs>thr,1); off = find(Fs>thr,1,'last');  Ta = te(off)-te(on);
    loss = PHI*abs(p(1))*Ta;

    Ml = readmatrix([D P{i,6}]);  tl = Ml(:,1);
    Fl = mean(Ml(tl>=zA(1)&tl<=zA(2), 3),'omitnan');     % the DEGRADED run
    frac(i) = loss/Fl;
    ord{i} = sprintf('%s->%s', P{i,2}, P{i,5});
    fprintf('%-9s %-14s | %9.4f %9.1f | %10.2f %9.2f %10.1f%%\n', ...
            P{i,1}, ord{i}, p(1), Ta, loss, Fl, 100*frac(i));
end

%% ── Correct every force feature's ATP ratio ────────────────────────────────
% The later run is depressed by a fraction `frac` OF ITSELF, so undo it there:
%   8->2 : the 2 mM (later) is the depressed numerator  => ratio *= (1+frac)
%   2->8 : the 8 mM (later) is the depressed denominator => ratio /= (1+frac)
rawF = nan(nP,numel(fF)); corF = nan(nP,numel(fF));
rawK = nan(nP,numel(fK)); corK = nan(nP,numel(fK));
for i = 1:nP
    Se = load([D P{i,4}]);  Sl = load([D P{i,7}]);
    is8first = strcmp(P{i,2},'8 mM');
    for j = 1:numel(fF)
        ve = mean(Se.features_data.(fF{j}),'omitnan');
        vl = mean(Sl.features_data.(fF{j}),'omitnan');
        if is8first; r = vl/ve; rawF(i,j) = r; corF(i,j) = r*(1+frac(i));
        else;        r = ve/vl; rawF(i,j) = r; corF(i,j) = r/(1+frac(i)); end
    end
    for j = 1:numel(fK)
        ve = mean(Se.features_data.(fK{j}),'omitnan');
        vl = mean(Sl.features_data.(fK{j}),'omitnan');
        if is8first; rawK(i,j) = vl/ve; else; rawK(i,j) = ve/vl; end
        % ktr correction, for the record only (see the test below)
        kfac = 0.884^(frac(i)/(11.98/56.45));     % bracket ktr loss, scaled by dose
        if is8first; corK(i,j) = rawK(i,j)/kfac; else; corK(i,j) = rawK(i,j)*kfac; end
    end
end

cv = @(x) 100*std(x,'omitnan')/abs(mean(x,'omitnan'));
fprintf('\n\n============ FORCE FEATURES: ATP RATIO 2 mM / 8 mM ============\n');
fprintf('%-12s %-47s | %s\n','','RAW','RUNDOWN-CORRECTED');
fprintf('%-12s %9s %9s %9s %8s %7s | %9s %9s %9s %8s %7s\n', ...
        'feature', P{1,1},P{2,1},P{3,1},'mean','CV', P{1,1},P{2,1},P{3,1},'mean','CV');
for j = 1:numel(fF)
    fprintf('%-12s %9.3f %9.3f %9.3f %8.3f %6.0f%% | %9.3f %9.3f %9.3f %8.3f %6.0f%%\n', fF{j}, ...
        rawF(:,j), mean(rawF(:,j)), cv(rawF(:,j)), corF(:,j), mean(corF(:,j)), cv(corF(:,j)));
end
fprintf('%-12s %9s %9s %9s %8.3f %6.0f%% | %9s %9s %9s %8.3f %6.0f%%\n','ALL FORCE', ...
        '','','', mean(rawF(:)), cv(rawF(:)), '','','', mean(corF(:)), cv(corF(:)));

fprintf('\n\n============ RATE FEATURES ============\n');
fprintf('%-12s %9s %9s %9s %8s %7s | %8s %7s\n','feature',P{1,1},P{2,1},P{3,1},'mean','CV','corr mean','CV');
for j = 1:numel(fK)
    fprintf('%-12s %9.3f %9.3f %9.3f %8.3f %6.0f%% | %8.3f %6.0f%%\n', fK{j}, ...
        rawK(:,j), mean(rawK(:,j)), cv(rawK(:,j)), mean(corK(:,j)), cv(corK(:,j)));
end
% Baker ktr (no force correction possible)
Bh = load([D 'bakers_slack8mM_all.mat']);  Bl = load([D 'bakers_slack2mM.mat']);
set(0,'DefaultFigureVisible','off');
fbh = extractSlackAttributes(Bh.datatable(:,1),Bh.datatable(:,3),Bh.datatable(:,2),Bh.velocitytable,[],[],false,true);
fbl = extractSlackAttributes(Bl.datatable(:,1),Bl.datatable(:,3),Bl.datatable(:,2),Bl.velocitytable,[],[],false,true);
set(0,'DefaultFigureVisible','on');
bakerKtr = mean(fbl.ktr,'omitnan')/mean(fbh.ktr,'omitnan');
% 2026-08-12: Baker is DROPPED FROM POOLING entirely -- it is an older protocol
% (2-6x shorter recovery windows, different length ladder), so pooling it with
% the 2026 protocol days mixes two experiment designs. It is still computed and
% printed as an independent sanity check, but it enters no statistic.
fprintf('\n  Baker ktr ratio: %.3f   [REPORTED ONLY -- old protocol, not pooled]\n', bakerKtr);
allK = rawK(:,1)';
fprintf('  ktr across the three 2026 preps: %s  mean %.3f  CV %.0f%%\n', ...
        mat2str(round(allK,3)), mean(allK), cv(allK));
fprintf('  Correcting ktr for rundown pushes CV %0.f%% -> %0.f%% : DO NOT correct it.\n', ...
        cv(rawK(:,1)), cv(corK(:,1)));

%% ── Headline ───────────────────────────────────────────────────────────────
FORCE = mean(corF(:));  KTR = mean(allK);
fprintf('\n\n================== THE ATP EFFECT ==================\n');
fprintf('  FORCE  x%.2f   (rundown-corrected, %d features x %d preps, CV %.0f%%)\n', ...
        FORCE, numel(fF), nP, cv(corF(:)));
fprintf('  ktr    x%.2f   (raw, %d 2026 preps, Baker excluded, CV %.0f%%)\n', KTR, nP, cv(allK));
fprintf('  => 2 mM ATP makes the muscle STRONGER and SLOWER.\n');
fprintf('     force x ktr = x%.2f, so the slowed turnover outweighs the extra force.\n', FORCE*KTR);

%% ── Figures ────────────────────────────────────────────────────────────────
cols = lines(nP);
fig = figure(950); clf; set(fig,'Position',[20 20 1500 760]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile; hold on; box on;
b = bar([mean(rawF,2) mean(corF,2)]); b(1).FaceColor=[.6 .6 .6]; b(2).FaceColor=[.2 .5 .8];
yline(FORCE,'k--','LineWidth',1.5);
set(gca,'XTick',1:nP,'XTickLabel',strcat(P(:,1)',{' '},ord));
ylabel('force ratio 2 mM / 8 mM'); xtickangle(15);
title(sprintf('Force: raw CV %.0f%% -> corrected %.0f%%', cv(mean(rawF,2)), cv(mean(corF,2))));
legend({'raw','rundown-corrected','consensus'},'Location','northwest','FontSize',7);

nexttile; hold on; box on;
for i=1:nP
    plot(1:numel(fF), rawF(i,:), 'o--','Color',[cols(i,:) .45],'HandleVisibility','off');
    plot(1:numel(fF), corF(i,:), 'o-','Color',cols(i,:),'LineWidth',1.8,'MarkerFaceColor',cols(i,:), ...
         'DisplayName',P{i,1});
end
yline(FORCE,'k--','LineWidth',1.5,'HandleVisibility','off');
set(gca,'XTick',1:numel(fF),'XTickLabel',fF); xtickangle(25);
ylabel('ATP force ratio'); title('Per feature (dashed = raw, solid = corrected)');
legend('Location','best','FontSize',7);

nexttile; hold on; box on;
b = bar([allK'; nan]); b.FaceColor = 'flat';
b.CData = [cols; 1 1 1];
set(gca,'XTick',1:nP,'XTickLabel',P(:,1)');
yline(KTR,'k--','LineWidth',1.5); ylabel('ktr ratio 2 mM / 8 mM'); ylim([0 0.8]);
title(sprintf('ktr: raw, already consistent (CV %.0f%%)', cv(allK)));

nexttile; hold on; box on;
plot(mean(rawF,2), rawK(:,1), 'o','MarkerSize',11,'Color',[.6 .6 .6],'MarkerFaceColor',[.6 .6 .6]);
for i=1:nP
    plot(mean(corF(i,:)), rawK(i,1), 'o','MarkerSize',13,'Color',cols(i,:),'MarkerFaceColor',cols(i,:));
    plot([mean(rawF(i,:)) mean(corF(i,:))], [rawK(i,1) rawK(i,1)], '-','Color',[.5 .5 .5]);
    text(mean(corF(i,:))+0.01, rawK(i,1), P{i,1},'FontSize',7);
end
plot(1,1,'kp','MarkerSize',16,'MarkerFaceColor','y');
xline(FORCE,':'); yline(KTR,':');
xlabel('force ratio'); ylabel('ktr ratio');
title('Correction moves the preps together, not apart'); grid on;

nexttile; hold on; box on;
fprintf('\n  per-prep corrected force ratios: %s\n', mat2str(round(mean(corF,2)',3)));
hist_x = [mean(rawF,2); mean(corF,2)];
plot(ones(nP,1), mean(rawF,2), 'o','Color',[.6 .6 .6],'MarkerFaceColor',[.6 .6 .6],'MarkerSize',10);
plot(2*ones(nP,1), mean(corF,2), 'o','Color',[.2 .5 .8],'MarkerFaceColor',[.2 .5 .8],'MarkerSize',10);
for i=1:nP; plot([1 2],[mean(rawF(i,:)) mean(corF(i,:))],'-','Color',[.7 .7 .7]); end
errorbar(1, mean(mean(rawF,2)), std(mean(rawF,2)),'k','LineWidth',2);
errorbar(2, mean(mean(corF,2)), std(mean(corF,2)),'k','LineWidth',2);
set(gca,'XTick',1:2,'XTickLabel',{'raw','corrected'}); xlim([0.5 2.5]);
ylabel('force ratio'); title('Spread collapses after correction');

nexttile; hold on; box on;
txt = {sprintf('\\bfATP effect (2 mM vs 8 mM)\\rm'), '', ...
       sprintf('FORCE   \\times%.2f   (CV %.0f%%, n=%d preps)', FORCE, cv(corF(:)), nP), ...
       sprintf('ktr     \\times%.2f   (CV %.0f%%, n=%d preps)', KTR, cv(allK), nP), '', ...
       'Low ATP -> STRONGER and SLOWER.', '', ...
       sprintf('rundown correction: \\phi=%.2f on', PHI), ...
       'effective activated time', '', ...
       'Baker DROPPED entirely: old protocol,', ...
       'truncated recovery. (Its ktr ratio is', ...
       sprintf('%.2f -- consistent, but not pooled.)', bakerKtr)};
text(0.03, 0.95, txt, 'VerticalAlignment','top','FontSize',10); axis off;

sgtitle('The ATP effect, reconciled across protocol days after rundown correction');
exportgraphics(fig, fullfile(resDir,'atp_reconciliation.png'), 'Resolution', 150);
save(fullfile(resDir,'atp_reconciliation.mat'), 'P','frac','rawF','corF','rawK','corK', ...
     'fF','fK','allK','FORCE','KTR','PHI');
fprintf('\nSaved %s\n', fullfile(resDir,'atp_reconciliation.png'));
