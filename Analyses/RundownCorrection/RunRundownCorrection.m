% RunRundownCorrection.m
% Build and calibrate the rundown correction, on EFFECTIVE ACTIVATED TIME.
%
% WHY NOT WALL-CLOCK TIME
% A skinned fibre degrades while it is generating force, not while it sits in
% relaxing solution. So the natural clock for rundown is accumulated ACTIVATED
% time, and the natural rate is the force decay each recording shows during its
% own activation (the within-recording slope, kPa/s, measured on the isometric
% plateaus that bracket the slack train).
%
% THE PREMISE, TESTED
% 03/27 recorded 8 mM twice: fresh (t=0) and again at the end of the session.
% Their within-run slopes are -0.4592 and -0.4618 kPa/s -- essentially IDENTICAL
% in absolute terms even though the second run is 17.5 % weaker. An exponential
% (force-proportional) rundown would have predicted -0.379. So the decay is
% LINEAR in effective time with a force-independent kPa/s rate, and the gap
% between two runs is naturally expressed as an EFFECTIVE TIME
%       tau_eff = dF / |slope|
% rather than as a real elapsed delay. (This is the T0 quantity that
% DataCuration/FitRundownCorrection.m already computes.)
%
% THE ONE FREE PARAMETER
% Integrating the within-run slope over the measured activation duration
% over-predicts the bracket loss ~2.3x, so most of the within-run decline is
% REVERSIBLE (product accumulation that washes out between activations) and only
% a fraction is permanent damage:
%       loss between runs = phi * SUM_i |slope_i| * T_act_i
% over the intervening force-generating activations. phi is calibrated two
% independent ways here -- from the 03/27 bracket, and from cross-prep
% consistency of the ATP effect -- and they are compared.
%
% Outputs -> results/ : rundown_correction.png, rundown_correction.mat
% See also: RunRundownMechanisms.m (what rundown IS, mechanically),
%           ../ATPEffectReconciliation (what this correction buys),
%           ../../DataCuration/FitRundownCorrection.m (the r(F,SL) surface).

clear; close all; clc;
cd(fileparts(which('RunRundownCorrection')));
addpath(genpath('../..'));
resDir = 'results';  if ~isfolder(resDir), mkdir(resDir); end
D = '../../data/';

%% ── Run registry ───────────────────────────────────────────────────────────
% {day, condition, merged file, activation index, elapsed s}
R = {
 '03/27 M','8 mM',    '03 27 2026 M/02_Merged_8mM_Active.txt',            1,   0
 '03/27 M','2 mM',    '03 27 2026 M/03_Merged_2mM_Active.txt',            2, 132
 '03/27 M','PNB+Mava','03 27 2026 M/06_Merged_8mM_Active_PNB_Mava.txt',   3, 348
 '03/27 M','8 mM rpt','03 27 2026 M/04_Merged_8mM_Active_repeat.txt',     4, 694
 '04/03 F','8 mM',    '04 03 2026 F/02_Merged_8mM_Active.txt',            1,   0
 '04/03 F','2 mM',    '04 03 2026 F/03_Merged_2mM_Active.txt',            2, 138
 '04/03 F','PNB+Mava','04 03 2026 F/04_Merged_8mM_Active_PNB_Mava.txt',   3, 288
 '04/10 M','2 mM',    '04 10 2026 Male 2-8/02_Merged_2mM_Active.txt',     1,   0
 '04/10 M','8 mM',    '04 10 2026 Male 2-8/03_Merged_8mM_Active.txt',     2, 308
 '04/10 M','PNB+Mava','04 10 2026 Male 2-8/04_Merged_8mM_Active_PNB_Mava.txt', 3, 687};

zA  = [71.50 72.35];   % isometric ML 1.0 (SL 2.0), before the staircase
z22 = [74.00 74.40];   % isometric ML 1.1 (SL 2.2), just before slack 1
zB  = [77.70 78.30];   % isometric ML 1.0, after the last restretch

n = size(R,1);
F20 = nan(n,1); F22 = nan(n,1); slope = nan(n,1); Tact = nan(n,1); isAct = false(n,1);
fprintf('%-9s %-9s %4s | %8s %8s | %11s %8s | %8s\n', ...
        'day','cond','act','F@SL2.0','F@SL2.2','slope kPa/s','%%/s','T_act s');
for i = 1:n
    M = readmatrix([D R{i,3}]);  t = M(:,1);  F = M(:,3);
    F20(i) = mean(F(t>=zA(1)  & t<=zA(2)),  'omitnan');
    F22(i) = mean(F(t>=z22(1) & t<=z22(2)), 'omitnan');
    m = (t>=zA(1)&t<=zA(2)) | (t>=zB(1)&t<=zB(2));
    p = polyfit(t(m), F(m), 1);  slope(i) = p(1);
    Fs = movmean(F,51);  thr = 0.2*prctile(Fs,99);
    on = find(Fs>thr,1); off = find(Fs>thr,1,'last');
    Tact(i) = t(off)-t(on);
    isAct(i) = ~contains(R{i,2},'PNB');           % PNB+Mava generates ~no force
    if ~isAct(i); Tact(i) = 0; end                % ...so it does no damage
    fprintf('%-9s %-9s %4d | %8.2f %8.2f | %11.4f %7.2f%% | %8.1f\n', ...
            R{i,1}, R{i,2}, R{i,4}, F20(i), F22(i), slope(i), 100*slope(i)/F20(i), Tact(i));
end
fprintf('\n  Activation lasts a very consistent ~%.1f s in every force-generating run.\n', ...
        mean(Tact(isAct)));
fprintf('  PNB+Mava is blebbistatin-suppressed (slope ~%.3f kPa/s) -> charged zero damage.\n', ...
        mean(slope(~isAct)));

%% ── 1. Is rundown linear or exponential in effective time? ─────────────────
i1 = 1; i4 = 4;                                   % 03/27 fresh and repeat
fprintf('\n\n============ 1. LINEAR vs EXPONENTIAL ============\n');
fprintf('  03/27 8 mM fresh : F = %.2f kPa, slope = %.4f kPa/s\n', F20(i1), slope(i1));
fprintf('  03/27 8 mM repeat: F = %.2f kPa, slope = %.4f kPa/s\n', F20(i4), slope(i4));
fprintf('  exponential (slope proportional to force) predicts slope_rpt = %.4f\n', ...
        slope(i1)*F20(i4)/F20(i1));
fprintf('  linear      (slope force-independent)      predicts slope_rpt = %.4f\n', slope(i1));
fprintf('  MEASURED %.4f  ->  LINEAR. Rundown is a constant kPa/s while activated.\n', slope(i4));
tau_eff = (F20(i1)-F20(i4))/abs(slope(i1));
fprintf('\n  => effective time between the two 8 mM runs: tau_eff = dF/|slope| = %.1f s\n', tau_eff);
fprintf('     (real elapsed %.0f s; actual activation in between %.1f s)\n', ...
        R{i4,5}-R{i1,5}, sum(Tact(2:3)));

%% ── 2. How much of the within-run decline is permanent? ────────────────────
predFull = sum(abs(slope(2:3)).*Tact(2:3)) + abs(slope(i1))*Tact(i1);
measLoss = F20(i1)-F20(i4);
phi20 = measLoss / predFull;
phi22 = (F22(i1)-F22(i4)) / (sum(abs(slope(2:3)).*Tact(2:3)) + abs(slope(i1))*Tact(i1));
fprintf('\n\n============ 2. REVERSIBLE vs PERMANENT ============\n');
fprintf('  integrating each run''s own slope over its activation predicts %.2f kPa of loss\n', predFull);
fprintf('  the bracket measures only %.2f kPa\n', measLoss);
fprintf('  => phi = permanent fraction = %.3f at SL 2.0  (%.3f at SL 2.2)\n', phi20, phi22);
fprintf('  i.e. ~%.0f%% of the within-run force decline RECOVERS between activations\n', 100*(1-phi20));
fprintf('     (product accumulation), ~%.0f%% is permanent damage.\n', 100*phi20);

%% ── 3. Correct each prep''s ATP ratio, and scan phi ─────────────────────────
% Damage accrued during the EARLIER of the two runs is what separates them.
pairs = {1,2,'03/27 M'; 5,6,'04/03 F'; 8,9,'04/10 M'};
Fsel  = {F20, F22};  zname = {'SL 2.0','SL 2.2'};

phis = linspace(0, 1.2, 241);
CV = nan(numel(phis), 2);  RAT = nan(numel(phis), 3, 2);
for z = 1:2
    Fz = Fsel{z};
    for j = 1:numel(phis)
        rr = nan(1,3);
        for q = 1:3
            ie = pairs{q,1}; il = pairs{q,2};                 % earlier, later
            loss = phis(j)*abs(slope(ie))*Tact(ie);
            if strcmp(R{ie,2},'8 mM')                          % 8 -> 2
                rr(q) = Fz(il) / (Fz(ie) - loss);
            else                                               % 2 -> 8 (reversed)
                rr(q) = Fz(ie) / (Fz(il) + loss);
            end
        end
        RAT(j,:,z) = rr;  CV(j,z) = 100*std(rr)/mean(rr);
    end
end
[~, jOpt] = min(CV(:,1));  phiOpt = phis(jOpt);
[~, jBr]  = min(abs(phis - phi20));

fprintf('\n\n============ 3. ATP RATIO vs phi ============\n');
for z = 1:2
    fprintf('  --- %s ---\n', zname{z});
    fprintf('  %-28s %8s %8s %8s %8s %7s\n','phi', pairs{1,3}, pairs{2,3}, pairs{3,3}, 'mean','CV');
    for j = [1 jBr jOpt numel(phis)]
        fprintf('  phi=%.2f %-20s %8.3f %8.3f %8.3f %8.3f %6.1f%%\n', phis(j), '', ...
                RAT(j,1,z), RAT(j,2,z), RAT(j,3,z), mean(RAT(j,:,z)), CV(j,z));
    end
end
fprintf('\n  phi from the 03/27 BRACKET (mechanical, independent) = %.3f\n', phi20);
fprintf('  phi minimising CROSS-PREP CV                         = %.3f\n', phiOpt);
fprintf('  Two unrelated calibrations land in the same region, and the corrected\n');
fprintf('  ATP effect is ~%.2f across that whole range (CV %.0f-%.0f%% vs %.0f%% raw).\n', ...
        mean([mean(RAT(jBr,:,1)) mean(RAT(jOpt,:,1))]), CV(jOpt,1), CV(jBr,1), CV(1,1));

%% ── 4. Variant: universal damage rate instead of slope-proportional ────────
% Alternative reading: every activation does the same damage per second (rate r),
% and the steeper 2 mM within-run slope is entirely reversible fatigue.
rUni = measLoss / sum(Tact(1:3));
rrU = nan(1,3);
for q = 1:3
    ie = pairs{q,1}; il = pairs{q,2};
    loss = rUni*Tact(ie);
    if strcmp(R{ie,2},'8 mM'); rrU(q) = F20(il)/(F20(ie)-loss);
    else;                      rrU(q) = F20(ie)/(F20(il)+loss); end
end
fprintf('\n\n============ 4. SENSITIVITY TO THE DAMAGE MODEL ============\n');
fprintf('  A) damage proportional to the run''s own slope (phi=%.2f): %s  mean %.3f  CV %.1f%%\n', ...
        phi20, mat2str(round(RAT(jBr,:,1),3)), mean(RAT(jBr,:,1)), CV(jBr,1));
fprintf('  B) universal damage rate r = %.3f kPa/s              : %s  mean %.3f  CV %.1f%%\n', ...
        rUni, mat2str(round(rrU,3)), mean(rrU), 100*std(rrU)/mean(rrU));
fprintf('  Both variants give ~1.3-1.4 and a large CV improvement, so the ATP number\n');
fprintf('  is robust to which of the two damage models is assumed.\n');

%% ── 5. Figures ─────────────────────────────────────────────────────────────
fig = figure(930); clf; set(fig,'Position',[20 20 1560 800]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
cols = lines(3); days = {'03/27 M','04/03 F','04/10 M'};

% (a) linear vs exponential
nexttile; hold on; box on;
bar([abs(slope(i1)) abs(slope(i4))]);
yline(abs(slope(i1)*F20(i4)/F20(i1)),'r--','LineWidth',2);
yline(abs(slope(i1)),'g--','LineWidth',2);
set(gca,'XTick',1:2,'XTickLabel',{sprintf('fresh (%.1f kPa)',F20(i1)), sprintf('repeat (%.1f kPa)',F20(i4))});
ylabel('|within-run slope| (kPa/s)'); ylim([0 0.55]);
title('Rundown is LINEAR, not exponential');
legend({'measured','exponential prediction','linear prediction'},'Location','south','FontSize',7);

% (b) within-run decline: predicted vs measured permanent loss
nexttile; hold on; box on;
bar([predFull measLoss]);
set(gca,'XTick',1:2,'XTickLabel',{'integrated slope','bracket (measured)'});
ylabel('force lost between the two 8 mM runs (kPa)');
title(sprintf('Only \\phi = %.2f of the decline is permanent', phi20));
text(1.5, measLoss*1.15, sprintf('%.0f%% recovers', 100*(1-phi20)), 'HorizontalAlignment','center','FontSize',9);

% (c) phi scan
nexttile; hold on; box on;
plot(phis, CV(:,1), 'k-', 'LineWidth', 1.8);
plot(phis, CV(:,2), '-', 'Color',[.5 .5 .5], 'LineWidth', 1.2);
xline(phi20,'g--','LineWidth',2); xline(phiOpt,'b--','LineWidth',2);
xlabel('\phi  (permanent fraction)'); ylabel('cross-prep CV of the ATP ratio (%)');
title('Two independent calibrations agree');
legend({'CV at SL 2.0','CV at SL 2.2','\phi from bracket','\phi from consistency'}, ...
       'Location','north','FontSize',7);

% (d) ATP ratio vs phi, per prep
nexttile; hold on; box on;
for q=1:3
    plot(phis, RAT(:,q,1), '-', 'Color', cols(q,:), 'LineWidth', 1.6, 'DisplayName', days{q});
end
xline(phi20,'g--','LineWidth',1.5,'HandleVisibility','off');
xline(phiOpt,'b--','LineWidth',1.5,'HandleVisibility','off');
xlabel('\phi'); ylabel('ATP force ratio  2 mM / 8 mM');
title('The three preps converge as \phi rises'); legend('Location','northwest','FontSize',7);

% (e) before / after
nexttile; hold on; box on;
b = bar([RAT(1,:,1); RAT(jBr,:,1)]', 'grouped');
b(1).FaceColor=[.6 .6 .6]; b(2).FaceColor=[.2 .5 .8];
yline(mean(RAT(jBr,:,1)),'k--','LineWidth',1.2);
set(gca,'XTick',1:3,'XTickLabel',days); ylabel('ATP force ratio');
title(sprintf('raw CV %.0f%%  ->  corrected CV %.0f%%', CV(1,1), CV(jBr,1)));
legend({'raw','corrected (\phi from bracket)','corrected mean'},'Location','northwest','FontSize',7);

% (f) the effective-time picture: every run placed on the damage axis
nexttile; hold on; box on;
fprintf('\n  effective-time placement (s, in each prep''s own 8 mM-slope units):\n');
for q = 1:3
    idx = find(strcmp(R(:,1), days{q}) & isAct);
    ref = abs(slope(idx(1)));
    te  = zeros(numel(idx),1);
    for kk = 2:numel(idx)
        te(kk) = te(kk-1) + phi20*abs(slope(idx(kk-1)))*Tact(idx(kk-1))/ref;
    end
    plot(te, F20(idx), '-', 'Color', [cols(q,:) 0.5], 'HandleVisibility','off');
    for kk = 1:numel(idx)
        if contains(R{idx(kk),2},'2 mM'); mk = 's'; else; mk = 'o'; end
        plot(te(kk), F20(idx(kk)), mk, 'Color',cols(q,:), ...
             'MarkerFaceColor',cols(q,:),'MarkerSize',9);
        text(te(kk)+0.5, F20(idx(kk)), R{idx(kk),2}, 'FontSize',7,'Color',cols(q,:));
        fprintf('    %-9s %-9s  te = %5.1f s   F = %5.2f kPa\n', ...
                days{q}, R{idx(kk),2}, te(kk), F20(idx(kk)));
    end
end
plot([0 tau_eff],[F20(i1) F20(i4)],'k--','LineWidth',1.5);
xlabel('effective activated time (s)'); ylabel('force at SL 2.0 (kPa)');
title('All runs on one damage axis (square = 2 mM, circle = 8 mM)');
xlim([-1 30]);

sgtitle('Rundown correction on effective activated time — 03/27 8 mM bracket calibrates it');
exportgraphics(fig, fullfile(resDir,'rundown_correction.png'), 'Resolution', 150);

save(fullfile(resDir,'rundown_correction.mat'), 'R','F20','F22','slope','Tact','isAct', ...
     'phi20','phi22','phiOpt','phis','CV','RAT','rUni','tau_eff','pairs');
fprintf('\nSaved %s\n', fullfile(resDir,'rundown_correction.png'));
