% RunDesignStrategy.m
% Which running order, how many repeats, and how to pool -- quantified.
%
% The order bias is NOT symmetric, and that changes the design advice.
%
%   8 -> 2 : the 8 mM run goes first and does the damage. It decays slowly
%            (~0.56 kPa/s) and the degraded run (2 mM) is the STRONG one, so the
%            same absolute damage is a small fraction of it.
%   2 -> 8 : the 2 mM run goes first. It decays ~1.95x faster (~1.09 kPa/s) AND
%            the degraded run (8 mM) is the WEAK one, so the same damage is a
%            much larger fraction.
%
% Both effects compound, so 2->8 carries roughly 3x the bias of 8->2. A
% "balanced design" therefore does NOT cancel the bias -- it averages a small
% negative bias against a large positive one.
%
% Outputs -> results/design_strategy.png, .mat

clear; close all; clc;
cd(fileparts(which('RunDesignStrategy')));
addpath(genpath('../..'));
resDir = 'results';  if ~isfolder(resDir), mkdir(resDir); end

%% ── Measured inputs (../RundownCorrection, ../ATPEffectReconciliation) ─────
PHI    = 0.581;              % permanent fraction, own-peak / T_dec footing
s8     = 0.560;  T8 = 10.4;  % 8 mM decay rate (kPa/s) and decay-phase duration
s2     = 1.092;  T2 = 11.6;  % 2 mM
F8     = 50;                 % typical 8 mM plateau force (kPa)
F2     = 70;                 % typical 2 mM plateau force
aTrue  = 1.344;              % rundown-corrected ATP force effect
sdPrep = 0.053;              % residual between-prep SD after correction

% ── The bias, derived explicitly ───────────────────────────────────────────
% DEFINITIONS
%   a       the TRUE ATP effect = F_2mM / F_8mM with BOTH measured at the SAME
%           fibre state. This is the quantity we want.
%   L       absolute damage (kPa) accrued during the EARLIER run's activation,
%           L = phi * |slope_earlier| * T_earlier.  Rundown is linear in kPa/s,
%           so the damage is ABSOLUTE, not fractional -- this is what makes the
%           bias asymmetric.
%   The LATER run is the degraded one:  F_later_measured = F_later_true - L.
%   f       = L / F_later_measured   (damage as a fraction of what was measured)
%
% BIAS = R_measured / a   -- the measured 2mM/8mM ratio against the true effect
%
%   8 -> 2 : R_meas = F2_meas / F8 ,  F2_true = F2_meas + L
%            a = (F2_meas + L)/F8 = R_meas*(1+f)   =>  BIAS = 1/(1+f)   < 1
%            the 2 mM numerator is the degraded one -> ratio UNDER-states a
%
%   2 -> 8 : R_meas = F2 / F8_meas ,  F8_true = F8_meas + L
%            a = F2/(F8_meas + L) = R_meas/(1+f)   =>  BIAS = (1+f)     > 1
%            the 8 mM denominator is the degraded one -> ratio OVER-states a
L82 = PHI*s8*T8;   f82 = L82/F2;   b82 = 1/(1+f82);
L28 = PHI*s2*T2;   f28 = L28/F8;   b28 = 1+f28;

fprintf('=== HOW THE BIAS ARISES, AND WHAT IS BIASED AGAINST WHAT ===\n');
fprintf('  BIAS = (measured 2mM/8mM force ratio) / (true ATP effect at fixed fibre state)\n\n');
fprintf('  8->2 : the 8 mM runs first and does the damage.\n');
fprintf('         L    = phi*|s8|*T8 = %.2f*%.3f*%.1f = %.2f kPa\n', PHI,s8,T8,L82);
fprintf('         the degraded run is the 2 mM one, F ~ %.0f kPa -> f = %.3f\n', F2, f82);
fprintf('         the DEGRADED run is the NUMERATOR -> ratio too SMALL\n');
fprintf('         BIAS = 1/(1+f) = %.3f   (%+.1f%%)\n\n', b82, 100*(b82-1));
fprintf('  2->8 : the 2 mM runs first and does the damage.\n');
fprintf('         L    = phi*|s2|*T2 = %.2f*%.3f*%.1f = %.2f kPa\n', PHI,s2,T2,L28);
fprintf('         the degraded run is the 8 mM one, F ~ %.0f kPa -> f = %.3f\n', F8, f28);
fprintf('         the DEGRADED run is the DENOMINATOR -> ratio too LARGE\n');
fprintf('         BIAS = (1+f)   = %.3f   (%+.1f%%)\n\n', b28, 100*(b28-1));
fprintf('  ASYMMETRY = %.1fx, from two compounding factors:\n', (b28-1)/(1-b82));
fprintf('    (i)  numerator : the 2 mM run decays x%.2f faster, so L is x%.2f larger\n', s2/s8, L28/L82);
fprintf('    (ii) denominator: it degrades the 8 mM run, which is x%.2f weaker,\n', F2/F8);
fprintf('         so the same kPa is a larger fraction\n');
fprintf('    net f: %.3f vs %.3f = x%.1f\n', f28, f82, f28/f82);
fprintf('\n  (F8=%.0f, F2=%.0f are representative plateau forces; the real bias is\n', F8, F2);
fprintf('   prep-specific and computed per prep in ../ATPEffectReconciliation.)\n');

%% ── Designs, with and without correction ──────────────────────────────────
% "correction" = the model-lesion correction, assumed to remove the bias with a
% relative error eCorr (it removed it to CV 4 % on 3 preps, so eCorr ~ 0.25 of
% the bias is a conservative residual).
eCorr = 0.25;
designs = { '6x  8->2',            6, 0
            '5x 8->2, 1x 2->8',    5, 1
            '4x 8->2, 2x 2->8',    4, 2
            '3x 8->2, 3x 2->8',    3, 3
            '2x 8->2, 4x 2->8',    2, 4
            '6x  2->8',            0, 6 };

fprintf('\n=== DESIGN COMPARISON (n = 6 preparations) ===\n');
fprintf('%-20s %12s %12s %12s %10s\n','design','bias RAW','bias CORR','SEM','can test?');
Res = nan(size(designs,1),3);
for i = 1:size(designs,1)
    n82 = designs{i,2};  n28 = designs{i,3};  n = n82+n28;
    mRaw  = (n82*aTrue*b82 + n28*aTrue*b28)/n;
    mCorr = (n82*aTrue*(1+eCorr*(b82-1)) + n28*aTrue*(1+eCorr*(b28-1)))/n;
    sem   = sdPrep/sqrt(n);
    Res(i,:) = [100*(mRaw/aTrue-1), 100*(mCorr/aTrue-1), sem];
    if n82>0 && n28>0; tst = 'YES'; else; tst = 'no'; end
    fprintf('%-20s %11.1f%% %11.1f%% %12.3f %10s\n', designs{i,1}, Res(i,1), Res(i,2), sem, tst);
end
fprintf('\n  "can test?" = both orders present, so a systematic error in the\n');
fprintf('  correction shows up as a difference between the two order groups.\n');
fprintf('\n  NOTE: a balanced design does NOT cancel the bias (%.1f%% at 3+3),\n', Res(4,1));
fprintf('  because the two biases are unequal. All-8->2 is the least biased RAW\n');
fprintf('  design (%.1f%%) but cannot detect a broken correction.\n', Res(1,1));

%% ── Value of a repeat ──────────────────────────────────────────────────────
fprintf('\n=== VALUE OF AN END-OF-SESSION REPEAT ===\n');
fprintf('  Without a repeat, phi must be ASSUMED (transferred from 03/27).\n');
fprintf('  With one, the damage is MEASURED for that preparation.\n\n');
fprintf('  8-2-8 : brackets the 2 mM run in 8 mM units. Directly gives the\n');
fprintf('          correction actually needed for the 8->2 comparison.  <-- most useful\n');
fprintf('  2-8-2 : brackets the 8 mM run in 2 mM units. Tests whether phi itself\n');
fprintf('          is ATP-dependent -- likely, since 2 mM decays x%.2f faster.\n', s2/s8);
fprintf('  Having ONE OF EACH tests the transferability of phi, which is currently\n');
fprintf('  the single largest assumption in the whole correction.\n');

%% ── Figure ─────────────────────────────────────────────────────────────────
fig = figure(1020); clf; set(fig,'Position',[30 30 1400 640]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile; hold on; box on;
b = bar([100*(1-b82) 100*(b28-1)]); b.FaceColor='flat'; b.CData=[0 0.45 0.74; 0.85 0.33 0.10];
set(gca,'XTick',1:2,'XTickLabel',{'8\rightarrow2','2\rightarrow8'});
ylabel('|bias| on the measured ATP ratio (%)');
title(sprintf('(a) 2\\rightarrow8 carries %.1fx the bias', (b28-1)/(1-b82)));

nexttile; hold on; box on;
b = bar([Res(:,1) Res(:,2)]); b(1).FaceColor=[.6 .6 .6]; b(2).FaceColor=[.2 .5 .8];
yline(0,'k-'); set(gca,'XTick',1:size(designs,1),'XTickLabel',designs(:,1)); xtickangle(30);
ylabel('bias in the pooled ATP effect (%)');
title('(b) Balancing does NOT cancel the bias');
legend({'uncorrected','corrected'},'Location','northwest','FontSize',7);

nexttile; hold on; box on;
nn = 3:9;
plot(nn, sdPrep./sqrt(nn), 'ko-','LineWidth',2,'MarkerFaceColor','k');
xlabel('number of preparations'); ylabel('SEM of the pooled ATP effect');
title('(c) Precision vs n (random error only)');
for k=[3 6 9]; text(k, sdPrep/sqrt(k)+0.001, sprintf('%.3f',sdPrep/sqrt(k)),'FontSize',8); end

nexttile; hold on; box on;
xq = 0:0.02:0.30;
plot(100*xq, 100*(1./(1+xq)-1), '-','Color',[0 0.45 0.74],'LineWidth',2);
plot(100*xq, 100*((1+xq)-1),    '-','Color',[0.85 0.33 0.10],'LineWidth',2);
plot(100*f82, 100*(b82-1),'o','Color',[0 0.45 0.74],'MarkerFaceColor',[0 0.45 0.74],'MarkerSize',10);
plot(100*f28, 100*(b28-1),'o','Color',[0.85 0.33 0.10],'MarkerFaceColor',[0.85 0.33 0.10],'MarkerSize',10);
yline(0,'k:'); xlabel('f = damage / force of the DEGRADED run (%)'); ylabel('bias (%)');
title('(d) Same f gives a bigger bias when it hits the denominator');
legend({'8\rightarrow2: degraded run is the numerator','2\rightarrow8: degraded run is the denominator'}, ...
       'Location','northwest','FontSize',7);

nexttile; hold on; box on;
txt = {'\bfRECOMMENDED DESIGN\rm','', ...
       'Majority \bf8\rightarrow2\rm (3x smaller bias)','', ...
       'At least \bfone 2\rightarrow8\rm, to detect a', ...
       'broken correction','', ...
       'At least \bfone 8-2-8 repeat\rm, to MEASURE', ...
       '\phi instead of assuming it','', ...
       'If a second repeat is possible, make it', ...
       '\bf2-8-2\rm to test whether \phi is ATP-dependent','', ...
       'For the 3 incoming datasets:', ...
       '  2x  8-2  (one with 8-2-8 repeat)', ...
       '  1x  2-8  (with 2-8-2 repeat if possible)'};
text(0.02,0.97,txt,'VerticalAlignment','top','FontSize',9.5); axis off;

nexttile; hold on; box on;
txt2 = {'\bfWHY NOT SIMPLY BALANCE?\rm','', ...
        sprintf('8\\rightarrow2 bias  %+.1f%%', 100*(b82-1)), ...
        sprintf('2\\rightarrow8 bias  %+.1f%%', 100*(b28-1)), '', ...
        'They are unequal, so 3+3 leaves', ...
        sprintf('%+.1f%% -- WORSE than all-8\\rightarrow2 (%+.1f%%).', Res(4,1), Res(1,1)), '', ...
        'Correction is needed either way.', ...
        'The role of the 2\rightarrow8 preps is to', ...
        'VALIDATE the correction, not to', ...
        'cancel the bias by averaging.'};
text(0.02,0.97,txt2,'VerticalAlignment','top','FontSize',9.5); axis off;

sgtitle('Design strategy — running order, repeats, and pooling');
exportgraphics(fig, fullfile(resDir,'design_strategy.png'), 'Resolution', 150);
save(fullfile(resDir,'design_strategy.mat'),'f82','f28','b82','b28','Res','designs','aTrue','sdPrep');
fprintf('\nSaved %s\n', fullfile(resDir,'design_strategy.png'));
