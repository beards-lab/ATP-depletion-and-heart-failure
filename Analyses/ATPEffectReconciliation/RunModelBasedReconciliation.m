% RunModelBasedReconciliation.m
% Does correcting with the MODEL LESION reconcile the reversed-order prep better
% than the empirical fractional correction?
%
% THE QUESTION
% 03/27 and 04/03 ran 8 mM -> 2 mM; 04/10 ran 2 mM -> 8 mM. The empirical
% correction (../RundownCorrection) fixes force to CV 6 % but cannot touch ktr.
% The model route promises more: the rundown lesion identified in
% RundownCorrection/RunModelRundownFit.m (kstiff x0.84 + kSE x0.65 at the full
% bracket dose) predicts BOTH a force change and a ktr change, so in principle it
% could correct both from one calibrated lesion.
%
% METHOD
% Assume the lesion grows linearly with damage dose:
%     lambda   = frac_damage(gap) / frac_damage(bracket)
%     kstiff  *= 1 - 0.16*lambda        (16 % of heads lost at the full dose)
%     kSE     *= 1 - 0.35*lambda        (35 % more series compliance)
% Simulate at each preparation's own lambda, read the predicted rundown factors
% for force and ktr, and divide them out of the measured ATP ratios.
%
% Outputs -> results/model_based_reconciliation.png, .mat

clear; close all; clc;
cd(fileparts(which('RunModelBasedReconciliation')));
addpath(genpath('../..'));
resDir = 'results';  if ~isfolder(resDir), mkdir(resDir); end
D = '../../data/';

% damage of the DEGRADED run, as a fraction of its own force (RunJointCorrection.m)
days   = {'03/27 M','04/03 F','04/10 M'};
fracG  = [0.029 0.064 0.141];        % gap damage
fracBr = 0.169;                      % bracket damage (11.59 kPa on F_pk 68.64)
lam    = fracG / fracBr;
order8first = [true true false];     % 04/10 is reversed

% measured ATP ratios (own-peak referenced force; raw ktr)
rawF = [1.269 1.233 1.649];
rawK = [0.537 0.593 0.492];

%% ── Simulate the lesion at each preparation's dose ─────────────────────────
p0 = getParams(loadParams('ThursdayNightFever'), [], true, false);
p0.PlotEachSeparately=0; p0.PlotFullscreen=0; p0.PlotFeatureFitting=0;
p0.ghostSave=''; p0.ghostLoad=''; p0.EvalFeatures=1; p0.RunSlackSegments='All';
p0.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';

fprintf('=== simulating the rundown lesion at each prep''s dose ===\n');
tic; [~,~,fb] = runSlackExperiment(p0);
Fb = mean(fb.A,'omitnan');  Kb = mean(fb.ktr,'omitnan');
fprintf('  baseline                       %5.1f s   A %.2f  ktr %.2f\n', toc, Fb, Kb);

rdF = nan(1,3); rdK = nan(1,3);
for d = 1:3
    p = p0;
    p.kstiff1 = p.kstiff1*(1-0.16*lam(d));
    p.kstiff2 = p.kstiff2*(1-0.16*lam(d));
    p.kSE     = p.kSE    *(1-0.35*lam(d));
    tic; [~,~,fm] = runSlackExperiment(p);
    rdF(d) = mean(fm.A,'omitnan')/Fb;  rdK(d) = mean(fm.ktr,'omitnan')/Kb;
    fprintf('  %-9s lambda %.2f  kstiff x%.3f kSE x%.3f | %5.1f s | force x%.3f  ktr x%.3f\n', ...
            days{d}, lam(d), 1-0.16*lam(d), 1-0.35*lam(d), toc, rdF(d), rdK(d));
end

%% ── Divide the predicted rundown out of the measured ratios ────────────────
corF = nan(1,3); corK = nan(1,3);
for d = 1:3
    if order8first(d)                 % later (degraded) run is the 2 mM numerator
        corF(d) = rawF(d)/rdF(d);   corK(d) = rawK(d)/rdK(d);
    else                              % later (degraded) run is the 8 mM denominator
        corF(d) = rawF(d)*rdF(d);   corK(d) = rawK(d)*rdK(d);
    end
end
cv = @(x) 100*std(x)/mean(x);

fprintf('\n=== FORCE ===\n');
fprintf('%-9s %10s %12s %12s\n','day','raw','model-corr','empirical');
emp = [1.306 1.311 1.445];            % ../RundownCorrection/RunJointCorrection.m
for d=1:3; fprintf('%-9s %10.3f %12.3f %12.3f\n', days{d}, rawF(d), corF(d), emp(d)); end
fprintf('%-9s %10.3f %12.3f %12.3f\n','mean',mean(rawF),mean(corF),mean(emp));
fprintf('%-9s %9.0f%% %11.0f%% %11.0f%%\n','CV',cv(rawF),cv(corF),cv(emp));

fprintf('\n=== ktr ===\n');
fprintf('%-9s %10s %12s\n','day','raw','model-corr');
for d=1:3; fprintf('%-9s %10.3f %12.3f\n', days{d}, rawK(d), corK(d)); end
fprintf('%-9s %10.3f %12.3f\n','mean',mean(rawK),mean(corK));
fprintf('%-9s %9.0f%% %11.0f%%\n','CV',cv(rawK),cv(corK));

fprintf('\n=== VERDICT ===\n');
if cv(corF) <= cv(rawF); fprintf('  FORCE: model correction helps (CV %.0f%% -> %.0f%%)\n', cv(rawF), cv(corF));
else;                    fprintf('  FORCE: model correction does NOT help (CV %.0f%% -> %.0f%%)\n', cv(rawF), cv(corF)); end
if cv(corK) <= cv(rawK); fprintf('  ktr  : model correction helps (CV %.0f%% -> %.0f%%)\n', cv(rawK), cv(corK));
else;                    fprintf('  ktr  : model correction does NOT help (CV %.0f%% -> %.0f%%)\n', cv(rawK), cv(corK)); end
fprintf('\n  The lesion''s force:ktr coupling is %.2f, essentially the same ratio the\n', ...
        (1-mean(rdK))/(1-mean(rdF)));
fprintf('  empirical scaling assumed -- so the model route inherits the same problem.\n');
fprintf('  04/10 has the LARGEST damage but the LOWEST raw ktr ratio, so any\n');
fprintf('  dose-proportional ktr correction moves it further from the others.\n');

%% ── Figure ─────────────────────────────────────────────────────────────────
fig = figure(1000); clf; set(fig,'Position',[30 30 1350 640]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
cols = lines(3);

nexttile; hold on; box on;
b = bar([rawF; corF; emp]'); b(1).FaceColor=[.6 .6 .6]; b(2).FaceColor=[.2 .5 .8]; b(3).FaceColor=[.4 .75 .4];
set(gca,'XTick',1:3,'XTickLabel',days); ylabel('ATP force ratio');
title(sprintf('FORCE  CV: raw %.0f%%, model %.0f%%, empirical %.0f%%', cv(rawF), cv(corF), cv(emp)));
legend({'raw','model-based','empirical'},'Location','northwest','FontSize',7);

nexttile; hold on; box on;
b = bar([rawK; corK]'); b(1).FaceColor=[.6 .6 .6]; b(2).FaceColor=[.85 .33 .1];
set(gca,'XTick',1:3,'XTickLabel',days); ylabel('ktr ratio'); ylim([0 0.8]);
title(sprintf('ktr  CV: raw %.0f%% -> model %.0f%%', cv(rawK), cv(corK)));
legend({'raw','model-based'},'Location','northwest','FontSize',7);

nexttile; hold on; box on;
plot(lam, rdF, 'o-','Color',[.2 .5 .8],'LineWidth',2,'MarkerFaceColor',[.2 .5 .8]);
plot(lam, rdK, 's-','Color',[.85 .33 .1],'LineWidth',2,'MarkerFaceColor',[.85 .33 .1]);
xlabel('\lambda  (dose / bracket dose)'); ylabel('predicted rundown factor');
title('What the lesion does at each dose');
legend({'force','ktr'},'Location','southwest','FontSize',7);
for d=1:3; text(lam(d), rdF(d)+0.012, days{d},'FontSize',7,'HorizontalAlignment','center'); end

nexttile; hold on; box on;
plot(fracG*100, rawK, 'ko','MarkerSize',10,'MarkerFaceColor','w');
for d=1:3
    plot(fracG(d)*100, rawK(d), 'o','Color',cols(d,:),'MarkerFaceColor',cols(d,:),'MarkerSize',12);
    text(fracG(d)*100+0.4, rawK(d), days{d},'FontSize',8);
end
xlabel('damage over the gap (%)'); ylabel('raw ktr ratio');
title('The obstacle: largest damage, lowest ktr ratio');

nexttile; hold on; box on;
plot(mean(rawF), mean(rawK),'ko','MarkerSize',12,'MarkerFaceColor',[.6 .6 .6]);
for d=1:3
    plot(rawF(d), rawK(d), 'o','Color',cols(d,:),'MarkerFaceColor','w','MarkerSize',10);
    plot(corF(d), corK(d), 'o','Color',cols(d,:),'MarkerFaceColor',cols(d,:),'MarkerSize',10);
    plot([rawF(d) corF(d)],[rawK(d) corK(d)],'-','Color',[cols(d,:) 0.5]);
end
xlabel('force ratio'); ylabel('ktr ratio'); grid on;
title('open = raw, filled = model-corrected');

nexttile; hold on; box on;
txt = {'\bfAnswer\rm','', ...
       sprintf('FORCE  model %.2f (CV %.0f%%)', mean(corF), cv(corF)), ...
       sprintf('       empirical %.2f (CV %.0f%%)', mean(emp), cv(emp)), ...
       '       -> same answer, no better', '', ...
       sprintf('ktr    CV %.0f%% -> %.0f%%  WORSE', cv(rawK), cv(corK)), ...
       '       -> lesion coupling ~0.7, same', ...
       '          problem as the empirical route', '', ...
       'The model route''s value is NOT better', ...
       'reconciliation. It is that the fibre', ...
       'state can be FITTED rather than the', ...
       'data corrected -- which matters once', ...
       'the ATP effect is itself in the model.'};
text(0.03,0.95,txt,'VerticalAlignment','top','FontSize',9.5); axis off;

sgtitle('Model-based vs empirical rundown correction — does it reconcile the reversed-order prep better?');
exportgraphics(fig, fullfile(resDir,'model_based_reconciliation.png'), 'Resolution', 150);
save(fullfile(resDir,'model_based_reconciliation.mat'),'lam','rdF','rdK','rawF','rawK','corF','corK','emp','days');
fprintf('\nSaved %s\n', fullfile(resDir,'model_based_reconciliation.png'));
