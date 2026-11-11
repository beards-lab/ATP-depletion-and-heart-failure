function plotNoiseBudget(B, resDir)
fig = figure('Position',[60 60 1400 760],'Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
lbl = arrayfun(@(b) sprintf('%s %s', b.cond, b.day), B, 'UniformOutput', false);
act = ~strcmp({B.cond},'PNB');

nexttile; hold on; box on
s = find(strcmp({B.type},'postRestretch'));
bar(1:numel(s), [[B(s).fracW]; [B(s).fracC]; [B(s).fracS]]'*100, 'stacked');
set(gca,'XTick',1:numel(s),'XTickLabel',lbl(s),'XTickLabelRotation',40);
ylabel('% of residual variance'); legend({'white','locked','unexplained'},'Box','off','Location','southoutside','Orientation','horizontal');
title('post-re-stretch: what the residual is');

nexttile; hold on; box on
s = find(strcmp({B.type},'postSlack'));
bar(1:numel(s), [[B(s).fracW]; [B(s).fracC]; [B(s).fracS]]'*100, 'stacked');
set(gca,'XTick',1:numel(s),'XTickLabel',lbl(s),'XTickLabelRotation',40);
ylabel('%'); title('post-slack');

nexttile; hold on; box on
% The 5-cycle average still carries sigma_w/sqrt(5) ~ 0.6 kPa of white noise,
% which is the same size as the locked part being looked for. A 5 ms movmean
% takes that down another 4x and leaves the locked structure visible; it is a
% display choice only -- every number in the table is from the raw average.
for i = 1:numel(B)
    if strcmp(B(i).type,'postSlack'); continue; end
    c = [0 .45 .74]; if strcmp(B(i).cond,'2mM'); c=[.85 .33 .10]; end
    if strcmp(B(i).cond,'PNB'); c=[.4 .4 .4]; end
    n = max(3, round(0.005*4000));
    plot(1e3*B(i).t, movmean(B(i).mR, n), '-', 'Color',[c 0.85], 'LineWidth',1.4);
end
yline(0,'k'); xlabel('ms after the re-stretch'); ylabel('kPa');
title('the LOCKED residual (5-cycle average, 5 ms smoothed)');

nexttile; hold on; box on
r2 = [B(act).r2]; rc = [B(act).r2ceil];
scatter(r2, rc, 40, 'filled'); plot([0 1],[0 1],'k--');
yline(0.999,'r-'); xlabel('R^2 achieved by the exponential'); ylabel('R^2 ceiling for a noise-free model');
title('the ceiling'); xlim([.6 1]); ylim([.6 1.02]);

nexttile; hold on; box on
sW = [B.sigW]; sD = [B.sigWdiff];
scatter(sW, sD, 40, 'filled'); mx = max([sW sD])*1.1; plot([0 mx],[0 mx],'k--');
xlabel('\sigma_w  (plateau, >60 Hz)'); ylabel('\sigma_w  (first difference)');
title('two independent noise estimators agree');

nexttile; hold on; box on
tot = [B(act).sigTot]; A = [B(act).A];
bar([mean(100*tot./A), 100*sqrt(mean([B(act).sigW].^2))/mean(A), ...
     100*sqrt(mean([B(act).sigC].^2))/mean(A), 100*sqrt(mean([B(act).sigS].^2))/mean(A)]);
set(gca,'XTick',1:4,'XTickLabel',{'total','white','locked','real'});
ylabel('RMS, % of the rise'); title('in units of the thing being measured');

exportgraphics(fig, fullfile(resDir,'noise_budget.png'), 'Resolution',140);
close(fig);
end
