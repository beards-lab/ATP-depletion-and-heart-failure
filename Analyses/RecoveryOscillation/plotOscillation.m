function plotOscillation(W, resDir)
%PLOTOSCILLATION  One page: the ripple, its spectrum, and its invariances.

cond = {W.cond};  type = {W.type};  f0 = [W.f0];
act  = ~strcmp(cond,'PNB');
col  = containers.Map({'8mM','2mM','PNB'}, {[0 .45 .74],[.85 .33 .10],[.4 .4 .4]});

fig = figure('Position',[60 60 1500 900],'Color','w');
tl  = tiledlayout(fig, 3, 3, 'Padding','compact','TileSpacing','compact');

% (1) one window, data + exponential
nexttile; hold on; box on
i = find(strcmp(cond,'8mM') & strcmp(type,'postRestretch'), 1);
plot(1e3*W(i).t, W(i).z, '-', 'Color',[.6 .6 .6]);
plot(1e3*W(i).t, W(i).pred, '-', 'Color',[0 .45 .74], 'LineWidth',1.8);
xlabel('ms'); ylabel('F - B (kPa)');
title(sprintf('%s %s %s c%d:  R^2 = %.3f', W(i).day, W(i).cond, W(i).type, W(i).cyc, W(i).r2));
legend({'data','single exponential'},'Location','southeast','Box','off');

% (2) the residual
nexttile; hold on; box on
plot(1e3*W(i).t, W(i).res, '-', 'Color',[.85 .33 .10]);
yline(0,'k-'); xlabel('ms'); ylabel('residual (kPa)');
title(sprintf('what is left over:  %.1f Hz, %.2f kPa', W(i).f0, W(i).pkAmp));

% (3) residual spectra, all active windows, normalised
nexttile; hold on; box on
for j = find(act)
    if isempty(W(j).pxx); continue; end
    plot(W(j).fax, W(j).pxx/max(W(j).pxx), '-', 'Color',[col(W(j).cond) 0.25]);
end
xlim([0 200]); set(gca,'YScale','log'); ylim([1e-4 1.5]);
xlabel('Hz'); ylabel('PSD (norm)'); title('S1: every residual, same peak');
xline(median(f0(act),'omitnan'),'k--','LineWidth',1.4);

% (4) S1/S2 f0 by group
nexttile; hold on; box on
grp = strcat(cond,{' '},type); u = unique(grp,'stable');
for i2 = 1:numel(u)
    s = strcmp(grp,u{i2}) & isfinite(f0);
    c = col(strtok(u{i2}));
    scatter(i2 + 0.12*randn(1,nnz(s)), f0(s), 26, c, 'filled', 'MarkerFaceAlpha',.6);
    plot(i2 + [-.25 .25], mean(f0(s))*[1 1], '-', 'Color',c, 'LineWidth',2.5);
end
set(gca,'XTick',1:numel(u),'XTickLabel',strrep(u,' ','\newline'),'XTickLabelRotation',0);
ylabel('f_0 (Hz)'); title('S1+S2: one frequency, ATP-blind'); ylim([0 max(f0)*1.15]);

% (5) S3 PNB
nexttile; hold on; box on
for c = {'8mM','2mM','PNB'}
    s = strcmp(cond,c{1}) & isfinite(f0);
    if ~any(s); continue; end
    histogram(f0(s), 0:5:150, 'FaceColor',col(c{1}), 'FaceAlpha',.55, 'EdgeColor','none');
end
xlabel('f_0 (Hz)'); ylabel('windows'); legend({'8 mM','2 mM','PNB'},'Box','off');
title('S3: cross-bridges off, ripple stays');

% (6) S4 amplitude vs excitation
nexttile; hold on; box on
A = [W.A]; pk = [W.pkAmp]; s = act & isfinite(A) & isfinite(pk);
scatter(A(s), pk(s), 30, cell2mat(cellfun(@(c) col(c), cond(s)','UniformOutput',false)), 'filled');
xlabel('rise amplitude A (kPa)'); ylabel('ripple amplitude (kPa)');
if nnz(s)>3
    pf = polyfit(A(s),pk(s),1); xl = xlim; plot(xl, polyval(pf,xl),'k--');
    title(sprintf('S4: r = %+.2f, slope %.3f', corr(A(s)',pk(s)'), pf(1)));
end

% (7) S5 length channel
nexttile; hold on; box on
f0L = [W.f0L]; sL = isfinite(f0L)&isfinite(f0);
scatter(f0(sL), f0L(sL), 30, [.2 .2 .2], 'filled', 'MarkerFaceAlpha',.5);
xl = [0 max(f0)*1.1]; plot(xl,xl,'k--'); xlim(xl); ylim(xl);
xlabel('force f_0 (Hz)'); ylabel('length f_0 (Hz)'); title('S5: does the motor ring too?');

% (8) R2 ceiling
nexttile; hold on; box on
r2 = [W.r2];
for tp = {'postSlack','postRestretch'}
    s = act & strcmp(type,tp{1}) & isfinite(r2);
    histogram(r2(s), 0:0.025:1, 'FaceAlpha',.55, 'EdgeColor','none');
end
xline(0.999,'r-','LineWidth',1.5);
xlabel('R^2 of a pure exponential'); ylabel('windows');
legend({'post-slack','post-re-stretch','model'},'Box','off','Location','northwest');
title('the ceiling the model is scored against');

% (9) relative ripple size
nexttile; hold on; box on
rel = [W.relRMS];
for c = {'8mM','2mM','PNB'}
    s = strcmp(cond,c{1}) & isfinite(rel);
    if ~any(s); continue; end
    scatter(find(strcmp({'8mM','2mM','PNB'},c{1}))+0.12*randn(1,nnz(s)), 100*rel(s), ...
            26, col(c{1}), 'filled','MarkerFaceAlpha',.6);
end
set(gca,'XTick',1:3,'XTickLabel',{'8 mM','2 mM','PNB'});
ylabel('residual RMS (% of rise)'); title('how big is it');

title(tl, 'The recovery ripple: rig resonance, not kinetics', 'FontWeight','bold');
exportgraphics(fig, fullfile(resDir,'oscillation_id.png'), 'Resolution',140);
close(fig);
end
