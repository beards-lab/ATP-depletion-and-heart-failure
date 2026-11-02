% PlotSteadyStateSD.m -- steady (post-re-stretch isometric) force per slack
% segment: every preparation, the pooled mean, and BOTH spreads -- the raw
% between-preparation SD and the moderated SD the objective actually uses --
% against the model at 8 mM and at 2 mM (k2-only best, g = 0.54).
%
% Why: the single-prep 03/27 trace (Fss22 79.8 kPa) sits ~15 % above the pooled
% level (Fref 69.35 kPa), so trace plots make the model's plateau look low. The
% objective scores the pooled level; this plot shows what it is scored against.
rootS = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(rootS); addpath(genpath(rootS));
resS = fullfile(rootS,'Analyses','LowATP_Refit9K2Sweep','results');
evalc('T = buildPooledTarget(struct());');
k = find(strcmp(T.name,'steady')); Fref = T.scale.mean; nS = size(T.m8,2); sg = 1:nS;
Y8 = squeeze(T.Y(:,1,k,:))*Fref; Y2 = squeeze(T.Y(:,2,k,:))*Fref;     % prep x seg, kPa
m8 = T.m8(k,:)*Fref; m2 = T.m2(k,:)*Fref;
S = load(fullfile(resS,'k2sweep.mat'), 'RK2', 'fm8K2'); ig = find(abs([S.RK2.g]-0.54) < 1e-9, 1);
mod8 = double(S.fm8K2.steady); mod2 = double(S.RK2(ig).fm.steady);
navy = [0 39 76]/255; red = [0.75 0.2 0.12]; mk = {'o','s','^'};
f = figure('Color','w','Position',[50 50 1500 700],'Visible','off'); hold on; box on; grid on;
for c = 1:2
    if c == 1; m = m8; sr = T.s8raw(k); su = T.s8(k); col = navy; Y = Y8; dx = -0.12;
    else;      m = m2; sr = T.s2raw(k); su = T.s2(k); col = red;  Y = Y2; dx = 0.12; end
    patch([sg fliplr(sg)]+dx, [m.*(1+su) fliplr(m.*(1-su))], col, 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    patch([sg fliplr(sg)]+dx, [m.*(1+sr) fliplr(m.*(1-sr))], col, 'FaceAlpha', 0.30, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(sg+dx, m, '-', 'Color', col, 'LineWidth', 3, 'HandleVisibility', 'off');
    for p = 1:size(Y,1)
        plot(sg+dx, Y(p,:), mk{p}, 'Color', col, 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end
plot(sg-0.12, mod8(1:nS), 'd--', 'Color', navy, 'MarkerFaceColor', navy, 'MarkerSize', 11, 'LineWidth', 2);
plot(sg+0.12, mod2(1:nS), 'd--', 'Color', red, 'MarkerFaceColor', red, 'MarkerSize', 11, 'LineWidth', 2);
% legend proxies
h = gobjects(0);
h(end+1) = plot(nan, nan, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 3);
h(end+1) = patch(nan, nan, [0.3 0.3 0.3], 'FaceAlpha', 0.45, 'EdgeColor', 'none');
h(end+1) = patch(nan, nan, [0.3 0.3 0.3], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
for p = 1:3; h(end+1) = plot(nan, nan, mk{p}, 'Color', [0.3 0.3 0.3], 'MarkerSize', 8, 'LineWidth', 1.5); end %#ok<SAGROW>
h(end+1) = plot(nan, nan, 'd--', 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 10, 'LineWidth', 2);
legend(h, {'pooled mean', sprintf('\\pm raw between-prep SD (8 mM %.0f %%, 2 mM %.0f %%)', 100*T.s8raw(k), 100*T.s2raw(k)), ...
           sprintf('\\pm SD used in the objective (8 mM %.0f %%, 2 mM %.0f %%)', 100*T.s8(k), 100*T.s2(k)), ...
           ['prep ' T.scale.days{1}], ['prep ' T.scale.days{2}], ['prep ' T.scale.days{3}], 'model (8 mM refit9; 2 mM R_2\times0.54)'}, ...
       'Location', 'eastoutside', 'FontSize', 12);
text(0.6, 45, 'blue = 8 mM, red = 2 mM', 'FontSize', 13, 'Color', [0.3 0.3 0.3]);
set(gca, 'XTick', sg, 'FontSize', 14, 'XLim', [0.5 nS+0.5]); xlabel('slack–restretch cycle'); ylabel('steady force after re-stretch (kPa)');
title(sprintf('Steady force: each prep normalised to its own pre-slack isometric, rescaled to the pooled level F_{ref} = %.1f kPa', Fref), 'FontSize', 14);
exportgraphics(f, fullfile(resS,'steady_state_sd.png'), 'Resolution', 200); close(f);
fprintf('8 mM: data %s | model %s\n2 mM: data %s | model %s\n', mat2str(m8,3), mat2str(mod8(1:nS),3), mat2str(m2,3), mat2str(mod2(1:nS),3));
fprintf('model-data in raw SDs, 8 mM: %s\n', mat2str((mod8(1:nS)-m8)./(m8*T.s8raw(k)),2));
fprintf('model-data in raw SDs, 2 mM: %s\n', mat2str((mod2(1:nS)-m2)./(m2*T.s2raw(k)),2));
fprintf('03/27 Fss22 %.1f kPa vs Fref %.1f (x%.3f)\n', T.scale.perPrep(1), Fref, T.scale.perPrep(1)/Fref);
