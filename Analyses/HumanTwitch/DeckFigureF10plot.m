%% DeckFigures2.m -- v3 deck figures: F6b paced-only Ca span, F9 two cores / same twitch (states, DRX/SRX, ATP),
%% F10 re-plot from saved Q
Q = load(fullfile(out, 'F10_compliance.mat'));  Q = Q.Q;
ik = 1:4;  is = [5 1 6 7];     % kSE scan; SL scan 1.9, 2.0 (= the kSE x1 run), 2.1, 2.2
cK = [KB; BL; OR; RD];  cS = [GR; KB; BL; PU];
fig = figure('Color', 'w', 'Position', [20 20 1600 900], 'Visible', 'off');  tl_ = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile; hold(ax, 'on');  plot(ax, 1e3*tb, Fd/max(Fd), '.', 'Color', GY, 'MarkerSize', 9, 'DisplayName', 'EHT data');
for i = 1:4, r = ik(i); plot(ax, 1e3*tfine, Q(r).fm/max(Q(r).fm), '-', 'Color', cK(i, :), 'DisplayName', sprintf('%s: rms %.1f %%', Q(r).lab, Q(r).rms)); end
title(ax, 'A  series stiffness: twitch shape, no refit', 'FontWeight', 'normal');  legend(ax, 'Box', 'off', 'FontSize', 10);  grid(ax, 'on');  ylabel(ax, 'F / peak');
ax = nexttile; hold(ax, 'on');
for i = 1:4, r = ik(i); plot(ax, 1e3*tfine, -1e3*Q(r).LSE, '-', 'Color', cK(i, :)); end
title(ax, 'B  contractile-element shortening (= series stretch)', 'FontWeight', 'normal');  ylabel(ax, 'shortening (nm per half-sarcomere)');  grid(ax, 'on');
ax = nexttile; hold(ax, 'on');
plot(ax, [Q(ik).kSE], [Q(ik).dF_Ca15], 'o-', 'Color', BL, 'MarkerFaceColor', BL, 'DisplayName', 'peak force for Ca amplitude x1.5 (%)');
plot(ax, [Q(ik).kSE], [Q(ik).dF_mava], 's-', 'Color', RD, 'MarkerFaceColor', RD, 'DisplayName', 'peak force for attachment x0.5 (%)');
plot(ax, [Q(ik).kSE], 100*([Q(ik).Fdev]/Q(1).Fdev - 1), 'd-', 'Color', KB, 'MarkerFaceColor', KB, 'DisplayName', 'peak force itself vs kSE x1 (%)');
set(ax, 'XScale', 'log');  title(ax, 'C  predictions vs series stiffness', 'FontWeight', 'normal');  xlabel(ax, 'kSE multiplier');  ylabel(ax, '% change');  grid(ax, 'on');  legend(ax, 'Box', 'off', 'FontSize', 10, 'Location', 'southeast');
ax = nexttile; hold(ax, 'on');  plot(ax, 1e3*tb, Fd/max(Fd), '.', 'Color', GY, 'MarkerSize', 9, 'DisplayName', 'EHT data');
for i = 1:4, r = is(i); plot(ax, 1e3*tfine, Q(r).fm/max(Q(r).fm), '-', 'Color', cS(i, :), 'DisplayName', sprintf('%s: rms %.1f %%', Q(r).lab, Q(r).rms)); end
title(ax, 'D  resting sarcomere length: twitch shape, no refit', 'FontWeight', 'normal');  legend(ax, 'Box', 'off', 'FontSize', 10);  grid(ax, 'on');  xlabel(ax, 'ms from stimulus');  ylabel(ax, 'F / peak');
ax = nexttile; hold(ax, 'on');
b = bar(ax, 1:4, [Q(is).Fdev], 0.6, 'FaceColor', 'flat');  b.CData = cS;  set(ax, 'XTick', 1:4, 'XTickLabel', {'1.9', '2.0', '2.1', '2.2 \mum'});
title(ax, 'E  peak twitch force vs sarcomere length', 'FontWeight', 'normal');  ylabel(ax, 'kPa');  grid(ax, 'on');
ax = nexttile; hold(ax, 'on');
plot(ax, [Q(is).SL0], [Q(is).dF_Ca15], 'o-', 'Color', BL, 'MarkerFaceColor', BL, 'DisplayName', 'Ca amplitude x1.5');
plot(ax, [Q(is).SL0], [Q(is).dF_mava], 's-', 'Color', RD, 'MarkerFaceColor', RD, 'DisplayName', 'attachment x0.5');
title(ax, 'F  predictions vs sarcomere length', 'FontWeight', 'normal');  xlabel(ax, 'SL0 (\mum)');  ylabel(ax, '% change of peak force');  grid(ax, 'on');  legend(ax, 'Box', 'off', 'FontSize', 10, 'Location', 'best');
exportgraphics(fig, fullfile(out, 'F10_compliance_SL.png'), 'Resolution', 140);  close(fig);
fprintf('deck figures 2 written\n');
