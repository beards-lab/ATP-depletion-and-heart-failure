%% DeckFigures2.m -- v3 deck figures: F6b paced-only Ca span, F9 two cores / same twitch (states, DRX/SRX, ATP),
%% F10 series compliance (kSE) and sarcomere length: how much they matter for the twitch and for predictions
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');  out = fullfile(res, 'deck');
set(groot, 'defaultAxesFontSize', 13, 'defaultAxesFontName', 'Arial', 'defaultLineLineWidth', 2);
L = load(fullfile(here, '..', 'CalciumTwitch', 'results', 'eht_controls.mat'));  D = L.D;
sel = D.tb >= 0;  tb = D.tb(sel);  Fd = D.meanF(sel);  sdF = D.sdF(sel);  tfine = (0:0.0025:1)';
KB = [0 0 0];  RD = [0.80 0.15 0.10];  BL = [0 0.35 0.70];  GY = [0.55 0.55 0.55];  OR = [0.9 0.5 0.1];  GR = [0.1 0.6 0.3];  PU = [0.49 0.18 0.56];
dataBand = @(ax) fill(ax, 1e3*[tb; flipud(tb)], [Fd - sdF; flipud(Fd + sdF)], GY, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

%% F9 -- two cores, same twitch ----------------------------------------------------------------------------
A = load(fullfile(res, 't1_isoT_fit.mat'));  B = load(fullfile(res, 'tB__r_net__8_fit.mat'));
cases = {A.R, B.R};  nm = {'humanized core: cycle x0.45 of mouse', 'fast core: cycle x10 of mouse, SRX exchange x0.22'};  cl = [RD; PU];
fig = figure('Color', 'w', 'Position', [20 20 1600 900], 'Visible', 'off');  tl_ = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile; hold(ax, 'on');  hb = dataBand(ax);  hb.HandleVisibility = 'off';  plot(ax, 1e3*tb, Fd, '.', 'Color', KB, 'MarkerSize', 10, 'DisplayName', 'EHT data');
for k = 1:2, plot(ax, 1e3*cases{k}.tfine, cases{k}.sigma*cases{k}.fm, '-', 'Color', cl(k, :), 'DisplayName', sprintf('%s (rms %.2f %%)', nm{k}, 100*cases{k}.rms)); end
title(ax, 'A  twitch: both cores fit the data', 'FontWeight', 'normal');  ylabel(ax, 'active force (\muN)');  grid(ax, 'on');  legend(ax, 'Box', 'off', 'FontSize', 10, 'Location', 'northeast');
ax = nexttile; hold(ax, 'on');
for k = 1:2, plot(ax, 1e3*cases{k}.tfine, cases{k}.S.gate, '-', 'Color', cl(k, :)); end
plot(ax, 1e3*cases{1}.tfine, (cases{1}.S.Ca - min(cases{1}.S.Ca))/(max(cases{1}.S.Ca) - min(cases{1}.S.Ca)), ':', 'Color', BL);
title(ax, 'B  thin-filament gate (dotted: the fixed Ca, normalised)', 'FontWeight', 'normal');  ylabel(ax, 'fraction of sites available');  grid(ax, 'on');  ylim(ax, [0 1]);
ax = nexttile; hold(ax, 'on');
for k = 1:2, plot(ax, 1e3*cases{k}.tfine, cases{k}.S.bound, '-', 'Color', cl(k, :)); end
title(ax, 'C  attached heads (A1 + A2)', 'FontWeight', 'normal');  ylabel(ax, 'fraction of heads');  grid(ax, 'on');
ax = nexttile; hold(ax, 'on');
for k = 1:2, plot(ax, 1e3*cases{k}.tfine, cases{k}.S.PD + cases{k}.S.PT, '-', 'Color', cl(k, :)); end
title(ax, 'D  DRX detached heads (PT + PD, ready to attach)', 'FontWeight', 'normal');  ylabel(ax, 'fraction of heads');  xlabel(ax, 'ms from stimulus');  grid(ax, 'on');  ylim(ax, [0 1]);
ax = nexttile; hold(ax, 'on');
for k = 1:2, plot(ax, 1e3*cases{k}.tfine, cases{k}.S.SRX, '-', 'Color', cl(k, :)); end
title(ax, 'E  SRX heads (SR + SRD, parked)', 'FontWeight', 'normal');  ylabel(ax, 'fraction of heads');  xlabel(ax, 'ms from stimulus');  grid(ax, 'on');  ylim(ax, [0 1]);
ax = nexttile; hold(ax, 'on');
for k = 1:2, plot(ax, 1e3*cases{k}.tfine, cases{k}.S.R2T, '-', 'Color', cl(k, :), 'DisplayName', sprintf('%.1f ATP per head per beat', trapz(cases{k}.tfine, cases{k}.S.R2T))); end
set(ax, 'YScale', 'log');  title(ax, 'F  ATP turnover (per head per s, log scale)', 'FontWeight', 'normal');  xlabel(ax, 'ms from stimulus');  grid(ax, 'on');  legend(ax, 'Box', 'off', 'FontSize', 11, 'Location', 'northeast');
exportgraphics(fig, fullfile(out, 'F9_two_cores.png'), 'Resolution', 140);  close(fig);
for k = 1:2, fprintf('F9 %s: rms %.2f, sigma %.1f, Fsat_dev %.1f, bound pk %.2f, SRX %.2f->%.2f, ATP/beat %.1f\n', nm{k}, 100*cases{k}.rms, cases{k}.sigma, cases{k}.Fsat_dev, max(cases{k}.S.bound), cases{k}.S.SRX(1), min(cases{k}.S.SRX), trapz(cases{k}.tfine, cases{k}.S.R2T)); end


disp('F9 done')
