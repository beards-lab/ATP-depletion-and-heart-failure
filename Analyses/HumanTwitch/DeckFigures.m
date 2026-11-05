%% DeckFigures.m -- clean, single-purpose figures for the storyline deck (results/deck/*.png)
% F2 failure: mouse core + rat (Rice 2008) Ca + literature thin filament vs the EHT twitch
% F3 fit C: Ca fitted, core fixed -> match; the fitted Ca against the rat one
% F4 mouse vs human: FV curves and ktr / ATP bars of the mouse core, the literature human core (rig T), the same at 37 C, the tier-1 core
% F5 tier 1: hiPSC Ca fixed, humanizing directions only -> match
% F6 literature Ca variability (Ca only, no force overlay)
% F7 consequence: identified core factor vs assumed Ca time to peak
% F8 case B: profile over the core speed and the spread of ktr / ATP per beat
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');  out = fullfile(res, 'deck');
if ~exist(out, 'dir'), mkdir(out); end
set(groot, 'defaultAxesFontSize', 14, 'defaultAxesFontName', 'Arial', 'defaultLineLineWidth', 2);
L = load(fullfile(here, '..', 'CalciumTwitch', 'results', 'eht_controls.mat'));  D = L.D;
sel = D.tb >= 0;  tb = D.tb(sel);  Fd = D.meanF(sel);  sdF = D.sdF(sel);  tfine = (0:0.0025:1)';
KB = [0 0 0];  RD = [0.80 0.15 0.10];  BL = [0 0.35 0.70];  GY = [0.55 0.55 0.55];  OR = [0.9 0.5 0.1];  GR = [0.1 0.6 0.3];
dataBand = @(ax) fill(ax, 1e3*[tb; flipud(tb)], [Fd - sdF; flipud(Fd + sdF)], GY, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
exampleBase = loadParams(fullfile(here, '..', '..', 'params', 'ExampleHighATP.m'));  exampleBase.SL0 = 2.0;
FC = load(fullfile(here, '..', 'CalciumTwitch', 'results', 'eht_fitC.mat'));  thC = FC.thC;
if isempty(gcp('nocreate')), parpool('local', 6); end

%% F2 -- total failure --------------------------------------------------------------------------
thRat = struct('CaKind', 'rice2008', 'CaOpt', struct(), 'tf_k_on', 101.185, 'tf_k_off', 723.852, 'tf_K_coop', 9.6846, 'pre', 2, 'nb', 2);
[f2, S2] = simHumanTwitch(exampleBase, thRat, tfine);
f2d = interp1(tfine, f2, tb);  sg2 = (f2d'*Fd)/(f2d'*f2d);  r2 = (sg2*f2d - Fd)/max(Fd);
m2 = twitchShape(tfine, f2);  mD = twitchShape(tb, Fd);
fig = figure('Color', 'w', 'Position', [20 20 1500 560], 'Visible', 'off');  tl_ = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile; hold(ax, 'on');  plot(ax, 1e3*tfine, S2.Ca, '-', 'Color', BL);
title(ax, 'Ca transient: rat ventricle, Rice 2008 (MyoSim-type input)', 'FontWeight', 'normal');  ylabel(ax, '[Ca^{2+}]_i (\muM)');  xlabel(ax, 'ms from stimulus');  grid(ax, 'on');  xlim(ax, [0 1000]);
text(ax, 0.97, 0.9, sprintf('peak at %.0f ms', 1e3*tfine(find(S2.Ca == max(S2.Ca), 1))), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 13);
ax = nexttile; hold(ax, 'on');  dataBand(ax);  plot(ax, 1e3*tb, Fd, '.', 'Color', KB, 'MarkerSize', 12);  plot(ax, 1e3*tfine, sg2*f2, '-', 'Color', RD);
title(ax, 'Twitch: EHT data (black) vs mouse model with that Ca (red, amplitude scaled)', 'FontWeight', 'normal');  ylabel(ax, 'active force (\muN)');  xlabel(ax, 'ms from stimulus');  grid(ax, 'on');  xlim(ax, [0 1000]);
text(ax, 0.97, 0.9, sprintf('model: peak %.0f ms, RT50 %.0f ms\ndata: peak %.0f ms, RT50 %.0f ms\nrms %.0f %% of peak', 1e3*m2.TTP, 1e3*m2.RT50, 1e3*mD.TTP, 1e3*mD.RT50, 100*sqrt(mean(r2.^2))), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 13);
exportgraphics(fig, fullfile(out, 'F2_failure.png'), 'Resolution', 150);  close(fig);
fprintf('F2: rms %.1f %%, TTP %.0f, RT50 %.0f, Fdev %.1f kPa (sat ~ %.0f kPa)\n', 100*sqrt(mean(r2.^2)), 1e3*m2.TTP, 1e3*m2.RT50, max(f2), 70);

%% F3 -- fit C: Ca fitted, core fixed ---------------------------------------------------------------
caC = struct('tau1', thC.tau1, 'tau2', thC.tau2 + 1e-3*(abs(thC.tau1 - thC.tau2) < 1e-3), 't0', 0, 'Cadia', thC.Cadia, 'Camax', thC.Camax);
thFitC = struct('CaKind', 'rice2008', 'CaOpt', caC, 'tf_k_on', thC.tf_k_on, 'tf_k_off', thC.tf_k_off, 'tf_K_coop', thC.tf_K_coop, 'pre', 2, 'nb', 2, 'shift', thC.shift);
[f3, S3] = simHumanTwitch(exampleBase, thFitC, tfine);
f3d = interp1(tfine, f3, tb);  sg3 = (f3d'*Fd)/(f3d'*f3d);  r3 = (sg3*f3d - Fd)/max(Fd);
fig = figure('Color', 'w', 'Position', [20 20 1500 560], 'Visible', 'off');  tl_ = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile; hold(ax, 'on');
plot(ax, 1e3*tfine, S2.Ca/max(S2.Ca), '-', 'Color', BL);  plot(ax, 1e3*tfine, (S3.Ca - min(S3.Ca))/(max(S3.Ca) - min(S3.Ca)), '-', 'Color', RD);
legend(ax, {'rat, Rice 2008 (peak 48 ms)', sprintf('fitted to the EHT twitch (peak %.0f ms)', 1e3*tfine(find(S3.Ca == max(S3.Ca), 1)))}, 'Box', 'off', 'Location', 'northeast');
title(ax, 'Ca transient (normalised): the one the twitch demands vs the rat one', 'FontWeight', 'normal');  ylabel(ax, 'Ca (normalised)');  xlabel(ax, 'ms from stimulus');  grid(ax, 'on');  xlim(ax, [0 990]);
ax = nexttile; hold(ax, 'on');  dataBand(ax);  plot(ax, 1e3*tb, Fd, '.', 'Color', KB, 'MarkerSize', 12);  plot(ax, 1e3*tfine, sg3*f3, '-', 'Color', RD);
title(ax, sprintf('Twitch: data vs model with the fitted Ca (rms %.1f %% of peak)', 100*sqrt(mean(r3.^2))), 'FontWeight', 'normal');  ylabel(ax, 'active force (\muN)');  xlabel(ax, 'ms from stimulus');  grid(ax, 'on');  xlim(ax, [0 1000]);
exportgraphics(fig, fullfile(out, 'F3_fitC.png'), 'Resolution', 150);  close(fig);
fprintf('F3: rms %.2f %%\n', 100*sqrt(mean(r3.^2)));

%% F4 -- mouse vs human: FV and ktr / ATP ------------------------------------------------------------
c0 = tierConfigs('t0_v0');  mouse = c0.base;  humanRig = mouse;
for f = fieldnames(c0.fixedScale)', humanRig.(f{1}) = humanRig.(f{1})*c0.fixedScale.(f{1}); end
human37 = humanRig;  for f = {'k2', 'k_2', 'kah', 'kamh', 'ka', 'kd', 'k1', 'k_1'}, human37.(f{1}) = human37.(f{1})*5; end
T1 = load(fullfile(res, 't1_isoT_fit.mat'));  tier1 = buildModel(T1.cfg, T1.R.x);
cores = {mouse, humanRig, human37, tier1};  names = {'mouse (fitted, rig 20-25 C)', 'human beta-MHC, literature ratios, rig T', 'human beta-MHC at 37 C (Q10 3, x5)', 'tier-1 core fitted to the EHT twitch'};
Rp = cell(1, 4);
parfor k = 1:4, Rp{k} = humanProbe(cores{k}, struct('Tiso', 4, 'Tktr', 2, 'vels', [0.1 0.25 0.5 0.75 1 1.5 2 3])); end
cc = [KB; BL; OR; RD];
fig = figure('Color', 'w', 'Position', [20 20 1500 560], 'Visible', 'off');  tl_ = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile([1 1]); hold(ax, 'on');
for k = 1:4, plot(ax, Rp{k}.FV_v, Rp{k}.FV_f, 'o-', 'Color', cc(k, :), 'MarkerFaceColor', cc(k, :), 'DisplayName', names{k}); end
title(ax, 'Force-velocity (active force, SL 2.2 \mum, 8 mM ATP)', 'FontWeight', 'normal');  xlabel(ax, 'shortening velocity (ML/s)');  ylabel(ax, 'force (kPa)');  grid(ax, 'on');  ylim(ax, [0 Inf]);
legend(ax, 'Box', 'off', 'Location', 'northeast', 'FontSize', 11);
q = {'F0', 'ktr', 'R2T'};  ttl = {'isometric force F0 (kPa)', 'k_{tr} (s^{-1})', 'ATP use per head (s^{-1})'};
for i = 2:3
    ax = nexttile; hold(ax, 'on');
    v = cellfun(@(r) r.(q{i}), Rp);  b = bar(ax, 1:4, v, 0.6, 'FaceColor', 'flat');  b.CData = cc;
    for k = 1:4, text(ax, k, v(k), sprintf(' %.3g', v(k)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 13); end
    set(ax, 'XTick', 1:4, 'XTickLabel', {'mouse', 'human, rig T', 'human, 37 C', 'tier-1 fit'});  title(ax, ttl{i}, 'FontWeight', 'normal');  grid(ax, 'on');  ylim(ax, [0 1.15*max(v)]);
end
exportgraphics(fig, fullfile(out, 'F4_mouse_vs_human.png'), 'Resolution', 150);  close(fig);
for k = 1:4, fprintf('F4 %-42s F0 %.0f ktr %.1f R2T %.1f Vmax %.2f a/F0 %.2f\n', names{k}, Rp{k}.F0, Rp{k}.ktr, Rp{k}.R2T, Rp{k}.Vmax_hill, Rp{k}.a_over_F0); end
save(fullfile(out, 'F4_probes.mat'), 'Rp', 'names');

%% F5 -- tier 1: hiPSC Ca fixed, humanizing directions only ------------------------------------------
R1 = T1.R;
fig = figure('Color', 'w', 'Position', [20 20 1500 560], 'Visible', 'off');  tl_ = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile; hold(ax, 'on');  plot(ax, 1e3*R1.tfine, R1.S.Ca, '-', 'Color', BL);
title(ax, 'Ca transient FIXED: human iPSC tissue, Fura-2, 37 C, 1 Hz (de Lange 2021)', 'FontWeight', 'normal');  ylabel(ax, '[Ca^{2+}]_i (\muM)');  xlabel(ax, 'ms from stimulus');  grid(ax, 'on');  xlim(ax, [0 1000]);
text(ax, 0.97, 0.9, 'peak 93 ms, peak\rightarrow50 % 186 ms, peak\rightarrow90 % 448 ms', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 13);
ax = nexttile; hold(ax, 'on');  dataBand(ax);  plot(ax, 1e3*tb, Fd, '.', 'Color', KB, 'MarkerSize', 12);  plot(ax, 1e3*R1.tfine, R1.sigma*R1.fm, '-', 'Color', RD);
title(ax, sprintf('Twitch: data vs humanized mouse model (rms %.2f %% of peak)', 100*R1.rms), 'FontWeight', 'normal');  ylabel(ax, 'active force (\muN)');  xlabel(ax, 'ms from stimulus');  grid(ax, 'on');  xlim(ax, [0 1000]);
text(ax, 0.97, 0.9, {'moved: ADP release, hydrolysis,', 'attachment, weak detachment', 'together \times0.45 (isoform \times temperature)', 'thin filament fitted (never calibrated)'}, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 13);
exportgraphics(fig, fullfile(out, 'F5_tier1.png'), 'Resolution', 150);  close(fig);

%% F6 -- literature Ca variability (Ca only) ----------------------------------------------------------
V = load(fullfile(res, 'stage1b_ca_variability.mat'));  Sv = V.S;
fig = figure('Color', 'w', 'Position', [20 20 1500 600], 'Visible', 'off');  ax = axes(fig); hold(ax, 'on');
tt = (0:0.001:1)';
cols = containers.Map({'0', '1', '2', '3', '4'}, {KB, RD, BL, GR, [0.5 0.2 0.6]});
for k = 1:numel(Sv)
    c = Sv(k).Ca;  lw = 2 + 2*(Sv(k).cls == 0);  ls = '-';  if Sv(k).cls == 3, ls = '--'; end
    plot(ax, 1e3*tt, (c - min(c))/(max(c) - min(c)), ls, 'Color', cols(num2str(Sv(k).cls)), 'LineWidth', lw, 'DisplayName', sprintf('%s  (peak %.0f ms)', Sv(k).name, 1e3*Sv(k).m.ttp));
end
legend(ax, 'Box', 'off', 'FontSize', 12, 'Location', 'northeast');  grid(ax, 'on');  xlim(ax, [0 900]);
title(ax, '"Human iPSC-CM Ca transient" in the literature: each curve reproduces one source''s metrics', 'FontWeight', 'normal');  xlabel(ax, 'ms from stimulus');  ylabel(ax, 'Ca (normalised)');
exportgraphics(fig, fullfile(out, 'F6_ca_variability.png'), 'Resolution', 150);  close(fig);

%% F7 -- consequence: identified core factor vs assumed Ca shape -------------------------------------
vars = {'t1_caFast', 't1_isoT', 't1_caSlow'};  lab = {'fast literature transient', 'de Lange 2021 (used)', 'slow literature transient'};
netf = zeros(1, 3);  ttpv = netf;  rmsv = netf;
for k = 1:3
    F = load(fullfile(res, [vars{k} '_fit.mat']));  m = caTransientMetrics('risedecay', F.cfg.ca.opt, 1);
    netf(k) = 10^(F.R.x(1) + F.R.x(2));  ttpv(k) = 1e3*m.ttp;  rmsv(k) = 100*F.R.rms;
end
fig = figure('Color', 'w', 'Position', [20 20 1500 560], 'Visible', 'off');  tl_ = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
cc3 = [GR; KB; OR];
ax = nexttile; hold(ax, 'on');
for k = 1:3, F = load(fullfile(res, [vars{k} '_fit.mat'])); plot(ax, 1e3*F.R.tfine, F.R.S.Ca, '-', 'Color', cc3(k, :), 'DisplayName', sprintf('%s: peak %.0f ms', lab{k}, ttpv(k))); end
legend(ax, 'Box', 'off', 'Location', 'northeast');  grid(ax, 'on');  xlim(ax, [0 1000]);
title(ax, 'Three assumed Ca transients, all within the literature range', 'FontWeight', 'normal');  xlabel(ax, 'ms from stimulus');  ylabel(ax, '[Ca^{2+}]_i (\muM)');
ax = nexttile; hold(ax, 'on');
b = bar(ax, 1:3, netf, 0.55, 'FaceColor', 'flat');  b.CData = cc3;  set(ax, 'YScale', 'log');  ylim(ax, [0.05 10]);
for k = 1:3, text(ax, k, netf(k)*1.25, sprintf('\\times%.2f\n(fit %.2f %%)', netf(k), rmsv(k)), 'HorizontalAlignment', 'center', 'FontSize', 13); end
yline(ax, 1, ':', 'mouse', 'Color', GY, 'LabelHorizontalAlignment', 'left');  yline(ax, 0.1, ':', 'literature human/mouse isoform ratio', 'Color', GY, 'LabelHorizontalAlignment', 'left');
set(ax, 'XTick', 1:3, 'XTickLabel', {'fast', 'de Lange', 'slow'});  grid(ax, 'on');
title(ax, 'Cross-bridge rates the twitch then asks for', 'FontWeight', 'normal');  ylabel(ax, 'rate multiplier vs the mouse');
exportgraphics(fig, fullfile(out, 'F7_consequence.png'), 'Resolution', 150);  close(fig);

%% F8 -- case B: profile and prediction spread --------------------------------------------------------
PB = load(fullfile(res, 'predict_caseB.mat'));  P = PB.P;  n = 96;
m = strcmp({P.prof}, 'r_net');  v = [P(m).profVal];  rr = 100*sqrt([P(m).ssr]/n);  kt = [P(m).ktr];  at = [P(m).ATPbeat];  [v, o] = sort(v);  rr = rr(o);  kt = kt(o);  at = at(o);
fig = figure('Color', 'w', 'Position', [20 20 1500 560], 'Visible', 'off');  tl_ = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile; hold(ax, 'on');  plot(ax, v, rr, 'o-', 'Color', KB, 'MarkerFaceColor', KB);  set(ax, 'XScale', 'log');  grid(ax, 'on');
yline(ax, 0.35, '--', 'data noise (SEM)', 'Color', GY);  xline(ax, 1, ':', 'mouse', 'Color', GY);
title(ax, 'Ca transient known: how well does the twitch fit vs the cross-bridge speed?', 'FontWeight', 'normal');  xlabel(ax, 'cross-bridge cycle speed (multiplier on the mouse), held fixed');  ylabel(ax, 'rms misfit (% of peak)');  ylim(ax, [0 3]);
ax = nexttile; hold(ax, 'on');
yyaxis(ax, 'left');  plot(ax, v, kt, 's-', 'Color', BL, 'MarkerFaceColor', BL);  ylabel(ax, 'predicted k_{tr} (s^{-1})');  ax.YColor = BL;
yyaxis(ax, 'right');  plot(ax, v, at, 'd-', 'Color', RD, 'MarkerFaceColor', RD);  ylabel(ax, 'predicted ATP per head per beat');  ax.YColor = RD;
set(ax, 'XScale', 'log');  grid(ax, 'on');  xlabel(ax, 'cross-bridge cycle speed (multiplier on the mouse)');
title(ax, 'Equally good fits (speed \geq 0.8) predict 4\times different k_{tr} and 16\times different ATP cost', 'FontWeight', 'normal');
exportgraphics(fig, fullfile(out, 'F8_caseB.png'), 'Resolution', 150);  close(fig);
fprintf('deck figures written to %s\n', out);
