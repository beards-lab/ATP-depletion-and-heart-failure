%% DeckFigures2.m -- v3 deck figures: F6b paced-only Ca span, F9 two cores / same twitch (states, DRX/SRX, ATP),
%% F10 series compliance (kSE) and sarcomere length: how much they matter for the twitch and for predictions
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');  out = fullfile(res, 'deck');
set(groot, 'defaultAxesFontSize', 13, 'defaultAxesFontName', 'Arial', 'defaultLineLineWidth', 2);
L = load(fullfile(here, '..', 'CalciumTwitch', 'results', 'eht_controls.mat'));  D = L.D;
sel = D.tb >= 0;  tb = D.tb(sel);  Fd = D.meanF(sel);  sdF = D.sdF(sel);  tfine = (0:0.0025:1)';
KB = [0 0 0];  RD = [0.80 0.15 0.10];  BL = [0 0.35 0.70];  GY = [0.55 0.55 0.55];  OR = [0.9 0.5 0.1];  GR = [0.1 0.6 0.3];  PU = [0.49 0.18 0.56];
dataBand = @(ax) fill(ax, 1e3*[tb; flipud(tb)], [Fd - sdF; flipud(Fd + sdF)], GY, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
if isempty(gcp('nocreate')), parpool('local', 6); end

%% F6b -- paced Ca transients only --------------------------------------------------------------------
V = load(fullfile(res, 'stage1b_ca_variability.mat'));  Sv = V.S;
keep = find(contains({Sv.name}, 'de Lange') | contains({Sv.name}, 'Saleem') | contains({Sv.name}, 'adult human'));
tt = (0:0.001:1)';  cc = [RD; OR; PU];
fig = figure('Color', 'w', 'Position', [20 20 1500 560], 'Visible', 'off');  ax = axes(fig); hold(ax, 'on');
for i = 1:numel(keep)
    k = keep(i);  c = Sv(k).Ca;
    plot(ax, 1e3*tt, (c - min(c))/(max(c) - min(c)), '-', 'Color', cc(i, :), 'DisplayName', sprintf('%s   (peak %.0f ms, peak\\rightarrow50 %% %.0f ms)', Sv(k).name, 1e3*Sv(k).m.ttp, 1e3*Sv(k).m.dec50));
end
legend(ax, 'Box', 'off', 'FontSize', 13, 'Location', 'northeast');  grid(ax, 'on');  xlim(ax, [0 900]);
title(ax, 'Paced Ca transients at 37 C in the literature (normalised; each curve reproduces the source''s reported metrics)', 'FontWeight', 'normal');
xlabel(ax, 'ms from stimulus');  ylabel(ax, 'Ca (normalised)');
exportgraphics(fig, fullfile(out, 'F6b_ca_paced.png'), 'Resolution', 150);  close(fig);

%% F9 -- two cores, same twitch ----------------------------------------------------------------------------
A = load(fullfile(res, 't1_isoT_fit.mat'));  B = load(fullfile(res, 'tB__r_net__8_fit.mat'));
cases = {A.R, B.R};  nm = {'humanized core: cycle x0.45 of mouse', 'fast core: cycle x10 of mouse, SRX exchange x0.22'};  cl = [RD; PU];
fig = figure('Color', 'w', 'Position', [20 20 1600 900], 'Visible', 'off');  tl_ = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile; hold(ax, 'on');  dataBand(ax);  plot(ax, 1e3*tb, Fd, '.', 'Color', KB, 'MarkerSize', 10);
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

%% F10 -- series compliance and sarcomere length -----------------------------------------------------------
[base1, th1] = buildModel(A.cfg, A.R.x);
kf = [1 0.5 0.25 0.1];  sl = [1.9 2.0 2.1 2.2];
runs = [arrayfun(@(f) struct('lab', sprintf('kSE x%.2g', f), 'kSE', f, 'SL0', 2.0), kf), arrayfun(@(z) struct('lab', sprintf('SL %.1f um', z), 'kSE', 1, 'SL0', z), sl([1 3 4]))];
Q = cell(numel(runs), 1);
parfor r = 1:numel(runs)
    b = base1;  b.kSE = base1.kSE*runs(r).kSE;  b.SL0 = runs(r).SL0;
    [fm, S] = simHumanTwitch(b, th1, tfine);
    fmd = interp1(tfine, fm, tb);  sg = (fmd'*Fd)/(fmd'*fmd);  q = struct('lab', runs(r).lab, 'kSE', runs(r).kSE, 'SL0', runs(r).SL0);
    q.rms = 100*sqrt(mean(((sg*fmd - Fd)/max(Fd)).^2));  q.Fdev = max(fm);  q.sigma = sg;
    q.dSL = min(S.SL) - S.SL(1);  q.SLdia = S.SL(1);  q.fm = fm;  q.SL = S.SL;  q.ATP = trapz(tfine, S.R2T);  q.LSEmax = max(S.LSE);
    o = th1.CaOpt;  thc = th1;  thc.CaOpt.Camax = o.Cadia + 1.5*(o.Camax - o.Cadia);  q.dF_Ca15 = 100*(max(simHumanTwitch(b, thc, tb))/q.Fdev - 1);
    bm = b;  bm.ka = 0.5*bm.ka;  q.dF_mava = 100*(max(simHumanTwitch(bm, th1, tb))/q.Fdev - 1);
    m = twitchShape(tfine, fm);  q.TTP = 1e3*m.TTP;  q.RT50 = 1e3*m.RT50;
    Q{r} = q;
end
Q = [Q{:}];
for r = 1:numel(Q), fprintf('F10 %-12s rms %5.2f %%  Fdev %5.1f kPa  sigma %5.1f  dSL %+.3f um (%.1f %%)  LSE max %.3f  TTP %.0f RT50 %.0f  dF Ca1.5 %+.1f %%  dF mava %+.1f %%  ATP %.2f\n', ...
    Q(r).lab, Q(r).rms, Q(r).Fdev, Q(r).sigma, Q(r).dSL, 100*Q(r).dSL/Q(r).SLdia, Q(r).LSEmax, Q(r).TTP, Q(r).RT50, Q(r).dF_Ca15, Q(r).dF_mava, Q(r).ATP); end
save(fullfile(out, 'F10_compliance.mat'), 'Q');
ik = 1:4;  is = [2 5 6 7];     % kSE scan; SL scan (SL 2.0 is the kSE x1 run)
cK = [KB; BL; OR; RD];  cS = [GR; KB; BL; PU];
fig = figure('Color', 'w', 'Position', [20 20 1600 900], 'Visible', 'off');  tl_ = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile; hold(ax, 'on');  plot(ax, 1e3*tb, Fd/max(Fd), '.', 'Color', GY, 'MarkerSize', 9);
for i = 1:4, r = ik(i); plot(ax, 1e3*tfine, Q(r).fm/max(Q(r).fm), '-', 'Color', cK(i, :), 'DisplayName', sprintf('%s: rms %.1f %%', Q(r).lab, Q(r).rms)); end
title(ax, 'A  series stiffness: twitch shape (normalised), no refit', 'FontWeight', 'normal');  legend(ax, 'Box', 'off', 'FontSize', 10);  grid(ax, 'on');  ylabel(ax, 'F / peak');
ax = nexttile; hold(ax, 'on');
for i = 1:4, r = ik(i); plot(ax, 1e3*tfine, 1e3*(Q(r).SL - Q(r).SLdia), '-', 'Color', cK(i, :)); end
title(ax, 'B  sarcomere shortening during the twitch (model, fixed ends)', 'FontWeight', 'normal');  ylabel(ax, '\DeltaSL (nm)');  grid(ax, 'on');
ax = nexttile; hold(ax, 'on');
plot(ax, [Q(ik).kSE], [Q(ik).dF_Ca15], 'o-', 'Color', BL, 'MarkerFaceColor', BL, 'DisplayName', 'peak force for Ca amplitude x1.5 (%)');
plot(ax, [Q(ik).kSE], [Q(ik).dF_mava], 's-', 'Color', RD, 'MarkerFaceColor', RD, 'DisplayName', 'peak force for attachment x0.5 (%)');
plot(ax, [Q(ik).kSE], 100*([Q(ik).Fdev]/Q(1).Fdev - 1), 'd-', 'Color', KB, 'MarkerFaceColor', KB, 'DisplayName', 'peak force itself vs kSE x1 (%)');
set(ax, 'XScale', 'log');  title(ax, 'C  predictions vs series stiffness', 'FontWeight', 'normal');  xlabel(ax, 'kSE multiplier');  ylabel(ax, '% change');  grid(ax, 'on');  legend(ax, 'Box', 'off', 'FontSize', 10, 'Location', 'best');
ax = nexttile; hold(ax, 'on');  plot(ax, 1e3*tb, Fd/max(Fd), '.', 'Color', GY, 'MarkerSize', 9);
for i = 1:4, r = is(i); plot(ax, 1e3*tfine, Q(r).fm/max(Q(r).fm), '-', 'Color', cS(i, :), 'DisplayName', sprintf('%s: rms %.1f %%', Q(r).lab, Q(r).rms)); end
title(ax, 'D  sarcomere length: twitch shape (normalised), no refit', 'FontWeight', 'normal');  legend(ax, 'Box', 'off', 'FontSize', 10);  grid(ax, 'on');  xlabel(ax, 'ms from stimulus');  ylabel(ax, 'F / peak');
ax = nexttile; hold(ax, 'on');
b = bar(ax, 1:4, [Q(is).Fdev], 0.6, 'FaceColor', 'flat');  b.CData = cS;  set(ax, 'XTick', 1:4, 'XTickLabel', {Q(is).lab});
title(ax, 'E  peak twitch force (model kPa) vs sarcomere length', 'FontWeight', 'normal');  ylabel(ax, 'kPa');  grid(ax, 'on');
ax = nexttile; hold(ax, 'on');
plot(ax, [Q(is).SL0], [Q(is).dF_Ca15], 'o-', 'Color', BL, 'MarkerFaceColor', BL, 'DisplayName', 'Ca amplitude x1.5');
plot(ax, [Q(is).SL0], [Q(is).dF_mava], 's-', 'Color', RD, 'MarkerFaceColor', RD, 'DisplayName', 'attachment x0.5');
title(ax, 'F  predictions vs sarcomere length', 'FontWeight', 'normal');  xlabel(ax, 'SL0 (\mum)');  ylabel(ax, '% change of peak force');  grid(ax, 'on');  legend(ax, 'Box', 'off', 'FontSize', 10, 'Location', 'best');
exportgraphics(fig, fullfile(out, 'F10_compliance_SL.png'), 'Resolution', 140);  close(fig);
fprintf('deck figures 2 written\n');
