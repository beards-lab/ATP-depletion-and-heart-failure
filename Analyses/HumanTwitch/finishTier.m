function R = finishTier(cfg, x, tb, Fd, D, outf, logf)
%FINISHTIER  Diagnostics of a tier at parameter vector X: fine twitch, metrics, states, ATP,
%   saturating-Ca steady force (activation), parameter table; saves OUTF and a png beside it.
here = fileparts(mfilename('fullpath'));
tfine = (0:0.0025:1)';
[base, th] = buildModel(cfg, x);
[fm, S] = simHumanTwitch(base, th, tfine);
fmd = interp1(tfine, fm, tb);  sig = max((fmd'*Fd)/(fmd'*fmd), 0);
r = (sig*fmd - Fd)/max(Fd);
R.name = cfg.name;  R.x = x;  R.par = cfg.par;  R.ssr = sum(r.^2);  R.rms = sqrt(mean(r.^2));
R.sigma = sig;  R.fm = fm;  R.tfine = tfine;  R.S = S;  R.r = r;  R.tb = tb;  R.Fd = Fd;
R.m = twitchShape(tfine, fm);  R.mD = twitchShape(tb, Fd);
R.ATP = trapz(tfine, S.hyd);  R.ATP_R2T = trapz(tfine, S.R2T);
[~, iCa] = max(S.Ca);  [~, iF] = max(fm);  R.lagCaF = tfine(iF) - tfine(iCa);
% steady saturating force with the same core and thin filament (pCa 4.51), 6 s (slow cores)
P = base;  P.CaModel = th.CaModel;  P.CaInput = 'const';  P.pCa = 4.51;
if ~isfield(th, 'CaModel'), P.CaModel = 'twoState'; end
for f = ["tf_k_on" "tf_k_off" "tf_K_coop"], if isfield(th, f), P.(f) = th.(f); end, end
P = getParams(P, [], true, false);  P.Vums = 0;  P.v = 0;
[~, y] = ode15s(@dPUdT_CombinedTransitions, [0 3 6], P.PU0(:), odeset('AbsTol', 1e-4, 'RelTol', 1e-3, 'MaxStep', 2e-3), P);
[~, o, rr] = dPUdT_CombinedTransitions(6, y(end, :)', P);
R.Fsat = o(1);  R.Fsat_dev = o(1) - S.Fdia;  R.act = max(fm)/R.Fsat_dev;  R.hydSat = (rr(1) - rr(2)) + (rr(13) - rr(14));
% parameter table (value, multiplier vs seed where applicable)
tab = cell(numel(cfg.par), 4);
for i = 1:numel(cfg.par)
    p = cfg.par(i);  v = x(i);  if p.log, v = 10^v; end
    tab(i, :) = {p.name, p.type, v, strjoin(cellstr(string(p.fields)), ',')};
end
R.tab = tab;
% key resolved core values
R.core = struct();
for f = ["ka" "kd" "k1" "k_1" "k2" "k_2" "kah" "kamh" "ksr0" "kmsrd" "ksrd" "kSE" "k_pas" "kstiff1" "kstiff2" "A0" "SL0"]
    if isfield(S.P, f), R.core.(f) = S.P.(f); end
end
R.tf = struct('k_on', S.P.tf_k_on, 'k_off', S.P.tf_k_off, 'K_coop', S.P.tf_K_coop);
R.shift = 0;  if isfield(th, 'shift'), R.shift = th.shift; end
save(outf, 'R', 'cfg');
msg = sprintf(['RESULT %s: rms %.3f %% | sigma %.3f uN/kPa | Fdev %.2f kPa | TTP %.0f/%.0f  10->pk %.0f/%.0f  RT50 %.0f/%.0f  RT90 %.0f/%.0f  W50 %.0f/%.0f ms (model/data) | ' ...
    'Ca->F lag %.0f ms | act %.0f %% of sat (Fsat_dev %.1f kPa) | gate pk %.3f bound pk %.3f SRX %.3f->%.3f | ATP net %.3f, R2T %.3f /head/beat | drift %.2g | shift %+.1f ms\n'], ...
    cfg.name, 100*R.rms, sig, max(fm), 1e3*[R.m.TTP R.mD.TTP R.m.C10toPk R.mD.C10toPk R.m.RT50 R.mD.RT50 R.m.RT90 R.mD.RT90 R.m.W50 R.mD.W50], ...
    1e3*R.lagCaF, 100*R.act, R.Fsat_dev, max(S.gate), max(S.bound), S.SRX(1), min(S.SRX), R.ATP, R.ATP_R2T, S.drift, 1e3*R.shift);
for i = 1:size(tab, 1), msg = [msg sprintf('   %-10s %-6s %10.4g   [%s]\n', tab{i, 1}, tab{i, 2}, tab{i, 3}, tab{i, 4})]; end %#ok<AGROW>
fid = fopen(logf, 'a');  fprintf(fid, '%s', msg);  fclose(fid);  fprintf('%s', msg);
% figure
try
    fig = figure('Color', 'w', 'Position', [30 30 1400 800], 'Visible', 'off');
    tl_ = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl_, sprintf('%s: rms %.2f %%, sigma %.2f uN/kPa, activation %.0f %%', strrep(cfg.name, '_', '\_'), 100*R.rms, sig, 100*R.act), 'FontWeight', 'normal');
    ax = nexttile([1 2]); hold(ax, 'on');
    fill(ax, 1e3*[tb; flipud(tb)], [Fd - D.sdF(D.tb >= 0); flipud(Fd + D.sdF(D.tb >= 0))], [.6 .6 .6], 'FaceAlpha', .2, 'EdgeColor', 'none');
    plot(ax, 1e3*tb, Fd, 'k.', 'MarkerSize', 9);  plot(ax, 1e3*tfine, sig*fm, 'r-', 'LineWidth', 1.8);
    title(ax, 'A  data (mean \pm SD, 15 wells) and fit', 'FontWeight', 'normal');  ylabel(ax, 'active force (\muN)');  grid(ax, 'on');
    ax = nexttile; hold(ax, 'on');  plot(ax, 1e3*tb, 100*r, 'r-o', 'MarkerSize', 3);  yline(ax, 0);
    title(ax, 'B  residual (% of peak)', 'FontWeight', 'normal');  grid(ax, 'on');
    ax = nexttile; hold(ax, 'on');
    plot(ax, 1e3*tfine, (S.Ca - min(S.Ca))/(max(S.Ca) - min(S.Ca)), 'b-', 'LineWidth', 1.5);  plot(ax, 1e3*tfine, 1 - S.N, 'm-');  plot(ax, 1e3*tfine, S.gate, 'r--');  plot(ax, 1e3*tfine, fm/max(fm), 'k:', 'LineWidth', 1.5);
    legend(ax, {'Ca (norm.)', 'permissive 1-N', 'gate', 'force/peak'}, 'Box', 'off', 'FontSize', 8);  title(ax, 'C  Ca -> filament -> force', 'FontWeight', 'normal');  grid(ax, 'on');
    ax = nexttile; hold(ax, 'on');
    plot(ax, 1e3*tfine, S.PD, 'b-');  plot(ax, 1e3*tfine, S.bound, 'r-');  plot(ax, 1e3*tfine, S.SRX, '-', 'Color', [.9 .6 0]);  plot(ax, 1e3*tfine, S.PT, 'g-');
    legend(ax, {'PD', 'attached', 'SRX', 'PT'}, 'Box', 'off', 'FontSize', 8);  title(ax, 'D  head pools', 'FontWeight', 'normal');  grid(ax, 'on');
    ax = nexttile; hold(ax, 'on');  area(ax, 1e3*tfine, S.hyd, 'FaceColor', [.2 .7 .4], 'FaceAlpha', .3);  plot(ax, 1e3*tfine, S.R2T, 'k-');
    title(ax, sprintf('E  ATP: net %.2f, R2T %.2f /head/beat; sat %.2f /s', R.ATP, R.ATP_R2T, R.hydSat), 'FontWeight', 'normal');  grid(ax, 'on');
    for a = tl_.Children', if isa(a, 'matlab.graphics.axis.Axes'), xlabel(a, 'ms from stimulus'); end, end
    exportgraphics(fig, strrep(outf, '_fit.mat', '.png'), 'Resolution', 120);  close(fig);
catch e
    fprintf('figure failed: %s\n', e.message);
end
end
