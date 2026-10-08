%% Stage1_CaTransient.m -- the fixed human-iPSC Ca transient, its metrics, and zero-tuning baselines
% 1) the prescribed hiPSC-EHT transient (literature_ca_transient.md, de Lange 2021 anchor) against
%    the adult-human form and the transient the EHT fit C had fitted freely (rice2008, 112/108 ms);
% 2) what the existing models do with the new transient and NO retuning:
%    (a) mouse incumbent (robust1_seed) + fit C thin filament;  (b) human_v0 ratios + fit C filament;
%    (c) sanity: mouse incumbent + fit C filament + fit C's own Ca (should be ~fit C, 0.8 % rms).
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');
L = load(fullfile(here, '..', 'CalciumTwitch', 'results', 'eht_controls.mat'));  D = L.D;
FC = load(fullfile(here, '..', 'CalciumTwitch', 'results', 'eht_fitC.mat'));  thC = FC.thC;
sel = D.tb >= 0;  tb = D.tb(sel);  Fd = D.meanF(sel);  tfine = (0:0.0025:1)';
cfg = tierConfigs('t0_v0');
tr = struct('name', {'hiPSC-EHT (prescribed, de Lange 2021 anchor)', 'adult human LV (literature form)', 'EHT fit C (freely fitted rice2008)'}, ...
    'kind', {'risedecay', 'risedecay', 'rice2008'}, ...
    'opt', {cfg.ca.opt, struct('Cadia', 0.10, 'Camax', 0.45, 'tau_r', 0.018, 'm', 1, 't50', 0.210, 'hdec', 2.8, 't0', 0.005, 'nsup', 3), ...
            struct('tau1', thC.tau1, 'tau2', thC.tau2 + 1e-3*(abs(thC.tau1 - thC.tau2) < 1e-3), 't0', 0, 'Cadia', thC.Cadia, 'Camax', thC.Camax)});
fprintf('%-46s %5s %5s %5s %6s %5s %5s %5s %5s %6s\n', 'transient', 'dia', 'peak', 'TTP', 't10-90', 'dec50', 'dec90', 'CTD50', 'CTD90', 'carry');
for k = 1:numel(tr)
    m = caTransientMetrics(tr(k).kind, tr(k).opt, 1);  tr(k).m = m;
    fprintf('%-46s %5.2f %5.2f %5.0f %6.0f %5.0f %5.0f %5.0f %5.0f %6.3f\n', tr(k).name, m.dia, m.peak, 1e3*m.ttp, 1e3*m.t10_90, 1e3*m.dec50, 1e3*m.dec90, 1e3*m.ctd50, 1e3*m.ctd90, m.carry);
end
%% zero-tuning baselines
thTF = struct('tf_k_on', thC.tf_k_on*(thC.Camax - thC.Cadia)/(cfg.ca.opt.Camax - cfg.ca.opt.Cadia), 'tf_k_off', thC.tf_k_off, 'tf_K_coop', thC.tf_K_coop, ...
    'CaKind', 'risedecay', 'CaOpt', cfg.ca.opt, 'pre', 2, 'nb', 3, 'shift', 0);   % k_on rescaled to the new Ca amplitude (k_on*dCa invariance)
baseM = cfg.base;
baseH = cfg.base;  for f = fieldnames(cfg.fixedScale)', baseH.(f{1}) = baseH.(f{1})*cfg.fixedScale.(f{1}); end
thCa = thC;  thCa.CaKind = 'rice2008';  thCa.CaOpt = tr(3).opt;  thCa.pre = 2;  thCa.nb = 3;  thCa.shift = thC.shift;
runs = {'mouse core + fit C filament + hiPSC Ca', baseM, thTF; 'human_v0 core + fit C filament + hiPSC Ca', baseH, thTF; 'sanity: mouse core + fit C filament + fit C Ca', baseM, thCa};
B = struct('name', {}, 'rms', {}, 'sigma', {}, 'm', {}, 'F', {}, 'S', {});
for k = 1:size(runs, 1)
    tic; [fm, S] = simHumanTwitch(runs{k, 2}, runs{k, 3}, tfine);  el = toc;
    fmd = interp1(tfine, fm, tb);  sg = max((fmd'*Fd)/(fmd'*fmd), 0);  r = (sg*fmd - Fd)/max(Fd);
    m = twitchShape(tfine, fm);  [~, iC] = max(S.Ca);  [~, iF] = max(fm);
    B(k) = struct('name', runs{k, 1}, 'rms', sqrt(mean(r.^2)), 'sigma', sg, 'm', m, 'F', fm, 'S', S);
    fprintf('%-48s rms %5.1f %%  sigma %6.2f  Fdev %5.1f kPa  TTP %3.0f 10->pk %3.0f RT50 %3.0f RT90 %3.0f ms  Ca->F lag %3.0f ms  gate pk %.2f  SRX %.2f->%.2f  drift %.1e  (%.0f s)\n', ...
        B(k).name, 100*B(k).rms, sg, max(fm), 1e3*[m.TTP m.C10toPk m.RT50 m.RT90], 1e3*(tfine(iF) - tfine(iC)), max(S.gate), S.SRX(1), min(S.SRX), S.drift, el);
end
mD = twitchShape(tb, Fd);
fprintf('%-48s %22s  TTP %3.0f 10->pk %3.0f RT50 %3.0f RT90 %3.0f ms\n', 'DATA', '', 1e3*[mD.TTP mD.C10toPk mD.RT50 mD.RT90]);
save(fullfile(res, 'stage1_ca.mat'), 'tr', 'B', 'tb', 'Fd', 'tfine');
%% figure
fig = figure('Color', 'w', 'Position', [30 30 1500 480], 'Visible', 'off');
tl_ = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
cols = [0 0.45 0.75; 0.85 0.33 0.1; 0.5 0.5 0.5];
ax = nexttile; hold(ax, 'on');
for k = 1:3, c = caTransient(tr(k).kind, tfine, setfield(tr(k).opt, 'T', 1)); plot(ax, 1e3*tfine, c, '-', 'Color', cols(k,:), 'LineWidth', 1.6, 'DisplayName', tr(k).name); end
title(ax, 'A  prescribed Ca transients (\muM)', 'FontWeight', 'normal');  legend(ax, 'Box', 'off', 'FontSize', 8);  grid(ax, 'on');  xlabel(ax, 'ms from stimulus');
ax = nexttile; hold(ax, 'on');
for k = 1:3, c = caTransient(tr(k).kind, tfine, setfield(tr(k).opt, 'T', 1)); plot(ax, 1e3*tfine, (c - min(c))/(max(c) - min(c)), '-', 'Color', cols(k,:), 'LineWidth', 1.6); end
plot(ax, 1e3*tb, Fd/max(Fd), 'k.', 'MarkerSize', 8, 'DisplayName', 'EHT force (norm.)');
title(ax, 'B  normalised shape; dots: the measured twitch', 'FontWeight', 'normal');  grid(ax, 'on');  xlabel(ax, 'ms from stimulus');
ax = nexttile; hold(ax, 'on');
plot(ax, 1e3*tb, Fd/max(Fd), 'k.', 'MarkerSize', 8, 'DisplayName', 'data');
for k = 1:numel(B), plot(ax, 1e3*tfine, B(k).F/max(B(k).F), '-', 'Color', cols(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%s (rms %.1f %%)', B(k).name, 100*B(k).rms)); end
title(ax, 'C  zero-tuning baselines (each normalised to its peak)', 'FontWeight', 'normal');  legend(ax, 'Box', 'off', 'FontSize', 7, 'Location', 'northeast');  grid(ax, 'on');  xlabel(ax, 'ms from stimulus');
exportgraphics(fig, fullfile(res, 'stage1_ca_transients.png'), 'Resolution', 130);  close(fig);
fprintf('stage 1 figure written\n');
