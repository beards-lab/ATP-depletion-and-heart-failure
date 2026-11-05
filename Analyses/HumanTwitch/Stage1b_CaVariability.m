%% Stage1b_CaVariability.m -- how variable is the "human iPSC Ca transient" in the literature, and what it costs
% Each source's reported shape metrics (literature_ca_transient.md) are reproduced by a 'risedecay' transient
% (tau_r, t50, hdec fitted to the metrics); the family is plotted, and the identified net core slowing of the
% tier-1 refits (t1_caFast / t1_isoT / t1_caSlow, and the caE indicator-error refits when present) is shown
% against the Ca time to peak. Output: results/stage1b_ca_variability.png, results/stage1b_ca_variability.mat
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');
% metrics in ms (from the stimulus/onset); fields used: ttp, t80rise, dec50, dec80, dec90, dt90_10, ctd90
S = struct('name', {}, 'cls', {}, 'met', {});
S(end+1) = struct('name', 'de Lange 2021, 3D hiPSC tissue, Fura-2, 37 C, 1 Hz', 'cls', 1, 'met', struct('ttp', 88, 'dec50', 185));
S(end+1) = struct('name', 'Saleem 2020, Hamburg EHT, GCaMP6f, 37 C, 1.5-2 Hz', 'cls', 1, 'met', struct('t80rise', 116, 'dec80', 173));
S(end+1) = struct('name', 'Paci 2018 exp., single hiPSC-CM, Fluo-4, 37 C, spontaneous', 'cls', 2, 'met', struct('ttp', 270, 'dt90_10', 410));
S(end+1) = struct('name', 'Paci 2020 model (hiPSC-CM, spontaneous)', 'cls', 3, 'met', struct('ttp', 184, 'dt90_10', 341));
S(end+1) = struct('name', 'Forouzandehmehr 2021 model (hiPSC-CM-CE)', 'cls', 3, 'met', struct('ttp', 163, 'dt90_10', 343));
S(end+1) = struct('name', 'Akwaboah 2021 model (hiPSC-CM, 2.3 Hz)', 'cls', 3, 'met', struct('ttp', 66, 'dt90_10', 127));
S(end+1) = struct('name', 'adult human LV (Coppini 2013 / ToR-ORd), 37 C, 1 Hz', 'cls', 4, 'met', struct('ttp', 60, 'ctd90', 460));
S(end+1) = struct('name', 'this campaign: prescribed hiPSC-EHT set', 'cls', 0, 'met', struct('ttp', 93, 'dec50', 186, 'dec90', 448));
tfine = (0:0.001:1)';
for k = 1:numel(S)
    f = fieldnames(S(k).met);  tgt = cellfun(@(q) S(k).met.(q), f);
    obj = @(z) sum((log(max(arrayfun(@(i) getm(z, f{i}), 1:numel(f)), 1e-3)) - log(tgt(:)')).^2) + 0.02*(log(z(3)/3))^2 + 10*max(0, log(abs(z(1))/0.15))^2;   % weak prior h ~ 3; tau_r <= 150 ms
    z = fminsearch(obj, [0.025 0.265 3], optimset('Display', 'off', 'TolX', 1e-4, 'MaxFunEvals', 600));
    S(k).opt = struct('Cadia', 0.15, 'Camax', 0.60, 'tau_r', abs(z(1)), 'm', 1, 't50', abs(z(2)), 'hdec', abs(z(3)), 't0', 0.005, 'nsup', 3);
    S(k).m = caTransientMetrics('risedecay', S(k).opt, 1);
    S(k).Ca = caTransient('risedecay', tfine, setfield(S(k).opt, 'T', 1));
    fprintf('%-62s TTP %4.0f  dec50 %4.0f  dec90 %4.0f ms  (tau_r %.0f t50 %.0f h %.1f)\n', S(k).name, 1e3*S(k).m.ttp, 1e3*S(k).m.dec50, 1e3*S(k).m.dec90, 1e3*S(k).opt.tau_r, 1e3*S(k).opt.t50, S(k).opt.hdec);
end
% identified net slowing vs the assumed Ca shape (tier-1 refits)
V = {'t1_caFast', 't1_isoT', 't1_caSlow', 'tB_caE1', 'tB_caE2', 'tB_caE3', 'tB_caE4', 'tB'};
pts = struct('name', {}, 'ttp', {}, 'dec50', {}, 'net', {}, 'rms', {}, 'cls', {});
for k = 1:numel(V)
    fn = fullfile(res, [V{k} '_fit.mat']);  if ~exist(fn, 'file'), continue; end
    F = load(fn);  m = caTransientMetrics('risedecay', F.cfg.ca.opt, 1);  nm = {F.cfg.par.name};
    if any(strcmp(nm, 'r_net')), net = 10^F.R.x(strcmp(nm, 'r_net')); else, net = 10^(F.R.x(strcmp(nm, 'r_iso')) + F.R.x(strcmp(nm, 'r_T'))); end
    pts(end+1) = struct('name', V{k}, 'ttp', 1e3*m.ttp, 'dec50', 1e3*m.dec50, 'net', net, 'rms', 100*F.R.rms, 'cls', startsWith(V{k}, 'tB')); %#ok<SAGROW>
end
save(fullfile(res, 'stage1b_ca_variability.mat'), 'S', 'pts');
%% figure
cols = [0 0 0; 0.85 0.33 0.1; 0 0.45 0.75; 0.47 0.67 0.19; 0.49 0.18 0.56];   % 0 ours, 1 paced 3D tissue 37C, 2 single cell, 3 models, 4 adult
fig = figure('Color', 'w', 'Position', [20 20 1500 520], 'Visible', 'off');
tl_ = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile([1 2]); hold(ax, 'on');
for k = 1:numel(S)
    c = S(k).Ca;  lw = 1.4 + 1.2*(S(k).cls == 0);  ls = '-';  if S(k).cls == 3, ls = '--'; end
    plot(ax, 1e3*tfine, (c - min(c))/(max(c) - min(c)), ls, 'Color', cols(S(k).cls + 1, :), 'LineWidth', lw, 'DisplayName', S(k).name);
end
D = load(fullfile(here, '..', 'CalciumTwitch', 'results', 'eht_controls.mat'));  D = D.D;  sel = D.tb >= 0;
plot(ax, 1e3*D.tb(sel), D.meanF(sel)/max(D.meanF(sel)), 'k.', 'MarkerSize', 7, 'DisplayName', 'measured EHT force (normalised)');
legend(ax, 'Box', 'off', 'FontSize', 8, 'Location', 'northeast');  grid(ax, 'on');  xlim(ax, [0 900]);
title(ax, 'A  "human iPSC-CM Ca transient" in the literature: each curve reproduces one source''s reported metrics (black: the set used here)', 'FontWeight', 'normal');
xlabel(ax, 'ms from stimulus');  ylabel(ax, 'normalised Ca');
ax = nexttile; hold(ax, 'on');
for k = 1:numel(pts)
    mk = 'o';  if pts(k).cls, mk = 's'; end
    plot(ax, pts(k).ttp, pts(k).net, mk, 'MarkerSize', 8, 'MarkerFaceColor', cols(1 + 2*pts(k).cls, :), 'Color', cols(1 + 2*pts(k).cls, :));
    text(ax, pts(k).ttp + 1.5, pts(k).net, sprintf('%s (%.2f %%)', strrep(pts(k).name, '_', '\_'), pts(k).rms), 'FontSize', 7);
end
yline(ax, 1, ':', 'mouse core', 'Color', [.5 .5 .5]);  yline(ax, 0.1, ':', 'literature isoform ratio', 'Color', [.5 .5 .5]);
set(ax, 'YScale', 'log');  grid(ax, 'on');  xlabel(ax, 'assumed Ca time to peak (ms)');  ylabel(ax, 'identified net multiplier on k2, k_2, kah, kamh, ka, kd');
title(ax, {'B  what the assumed shape does to the identified', 'human core (tier-1 refits; squares: \pm indicator error)'}, 'FontWeight', 'normal');
exportgraphics(fig, fullfile(res, 'stage1b_ca_variability.png'), 'Resolution', 130);  close(fig);
fprintf('figure written\n');

function v = getm(z, fld)
m = caTransientMetrics('risedecay', struct('Cadia', 0.15, 'Camax', 0.60, 'tau_r', abs(z(1)), 'm', 1, 't50', abs(z(2)), 'hdec', abs(z(3)), 't0', 0.005, 'nsup', 3), 1);
v = 1e3*m.(fld);
end
