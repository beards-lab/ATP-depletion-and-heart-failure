%% OverviewTiers.m -- table + figure over every finished tier (results/*_fit.mat)
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');
order = {'t0_v0', 't1_isoT', 't1_isoT_ex', 't1_caFast', 't1_caSlow', 't1_caHiLvl', 't2a_Pi', 't2b_SRX', 't2b_SRXc', 't2c_SE', 't2d_pas', 't2e_SL', 't2_all', 't2_allc', 's3_free', 's3_free_b', 's3_free_c'};
T = table();  Rs = struct([]);
for k = 1:numel(order)
    f = fullfile(res, [order{k} '_fit.mat']);  if ~exist(f, 'file'), continue; end
    F = load(f);  R = F.R;  S = R.S;  c = R.core;
    row = table(string(R.name), 100*R.rms, R.sigma, max(R.fm), 100*R.act, R.Fsat_dev, R.hydSat, R.ATP_R2T, S.SRX(1), min(S.SRX), max(S.gate), max(S.bound), ...
        1e3*R.m.TTP, 1e3*R.m.C10toPk, 1e3*R.m.RT50, 1e3*R.m.RT90, 1e3*R.lagCaF, 1e3*R.shift, ...
        c.ka, c.kd, c.k1, c.k2, c.kah, c.ksr0, c.kmsrd, c.kSE, R.tf.k_on, R.tf.k_off, R.tf.K_coop, ...
        'VariableNames', {'tier', 'rms_pct', 'sigma', 'Fdev_kPa', 'act_pct', 'Fsat_dev', 'ATPase_sat', 'ATP_per_beat', 'SRX_dia', 'SRX_min', 'gate_pk', 'bound_pk', ...
        'TTP', 'C10pk', 'RT50', 'RT90', 'lagCaF', 'shift_ms', 'ka', 'kd', 'k1', 'k2', 'kah', 'ksr0', 'kmsrd', 'kSE', 'k_on', 'k_off', 'K_coop'});
    T = [T; row]; %#ok<AGROW>
    if isempty(Rs), Rs = R; else, Rs(end+1) = R; end %#ok<SAGROW>
end
disp(T(:, 1:18));  disp(T(:, [1 19:end]));
writetable(T, fullfile(res, 'overview_tiers.csv'));
% figure: data + every fit (absolute, sigma applied), residuals, head pools of the headline tiers
D = load(fullfile(here, '..', 'CalciumTwitch', 'results', 'eht_controls.mat'));  D = D.D;  sel = D.tb >= 0;
fig = figure('Color', 'w', 'Position', [20 20 1500 900], 'Visible', 'off');
tl_ = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
cc = lines(numel(Rs));
ax = nexttile; hold(ax, 'on');  plot(ax, 1e3*D.tb(sel), D.meanF(sel), 'k.', 'MarkerSize', 9, 'DisplayName', 'data');
for k = 1:numel(Rs), plot(ax, 1e3*Rs(k).tfine, Rs(k).sigma*Rs(k).fm, '-', 'Color', cc(k,:), 'LineWidth', 1.2, 'DisplayName', sprintf('%s (%.2f %%)', strrep(Rs(k).name, '_', '\_'), 100*Rs(k).rms)); end
legend(ax, 'Box', 'off', 'FontSize', 7, 'Location', 'northeast');  title(ax, 'A  twitch: data and every tier', 'FontWeight', 'normal');  ylabel(ax, '\muN');  grid(ax, 'on');
ax = nexttile; hold(ax, 'on');
for k = 1:numel(Rs), plot(ax, 1e3*Rs(k).tb, 100*Rs(k).r, '-', 'Color', cc(k,:), 'LineWidth', 1.1); end
yline(ax, 0);  ylim(ax, [-4 4]);  title(ax, 'B  residuals (% of peak), same colours', 'FontWeight', 'normal');  grid(ax, 'on');
ax = nexttile; hold(ax, 'on');
for k = 1:numel(Rs), plot(ax, 1e3*Rs(k).tfine, Rs(k).S.SRX, '-', 'Color', cc(k,:), 'LineWidth', 1.1); end
title(ax, 'C  SRX fraction over the beat', 'FontWeight', 'normal');  grid(ax, 'on');
ax = nexttile; hold(ax, 'on');
for k = 1:numel(Rs), plot(ax, 1e3*Rs(k).tfine, Rs(k).S.R2T, '-', 'Color', cc(k,:), 'LineWidth', 1.1); end
set(ax, 'YScale', 'log');  title(ax, 'D  ATP use R2T (per head per s, log)', 'FontWeight', 'normal');  grid(ax, 'on');
for a = tl_.Children', if isa(a, 'matlab.graphics.axis.Axes'), xlabel(a, 'ms from stimulus'); end, end
exportgraphics(fig, fullfile(res, 'overview_tiers.png'), 'Resolution', 130);  close(fig);
fprintf('overview written\n');
