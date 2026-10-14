%% CompareTwitchFits.m -- the EHT twitch fitted with T-World's Ca: tw1_isoT (adult beat as granted)
% vs tw2_caH (its Ca handling freed), against the de Lange-transient fit t1_isoT.
% Reads the tier results Analyses/HumanTwitch/results/<tier>_fit.mat written by finishTier
% (R.fm model force on R.tfine, R.sigma uN/kPa, R.S.Ca/N/gate/bound diagnostics, R.Fd/R.tb data).
% Writes results/twitch_fit_tw1_tw2.png and results/twitch_fit_tw1_tw2.csv.

here = fileparts(mfilename('fullpath'));
res  = fullfile(here, 'results');  if ~exist(res, 'dir'), mkdir(res); end
fitd = fullfile(here, '..', 'HumanTwitch', 'results');
tiers = {'t1_isoT', 'tw1_isoT', 'tw2_caH'};
lab   = {'de Lange transient (t1\_isoT)', 'T-World as granted (tw1\_isoT)', 'T-World Ca handling freed (tw2\_caH)'};
col   = [0.45 0.45 0.45; 0 0.45 0.74; 0.85 0.33 0.1];
F = struct();
for i = 1:numel(tiers), S = load(fullfile(fitd, [tiers{i} '_fit.mat']));  F.(tiers{i}) = S.R; end
R0 = F.(tiers{1});  tb = R0.tb;  Fd = R0.Fd;

%% metrics table
rows = {};
for i = 1:numel(tiers)
    R = F.(tiers{i});  m = R.m;  S = R.S;
    x = R.x;  v = @(nm) pick(R, nm);
    rows(end+1, :) = {tiers{i}, 100*R.rms, R.sigma, max(R.fm), 1e3*m.TTP, 1e3*m.RT50, 1e3*m.RT90, 1e3*m.W50, 1e3*R.lagCaF, 100*R.act, ...
        min(S.Ca), max(S.Ca), 1e3*R.tfine(find(S.Ca == max(S.Ca), 1)), max(1 - S.N), max(S.gate), ...
        v('r_iso'), v('r_T'), v('r_iso')*v('r_T'), R.tf.k_on, R.tf.k_off, R.tf.K_coop, 1e3*R.shift, v('Jup_Multiplier'), v('ICaL_Multiplier'), v('INaCa_Multiplier')}; %#ok<SAGROW>
end
Tb = cell2table(rows, 'VariableNames', {'tier', 'rms_pct', 'sigma_uN_per_kPa', 'Fdev_kPa', 'TTP_ms', 'RT50_ms', 'RT90_ms', 'W50_ms', 'lagCaF_ms', 'act_pct_of_sat', ...
    'Ca_dia_uM', 'Ca_peak_uM', 'Ca_TTP_ms', 'permissive_max', 'gate_max', 'r_iso', 'r_T', 'r_iso_x_r_T', 'tf_k_on', 'tf_k_off', 'tf_K_coop', 'shift_ms', 'Jup_mult', 'ICaL_mult', 'INaCa_mult'});
dm = R0.mD;
fprintf('data: Fpk %.1f uN, TTP %.0f, RT50 %.0f, RT90 %.0f, W50 %.0f ms\n', dm.Fpk, 1e3*[dm.TTP dm.RT50 dm.RT90 dm.W50]);
disp(Tb);
writetable(Tb, fullfile(res, 'twitch_fit_tw1_tw2.csv'));

%% figure
fh = figure('Color', 'w', 'Position', [40 40 1400 820], 'Visible', 'off');
tl = tiledlayout(fh, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'EHT twitch fitted with T-World''s Ca transient vs the de Lange transient', 'FontWeight', 'normal');
% A force vs data
ax = nexttile([1 2]);  hold(ax, 'on');
plot(ax, 1e3*tb, Fd, 'k.', 'MarkerSize', 9, 'DisplayName', 'data (15 control wells, mean)');
for i = 1:numel(tiers), R = F.(tiers{i});  plot(ax, 1e3*R.tfine, R.sigma*R.fm, '-', 'Color', col(i, :), 'LineWidth', 1.6, 'DisplayName', sprintf('%s, rms %.2f %%', lab{i}, 100*R.rms)); end
legend(ax, 'Location', 'northeast', 'Box', 'off');  grid(ax, 'on');
xlabel(ax, 'ms after stimulus');  ylabel(ax, 'active force (\muN)');  title(ax, 'A  fit (sigma solved in closed form)', 'FontWeight', 'normal');
% B residuals
ax = nexttile;  hold(ax, 'on');
for i = 1:numel(tiers), R = F.(tiers{i});  plot(ax, 1e3*tb, 100*R.r, '-', 'Color', col(i, :), 'LineWidth', 1.2); end
yline(ax, 0);  grid(ax, 'on');  xlabel(ax, 'ms after stimulus');  ylabel(ax, '% of data peak');  title(ax, 'B  residual', 'FontWeight', 'normal');
% C Ca transients (absolute)
ax = nexttile;  hold(ax, 'on');
for i = 1:numel(tiers), R = F.(tiers{i});  plot(ax, 1e3*R.tfine, R.S.Ca, '-', 'Color', col(i, :), 'LineWidth', 1.5); end
grid(ax, 'on');  xlabel(ax, 'ms after stimulus');  ylabel(ax, '[Ca^{2+}]_i (\muM)');  title(ax, 'C  the Ca transient each fit used', 'FontWeight', 'normal');
% D permissive fraction and gate
ax = nexttile;  hold(ax, 'on');
for i = 1:numel(tiers), R = F.(tiers{i});  plot(ax, 1e3*R.tfine, 1 - R.S.N, '-', 'Color', col(i, :), 'LineWidth', 1.5);  plot(ax, 1e3*R.tfine, R.S.gate, '--', 'Color', col(i, :)); end
grid(ax, 'on');  xlabel(ax, 'ms after stimulus');  ylabel(ax, 'fraction');  title(ax, 'D  thin filament: permissive 1-N (solid), gate (dashed)', 'FontWeight', 'normal');
% E attached heads and force shape
ax = nexttile;  hold(ax, 'on');
for i = 1:numel(tiers), R = F.(tiers{i});  plot(ax, 1e3*R.tfine, R.S.bound, '-', 'Color', col(i, :), 'LineWidth', 1.5); end
grid(ax, 'on');  xlabel(ax, 'ms after stimulus');  ylabel(ax, 'attached fraction');  title(ax, 'E  attached heads (sigma absorbs the level)', 'FontWeight', 'normal');
exportgraphics(fh, fullfile(res, 'twitch_fit_tw1_tw2.png'), 'Resolution', 110);
fprintf('written %s\n', fullfile(res, 'twitch_fit_tw1_tw2.png'));

function v = pick(R, nm)
j = find(strcmp({R.par.name}, nm), 1);
if isempty(j), v = NaN; return; end
v = R.x(j);  if R.par(j).log, v = 10^v; end
end
