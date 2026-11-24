%% PlotTrisegHarness  Figures for the closed-loop harness (#97). Reads results/harness_<TAG>.mat.
here = fileparts(mfilename('fullpath'));  root = fileparts(fileparts(here));  addpath(fullfile(root, 'Organ'));
if ~exist('TAG', 'var'), TAG = 'kon3'; end
L = load(fullfile(here, 'results', ['harness_' TAG '.mat']));  R = L.R;  B = L.B;  H = L.H;  r = H.ref;
fd = fullfile(here, 'figures');  if ~exist(fd, 'dir'), mkdir(fd); end
T = H.T;  c = R.col;  k = R.t >= 9*T - 1e-9;  tk = 1e3*(R.t(k) - 9*T);  th = 1e3*r.t;
set(0, 'DefaultAxesFontSize', 10);  cw = lines(3);  wn = {'LV', 'SEP', 'RV'};

f = figure('Visible', 'off', 'Position', [50 50 1500 780]);  tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on;                                   % PV loops
plot(1e3*R.V_LV(k), R.out(k, c.P_LV), 'Color', cw(1, :), 'LineWidth', 1.6);  plot(1e3*R.V_RV(k), R.out(k, c.P_RV), 'Color', cw(3, :), 'LineWidth', 1.6);
plot(1e3*r.V(:, 1), r.P(:, 1), '--', 'Color', cw(1, :));  plot(1e3*r.V(:, 2), r.P(:, 4), '--', 'Color', cw(3, :));
xlabel('volume (\muL)'); ylabel('pressure (kPa)'); title('PV loops, last beat'); grid on;
legend('LV, cross-bridge walls', 'RV, cross-bridge walls', 'LV, host fit', 'RV, host fit', 'Location', 'best', 'FontSize', 8);
nexttile; hold on;                                   % pressures in time
plot(tk, R.out(k, c.P_LV), 'Color', cw(1, :), 'LineWidth', 1.4);  plot(tk, R.out(k, c.P_SA), 'k', 'LineWidth', 1.1);  plot(tk, R.out(k, c.P_RV), 'Color', cw(3, :), 'LineWidth', 1.4);
plot(th, r.P(:, 1), '--', 'Color', cw(1, :));  plot(th, r.P(:, 2), 'k--');  plot(th, r.P(:, 4), '--', 'Color', cw(3, :));
xlabel('t (ms)'); ylabel('pressure (kPa)'); title('Pressures (dashed: host fit)'); legend('P_{LV}', 'P_{SA}', 'P_{RV}', 'Location', 'best', 'FontSize', 8); grid on;
nexttile; hold on;                                   % sarcomere length
for w = 1:3, plot(tk, R.out(k, c.Ls(w)), 'Color', cw(w, :), 'LineWidth', 1.4, 'DisplayName', wn{w}); plot(th, r.Ls(:, w), '--', 'Color', cw(w, :), 'HandleVisibility', 'off'); end
xlabel('t (ms)'); ylabel('Ls (\mum)'); title('Sarcomere length (dashed: host fit)'); legend('Location', 'best', 'FontSize', 8); grid on;
nexttile; hold on;                                   % wall stress
for w = 1:3, plot(tk, R.out(k, c.sigma(w)), 'Color', cw(w, :), 'LineWidth', 1.4, 'DisplayName', wn{w}); plot(th, r.sigma(:, w), '--', 'Color', cw(w, :), 'HandleVisibility', 'off'); end
xlabel('t (ms)'); ylabel('wall stress (kPa)'); title('Wall stress (dashed: host fit)'); legend('Location', 'best', 'FontSize', 8); grid on;
nexttile; hold on;                                   % activation, LV
yyaxis left;  plot(tk, R.out(k, c.gate(1)), 'LineWidth', 1.3);  plot(tk, R.out(k, c.bound(1)), '--', 'LineWidth', 1.1);  ylabel('Ca gate, bound fraction (LV)');
yyaxis right; plot(tk, R.out(k, c.Ca), 'LineWidth', 1.1);  ylabel('Ca (\muM)');
xlabel('t (ms)'); title('Activation, LV'); legend('gate', 'bound', 'Ca', 'Location', 'best', 'FontSize', 8); grid on;
nexttile; hold on;                                   % approach to the periodic state
yyaxis left;  plot(B.beat, B.Pmax_LV, 'o-');  ylabel('LV peak pressure (kPa)');
yyaxis right; plot(B.beat, 1e3*B.EDV_LV, 's-');  plot(B.beat, 1e3*B.ESV_LV, 'v-');  ylabel('LV volume (\muL)');
xlabel('beat'); title('Approach to the periodic state'); legend('P_{max}', 'EDV', 'ESV', 'Location', 'best', 'FontSize', 8); grid on;
exportgraphics(f, fullfile(fd, ['harness_' TAG '.png']), 'Resolution', 130);  close(f);
disp('figure written');
