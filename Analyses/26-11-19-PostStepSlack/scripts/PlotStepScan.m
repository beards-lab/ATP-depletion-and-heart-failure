%% PlotStepScan  Figure for #77: the slack mechanism, the step scan, the levers and their cost.
here = fileparts(mfilename('fullpath'));  ana = fileparts(here);  root = fileparts(fileparts(ana));
addpath(genpath(fullfile(root, 'Model')));  addpath(genpath(fullfile(root, 'Auxiliary')));  addpath(fullfile(root, 'params'));  addpath(here);
L = load(fullfile(ana, 'results', 'step_scan.mat'));  Sc = load(fullfile(ana, 'results', 'score.mat'));
pP = basePars();  n0 = pP.ekSE;  n = 3;
o0 = simStepRamp(pP, 0.054, 4);
o3 = simStepRamp(pP, 0.054, 4, struct('set', struct('ekSE', n, 'kSE', L.A.Fss^(1/n)/(L.A.LSE0*n/n0))));
if ~isfolder(fullfile(ana, 'figures')), mkdir(fullfile(ana, 'figures')); end
fig = figure('Visible', 'off', 'Position', [50 50 1500 780], 'Color', 'w');
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on;
f = @(o) o.Force(:)/interp1(o.t(:), o.Force(:), -0.001);
plot(1e3*o0.t, f(o0), 'k', 'LineWidth', 1.4, 'DisplayName', 'incumbent'); plot(1e3*o3.t, f(o3), 'Color', [0.85 0.33 0.1], 'LineWidth', 1.4, 'DisplayName', 'series element with a toe (n = 3)');
yline(0.09, 'b--', 'DisplayName', 'data: lowest force in the ramp'); xlim([-1 18]); ylim([-0.05 1.05]);
xlabel('time from the step (ms)'); ylabel('force / isometric'); title('Step 0.054 L_o, ramp 4 ML/s'); legend('Box', 'off');
nexttile; hold on;
plot(1e3*o0.t, 1e3*o0.LSE, 'k', 'LineWidth', 1.4); plot(1e3*o3.t, 1e3*o3.LSE, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.4); yline(0, ':');
xlim([-1 18]); xlabel('time from the step (ms)'); ylabel('series-element extension (nm of SL)'); title('Below zero the preparation is slack');
nexttile; hold on;
plot(1e3*o0.t, -o0.Vxb, 'k', 'LineWidth', 1.4); yline(4*2.2, 'b--'); xlim([-1 18]); ylim([0 30]);
xlabel('time from the step (ms)'); ylabel('contractile shortening (\mum/s)'); title('Taking up the slack (dashed: the ramp)');
nexttile; hold on;
cm = lines(numel(L.vv));
for i = 1:numel(L.vv), plot(L.steps, L.A.slack(i, :), 'o-', 'Color', cm(i, :), 'LineWidth', 1.3, 'DisplayName', sprintf('%g ML/s', L.vv(i))); end
xline(L.A.LSE0/2, 'k:', 'DisplayName', 'series extension at isometric force');
xlabel('step (L_o)'); ylabel('slack time (ms)'); title('Model: slack begins near 0.03 L_o at every speed'); legend('Box', 'off', 'Location', 'northwest'); grid on;
nexttile; hold on;
ds = [0.021 0.042 0.043 0.054 0.065 0.066];  dm = [0.45 0.18 0.16 0.09 0.04 0.03];
plot(ds, dm, 'bs-', 'MarkerFaceColor', 'b', 'DisplayName', 'data, 8 mM (lowest force after the step)');
plot(L.steps, L.A.Fmin(1, :), 'ko-', 'DisplayName', 'model, ramp 1 ML/s'); plot(L.steps, L.A.Fmin(3, :), 'k^--', 'DisplayName', 'model, ramp 4 ML/s');
xlabel('step (L_o)'); ylabel('lowest force / isometric'); title('The data keep force after every step'); legend('Box', 'off'); grid on;
nexttile;
S = Sc.S;  bar(categorical({S.name}, {S.name}), [S.cost]); ylabel('8 mM objective E_t'); title('What the toe costs on the slack battery'); grid on;
text(1:numel(S), [S.cost], compose('%.2f', [S.cost]), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
title(tl, 'Post-step slack on the isovelocity cycle (#77): robust1\_seed, 8 mM');
exportgraphics(fig, fullfile(ana, 'figures', 'post_step_slack.png'), 'Resolution', 150);  close(fig);
