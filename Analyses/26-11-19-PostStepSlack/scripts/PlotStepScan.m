%% PlotStepScan  Figure for #77: the measured force time course against the model on three
% cycles, the slack mechanism, the step scan, the levers and their cost.
%
% Top row: force through the 1, 4 and 6 ML/s cycles. Data = the twelve 8 mM records of the
% Beard2022 set (mean and full range; time from the step onset; force above the slack baseline,
% as a fraction of the force in the hold before the step). Model = the same cycle with the
% median measured step. Needs Analyses/FV_Beard2022/results/beard2022_fv_raw.mat (git-ignored).
here = fileparts(mfilename('fullpath'));  ana = fileparts(here);  root = fileparts(fileparts(ana));
addpath(genpath(fullfile(root, 'Model')));  addpath(genpath(fullfile(root, 'Auxiliary')));  addpath(fullfile(root, 'params'));  addpath(here);
addpath(fullfile(root, 'Analyses', 'FV_Beard2022'));
L = load(fullfile(ana, 'results', 'step_scan.mat'));  Sc = load(fullfile(ana, 'results', 'score.mat'));
pP = basePars();  n0 = pP.ekSE;  n = 3;
toe = struct('set', struct('ekSE', n, 'kSE', L.A.Fss^(1/n)/(L.A.LSE0*n/n0)));
D = loadBeard2022FV();  D = D([D.atp] == 8);
vShow = [1 4 6];  tEnd = [70 18 13];                    % ML/s; panel length (ms)
orange = [0.85 0.33 0.1];  blue = [0.2 0.4 0.8];
fNorm = @(o) o.Force(:)/interp1(o.t(:), o.Force(:), -0.001);
if ~isfolder(fullfile(ana, 'figures')), mkdir(fullfile(ana, 'figures')); end
fig = figure('Visible', 'off', 'Position', [50 50 1900 800], 'Color', 'w');
tl = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
Tc = table();
for j = 1:numel(vShow)
    [tg, Fd, stp, tRamp] = dataCycle(D, vShow(j), tEnd(j));
    o0 = simStepRamp(pP, median(stp), vShow(j));  o3 = simStepRamp(pP, median(stp), vShow(j), toe);
    if vShow(j) == 4, o4 = o0; o4t = o3; end
    nexttile; hold on;
    fill([tg, fliplr(tg)], [min(Fd, [], 1), fliplr(max(Fd, [], 1))], blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', sprintf('data, 8 mM: range of %d preparations', size(Fd, 1)));
    plot(tg, mean(Fd, 1), 'Color', blue, 'LineWidth', 1.6, 'DisplayName', 'data: mean');
    plot(1e3*o0.t, fNorm(o0), 'k', 'LineWidth', 1.4, 'DisplayName', 'model: incumbent');
    plot(1e3*o3.t, fNorm(o3), 'Color', orange, 'LineWidth', 1.4, 'DisplayName', 'model: series element with a toe (n = 3)');
    yline(0, ':', 'HandleVisibility', 'off');
    xlim([-1 min(tEnd(j), median(tRamp))]); ylim([-0.05 1.05]);
    xlabel('time from the step (ms)'); ylabel('force / isometric');
    title(sprintf('Ramp %g ML/s after a step of %.3f L_o', vShow(j), median(stp)));
    if j == 1, legend('Box', 'off', 'Location', 'northeast'); end
    % numbers for the report: force at fixed times into the ramp
    tq = [2 5 10];  tq = tq(tq < min(tEnd(j), median(tRamp)));
    for q = tq
        Tc = [Tc; table(vShow(j), median(stp), q, mean(interp1(tg, Fd', q)), min(interp1(tg, Fd', q)), max(interp1(tg, Fd', q)), ...
            interp1(1e3*o0.t, fNorm(o0), q), interp1(1e3*o3.t, fNorm(o3), q), ...
            'VariableNames', {'v', 'step', 't_ms', 'data_mean', 'data_min', 'data_max', 'incumbent', 'toe'})]; %#ok<AGROW>
    end
end
disp(Tc);  writetable(Tc, fullfile(ana, 'results', 'timecourse_table.csv'));
nexttile; hold on;
plot(1e3*o4.t, 1e3*o4.LSE, 'k', 'LineWidth', 1.4); plot(1e3*o4t.t, 1e3*o4t.LSE, 'Color', orange, 'LineWidth', 1.4); yline(0, ':');
xlim([-1 18]); xlabel('time from the step (ms)'); ylabel('series-element extension (nm of SL)'); title('4 ML/s cycle: below zero the model is slack');
nexttile; hold on;
plot(1e3*o4.t, -o4.Vxb, 'k', 'LineWidth', 1.4); yline(4*2.2, 'b--'); xlim([-1 18]); ylim([0 30]);
xlabel('time from the step (ms)'); ylabel('contractile shortening (\mum/s)'); title('Taking up the slack (dashed: the ramp)');
nexttile; hold on;
cm = lines(numel(L.vv));
for i = 1:numel(L.vv), plot(L.steps, L.A.slack(i, :), 'o-', 'Color', cm(i, :), 'LineWidth', 1.3, 'DisplayName', sprintf('%g ML/s', L.vv(i))); end
xline(L.A.LSE0/2, 'k:', 'DisplayName', 'series extension at isometric force');
xlabel('step (L_o)'); ylabel('slack time (ms)'); title('Model: slack begins near 0.03 L_o at every speed'); legend('Box', 'off', 'Location', 'northwest'); grid on;
nexttile; hold on;
ds = [0.021 0.042 0.043 0.054 0.065 0.066];  dm = [0.45 0.18 0.16 0.09 0.04 0.03];
plot(ds, dm, 's-', 'Color', blue, 'MarkerFaceColor', blue, 'DisplayName', 'data, 8 mM (lowest force after the step)');
plot(L.steps, L.A.Fmin(1, :), 'ko-', 'DisplayName', 'model, ramp 1 ML/s'); plot(L.steps, L.A.Fmin(3, :), 'k^--', 'DisplayName', 'model, ramp 4 ML/s');
xlabel('step (L_o)'); ylabel('lowest force / isometric'); title('The data keep force after every step'); legend('Box', 'off'); grid on;
nexttile;
S = Sc.S;  bar(categorical({S.name}, {S.name}), [S.cost]); ylabel('8 mM objective E_t'); title('What the toe costs on the slack battery'); grid on;
text(1:numel(S), [S.cost], compose('%.2f', [S.cost]), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
title(tl, 'Post-step slack on the isovelocity cycle (#77): robust1\_seed, 8 mM');
exportgraphics(fig, fullfile(ana, 'figures', 'post_step_slack.png'), 'Resolution', 150);  close(fig);

function [tg, Fd, stp, tRamp] = dataCycle(D, v, tEnd)
% One row of Fd per record: force above the slack baseline over the force in the hold before
% the step, on the grid tg (ms from the step onset). stp = step size (Lo), tRamp = ramp length (ms).
vNom = [6 1 3 5 0.5 4 2];  c = find(vNom == v);          % recorded order of the cycles
tg = -2:0.05:tEnd;  Fd = nan(numel(D), numel(tg));  stp = nan(numel(D), 1);  tRamp = stp;
for k = 1:numel(D)
    t = D(k).t;  Lk = D(k).L;  F = D(k).F;
    R = readIsoV(t, Lk, F);  t0 = R.tStart(c);
    Fpre = median(F(t >= t0 - 20 & t < t0 - 1));
    Fd(k, :) = (interp1(t - t0, F, tg) - R.F0(c))/(Fpre - R.F0(c));
    stp(k) = median(Lk(t >= t0 - 20 & t < t0 - 1)) - mean(Lk(t >= t0 + 1.25 & t <= t0 + 1.75));
    tRamp(k) = R.tSlack(c) - t0;
end
end
