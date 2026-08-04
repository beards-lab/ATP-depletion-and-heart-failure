% PlotWindows.m — overlay data vs model post-restretch recovery windows,
% with the first-order fits drawn on top, to decide which estimator to adopt
% as a fit target. Requires results/baseline.mat (run Baseline.m first).

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

S = load(fullfile(here, 'results', 'baseline.mat'));

fig = figure(701); clf;
set(fig, 'Position', [60 60 1500 760], 'Color', 'w');
tl = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

for row = 1:2
    if row == 1; ty = 'postRestretch'; else; ty = 'postSlack'; end
    sd = S.res.(ty).data;  sm = S.res.(ty).model;
    for c = 1:5
        nexttile((row-1)*5 + c); hold on; box on;
        d = sd(c); m = sm(c);

        % traces, both referenced to their own baseline B, in F_iso units
        plot(d.t*1e3, d.F - d.B, 'k-',  'LineWidth', 1.0);
        plot(m.t*1e3, m.F - m.B, 'b-',  'LineWidth', 1.4);

        % A-fixed first-order fits
        td = linspace(0, max(d.t), 400);
        tm = linspace(0, max(m.t), 400);
        plot(td*1e3, d.A*(1-exp(-d.kfixA*max(td-d.t0A,0))), 'r--', 'LineWidth', 1.4);
        plot(tm*1e3, m.A*(1-exp(-m.kfixA*max(tm-m.t0A,0))), 'c--', 'LineWidth', 1.4);

        xlim([0 250]);
        title(sprintf('%s c%d\nk: D %.0f / M %.0f (%.1f\\times)', ...
              ty, c, d.kfixA, m.kfixA, m.kfixA/d.kfixA), 'FontSize', 9);
        if c == 1
            ylabel('F - F_{baseline}  (F_{iso})');
            legend({'data','model','fit data','fit model'}, 'Location','southeast','FontSize',7);
        end
        if row == 2; xlabel('t since window start (ms)'); end
    end
end
title(tl, sprintf('Recovery windows, %s — data vs model', S.SNAP), 'Interpreter','none');

exportgraphics(fig, fullfile(here, 'results', 'baseline_windows.png'), 'Resolution', 140);
fprintf('wrote results/baseline_windows.png\n');

%% How much of the data's post-restretch "slowness" is the slow tail?
% Refit the data over progressively shorter spans after the vall2 anchor.
spans = [0.030 0.050 0.080 0.120 0.200 Inf];
fprintf('\nDATA postRestretch kfixA vs fit span (ms):\n%-8s', 'cycle');
fprintf('%9.0f', spans*1e3); fprintf('\n');
for c = 1:5
    d = S.res.postRestretch.data(c);
    fprintf('%-8d', c);
    for s = spans
        mk = d.t <= d.t(1) + s;
        [~, k] = refit(d.t(mk), d.F(mk), d.B, d.A);
        fprintf('%9.1f', k);
    end
    fprintf('\n');
end
fprintf('%-8s', 'MODEL');fprintf('\n');
for c = 1:5
    m = S.res.postRestretch.model(c);
    fprintf('%-8d', c);
    for s = spans
        mk = m.t <= m.t(1) + s;
        [~, k] = refit(m.t(mk), m.F(mk), m.B, m.A);
        fprintf('%9.1f', k);
    end
    fprintf('\n');
end

function [A, k] = refit(t, F, B, A)
    tt = t - t(1); y = F - B;
    obj = @(q) sum((y - A.*(1-exp(-abs(q(1)).*max(tt-max(q(2),0),0)))).^2);
    q = fminsearch(obj, [45, 0.004], optimset('Display','off'));
    k = abs(q(1));
end
