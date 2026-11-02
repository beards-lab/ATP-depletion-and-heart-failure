% PlotSrxCompensation.m -- figures for RunSrxCompensation (results/compensation.mat).
%   comp_residuals.png  every feature vs the data (data-SD units) at pool x0.5 and x2:
%                       incumbent, uncompensated, compensated with set A and set B
%   comp_chi2.png       total chi2 of each case (the identifiability verdict)
%   comp_knobs.png      the knob multipliers each compensation chose
%   comp_traces.png     slack cycle 1 and ktr protocol: force and SRX pool
%   compensation.txt    the numbers behind the figures
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsC = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
C = load(fullfile(rsC, 'compensation.mat'));
[~, ~, labs] = srxResiduals(C.E0.fm, C.E0.fd); nz = numel(labs);
ink = [0.13 0.13 0.13]; grid_ = [0.88 0.88 0.88]; band = [0.93 0.93 0.93];
cA = [42 120 214]/255; cB = [27 175 122]/255; cU = [235 104 52]/255;      % slots 1, 3, 2
nT = numel(C.TGT);
get = @(t, s) C.R(find([C.R.target] == t & strcmp({C.R.set}, s), 1));

%% ---- residuals -----------------------------------------------------------------------
fg = figure('Color', 'w', 'Position', [20 20 1500 760], 'Visible', 'off');
tl = tiledlayout(fg, 1, nT, 'TileSpacing', 'compact', 'Padding', 'compact');
y = nz:-1:1;
for it = 1:nT
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
    patch(ax, [-1 1 1 -1], [0.4 0.4 nz + 0.6 nz + 0.6], band, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    for k = 1:nz; plot(ax, [-8 8], y(k)*[1 1], '-', 'Color', grid_, 'HandleVisibility', 'off'); end
    plot(ax, [0 0], [0.4 nz + 0.6], '-', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    t = C.TGT(it); a = get(t, 'A'); b = get(t, 'B'); u = C.U{it};
    S = {C.E0, 'incumbent', ink, 'o'; u, 'pool pinned, uncompensated', cU, 's'; ...
         a.E, 'compensated, set A (cross-bridge knobs)', cA, 'd'; b.E, 'compensated, set B (A + SRX kinetics)', cB, '^'};
    for s = 1:size(S, 1)
        z = max(-8, min(8, S{s, 1}.z));
        plot(ax, z, y + (s - 2.5) * 0.14, S{s, 4}, 'Color', S{s, 3}, 'MarkerFaceColor', S{s, 3}, 'MarkerSize', 7, ...
             'LineStyle', 'none', 'DisplayName', sprintf('%s  (\\chi^2 %.1f)', S{s, 2}, S{s, 1}.chi2));
    end
    set(ax, 'YTick', 1:nz, 'YTickLabel', flip(labs), 'XColor', ink, 'YColor', ink, 'FontSize', 11, ...
        'TickLength', [0 0]); xlim(ax, [-8 8]); ylim(ax, [0.4 nz + 0.6]);
    xlabel(ax, 'model - data, in data-SD units (clipped at \pm8)');
    if t == 1; tt = sprintf('CONTROL: incumbent pool (%.3f), refitted the same way', u.srx);
    else; tt = sprintf('SRX pool pinned at x%g (isometric %.3f vs %.3f)', t, u.srx, C.E0.srx); end
    title(ax, tt, 'Color', ink, 'FontWeight', 'normal');
    legend(ax, 'Location', 'southoutside', 'Box', 'off', 'NumColumns', 1);
    if it > 1; set(ax, 'YTickLabel', []); end
end
title(tl, 'Can other knobs hide a different SRX pool?  grey band = within one preparation-width of the data', 'Color', ink);
exportgraphics(fg, fullfile(rsC, 'comp_residuals.png'), 'Resolution', 110);

%% ---- chi2 (data-SD units) and the project objective -----------------------------------
V = zeros(nT, 3); O = zeros(nT, 3);
for it = 1:nT
    a_ = get(C.TGT(it), 'A'); b_ = get(C.TGT(it), 'B');
    V(it, :) = [C.U{it}.chi2, a_.E.chi2, b_.E.chi2]; O(it, :) = [C.U{it}.cost, a_.E.cost, b_.E.cost];
end
xl = arrayfun(@(t) sprintf('pool x%g', t), C.TGT, 'uni', 0); xl(C.TGT == 1) = {'x1 (control)'};
fg2 = figure('Color', 'w', 'Position', [20 20 1500 520], 'Visible', 'off');
tl2 = tiledlayout(fg2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for pnl = 1:2
    ax = nexttile(tl2); hold(ax, 'on'); grid(ax, 'on'); ax.GridColor = grid_; ax.GridAlpha = 1; ax.XGrid = 'off';
    if pnl == 1; W = V; ref = C.E0.chi2; else; W = O; ref = C.E0.cost; end
    bh = bar(ax, 1:nT, W, 0.8, 'EdgeColor', 'w', 'LineWidth', 2);
    bh(1).FaceColor = cU; bh(2).FaceColor = cA; bh(3).FaceColor = cB;
    for j = 1:3
        xt = bh(j).XEndPoints;
        for it = 1:nT; text(ax, xt(it), W(it, j), sprintf('%.1f', W(it, j)), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'Color', ink, 'FontSize', 10); end
    end
    yline(ax, ref, '-', sprintf('incumbent %.2f', ref), 'Color', ink, 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    set(ax, 'XTick', 1:nT, 'XTickLabel', xl, 'XColor', ink, 'YColor', ink, 'FontSize', 12);
    if pnl == 1
        yline(ax, ref + 4, '--', 'incumbent + 4 (\approx 2\sigma)', 'Color', ink, 'LabelHorizontalAlignment', 'left');
        ylabel(ax, '\chi^2 = \Sigma z^2 over 17 features (data-SD units)');
        title(ax, 'Distance from the data in noise units (what identifiability asks)', 'Color', ink, 'FontWeight', 'normal');
        legend(bh, {'uncompensated', 'set A: cross-bridge knobs', 'set B: + SRX kinetics'}, 'Location', 'northeast', 'Box', 'off');
    else
        set(ax, 'YScale', 'log'); ylabel(ax, 'project objective (weighted, floored; log scale)');
        title(ax, 'The project objective (not what the refit minimised)', 'Color', ink, 'FontWeight', 'normal');
    end
end
exportgraphics(fg2, fullfile(rsC, 'comp_chi2.png'), 'Resolution', 110);

%% ---- knobs -------------------------------------------------------------------------
nK = size(C.KN, 1); yk = nK:-1:1;
fg3 = figure('Color', 'w', 'Position', [20 20 900 620], 'Visible', 'off');
ax = axes(fg3); hold(ax, 'on'); box(ax, 'off');
for k = 1:nK; plot(ax, [1/5 5], yk(k)*[1 1], '-', 'Color', grid_, 'HandleVisibility', 'off'); end
plot(ax, [1 1], [0.4 nK + 0.6], '-', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
mk = {'o', 'd', 's'};   % x0.5 hollow circle, x1 (control) diamond, x2 filled square
for it = 1:nT
    for s = {'A', 'B'}
        r = get(C.TGT(it), s{1}); col = cA; if s{1} == 'B'; col = cB; end
        fc = col; if C.TGT(it) < 1; fc = 'w'; elseif C.TGT(it) == 1; fc = [0.8 0.8 0.8]; end
        plot(ax, exp(r.theta), yk + (it - (nT + 1)/2) * 0.22, mk{it}, 'Color', col, 'MarkerFaceColor', fc, 'MarkerSize', 8, ...
             'LineWidth', 1.5, 'LineStyle', 'none', 'DisplayName', sprintf('pool x%g, set %s', C.TGT(it), s{1}));
    end
end
set(ax, 'XScale', 'log', 'XTick', [1/4 1/2 1 2 4], 'XTickLabel', {'/4', '/2', '1', 'x2', 'x4'}, ...
    'YTick', 1:nK, 'YTickLabel', flip(C.KN(:, 1)), 'TickLabelInterpreter', 'none', 'TickLength', [0 0], 'XColor', ink, 'YColor', ink, 'FontSize', 11);
xlim(ax, [1/5 5]); ylim(ax, [0.4 nK + 0.6]); xlabel(ax, 'knob multiplier chosen by the compensation');
legend(ax, 'Location', 'eastoutside', 'Box', 'off'); title(ax, 'What the refit changed (hollow = pool x0.5, grey diamond = x1 control, filled = pool x2)', 'Color', ink, 'FontWeight', 'normal');
exportgraphics(fg3, fullfile(rsC, 'comp_knobs.png'), 'Resolution', 110);

%% ---- traces: incumbent vs best compensated, and uncompensated -------------------------
best = cell(1, nT);
for it = 1:nT
    a = get(C.TGT(it), 'A'); b = get(C.TGT(it), 'B'); if b.E.chi2 <= a.E.chi2; best{it} = b; else; best{it} = a; end
end
iTL = @(s) find(strcmp(C.TL, s), 1);
cT = [cA; ink; cU];                              % pool x0.5 blue, (x1 unused), x2 orange
fg4 = figure('Color', 'w', 'Position', [20 20 1600 860], 'Visible', 'off');
tl4 = tiledlayout(fg4, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ttl = {'Slack-restretch cycle 1: force (kPa)', 'Slack-restretch cycle 1: SRX pool', 'k_{tr} protocol: F / F_{pre}', 'k_{tr} protocol: SRX pool'};
for p = 1:4
    ax = nexttile(tl4); hold(ax, 'on'); grid(ax, 'on'); ax.GridColor = grid_; ax.GridAlpha = 1;
    tc = 2 + (mod(p - 1, 2) == 1); tcol = 2 + (p >= 3);           % trace column; TR column (2 slack, 3 ktr)
    nrmf = @(T) 1; if p == 3; nrmf = @(T) interp1(T(:, 1), T(:, 2), -1e-3); end
    T0 = C.TR{1, tcol}; plot(ax, T0(:, 1), T0(:, tc) / nrmf(T0), '-', 'Color', ink, 'LineWidth', 2.4, 'DisplayName', 'incumbent');
    for it = find(C.TGT ~= 1)
        Tu = C.TR{iTL(sprintf('pool x%g, uncompensated', C.TGT(it))), tcol};
        Tb = C.TR{iTL(sprintf('pool x%g, set %s', C.TGT(it), best{it}.set)), tcol};
        plot(ax, Tu(:, 1), Tu(:, tc) / nrmf(Tu), ':', 'Color', cT(it, :), 'LineWidth', 1.8, 'DisplayName', sprintf('pool x%g uncompensated', C.TGT(it)));
        plot(ax, Tb(:, 1), Tb(:, tc) / nrmf(Tb), '-', 'Color', cT(it, :), 'LineWidth', 1.8, 'DisplayName', sprintf('pool x%g compensated (set %s)', C.TGT(it), best{it}.set));
    end
    if p <= 2; xlim(ax, [-0.03 0.55]); xlabel(ax, 't (s), release at 0, restretch at ~0.30 s'); else; xlim(ax, [-0.005 0.075]); xlabel(ax, 't (s), release at 0'); end
    title(ax, ttl{p}, 'Color', ink, 'FontWeight', 'normal'); ax.XColor = ink; ax.YColor = ink;
    if p == 1; legend(ax, 'Location', 'southeast', 'Box', 'off'); end
end
exportgraphics(fg4, fullfile(rsC, 'comp_traces.png'), 'Resolution', 110);

%% ---- numbers -----------------------------------------------------------------------
fid = fopen(fullfile(rsC, 'compensation.txt'), 'w'); c = onCleanup(@() fclose(fid));
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('SRX pool pinned and compensated (robust1_seed, 8 mM). chi2 = sum z^2, 17 features, data-SD units\n\n');
out('incumbent: chi2 %.2f  SRX %.4f  F_iso %.2f  objective %.4f\n', C.E0.chi2, C.E0.srx, C.E0.F, C.E0.cost);
for it = 1:nT
    u = C.U{it}; out('\npool x%g: uncompensated chi2 %.2f  SRX %.4f  F %.2f  objective %.4f\n', C.TGT(it), u.chi2, u.srx, u.F, u.cost);
    for s = {'A', 'B'}
        r = get(C.TGT(it), s{1});
        out('  set %s: chi2 %.2f  SRX %.4f  F %.2f  objective %.4f  steps tried %d\n', s{1}, r.E.chi2, r.E.srx, r.E.F, r.E.cost, numel(r.hist) - 1);
        out('    knobs: %s\n', strjoin(arrayfun(@(k) sprintf('%s x%.2f', C.KN{k, 1}, exp(r.theta(k))), find(abs(r.theta) > 0.02)', 'uni', 0), ', '));
        [~, o] = sort(abs(r.E.z), 'descend', 'MissingPlacement', 'last');
        out('    largest residuals: %s\n', strjoin(arrayfun(@(k) sprintf('%s %+.1f', C.E0.names{k}, r.E.z(k)), o(1:4)', 'uni', 0), ', '));
    end
end
out('\nJacobian: dlog(SRX)/dlog(knob): %s\n', strjoin(arrayfun(@(k) sprintf('%s %+.2f', C.KN{k, 1}, C.Js(k)), 1:nK, 'uni', 0), ', '));
fprintf('DONE PlotSrxCompensation\n');
