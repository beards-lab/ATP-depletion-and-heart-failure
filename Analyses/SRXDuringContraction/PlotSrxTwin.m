% PlotSrxTwin.m -- figures for RunSrxTwin (results/twin.mat).
%   twin_fim.png       Fisher information: pool resolution per data set, preparations needed,
%                      and which knobs the pool trades against
%   twin_profile.png   +-20 % twin profile: chi2 per data block after each refit
%   twin_protocolB.png the release-hold-restretch curves of the pinned and refitted fibres
%   twin_mc.png        the "noised seed" test, linearised Monte Carlo
%   twin.txt           the numbers
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsT = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
W = load(fullfile(rsT, 'twin.mat'));
ink = [0.13 0.13 0.13]; grid_ = [0.88 0.88 0.88]; band = [0.92 0.92 0.92];
c1 = [42 120 214]/255; c2 = [235 104 52]/255; c3 = [27 175 122]/255; c4 = [237 161 0]/255;
iC2 = @(z) sum(z(isfinite(z)).^2);
fx = @(x) iFx(x);                                  % x/factor label, capped
% extra rows: all data, but the load sensitivities taken as KNOWN (column dropped)
iS2 = find(strcmp(W.KN(:, 1), 'sigma2')); iS1 = find(strcmp(W.KN(:, 1), 'sigma_srd1'));
Jall = W.CFG(end).J; CF = W.CFG;
for ex = {iS2, [iS2 iS1]}
    keep = setdiff(1:size(Jall, 2), ex{1}); [sd, ~] = iPoolSD(Jall(:, keep), W.Js(keep), 0); sp = iPoolSD(Jall(:, keep), W.Js(keep), 1/log(2)^2);
    if numel(ex{1}) == 1; nmx = 'A + B, \sigma_2 known'; else; nmx = 'A + B, \sigma_2 & \sigma_{srd1} known'; end
    CF(end+1) = W.CFG(end); CF(end).name = nmx; CF(end).sd = sd; CF(end).sd_prior = sp; %#ok<SAGROW>
end
W.CFG = CF;
cfgCol = [0.55 0.55 0.55; c4; c1; c3; 0.30 0.23 0.65; 0.85 0.30 0.45];
nC = numel(W.CFG);

%% ---- Fisher information ------------------------------------------------------------
fg = figure('Color', 'w', 'Position', [20 20 1650 520], 'Visible', 'off');
tl = tiledlayout(fg, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
patch(ax, [0.8 1.2 1.2 0.8], [0.4 0.4 nC + 0.6 nC + 0.6], band, 'EdgeColor', 'none', 'HandleVisibility', 'off');
for c = 1:nC
    y = nC + 1 - c; s = W.CFG(c).sd; sp = W.CFG(c).sd_prior;
    plot(ax, exp([-2 2]*sp), [y y] - 0.15, '-', 'Color', cfgCol(c, :)*0.5 + 0.5, 'LineWidth', 6, 'HandleVisibility', 'off');
    plot(ax, exp([-2 2]*min(s, 5)), [y y] + 0.15, '-', 'Color', cfgCol(c, :), 'LineWidth', 6, 'HandleVisibility', 'off');
    text(ax, 12, y, fx(exp(2*s)), 'Color', ink, 'FontSize', 11, 'VerticalAlignment', 'middle');
end
plot(ax, [1 1], [0.4 nC + 0.6], '-', 'Color', ink, 'HandleVisibility', 'off');
set(ax, 'XScale', 'log', 'XTick', [0.25 0.5 0.8 1 1.25 2 4], 'XTickLabel', {'/4', '/2', '0.8', '1', '1.25', 'x2', 'x4'}, ...
    'YTick', 1:nC, 'YTickLabel', flip({W.CFG.name}), 'TickLength', [0 0], 'XColor', ink, 'YColor', ink, 'FontSize', 11);
xlim(ax, [0.2 30]); ylim(ax, [0.4 nC + 0.6]);
xlabel(ax, 'pool, 95 % range (\pm2 SD), one preparation');
title(ax, {'Attainable pool precision (all 13 knobs free)', 'dark: data only   light: knobs known within \times2   grey band: \pm20 %'}, 'Color', ink, 'FontWeight', 'normal');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); ax.GridColor = grid_; ax.GridAlpha = 1;
n = 1:20;
for c = 1:nC
    plot(ax, n, 100*(exp(2*W.CFG(c).sd ./ sqrt(n)) - 1), '-', 'Color', cfgCol(c, :), 'LineWidth', 2.2, 'DisplayName', W.CFG(c).name);
end
yline(ax, 20, '--', '\pm20 %', 'Color', ink, 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
set(ax, 'YScale', 'log', 'XColor', ink, 'YColor', ink, 'FontSize', 11); ylim(ax, [5 1e4]);
xlabel(ax, 'number of preparations'); ylabel(ax, '95 % pool uncertainty (%)');
legend(ax, 'Location', 'northeast', 'Box', 'off'); title(ax, 'How many preparations to resolve \pm20 %', 'Color', ink, 'FontWeight', 'normal');
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off'); nK = size(W.KN, 1);
cb = W.CFG(1).corr; cB = W.CFG(3).corr; yk = nK:-1:1;
for k = 1:nK; plot(ax, [-1 1], yk(k)*[1 1], '-', 'Color', grid_, 'HandleVisibility', 'off'); end
plot(ax, [0 0], [0.4 nK + 0.6], '-', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
barh(ax, yk + 0.18, cb, 0.35, 'FaceColor', cfgCol(1, :), 'EdgeColor', 'none', 'DisplayName', 'battery');
barh(ax, yk - 0.18, cB, 0.35, 'FaceColor', cfgCol(3, :), 'EdgeColor', 'none', 'DisplayName', 'battery + B');
set(ax, 'YTick', 1:nK, 'YTickLabel', flip(W.KN(:, 1)), 'TickLabelInterpreter', 'none', 'TickLength', [0 0], ...
    'XColor', ink, 'YColor', ink, 'FontSize', 11); xlim(ax, [-1 1]); ylim(ax, [0.4 nK + 0.6]);
xlabel(ax, 'correlation of the knob estimate with the pool estimate');
legend(ax, 'Location', 'southeast', 'Box', 'off'); title(ax, 'What the pool trades against', 'Color', ink, 'FontWeight', 'normal');
exportgraphics(fg, fullfile(rsT, 'twin_fim.png'), 'Resolution', 110);

%% ---- twin profile --------------------------------------------------------------------
nTg = numel(W.TGT); cases = {'pinned, uncompensated', 'refit to battery only', 'refit to battery + B'};
X = zeros(nTg, 3, 3);                      % target x case x block (battery, B, A)
for it = 1:nTg
    r = W.PIN{it}.r; X(it, 1, :) = [iC2(r.bat), iC2(r.B), iC2(r.A)];
    a = W.PR([W.PR.target] == W.TGT(it) & strcmp({W.PR.fit}, 'bat')); b = W.PR([W.PR.target] == W.TGT(it) & strcmp({W.PR.fit}, 'batB'));
    X(it, 2, :) = a.chi2; X(it, 3, :) = b.chi2;
end
fg2 = figure('Color', 'w', 'Position', [20 20 1400 520], 'Visible', 'off');
tl2 = tiledlayout(fg2, 1, nTg, 'TileSpacing', 'compact', 'Padding', 'compact');
for it = 1:nTg
    ax = nexttile(tl2); hold(ax, 'on'); grid(ax, 'on'); ax.GridColor = grid_; ax.GridAlpha = 1; ax.XGrid = 'off';
    V = squeeze(X(it, :, :)); bh = bar(ax, 1:3, max(V, 1e-3), 0.8, 'EdgeColor', 'w', 'LineWidth', 2);
    bh(1).FaceColor = cfgCol(1, :); bh(2).FaceColor = c1; bh(3).FaceColor = c4;
    for j = 1:3
        xt = bh(j).XEndPoints;
        for k = 1:3
            if V(k, j) < 0.05; continue; end              % below the axis floor: no label
            text(ax, xt(k), V(k, j), sprintf('%.1f', V(k, j)), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'Color', ink, 'FontSize', 10);
        end
    end
    yline(ax, 4, '--', '\Delta\chi^2 = 4 (\approx 2\sigma)', 'Color', ink, 'LabelHorizontalAlignment', 'left');
    set(ax, 'YScale', 'log', 'XTick', 1:3, 'XTickLabel', cases, 'YTick', [0.01 0.1 1 10 100 1000], 'YTickLabel', {'<0.01', '0.1', '1', '10', '100', '1000'}, 'XColor', ink, 'YColor', ink, 'FontSize', 11); ylim(ax, [1e-2 1e3]);
    ylabel(ax, '\chi^2 vs the twin truth (one preparation, noise units)');
    title(ax, sprintf('true pool, fitted at x%g', W.TGT(it)), 'Color', ink, 'FontWeight', 'normal');
    if it == 1; legend(bh, {'battery (17 features)', 'protocol B: release-hold k63, F20', 'protocol A: activation t50, t90'}, 'Location', 'northeast', 'Box', 'off'); end
end
title(tl2, 'Twin test: can a \pm20 % pool be hidden? Truth = incumbent; the refit chooses all 13 knobs', 'Color', ink);
exportgraphics(fg2, fullfile(rsT, 'twin_profile.png'), 'Resolution', 110);

%% ---- protocol B curves ------------------------------------------------------------------
fg3 = figure('Color', 'w', 'Position', [20 20 1400 520], 'Visible', 'off');
tl3 = tiledlayout(fg3, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
T = 1e3 * W.TB; fld = {'k63', 'F20'}; yl = {'k_{63} after restretch (s^{-1})', 'redeveloped fraction at 20 ms'}; sdv = {W.sdP(1:numel(T)), W.sdP(numel(T)+1:end)};
tc = [c1; c2];
for p = 1:2
    ax = nexttile(tl3); hold(ax, 'on'); grid(ax, 'on'); ax.GridColor = grid_; ax.GridAlpha = 1;
    v0 = W.E0.P.(fld{p});
    patch(ax, [T fliplr(T)], [v0 + sdv{p}, fliplr(v0 - sdv{p})], band, 'EdgeColor', 'none', 'DisplayName', 'truth \pm1 SD (one prep)');
    plot(ax, T, v0, '-', 'Color', ink, 'LineWidth', 2.4, 'DisplayName', 'truth (pool x1)');
    for it = 1:nTg
        u = W.PIN{it}.E.P.(fld{p});
        a = W.PR([W.PR.target] == W.TGT(it) & strcmp({W.PR.fit}, 'bat')).E.P.(fld{p});
        b = W.PR([W.PR.target] == W.TGT(it) & strcmp({W.PR.fit}, 'batB')).E.P.(fld{p});
        plot(ax, T, u, ':o', 'Color', tc(it, :), 'LineWidth', 1.6, 'MarkerSize', 5, 'DisplayName', sprintf('x%g uncompensated', W.TGT(it)));
        plot(ax, T, a, '--s', 'Color', tc(it, :), 'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', sprintf('x%g refit to battery', W.TGT(it)));
        plot(ax, T, b, '-d', 'Color', tc(it, :), 'MarkerFaceColor', tc(it, :), 'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', sprintf('x%g refit to battery + B', W.TGT(it)));
    end
    set(ax, 'XScale', 'log', 'XTick', T, 'XColor', ink, 'YColor', ink, 'FontSize', 11);
    xlabel(ax, 'hold at short length T (ms)'); ylabel(ax, yl{p});
    if p == 1; legend(ax, 'Location', 'eastoutside', 'Box', 'off'); end
end
title(tl3, 'Protocol B (release 0.15 ML, hold T, restretch): does it expose a pool the battery refit hid?', 'Color', ink);
exportgraphics(fg3, fullfile(rsT, 'twin_protocolB.png'), 'Resolution', 110);

%% ---- noised-seed Monte Carlo ---------------------------------------------------------------
fg4 = figure('Color', 'w', 'Position', [20 20 1300 700], 'Visible', 'off');
tl4 = tiledlayout(fg4, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
cfgs = unique({W.MC.cfg}, 'stable'); prs = {'none', 'seed 5 %'};
ed = linspace(log(0.25), log(4), 121);
for ip = 1:2
    for ic = 1:numel(cfgs)
        ax = nexttile(tl4); hold(ax, 'on'); grid(ax, 'on'); ax.GridColor = grid_; ax.GridAlpha = 1;
        for it = 1:nTg
            m = W.MC(strcmp({W.MC.cfg}, cfgs{ic}) & strcmp({W.MC.prior}, prs{ip}) & abs([W.MC.truth] - log(W.TGT(it))) < 1e-9);
            e = max(ed(1), min(ed(end), m.est)); h = histcounts(e, ed, 'Normalization', 'pdf');
            stairs(ax, exp(ed), [h h(end)], '-', 'Color', tc(it, :), 'LineWidth', 1.8, 'DisplayName', ...
                   sprintf('truth x%g: estimate x%.2f, 95 %% %s', W.TGT(it), exp(mean(m.est)), fx(exp(2*std(m.est)))));
            xline(ax, W.TGT(it), '--', 'Color', tc(it, :), 'HandleVisibility', 'off');
        end
        set(ax, 'XScale', 'log', 'XTick', [0.25 0.5 0.8 1 1.2 2 4], 'XColor', ink, 'YColor', ink, 'FontSize', 11); xlim(ax, [0.25 4]);
        xlabel(ax, 'recovered pool (x incumbent)'); ylabel(ax, 'density');
        if ip == 1; pt = 'data only (seed irrelevant)'; else; pt = 'optimiser stays near a 5 % noised seed'; end
        title(ax, sprintf('%s  |  %s', cfgs{ic}, pt), 'Color', ink, 'FontWeight', 'normal');
        legend(ax, 'Location', 'northwest', 'Box', 'off');
    end
end
title(tl4, 'The proposed noised-seed test, linearised (4000 draws, one preparation of measurement noise)', 'Color', ink);
exportgraphics(fg4, fullfile(rsT, 'twin_mc.png'), 'Resolution', 110);

%% ---- numbers ----------------------------------------------------------------------------
fid = fopen(fullfile(rsT, 'twin.txt'), 'w'); cl = onCleanup(@() fclose(fid));
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('Twin identifiability of the SRX pool (truth robust1_seed, pool %.4f), one preparation of noise\n\n', W.E0.srx);
out('Fisher information, SD of log(pool), all 13 knobs free:\n');
for c = 1:nC; out('  %-16s SD %.3f -> 95 %% x/%.2f | with x2 prior x/%.2f | preps for +-20 %%: %.0f\n', W.CFG(c).name, ...
        W.CFG(c).sd, exp(2*W.CFG(c).sd), exp(2*W.CFG(c).sd_prior), ceil((2*W.CFG(c).sd/log(1.2))^2)); end
out('\nTwin profile (chi2 battery / B / A):\n');
for it = 1:nTg
    out('  x%g  uncompensated %s\n', W.TGT(it), mat2str(squeeze(X(it, 1, :))', 3));
    out('       refit battery %s   refit battery+B %s\n', mat2str(squeeze(X(it, 2, :))', 3), mat2str(squeeze(X(it, 3, :))', 3));
end
out('\nMonte Carlo (estimate, 95 %% factor):\n');
for m = W.MC; out('  %-12s %-9s truth x%.2f: x%.2f  x/%.2f\n', m.cfg, m.prior, exp(m.truth), exp(mean(m.est)), exp(2*std(m.est))); end
fprintf('DONE PlotSrxTwin\n');

function s = iFx(x)
    if x > 1000; s = '> \times/1000'; else; s = sprintf('\\times/%.2f', x); end
end

function [sd, corr] = iPoolSD(J, c, prior)
    F = J' * J + prior * eye(size(J, 2)) + 1e-10 * eye(size(J, 2));
    S = inv(F); sd = sqrt(c * S * c'); corr = (S * c') ./ sqrt(diag(S) * (c * S * c'));
end
