% PlotSrxDrug.m -- figures for RunSrxDrug (results/drug.mat).
%   drug_srx.png      model-predicted ACTIVATED SRX, control -> drug, for fibres that are
%                     indistinguishable at control (does the prediction depend on which one?)
%   drug_effects.png  the measurable drug effect on key features, per fibre
%   drug_discrim.png  can the drug effect tell those fibres apart? (chi2 of the difference)
%   drug.txt          the numbers
% Noise of a drug EFFECT measured within one preparation = sqrt(2) x one-prep SD.
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsD = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
S = load(fullfile(rsD, 'drug.mat')); FB = S.FB; DR = S.DR; R = S.R; nF = numel(FB); nD = numel(DR);
ink = [0.13 0.13 0.13]; grid_ = [0.88 0.88 0.88];
col = [42 120 214; 235 104 52; 27 175 122; 237 161 0; 232 123 164; 74 58 167]/255;   % slots 1-5, 7
get = @(i, d) R([R.fibre] == i & [R.drug] == d).E;

% feature vector of one run and its one-prep SD
[~, names, labs] = srxResiduals(get(1, 1).fm, get(1, 1).fd);
[~, ~, ~, ~, ~, sdb] = srxResiduals(get(1, 1).fm, get(1, 1).fd);
vec = @(E) iVec(E);
P0 = get(1, 1).P; nT = numel(P0.k63);
sdAll = [sdb; 0.08 * P0.k63(:); 0.05 * P0.F20(:); 0.02 * get(1, 1).F];         % + isometric force 2 %
labAll = [labs; arrayfun(@(t) sprintf('B k63, T=%g ms', 1e3*t), S.TB(:), 'uni', 0); ...
          arrayfun(@(t) sprintf('B F20, T=%g ms', 1e3*t), S.TB(:), 'uni', 0); {'isometric force'}];
sdDiff = sqrt(2) * sdAll;

%% ---- numbers ------------------------------------------------------------------------
EF = zeros(numel(sdAll), nF, nD); SRX = zeros(nF, nD); FIS = zeros(nF, nD);
for i = 1:nF
    v0 = vec(get(i, 1));
    for d = 1:nD
        E = get(i, d); EF(:, i, d) = vec(E) - v0; SRX(i, d) = E.srx; FIS(i, d) = E.F;
    end
end
grp = {FB.group}; ref = struct('twin', find(strcmp({FB.name}, 'twin truth x1')), 'real', find(strcmp({FB.name}, 'real refit x1 (control)')));
fid = fopen(fullfile(rsD, 'drug.txt'), 'w'); c = onCleanup(@() fclose(fid));
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('Drug test: fibres identical at control, different activated SRX pools (8 mM)\n\n');
out('%-26s %9s | %s\n', 'fibre', 'SRX ctrl', strjoin(arrayfun(@(d) sprintf('%-28s', DR(d).name), 2:nD, 'uni', 0), ''));
for i = 1:nF
    out('%-26s %9.4f |', FB(i).name, SRX(i, 1));
    for d = 2:nD; out('  SRX %.4f (x%.2f) F %+5.1f%%    ', SRX(i, d), SRX(i, d)/SRX(i, 1), 100*(FIS(i, d)/FIS(i, 1) - 1)); end
    out('\n');
end
out('\nDiscrimination: chi2 of (drug effect of fibre - drug effect of its group reference), within-prep SD\n');
CH = nan(nF, nD);
for i = 1:nF
    r0 = ref.(grp{i}); if i == r0; continue; end
    for d = 2:nD
        z = (EF(:, i, d) - EF(:, r0, d)) ./ sdDiff; CH(i, d) = sum(z(isfinite(z)).^2);
        [~, o] = sort(abs(z), 'descend', 'MissingPlacement', 'last');
        out('  %-24s vs %-24s %-14s chi2 %6.2f | top: %s\n', FB(i).name, FB(r0).name, DR(d).name, CH(i, d), ...
            strjoin(arrayfun(@(k) sprintf('%s %+.1f', labAll{k}, z(k)), o(1:3)', 'uni', 0), ', '));
    end
end
% control-level check: are the fibres really indistinguishable at control?
out('\nControl check, chi2 vs group reference (one-prep SD): ');
for i = 1:nF
    r0 = ref.(grp{i}); if i == r0; continue; end
    z = (vec(get(i, 1)) - vec(get(r0, 1))) ./ sdAll; out('%s %.2f; ', FB(i).name, sum(z(isfinite(z)).^2));
end
out('\n');

%% ---- figure 1: predicted activated SRX ------------------------------------------------
fg = figure('Color', 'w', 'Position', [20 20 1300 520], 'Visible', 'off');
tl = tiledlayout(fg, 1, nD - 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for d = 2:nD
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); ax.GridColor = grid_; ax.GridAlpha = 1;
    for i = 1:nF
        plot(ax, [1 2], [SRX(i, 1) SRX(i, d)], '-o', 'Color', col(i, :), 'MarkerFaceColor', col(i, :), 'LineWidth', 2, ...
             'MarkerSize', 8, 'DisplayName', sprintf('%s (x%.2f)', FB(i).name, SRX(i, d)/SRX(i, 1)));
    end
    set(ax, 'XTick', [1 2], 'XTickLabel', {'control', DR(d).name}, 'XColor', ink, 'YColor', ink, 'FontSize', 12); xlim(ax, [0.7 2.3]);
    ylabel(ax, 'activated SRX (isometric, SL 2.2)'); title(ax, sprintf('model-predicted activated SRX: %s', DR(d).name), 'Color', ink, 'FontWeight', 'normal');
    legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 10);
end
title(tl, 'Fibres that fit the control data equally well: what does each predict for the drug?', 'Color', ink);
exportgraphics(fg, fullfile(rsD, 'drug_srx.png'), 'Resolution', 110);

%% ---- figure 2: measurable drug effects -------------------------------------------------
key = [numel(sdAll), find(strcmp(names, 'ktrfit')), find(strcmp(names, 'ktrP_k63')), find(strcmp(names, 'ktrP_F20')), ...
       find(strcmp(names, 'FVp4')), find(strcmp(names, 'steady')), numel(sdb) + 2, numel(sdb) + 4];
fg2 = figure('Color', 'w', 'Position', [20 20 1500 560], 'Visible', 'off');
tl2 = tiledlayout(fg2, 1, nD - 1, 'TileSpacing', 'compact', 'Padding', 'compact');
yk = numel(key):-1:1;
for d = 2:nD
    ax = nexttile(tl2); hold(ax, 'on'); box(ax, 'off');
    for k = 1:numel(key); plot(ax, [-10 10], yk(k)*[1 1], '-', 'Color', grid_, 'HandleVisibility', 'off'); end
    plot(ax, [0 0], [0.4 numel(key) + 0.6], '-', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    for i = 1:nF
        z = EF(key, i, d) ./ sdDiff(key);
        plot(ax, z, yk + (i - (nF + 1)/2) * 0.1, 'o', 'Color', col(i, :), 'MarkerFaceColor', col(i, :), 'MarkerSize', 7, ...
             'LineStyle', 'none', 'DisplayName', FB(i).name);
    end
    set(ax, 'YTick', 1:numel(key), 'YTickLabel', flip(labAll(key)), 'TickLength', [0 0], 'XColor', ink, 'YColor', ink, 'FontSize', 11);
    ylim(ax, [0.4 numel(key) + 0.6]); xlabel(ax, 'drug effect (drug - control), in within-prep SD units');
    title(ax, DR(d).name, 'Color', ink, 'FontWeight', 'normal');
    if d == 2; legend(ax, 'Location', 'southoutside', 'Box', 'off', 'NumColumns', 2); end
end
title(tl2, 'What the drug would do to measurable features: do the fibres disagree?', 'Color', ink);
exportgraphics(fg2, fullfile(rsD, 'drug_effects.png'), 'Resolution', 110);

%% ---- figure 3: discrimination -------------------------------------------------------------
fg3 = figure('Color', 'w', 'Position', [20 20 900 480], 'Visible', 'off');
ax = axes(fg3); hold(ax, 'on'); grid(ax, 'on'); ax.GridColor = grid_; ax.GridAlpha = 1; ax.XGrid = 'off';
ii = find(~isnan(CH(:, 2))); V = CH(ii, 2:nD);
bh = bar(ax, 1:numel(ii), max(V, 1e-3), 0.8, 'EdgeColor', 'w', 'LineWidth', 2);
bh(1).FaceColor = col(1, :); if numel(bh) > 1; bh(2).FaceColor = col(2, :); end
for j = 1:numel(bh)
    xt = bh(j).XEndPoints;
    for k = 1:numel(ii); text(ax, xt(k), max(V(k, j), 1e-3), sprintf('%.1f', V(k, j)), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Color', ink, 'FontSize', 10); end
end
yline(ax, 4, '--', '\chi^2 = 4 (\approx 2\sigma)', 'Color', ink, 'LabelHorizontalAlignment', 'left');
xl = arrayfun(@(i) sprintf('%s vs %s', erase(FB(i).name, {'twin ', 'real '}), erase(FB(ref.(grp{i})).name, {'twin ', 'real ', ' (control)'})), ii, 'uni', 0);
set(ax, 'YScale', 'log', 'XTick', 1:numel(ii), 'XTickLabel', xl, 'XColor', ink, 'YColor', ink, 'FontSize', 11); ylim(ax, [1e-2 1e3]);
ylabel(ax, '\chi^2 of the drug-effect difference (one prep)');
legend(bh, {DR(2:end).name}, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
exportgraphics(fg3, fullfile(rsD, 'drug_discrim.png'), 'Resolution', 110);
fprintf('DONE PlotSrxDrug\n');

function v = iVec(E)
    [~, ~, ~, m] = srxResiduals(E.fm, E.fd);
    v = [m; E.P.k63(:); E.P.F20(:); E.F];
end
