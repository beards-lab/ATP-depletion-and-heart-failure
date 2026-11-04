% SummarizeSrxRateKicks.m -- what the SRX pool does to the features (RunSrxRateKicks).
% SRX is read as ONE pool:
%   SRX+  inflow  x f   (ksr0 and ksrd together: PT -> SRX.ATP, PD -> SRX.ADP.Pi)
%   SRX-  outflow x f   (kmsrd: SRX.ADP.Pi -> PD)
% f = 1/11, 1/3, 3, 11 (x11 = +1000 %). UNCOMPENSATED: force is free to move.
% Each shift is given in % of the incumbent and in data-SD units
% (shift / (|data| * relative between-prep SD)).
% Output: results/kicks.txt, results/kicks_summary.json, results/kicks_dose.png,
%         results/kicks_heatmap.png, results/kicks_traces.png, results/kicks_fv.png
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsU = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
S = load(fullfile(rsU, 'kicks.mat')); RK = S.RK; nm = {RK.name};
PT8 = pooledTarget8mM(); KT8 = pooledKtrTarget(8);
ix = @(n) find(strcmp(nm, n), 1);
dose = [1/11 1/3 1 3 11];
serP = arrayfun(@(k) ix(k{1}), {'in /11', 'in /3', 'ref', 'in x3', 'in x11'});
serM = arrayfun(@(k) ix(k{1}), {'out /11', 'out /3', 'ref', 'out x3', 'out x11'});

% ---- features: {field, label, segments, source P(ooled slack)/K(tr protocol)/I(so)/F(-V)}
FL = { 'F_iso',    'isometric force (kPa)',       [], 'I'
       'SRX_iso',  'isometric SRX fraction',      [], 'I'
       'steady',   'steady force (kPa)',          1:5, 'P'
       'A',        'redeveloped force A (kPa)',   1:5, 'P'
       'peak1_y',  're-stretch peak (kPa)',       1:5, 'P'
       'vall_y',   'valley (kPa)',                1:5, 'P'
       'ktrfit',   'slack k_{tr} fit (s^{-1})',   1:5, 'P'
       'ktr80',    'slack k_{tr}80 (s^{-1})',     1:5, 'P'
       't0',       'slack delay t_0 (s)',         1:5, 'P'
       'rsK80',    'post-restretch rate (s^{-1})',1:4, 'P'
       'vall_t',   'valley time (s)',             1:5, 'P'
       'ktrP_k80', 'k_{tr} protocol k80 (s^{-1})',1,   'K'
       'ktrP_k63', 'k_{tr} protocol k63 (s^{-1})',1,   'K'
       'ktrP_F20', 'k_{tr} protocol F(20 ms)/F',  1,   'K'
       'FVp2',     'F-V power, 2 ML/s',           3,   'F'
       'FVp4',     'F-V power, 4 ML/s',           4,   'F' };
nF = size(FL, 1); nV = numel(RK);
M = nan(nV, nF); D = nan(1, nF); SD = nan(1, nF);
for i = 1:nV
    fm = RK(i).fm; if ~isstruct(fm); continue; end
    for k = 1:nF
        f = FL{k, 1}; sg = FL{k, 3};
        switch FL{k, 4}
            case 'I'
                if strcmp(f, 'F_iso'); M(i, k) = RK(i).iso(1); else; M(i, k) = RK(i).iso(2); end
            case {'P', 'K'}
                if isfield(fm, f); x = fm.(f); M(i, k) = mean(x(sg(sg <= numel(x)))); end
            case 'F'
                if isfield(fm, 'FV_normpowerAvg'); M(i, k) = fm.FV_normpowerAvg(sg); end
        end
        if i == 1
            switch FL{k, 4}
                case 'P'; j = find(strcmp(PT8.name, f), 1); D(k) = mean(PT8.val(j, sg)); SD(k) = abs(D(k))*PT8.s(j);
                case 'K'; j = find(strcmp(KT8.name, f), 1);
                    if ~isempty(j); D(k) = KT8.val(j); SD(k) = abs(D(k))*KT8.s(j); end
                case 'F'; D(k) = RK(i).fd.FV_normpowerAvg(sg);
                    SD(k) = max(sqrt(RK(i).fd.FV_normpowerVar(sg)), 0.136*abs(D(k)));
            end
        end
    end
end
iref = ix('ref');
P = 100 * (M ./ M(iref, :) - 1); Z = (M - M(iref, :)) ./ SD;

%% ---- text ----------------------------------------------------------------------
fid = fopen(fullfile(rsU, 'kicks.txt'), 'w'); c = onCleanup(@() fclose(fid));
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('SRX as one pool, UNCOMPENSATED kicks on robust1_seed (8 mM).\n');
out('SRX+ = inflow (ksr0 & ksrd) x f;  SRX- = outflow (kmsrd) x f.  value  [shift in data-SD units]\n\n');
for sr = 1:2
    if sr == 1; ser = serP; t = 'SRX+ (inflow)'; else; ser = serM; t = 'SRX- (outflow)'; end
    out('%s\n%-30s %9s %9s | %s\n', t, 'feature', 'data', 'SD', strjoin(arrayfun(@(d) sprintf('%17s', iDose(d)), dose, 'uni', 0), ''));
    for k = 1:nF
        out('%-30s %9.4g %9.3g |', regexprep(FL{k, 2}, '[{}\\^_]', ''), D(k), SD(k));
        for j = 1:numel(ser)
            if j == 3; out(' %9.4g        ', M(ser(j), k));
            else; out(' %9.4g [%+5.1f]', M(ser(j), k), Z(ser(j), k)); end
        end
        out('\n');
    end
    out('\n');
end
out('Side runs (single rates / hydrolysis inside SRX), shift in data-SD units:\n');
for n = {'ksr0 x11', 'ksrd x11', 'in+out x11', 'ksr2srd x11', 'ksrd2sr x11'}
    i = ix(n{1}); if isempty(i); continue; end
    out('%-12s SRX_iso %.3f F_iso %.1f | %s\n', n{1}, M(i, 2), M(i, 1), strjoin(arrayfun(@(k) sprintf('%s %+.1f', FL{k, 1}, Z(i, k)), 3:nF, 'uni', 0), ', '));
end

%% ---- palette (reference instance: categorical slots, blue<->red diverging) --------
ink = [0.13 0.13 0.13]; grid_ = [0.87 0.87 0.87]; dataCol = [0.70 0.70 0.70];
cP = [42 120 214]/255; cM = [235 104 52]/255;                 % slot 1 blue, slot 2 orange
slot = [42 120 214; 235 104 52; 27 175 122; 237 161 0]/255;

%% ---- dose-response small multiples ------------------------------------------------
show = 2:nF;
fg = figure('Color', 'w', 'Position', [20 20 1700 980], 'Visible', 'off');
tl = tiledlayout(fg, 3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = show
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); ax.GridColor = grid_; ax.GridAlpha = 1;
    if isfinite(D(k))
        patch(ax, [dose(1) dose(end) dose(end) dose(1)], D(k) + SD(k)*[-1 -1 1 1], dataCol, ...
              'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax, dose([1 end]), D(k)*[1 1], '-', 'Color', dataCol, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    plot(ax, dose, M(serP, k), '-o', 'Color', cP, 'MarkerFaceColor', cP, 'MarkerSize', 7, 'LineWidth', 2, 'DisplayName', 'SRX+ (inflow)');
    plot(ax, dose, M(serM, k), '-o', 'Color', cM, 'MarkerFaceColor', cM, 'MarkerSize', 7, 'LineWidth', 2, 'DisplayName', 'SRX- (outflow)');
    plot(ax, 1, M(iref, k), 'o', 'Color', ink, 'MarkerFaceColor', ink, 'MarkerSize', 8, 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log', 'XTick', dose, 'XTickLabel', {'/11', '/3', '1', 'x3', 'x11'}, 'XColor', ink, 'YColor', ink, 'FontSize', 10);
    xlim(ax, [dose(1)/1.3 dose(end)*1.3]); title(ax, FL{k, 2}, 'FontWeight', 'normal', 'Color', ink);
    if k == show(1); legend(ax, 'Location', 'best', 'Box', 'off'); end
end
title(tl, 'SRX as one pool: inflow (SRX+) or outflow (SRX-) scaled, uncompensated  |  grey band = pooled data \pm SD', 'Color', ink);
exportgraphics(fg, fullfile(rsU, 'kicks_dose.png'), 'Resolution', 110);

%% ---- heat-map of the dose series (data-SD units) ----------------------------------
rows = [serP([1 2 4 5]), serM([1 2 4 5])];
rlab = {'SRX+ /11', 'SRX+ /3', 'SRX+ x3', 'SRX+ x11', 'SRX- /11', 'SRX- /3', 'SRX- x3', 'SRX- x11'};
cols = find(isfinite(SD));
blue = cP; red = [227 73 72]/255; mid = [240 239 236]/255;
tt = linspace(-1, 1, 129)'; cmap = (tt < 0).*((1+tt).*mid + (-tt).*blue) + (tt >= 0).*((1-tt).*mid + tt.*red);
Zc = max(-6, min(6, Z(rows, cols))); Zc(isnan(Z(rows, cols))) = NaN;   % min/max drop NaN
fg2 = figure('Color', 'w', 'Position', [20 20 1500 520], 'Visible', 'off');
ax = axes(fg2); imagesc(ax, Zc, [-6 6]); colormap(ax, cmap);
set(ax, 'XTick', 1:numel(cols), 'XTickLabel', FL(cols, 2), 'XTickLabelRotation', 30, 'YTick', 1:numel(rows), ...
    'YTickLabel', rlab, 'TickLength', [0 0], 'FontSize', 11, 'XColor', ink, 'YColor', ink);
hold(ax, 'on');
[ri, ki] = find(isnan(Zc));                       % feature not extractable: neutral cell
for n = 1:numel(ri); patch(ax, ki(n) + [-.5 .5 .5 -.5], ri(n) + [-.5 -.5 .5 .5], [0.97 0.97 0.97], 'EdgeColor', 'none'); end
for i = 1:numel(rows)
    for k = 1:numel(cols)
        v = Z(rows(i), cols(k)); if isnan(v); s = '–'; else; s = sprintf('%+.1f', v); end
        text(ax, k, i, s, 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', ink);
    end
end
for k = 0.5:1:numel(cols) + 0.5; plot(ax, [k k], [0.5 numel(rows) + 0.5], 'w-', 'LineWidth', 2); end
for i = 0.5:1:numel(rows) + 0.5; plot(ax, [0.5 numel(cols) + 0.5], [i i], 'w-', 'LineWidth', 2); end
plot(ax, [0.5 numel(cols) + 0.5], [4.5 4.5], '-', 'Color', ink, 'LineWidth', 1.5);
cb = colorbar(ax); cb.Label.String = 'shift vs incumbent, data-SD units (clipped \pm6)'; cb.Color = ink;
title(ax, 'How far each feature moves when the SRX inflow or outflow is scaled (uncompensated)', 'Color', ink);
exportgraphics(fg2, fullfile(rsU, 'kicks_heatmap.png'), 'Resolution', 120);

%% ---- traces: slack cycle 1 and the ktr protocol ------------------------------------
tsel = {'in x11', 'in /11', 'out x11', 'out /11'}; tlab = {'SRX+ x11', 'SRX+ /11', 'SRX- x11', 'SRX- /11'};
fg3 = figure('Color', 'w', 'Position', [20 20 1600 860], 'Visible', 'off');
tl3 = tiledlayout(fg3, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:4
    ax = nexttile(tl3); hold(ax, 'on'); grid(ax, 'on'); ax.GridColor = grid_; ax.GridAlpha = 1;
    isK = p >= 3; col = 2 + mod(p + 1, 2);                       % p=1,3 force; p=2,4 SRX
    if isK; fld = 'traceK'; else; fld = 'trace'; end
    tr0 = RK(iref).(fld); nrm = 1;
    if isK && col == 2; nrm = interp1(tr0(:, 1), tr0(:, 2), -1e-3); end
    if col == 2
        if isK; plot(ax, S.dataK(:, 1), S.dataK(:, 3), '-', 'Color', dataCol, 'LineWidth', 1.2, 'DisplayName', 'data 03/27');
        else;   plot(ax, S.dataC(:, 1), S.dataC(:, 3), '-', 'Color', dataCol, 'LineWidth', 1.2, 'DisplayName', 'data 03/27'); end
    end
    plot(ax, tr0(:, 1), tr0(:, col) / nrm, '-', 'Color', ink, 'LineWidth', 2.2, 'DisplayName', 'incumbent');
    for s = 1:numel(tsel)
        tr = RK(ix(tsel{s})).(fld); n2 = 1;
        if isK && col == 2; n2 = interp1(tr(:, 1), tr(:, 2), -1e-3); end
        plot(ax, tr(:, 1), tr(:, col) / n2, '-', 'Color', slot(s, :), 'LineWidth', 1.8, 'DisplayName', tlab{s});
    end
    ax.XColor = ink; ax.YColor = ink;
    if isK; xlim(ax, [-0.005 S.vtK(end, 1)]); xlabel(ax, 't (s), release at 0');
    else;   xlim(ax, [-0.03 S.vtC(end, 1)]); xlabel(ax, 't (s), release at 0, restretch at ~0.30 s'); end
    if col == 2 && isK; ylabel(ax, 'F / F_{pre}'); title(ax, 'k_{tr} protocol: force (each normalised to its own pre-release force)', 'Color', ink);
    elseif col == 2; ylabel(ax, 'Force (kPa)'); title(ax, 'Slack-restretch cycle 1: force', 'Color', ink);
    elseif isK; ylabel(ax, 'SRX'); title(ax, 'k_{tr} protocol: SRX pool', 'Color', ink);
    else; ylabel(ax, 'SRX'); title(ax, 'Slack-restretch cycle 1: SRX pool', 'Color', ink); end
    if p == 1; legend(ax, 'Location', 'southeast', 'Box', 'off'); end
end
exportgraphics(fg3, fullfile(rsU, 'kicks_traces.png'), 'Resolution', 110);

%% ---- F-V power --------------------------------------------------------------------
v = [0.5 1 2 4]; Dfv = RK(iref).fd.FV_normpowerAvg(1:4); SDfv = max(sqrt(RK(iref).fd.FV_normpowerVar(1:4)), 0.136*abs(Dfv));
fg4 = figure('Color', 'w', 'Position', [20 20 820 520], 'Visible', 'off');
ax = axes(fg4); hold(ax, 'on'); grid(ax, 'on'); ax.GridColor = grid_; ax.GridAlpha = 1;
errorbar(ax, v, Dfv, SDfv, 'o', 'Color', [0.55 0.55 0.55], 'MarkerFaceColor', dataCol, 'MarkerSize', 8, 'LineWidth', 1.2, 'DisplayName', 'data (Baker) \pm SD');
plot(ax, v, RK(iref).fm.FV_normpowerAvg(1:4), '-o', 'Color', ink, 'MarkerFaceColor', ink, 'LineWidth', 2.2, 'MarkerSize', 8, 'DisplayName', 'incumbent');
for s = 1:numel(tsel)
    plot(ax, v, RK(ix(tsel{s})).fm.FV_normpowerAvg(1:4), '-o', 'Color', slot(s, :), 'MarkerFaceColor', slot(s, :), 'LineWidth', 1.8, 'MarkerSize', 8, 'DisplayName', tlab{s});
end
set(ax, 'XScale', 'log', 'XTick', v, 'XColor', ink, 'YColor', ink); xlim(ax, [0.4 5]);
xlabel(ax, 'shortening velocity (ML/s)'); ylabel(ax, 'normalised power (P/P_0)\cdot v');
legend(ax, 'Location', 'northwest', 'Box', 'off'); title(ax, 'Force-velocity power under the SRX kicks', 'Color', ink);
exportgraphics(fg4, fullfile(rsU, 'kicks_fv.png'), 'Resolution', 110);

J = struct('kicks', {nm}, 'features', {FL(:, 1)'}, 'data', D, 'sd', SD, 'model', M, 'pct', P, 'z_shift', Z, ...
    'dose', dose, 'series_inflow', {nm(serP)}, 'series_outflow', {nm(serM)});
fj = fopen(fullfile(rsU, 'kicks_summary.json'), 'w'); fprintf(fj, '%s', jsonencode(J, 'PrettyPrint', true)); fclose(fj);
fprintf('DONE SummarizeSrxRateKicks\n');

function s = iDose(d)
    if d == 1; s = 'incumbent'; elseif d < 1; s = sprintf('/%g', round(1/d)); else; s = sprintf('x%g', d); end
end
