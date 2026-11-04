% SummarizeSrxFamily.m -- tables and figures from RunSrxFamily (reads the per-variant
% results/family_run_<k>.mat, so it also works on a partial run).
% Output: results/family.txt, results/family_summary.json, results/family_cost.png
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsS = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
dd = dir(fullfile(rsS, 'family_run_*.mat'));
[~, o] = sort(cellfun(@(s) sscanf(s, 'family_run_%d.mat'), {dd.name})); dd = dd(o);
RS = struct([]);
for i = 1:numel(dd); L = load(fullfile(rsS, dd(i).name)); RS = [RS, L.R1]; end %#ok<AGROW>
p0 = getParams(loadParams('robust1_seed'), [], true, false); fnS = p0.fn;
prim = regexprep(cellfun(@(s) strtok(s, '|'), fnS, 'uni', 0), '\[.*$', '');
prim{3} = 'bounds(XTOR,SRX,att,PT)';
iref = find(strcmp({RS.fam}, 'ref'), 1);

gv = @(f, nm, k) iGet(f, nm, k);
lbl = @(r) sprintf('%s x%g', r.fam, r.x);
fid = fopen(fullfile(rsS, 'family.txt'), 'w'); c = onCleanup(@() fclose(fid));
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('SRX family on robust1_seed, ka-compensated to F_iso(2.2) (8 mM, pooled target)\n');
out('data cost = objective without the SRX_ss literature bound\n\n');
out('%-10s %7s %8s %8s %8s %8s | %7s %7s %7s %7s %7s %7s\n', 'variant', 'ka x', 'SRX_ss', 'att', ...
    'cost', 'dataC', 'steady', 'ktrfit', 'rsK80', 't0 ms', 'ktrP80', 'XTOR');
for i = 1:numel(RS)
    r = RS(i); fm = r.fm;
    out('%-10s %7.3f %8.4f %8.4f %8.4f %8.4f | %7.2f %7.2f %7.2f %7.1f %7.2f %7.2f\n', lbl(r), r.ka/p0.ka, ...
        gv(fm, 'SRX_ss', 1), gv(fm, 'attached_ss', 1), r.total, r.dataTotal, gv(fm, 'steady', 1), ...
        mean(gv(fm, 'ktrfit', 1:5)), mean(gv(fm, 'rsK80', 1:4)), 1e3*mean(gv(fm, 't0', 1:5)), ...
        gv(fm, 'ktrP_k80', 1), gv(fm, 'XTOR', 1));
end
if ~isempty(iref)
    out('\nLargest term changes vs ref (data-cost terms, chi-2 units)\n');
    cref = RS(iref).costNoSRX(:)';
    for i = 1:numel(RS)
        if i == iref || ~isfinite(RS(i).dataTotal); continue; end
        dc = RS(i).costNoSRX(:)' - cref; [~, j] = sort(abs(dc), 'descend'); j = j(1:min(4, end));
        out('%-10s d=%+.3f : %s\n', lbl(RS(i)), sum(dc), strjoin(arrayfun(@(k) sprintf('%s %+.3f', prim{k}, dc(k)), j, 'uni', 0), ', '));
    end
end

%% objective-independent distance: z = (model - data) / SD, SD from the data
% mechanics: pooled 8 mM target, relative between-preparation SD (segment means);
% ktr protocol: pooledKtrTarget(8); F-V: normalised power vs the Baker sets' mean,
% SD = max(across-set SD, 13.6 % of the value -- DataPooling s6 reproducibility).
PT8 = pooledTarget8mM(); KT8 = pooledKtrTarget(8);
ZF = {'A', 1:5; 'steady', 1:5; 'peak1_y', 1:5; 'vall_y', 1:5; 'ktrfit', 1:5; 'rsK80', 1:4; 't0', 1:5};
zN = [ZF(:, 1)', {'ktrP_k80'}, arrayfun(@(k) sprintf('FVpow%d', k), 1:4, 'uni', 0)];
Z = nan(numel(RS), numel(zN));
for i = 1:numel(RS)
    fm = RS(i).fm; fd = RS(i).fd;
    for k = 1:size(ZF, 1)
        j = find(strcmp(PT8.name, ZF{k, 1}), 1); d = mean(PT8.val(j, ZF{k, 2}));
        Z(i, k) = (mean(gv(fm, ZF{k, 1}, ZF{k, 2})) - d) / (abs(d) * PT8.s(j));
    end
    j = find(strcmp(KT8.name, 'ktrP_k80'), 1);
    Z(i, size(ZF, 1) + 1) = (gv(fm, 'ktrP_k80', 1) - KT8.val(j)) / (abs(KT8.val(j)) * KT8.s(j));
    if isfield(fm, 'FV_normpowerAvg') && isfield(fd, 'FV_normpowerAvg')
        dv = fd.FV_normpowerAvg(:); sd = max(sqrt(fd.FV_normpowerVar(:)), 0.136 * abs(dv));
        mv = fm.FV_normpowerAvg(:); n = min([numel(mv), numel(dv), 4]);
        Z(i, size(ZF, 1) + 1 + (1:n)) = (mv(1:n) - dv(1:n)) ./ sd(1:n);
    end
end
iMech = 1:size(ZF, 1) + 1; iFV = size(ZF, 1) + 1 + (1:4);
out('\nDistance from the data in data-SD units (z); chi2 = sum z^2 per block\n');
out('%-10s %s | %8s %8s\n', 'variant', sprintf('%9s', zN{:}), 'chi2mech', 'chi2FV');
for i = 1:numel(RS)
    out('%-10s %s | %8.2f %8.2f\n', lbl(RS(i)), sprintf('%9.2f', Z(i, :)), sum(Z(i, iMech).^2, 'omitnan'), sum(Z(i, iFV).^2, 'omitnan'));
end
if ~isempty(iref)
    out('\nchange vs ref: dchi2 mechanics / F-V\n');
    for i = 1:numel(RS)
        out('%-10s SRX_ss %.3f  dchi2mech %+7.2f  dchi2FV %+7.2f\n', lbl(RS(i)), gv(RS(i).fm, 'SRX_ss', 1), ...
            sum(Z(i, iMech).^2, 'omitnan') - sum(Z(iref, iMech).^2, 'omitnan'), ...
            sum(Z(i, iFV).^2, 'omitnan') - sum(Z(iref, iFV).^2, 'omitnan'));
    end
end
if isfield(RS(iref).fm, 'FV_normpowerAvg')
    out('\nF-V normalised power at the 4 velocities (data | ref):  %s | %s\n', ...
        mat2str(RS(iref).fd.FV_normpowerAvg(:)', 3), mat2str(RS(iref).fm.FV_normpowerAvg(:)', 3));
    if isfield(RS(iref).fd, 'FV_powerv'); out('velocities (ML/s): %s\n', mat2str(RS(iref).fd.FV_powerv(:)', 3)); end
    for i = 1:numel(RS); out('  %-10s %s\n', lbl(RS(i)), mat2str(RS(i).fm.FV_normpowerAvg(:)', 3)); end
end
chi2 = struct('mech', num2cell(sum(Z(:, iMech).^2, 2, 'omitnan')'), 'FV', num2cell(sum(Z(:, iFV).^2, 2, 'omitnan')'));

%% chi2 vs SRX figure (data-SD units; mechanics and F-V apart)
fg = figure('Color', 'w', 'Position', [40 40 1400 480], 'Visible', 'off');
tiledlayout(fg, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; box on; grid on;
col = struct('L', [0.85 0.33 0.10], 'S', [0.00 0.45 0.74], 'M', [0.47 0.67 0.19], 'MA', [0.49 0.18 0.56]);
hl = []; lg = {};
for fam = {'L', 'S'}
    k = [iref, find(strcmp({RS.fam}, fam{1}))];
    x = arrayfun(@(r) gv(r.fm, 'SRX_ss', 1), RS(k)); [x, o] = sort(x); k = k(o);
    hl(end+1) = plot(x, [chi2(k).mech], '-o', 'Color', col.(fam{1}), 'MarkerFaceColor', col.(fam{1}), 'LineWidth', 1.8); %#ok<SAGROW>
    plot(x, [chi2(k).FV], '--s', 'Color', col.(fam{1}), 'LineWidth', 1.2);
    lg{end+1} = fam{1}; %#ok<SAGROW>
end
plot(gv(RS(iref).fm, 'SRX_ss', 1), chi2(iref).mech, 'kp', 'MarkerSize', 16, 'MarkerFaceColor', 'y');
set(gca, 'YScale', 'log'); xlabel('SRX fraction before the first slack (SRX\_ss)'); ylabel('\chi^2 in data-SD units');
legend(hl, {'L: recruitment \times\lambda', 'S: load sensitivity / s'}, 'Location', 'northwest');
title('solid: slack/ktr mechanics (8 features)   dashed: F-V power (4 points)');
nexttile; hold on; box on; grid on;
for fam = {'M', 'MA'}
    k = [iref, find(strcmp({RS.fam}, fam{1}))]; if numel(k) < 2; continue; end
    x = [RS(k).x]; [x, o] = sort(x); k = k(o);
    plot(x, [chi2(k).mech], '-o', 'Color', col.(fam{1}), 'MarkerFaceColor', col.(fam{1}), 'LineWidth', 1.8);
    plot(x, [chi2(k).FV], '--s', 'Color', col.(fam{1}), 'LineWidth', 1.2);
end
set(gca, 'XScale', 'log', 'YScale', 'log'); xline(1.1e-4, 'b--', 'mant-like');
xlabel('\mu : scale on SRX exchange rates'); ylabel('\chi^2 in data-SD units');
title('M: entry/escape only (pool shrinks)   MA: all SRX rates (pool kept)');
exportgraphics(fg, fullfile(rsS, 'family_cost.png'), 'Resolution', 110);

J = struct('variant', arrayfun(lbl, RS, 'uni', 0), 'ka_factor', num2cell([RS.ka]/p0.ka), ...
    'SRX_ss', arrayfun(@(r) gv(r.fm, 'SRX_ss', 1), RS, 'uni', 0), 'cost', {RS.total}, 'data_cost', {RS.dataTotal}, ...
    'chi2_mech', {chi2.mech}, 'chi2_FV', {chi2.FV});
fj = fopen(fullfile(rsS, 'family_summary.json'), 'w'); fprintf(fj, '%s', jsonencode(J, 'PrettyPrint', true)); fclose(fj);

function v = iGet(fm, nm, k)
    v = NaN(size(k));
    if isstruct(fm) && isfield(fm, nm); x = fm.(nm); k = k(k <= numel(x)); v = x(k); if isempty(v); v = NaN; end; end
end
