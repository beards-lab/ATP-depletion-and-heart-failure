% SummarizeKsr0Grid.m -- table, JSON and figure from results/ksr0grid.mat.
rootS = 'C:\home\git\ATP-depletion-and-heart-failure';
resS  = fullfile(rootS,'Analyses','LowATP_Refit9K2Sweep','results');
if ~exist('GR','var'); load(fullfile(resS,'ksr0grid.mat')); end          %#ok<LOAD>
KEYG = {'steady','peak1_y','ktrfit','ktrP_k80','rsK80','t0','rsTfitE','vall_tq','FV2'};
PT8 = pooledTarget8mM(); KT8 = pooledKtrTarget(8); KT2 = pooledKtrTarget(2);
dR = nan(1,numel(KEYG));
for i = 1:numel(KEYG)
    j = find(strcmp(PT8.name, KEYG{i}), 1); if ~isempty(j); dR(i) = PT8.r(j); end
    j = find(strcmp(KT8.name, KEYG{i}), 1); if ~isempty(j); dR(i) = KT2.val(j)/KT8.val(j); end
end
dR(end) = FV2L.r;
mv = @(st, nm) mean(double(st.(nm)(isfinite(st.(nm)))));
rat = @(R) [arrayfun(@(i) mv(R.fm, KEYG{i})/mv(fm8L, KEYG{i}), 1:numel(KEYG)-1), R.fv2/fv8L];

SR = reshape([GR.score_rel], numel(FGRID), numel(GGRID))';     % rows g, cols f
SP = reshape([GR.score_plain], numel(FGRID), numel(GGRID))';
[~, ib] = min([GR.score_rel]); B = GR(ib);
i0 = find([GR.g] == GGRID(1) & [GR.f] == 1, 1);                % k2-only reference
fid = fopen(fullfile(resS,'ksr0grid.txt'),'w');
fprintf(fid, 'k2 gate x ksr0 factor at 2 mM (%s). Relevance-weighted 2 mM chi-2 (plain in brackets)\n\n', CORE);
fprintf(fid, '%8s', 'g \\ f'); hdr = arrayfun(@(f) sprintf('ksr0 x%.1f', f), FGRID, 'uni', 0); fprintf(fid, '%16s', hdr{:}); fprintf(fid, '\n');
for a = 1:numel(GGRID)
    fprintf(fid, '%8.2f', GGRID(a));
    for b = 1:numel(FGRID); fprintf(fid, '%9.3f (%5.2f)', SR(a,b), SP(a,b)); end
    fprintf(fid, '\n');
end
fprintf(fid, '\nbest: g = %.2f, ksr0 x%.1f  rel %.4f (k2-only at g %.2f: %.4f)\n\n', B.g, B.f, B.score_rel, GR(i0).g, GR(i0).score_rel);
r0 = rat(GR(i0)); rB = rat(B);
fprintf(fid, '%-9s %8s %8s %8s\n', 'feature', 'k2 only', 'best', 'data');
for i = 1:numel(KEYG); fprintf(fid, '%-9s %8.3f %8.3f %8.3f\n', KEYG{i}, r0(i), rB(i), dR(i)); end
fclose(fid); type(fullfile(resS,'ksr0grid.txt'));
J = struct('g', GGRID, 'f', FGRID, 'score_rel', SR, 'score_plain', SP, 'best_g', B.g, 'best_f', B.f, ...
           'best_rel', B.score_rel, 'ref_rel', GR(i0).score_rel, 'features', {KEYG}, 'ratio_ref', r0, ...
           'ratio_best', rB, 'data', dR);
fw = fopen(fullfile(resS,'ksr0grid_summary.json'),'w'); fwrite(fw, jsonencode(J, 'PrettyPrint', true)); fclose(fw);

navy = [0 39 76]/255; red = [0.75 0.2 0.12]; gold = [1 0.8 0.02];
f = figure('Color','w','Position',[50 50 1700 700],'Visible','off');
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
nexttile; imagesc(SR); colormap(gca, flipud(parula)); cb = colorbar; cb.Label.String = 'relevance-weighted 2 mM \chi^2';
set(gca, 'XTick', 1:numel(FGRID), 'XTickLabel', arrayfun(@(x) sprintf('x%.1f', x), FGRID, 'uni', 0), ...
    'YTick', 1:numel(GGRID), 'YTickLabel', arrayfun(@(x) sprintf('%.2f', x), GGRID, 'uni', 0), 'FontSize', 14, 'TickLabelInterpreter', 'none');
xlabel('ksr0 factor at 2 mM (SRX entry)'); ylabel('g = R_2(2 mM)/R_2(8 mM)');
for a = 1:numel(GGRID); for b = 1:numel(FGRID)
    text(b, a, sprintf('%.2f', SR(a,b)), 'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold', 'Color', 'k');
end; end
title('k2 gate \times SRX-entry gate', 'FontSize', 16);
nexttile; hold on; box on; grid on;
xb = 1:numel(KEYG);
bar(xb - 0.27, dR, 0.26, 'FaceColor', navy); bar(xb, r0, 0.26, 'FaceColor', [0.6 0.6 0.6]); bar(xb + 0.27, rB, 0.26, 'FaceColor', gold);
yline(1, 'r--', 'LineWidth', 1.5);
set(gca, 'XTick', xb, 'XTickLabel', KEYG, 'FontSize', 13, 'TickLabelInterpreter', 'none'); xtickangle(35);
ylabel('2 mM / 8 mM ratio');
legend({'data', sprintf('k2 only (g %.2f)', GR(i0).g), sprintf('best (g %.2f, ksr0 x%.1f)', B.g, B.f), 'no effect'}, 'Location', 'northwest', 'FontSize', 12);
title('Ratios: k2 only vs k2 + SRX entry', 'FontSize', 16);
exportgraphics(f, fullfile(resS,'ksr0grid.png'), 'Resolution', 200); close(f);
