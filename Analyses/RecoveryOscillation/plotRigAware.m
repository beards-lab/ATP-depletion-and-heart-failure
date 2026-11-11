function plotRigAware(C, H, resDir, outstem)
%PLOTRIGAWARE  One page: the four descents, the ranking, and what moved.
%
%   OUTSTEM names the file, so two campaigns do not overwrite each other's
%   figure. Defaults to the original name for older callers.

if nargin < 4 || isempty(outstem); outstem = 'rigaware_report'; end

fig = figure('Position',[60 60 1400 620],'Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on
for i = 1:numel(H)
    plot(0:numel(H(i).cost), [C(1).cost H(i).cost], '-o', ...
         'LineWidth',1.6, 'MarkerSize',4);
end
yline(C(1).cost,'k--');
xlabel('improved round'); ylabel('objective');
if ~isempty(H); legend([{H.tag}], 'Box','off','Interpreter','none','Location','northeast'); end
title('the four descents');

nexttile; hold on; box on
nmC = {C.name}; v = [C.cost];
barh(1:numel(v), v); set(gca,'YTick',1:numel(v),'YTickLabel',nmC,'TickLabelInterpreter','none');
xline(C(1).cost,'k--'); xlabel('re-evaluated objective');
title('re-evaluated from disk, not from the header');
set(gca,'YDir','reverse');

nexttile; hold on; box on
if numel(C) > 1
    [~, iB] = min([C(2:end).cost]); iB = iB + 1;
    d = C(iB).featCost - C(1).featCost;
    keep = abs(d) > 1e-4;
    [ds, o] = sort(d(keep), 'ascend');
    nm = C(1).featNames(keep); nm = nm(o);
    b = barh(1:numel(ds), ds); b.FaceColor = 'flat';
    b.CData(ds > 0, :) = repmat([.85 .33 .10], nnz(ds > 0), 1);
    b.CData(ds <= 0, :) = repmat([0 .45 .74], nnz(ds <= 0), 1);
    set(gca,'YTick',1:numel(ds),'YTickLabel',nm,'TickLabelInterpreter','none','YDir','reverse');
    xline(0,'k-'); xlabel('cost change (blue = improved)');
    title(sprintf('%s: where the drop came from', C(iB).name), 'Interpreter','none');
end

exportgraphics(fig, fullfile(resDir,[outstem '.png']), 'Resolution',140);
close(fig);
end
