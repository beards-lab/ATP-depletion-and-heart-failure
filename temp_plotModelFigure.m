load dan/pcagraphingenv;
% load dan/relaxgraphingenv;
f = figure(18);
tic
aspect = 2;
% normal size of 2-col figure on page is 7.2 inches
% matlab's pixel is 1/96 of an inch
f.Position = [300 200 7.2*96 7.2*96/aspect];

clf;
% colors = lines(max(rampSet)+1); colors(1:end-1, :) = colors(2:end, :);
colors = gray(5);
% colors(1:end-1, :) = colors(2:end, :);
fs = 12;
tiledlayout(3, 4, 'TileSpacing', 'compact')
% prepare in advance so that it wont draw over my inset
sp = nexttile(1, [2 4]);hold on;

ym = 0;
for j = max(rampSet):-1:1
    if isempty(Force{j})
        continue;
    end
    hd{j} = errorbar(datatables{j}.Time-2,datatables{j}.F,datatables{j}.SD, '-', LineWidth=2, Color=colors(3, :), CapSize=0);
    set([hd{j}.Bar, hd{j}.Line], 'ColorType', 'truecoloralpha', 'ColorData', [hd{j}.Line.ColorData(1:3); 255*0.4])
    hm{j} = semilogx(Time{j},Force{j},'-', 'linewidth',2, 'Color', 'k'); 
    hph{j} = plot(nan, nan, 'x',Color=[1 1 1]); % just a placeholder
    % set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [h.Cap.EdgeColorData(1:3); 255*alpha])
    ym = max([ym Force{j}, datatables{j}.F']);
    set(sp, 'TickLength', [0.0125 0.05]);
    set(sp, 'TickLabelInterpreter', 'latex');
    ylabel('$T$ (kPa)', Interpreter='latex')
    
    if rds(j) > 1
        tp = [rds(j)*1.4, max(Force{j})*1.2];
    else
        tp = [rds(j)*0.95, max(Force{j})*0.95];
    end
    text(tp(1), tp(2), ...
        sprintf('$t_r$ = %g s', rds(j)),...
        'horizontalAlignment', 'right', VerticalAlignment='bottom', Interpreter='latex');    
end
xlim([1e-2, 160])
ylim([0 ceil((ym)/5)*5])
yl = ylim();
sp.XScale='log';
% hl = legend([hph{1} hph{2} hph{3} hph{4} hd{4} hd{3} hd{2} hd{1} hm{4} hm{3} hm{2} hm{1}], ...
%     't_r 0.1 s:', '  t_r 1 s:',' t_r 10 s:','t_r 100 s:',...
%     'Data', 'Data', 'Data', 'Data', ...    
%     'Model', 'Model', 'Model', 'Model',...
%     NumColumns=3);
hl = legend([hd{1} hm{1}], ...
    'Data',...    
    'Model',...
    NumColumns=3);

hl.ItemTokenSize = [30, 20];
hl.Box = 'off';
% title('Model fit', Interpreter='latex')

    % 'MData Ramp-up 0.1 s', 'MData Ramp-up 0.1 s', 'MData Ramp-up 0.1 s', 'MData Ramp-up 0.1 s'...
    % )

% nexttile;
% semilogx(PeakData(:, 1), PeakData(:, 2), 'ko', LineWidth=2);hold on;
% semilogx(PeakData(:, 1), PeakModel, 'x', 'MarkerEdgeColor', [1 1 1]*0.5, LineWidth=2, MarkerSize=8);
% axis([1e-1 1e2 0 ym])
% semilogx(PeakData(:, 1), PeakModel, '--', Color=[1 1 1]*0.5, LineWidth=1);
% hl = legend('Data', 'Model', 'Location', 'best')
% title(hl, 'Peaks')
% ylabel('Tension (kPa)', Interpreter='latex')

for j = max(rampSet):-1:1
    sp = nexttile;hold on;
    h = errorbar(datatables{j}.Time-2,datatables{j}.F,datatables{j}.SD, '-', LineWidth=2, Color=colors(3, :), CapSize=0);
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255])
    plot(Time{j},Force{j},'-', 'linewidth',2, 'Color', 'k'); 
    xlim([0 min(rds(j)*3, 160)]);
    xl = xlim;
    ylim(yl)
    yticks([0 20 40]);
    xticks([0 rds(j), rds(j)*2])
    if j == 4
        ylabel('$T$ (kPa)', Interpreter='latex')
    else
        yticklabels([]);
    end
    xlabel('$t$ (s)', Interpreter='latex', HorizontalAlignment='left');
    set(sp, 'TickLength', [0.05 0.05]);
    set(sp, 'TickLabelInterpreter', 'latex')
    tit = text(xl(2), min(yl(2)-5, max(Force{j})), sprintf('$t_r$ = %g', rds(j)), 'Interpreter', 'latex', ...
        HorizontalAlignment='right', VerticalAlignment='bottom');
    % tit.Position = tit.Position + [0 max(Force{j}) 0];

end    
fontsize(12, 'points');
toc
% exportgraphics(f,'../Figures/FigModelFitpCa4.png','Resolution',150)
% exportgraphics(f,'Figures/FigModelFitRelaxed.png','Resolution',150)
%% Model states