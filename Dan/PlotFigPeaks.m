% Plots comparison of peak heights for all ramps and for intermediate
% calcium levels

pCas = [4.4, 5.5, 5.75, 6, 11];

pCaPeaksModel = nan(length(pCas), 4);
pCaPeaksData = nan(length(pCas), 4);

for i_pCa = 1:length(pCas)
    pCa = pCas(i_pCa);
    model_data = load(sprintf('..\\pca%gmodeldata.mat', pCa));
    
    for j_rds = 1:length(model_data.Farr)
        if isempty(model_data.Farr{j_rds})
            continue;
        end
        
        pCaPeaksModel(i_pCa, j_rds) = max(model_data.Farr{j_rds});
        if isinf(pCa) || pCa >= 10
          %   % hack - the no-Ca noPNB experiments had higher ramps
          %   datatable = readtable(['..\Data\bakers_passiveStretch_' num2str(rds(i_rd)*1000) 'ms.csv']);
          %   datatable.Properties.VariableNames = {'Time'  'ML'  'F'  'SL'};
          %   datatables{i_rd} = datatable;
          % elseif isnan(pCa)
              % newest format of experiments    
            datatables{j_rds} = readtable(['..\Data\AvgRelaxed_' num2str(rds(j_rds)) 's.csv']);
        else
            datatables{j_rds} = readtable(['..\Data\AvgpCa' num2str(pCa) '_' num2str(rds(j_rds)) 's.csv']);
        end
        [m, i_m] = max(datatables{j_rds}.F);
        pCaPeaksData(i_pCa, j_rds) = m;
        pCaPeaksData_SD(i_pCa, j_rds) = datatables{j_rds}.SD(i_m);
    end
end
%%
f = figure(232);clf;
aspect = 2;
f.Position = [300 200 7.2*96 7.2*96/aspect];
tiledlayout('flow');
lw = 1;
ms = 6;
tl = nexttile;fontsize(12, 'points');hold on;

% for i = [1 5]
i = 1;
ple = errorbar(rds, pCaPeaksData(i, :), pCaPeaksData_SD(i, :), 'ok', LineWidth=lw);
pl1 = semilogx(rds, pCaPeaksModel(i, :), 'ks-.', 'MarkerFaceColor','k', LineWidth=lw, MarkerSize = ms)
i = 5;
errorbar(rds, pCaPeaksData(i, :), pCaPeaksData_SD(i, :), 'ok', LineWidth=lw);
pl2 = semilogx(rds, pCaPeaksModel(i, :), 'ks--', 'MarkerFaceColor','k', LineWidth=lw, MarkerSize = ms)

% end
tl.XScale = 'log';
set(gca, 'Xtick',fliplr(rds));
ylabel('$\Theta$ (kPa)', Interpreter='latex');
ylim([0, 50]);
xlabel('$t$ (s)', Interpreter='latex');
% legend([ple, pl1, pl2], 'Data', 'pCa 4.4 (model)', 'Relaxed (model)')
legend([pl1, pl2], 'pCa 4.4 (model)', 'Relaxed (model)')

nexttile;fontsize(12, 'points');hold on;
ylabel('$\Theta$ (kPa)', Interpreter='latex');
xlabel('pCa', Interpreter='latex')
ylim([0, 50]);
fontsize(12, 'points')
pCasP = pCas;
pCasP(end) = 8;

for i_rds = [4]
    ple = errorbar(-pCasP, pCaPeaksData(:, i_rds),pCaPeaksData_SD(:, i_rds), 'ok', LineWidth=lw);
    plot(-pCasP(1:end-1), pCaPeaksModel(1:end-1, i_rds), 'ks:', 'MarkerFaceColor','k',LineWidth=2, MarkerSize = ms)
    pl1 = plot(-pCasP(end), pCaPeaksModel(end), 'ks:', 'MarkerFaceColor','k',LineWidth=2, MarkerSize = ms);
end
legend([ple pl1], 'Data', 'Ramp 0.1s (model)', 'Location', 'southwest')
set(gca, 'XTick',-[8:-1:4]);
set(gca, 'Xticklabels', ["11" "" "6" "5" "4"]);
xlim([-8.5 -4])

exportgraphics(f,'../Figures/FigModelPeaks.png','Resolution',150)
exportgraphics(f,'../Figures/FigModelPeaks.eps')

