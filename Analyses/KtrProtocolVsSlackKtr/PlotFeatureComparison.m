% PlotFeatureComparison.m -- every feature, incumbent vs fast gate, plus a
% paired bar chart of the per-feature weighted costs.
%
% SCORING. The objective is the paramset's own 24 features. The ktr protocol is
% deliberately NOT in it (RunKtr = 0, no ktrP_* in fn) so that it stays a
% held-out validation. It is still SIMULATED here (RunKtr forced to 1) purely so
% its redevelopment rate can be drawn alongside the scored features -- it appears
% as the last tile with NO cost in its title, which is how you can tell it is
% held out. It contributes nothing to either E_t.
%
% Outputs: results/features_all.png, results/features_bars.png
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');

best=getParams(loadParams('fastgate_best'),[],true,false);
best.PlotEachSeparately=0; best.ghostSave=''; best.ghostLoad='';
best.EvalFeatures=1; best.BreakOnODEUnstable=0; best.SaveBest=0;
FN_SCORED = best.fn;                      % the 24 that pay
assert(~any(contains(FN_SCORED,'ktrP_')), 'ktr protocol must stay out of the objective');
FN_PLOT   = [FN_SCORED, {'ktrP_k63'}];    % + the held-out one, for the picture only

TAU  = [0.0266487, 2e-4];
NAME = {'incumbent tau_reg = 26.6 ms', 'fast gate tau_reg = 0.2 ms'};
FM = cell(1,2); COST = cell(1,2); FD = [];
for v = 1:2
    params0 = best; params0.tau_reg = TAU(v);
    params0.RunKtr = 1;                   % simulate it; it is still unscored
    params0 = getParams(params0, params0.g, true);
    LoadData;                             %#ok<*NASGU>  ATP_c / Data_ATP / Ktr_mean
    RunBakersExp;                         % fills features_model / features_data / E
    FM{v} = features_model;
    if v == 1; FD = features_data; end
    COST{v} = evalFeatureCost(features_data, features_model, FN_SCORED);
    fprintf('%-30s Et = %.4f  (24 scored features; ktrP_k63 held out at %.1f /s)\n', ...
            NAME{v}, sum(COST{v}), features_model.ktrP_k63);
end

%% ---------- 1. every feature, three series ----------
f1 = figure(90001); clf; set(f1,'Position',[20 20 1700 950],'Color','w');
cmap = containers.Map(FN_SCORED, num2cell(COST{2}));   % titles carry the FAST gate cost
                                                       % ktrP_k63 absent -> no cost shown

% Data FIRST: plotSingleBoundaryTile takes the LAST dataset that has the field,
% so the boundary / group tiles (ktr_rmse, rsR2, XTOR, k2) must end on a MODEL
% struct or they render the DATA's values instead. Occlusion of the data markers
% is then fixed by restacking below, not by reordering.
plotMultipleFeatures({FD, FM{1}, FM{2}}, ...
    {'data', 'incumbent', 'fast gate'}, ...
    [0 0 0; 0.85 0.33 0.10; 0 0.45 0.74], ...
    {'o','x','s'}, FN_PLOT, {':','-','-'}, [false false true], cmap);

% Bring every 'data' series to the front and enlarge it: the two model variants
% sit on top of the data in most tiles, which is exactly where you need to see it.
for ax = findobj(f1,'Type','axes')'
    ln = findobj(ax,'Type','line');
    if isempty(ln); continue; end
    ln = ln(strcmp(get(ln,'DisplayName'), 'data'));
    if isempty(ln); continue; end
    set(ln,'MarkerSize',8,'LineWidth',1.2);
    ax.Children = [ln(:); setdiff(ax.Children, ln, 'stable')];
end
% mark the held-out tile
axall = findobj(f1,'Type','axes');
for ax = axall'
    if contains(ax.Title.String, 'ktrP_k63')
        ax.Title.String = [ax.Title.String '  HELD OUT'];
        ax.Title.Color  = [0.85 0.33 0.10];
    end
end

sgtitle(['data (black o) vs incumbent 26.6 ms (orange x) vs fast gate 0.2 ms (blue s)' newline ...
         'titles carry the FAST gate cost; boundary/group tiles show the fast gate only; ' ...
         'ktrP_k63 is NOT scored'], 'FontSize', 11, 'Interpreter', 'none');
exportgraphics(f1, fullfile(resdir,'features_all.png'), 'Resolution', 110);

%% ---------- 2. paired bar chart (scored features only) ----------
lab = cellfun(@(s) regexprep(s,'\|.*$',''), FN_SCORED, 'UniformOutput', false);
c1 = COST{1}(:); c2 = COST{2}(:); d = c2 - c1;
[~,ord]  = sort(max(c1,c2), 'descend');
[~,ordd] = sort(d, 'descend');

f2 = figure(90002); clf; set(f2,'Position',[40 40 1250 760],'Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on; box on
y = 1:numel(ord);
b = barh(y, [c1(ord) c2(ord)], 0.82);
b(1).FaceColor = [0.85 0.33 0.10]; b(1).EdgeColor = 'none';
b(2).FaceColor = [0 0.45 0.74];    b(2).EdgeColor = 'none';
set(gca,'YTick',y,'YTickLabel',lab(ord),'YDir','reverse','FontSize',9, ...
        'YLim',[0.5 numel(ord)+0.5],'TickLabelInterpreter','none');
xlabel('weighted feature cost'); grid on
legend({sprintf('incumbent  (E_t = %.4f)', sum(c1)), sprintf('fast gate  (E_t = %.4f)', sum(c2))}, ...
       'Location','southeast','Box','off');
title('per-feature cost, sorted by the larger of the two');

nexttile; hold on; box on
for i = 1:numel(ordd)
    col = [0.80 0.20 0.20]; if d(ordd(i)) < 0; col = [0.20 0.60 0.30]; end
    barh(i, d(ordd(i)), 0.82, 'FaceColor', col, 'EdgeColor', 'none');
end
xline(0,'k-');
set(gca,'YTick',1:numel(ordd),'YTickLabel',lab(ordd),'YDir','reverse','FontSize',9, ...
        'YLim',[0.5 numel(ordd)+0.5],'TickLabelInterpreter','none');
xlabel('\Delta cost   (fast gate - incumbent)'); grid on
title(sprintf('what the fix trades   (net %+.4f)', sum(d)));
xl = xlim; xlim([xl(1)-0.05*diff(xl) xl(2)+0.02*diff(xl)]);
text(0.03,0.56,'red   = worse with the fast gate','Units','normalized','FontSize',9,'Color',[0.80 0.20 0.20]);
text(0.03,0.52,'green = better','Units','normalized','FontSize',9,'Color',[0.20 0.60 0.30]);
text(0.03,0.44,sprintf('held out, unscored: ktrP\\_k63  %.1f -> %.1f /s  (data 49.4)', ...
     FM{1}.ktrP_k63, FM{2}.ktrP_k63),'Units','normalized','FontSize',9,'Color',[0.3 0.3 0.3]);
exportgraphics(f2, fullfile(resdir,'features_bars.png'), 'Resolution', 130);

%% ---------- console table ----------
fprintf('\n%-26s %10s %10s %10s\n','feature','incumbent','fast','delta');
for i = ordd(:)'
    fprintf('%-26s %10.4f %10.4f %+10.4f\n', lab{i}, c1(i), c2(i), d(i));
end
fprintf('%-26s %10.4f %10.4f %+10.4f\n','TOTAL Et (24 scored)', sum(c1), sum(c2), sum(d));
fprintf('%-26s %10.1f %10.1f            <- HELD OUT, data 49.4 /s\n','ktrP_k63 [/s]', ...
        FM{1}.ktrP_k63, FM{2}.ktrP_k63);
save(fullfile(resdir,'feature_comparison.mat'),'FN_SCORED','lab','c1','c2','FM','FD');
fprintf('\nwrote results/features_all.png and results/features_bars.png\n');
