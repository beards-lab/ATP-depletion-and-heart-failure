%% PlotDrugLevers.m -- dose-response of each lever (force) and its shape signature (dTTP, dRT50), both cores.
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');
cores = {'t1_isoT', 's3_free'};  cols = [0 0.45 0.75; 0.85 0.33 0.1];
T = [];
for c = 1:numel(cores)
    f = fullfile(res, ['levers_' cores{c} '.mat']);  if exist(f, 'file'), S = load(f);  T = [T; S.T]; end
end
lev = unique(T.lever(~strcmp(T.lever, 'none')), 'stable');
% observed drug effects (drug - matched control, % points) to overlay as horizontal bands
obs = struct('name', {'mava 0.2/0.3/1 uM', 'iso 0.1-1 uM', 'diltiazem 3/5 uM', 'clozapine 7/10 uM 24 h', 'sunitinib 6/10 uM', 'dox 3/5 uM 48 h', 'carfilzomib 0.2-0.4 uM 48 h'}, ...
             'v', {[-33 -71 -91], [17 16 15], [-26 -31], [-21 -52], [-25 -52], [-30 -86], [-65 -77 -79]});
fig = figure('Color', 'w', 'Position', [20 20 1600 1000], 'Visible', 'off');
tl = tiledlayout(fig, 3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(lev)
    ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');
    for c = 1:numel(cores)
        m = strcmp(T.lever, lev{i}) & strcmp(T.core, cores{c});  if ~any(m), continue; end
        [x, o] = sort(T.f(m));  y = T.dF(m);  y = y(o);
        plot(ax, x, y, 'o-', 'Color', cols(c, :), 'LineWidth', 1.4, 'MarkerFaceColor', cols(c, :), 'DisplayName', cores{c});
    end
    yline(ax, 0, ':k');
    if ~strcmp(lev{i}, 'MgATP'), set(ax, 'XScale', 'log'); xline(ax, 1, ':k'); end
    title(ax, strrep(lev{i}, '_', '\_'), 'FontWeight', 'normal');  xlabel(ax, 'factor (MgATP: mM)');  ylabel(ax, '\DeltaF_{peak} (%)');
    ylim(ax, [-100 50]);
    if i == 1, legend(ax, 'Location', 'southeast'); end
end
% shape signature panel: dRT50 vs dF for every lever (realistic core), to show which levers are separable by kinetics
ax = nexttile(tl, [1 2]);  hold(ax, 'on');  grid(ax, 'on');
mk = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', 'x', '+', '*', '.'};
m0 = strcmp(T.core, 't1_isoT');
for i = 1:numel(lev)
    m = m0 & strcmp(T.lever, lev{i});  if ~any(m), continue; end
    [x, o] = sort(T.dF(m));  y = T.dRT50(m);  y = y(o);  z = T.dTTP(m);  z = z(o);
    plot(ax, x, y, ['-' mk{i}], 'LineWidth', 1.2, 'MarkerSize', 6, 'DisplayName', lev{i});
end
xline(ax, 0, ':k');  yline(ax, 0, ':k');
xlabel(ax, '\DeltaF_{peak} (%)');  ylabel(ax, '\DeltaRT50 (ms)');  title(ax, 'shape signature (t1\_isoT): relaxation change per lever', 'FontWeight', 'normal');
legend(ax, 'Location', 'eastoutside', 'NumColumns', 2);
exportgraphics(fig, fullfile(res, 'levers_doseresponse.png'), 'Resolution', 120);  close(fig);
fprintf('figure written: %s\n', fullfile(res, 'levers_doseresponse.png'));
