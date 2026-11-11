% ShowAllPlots.m
% Regenerate every plot from the current bounded fit (params/params_reseeded.m,
% iter 14) and export each open figure full-resolution to Figures/fit_plots/.
%
% Plot provenance (which code draws what):
%   fig 1      Force-velocity (top) + Slack (bottom), data vs sim
%              -> runFVExperiment.m + runSlackExperiment.m, called by RunBakersExp.m
%   fig 80085  Feature dashboard (per-segment ktr/A/peak1/vall/steady tiles +
%              "Cost culprits") -> Auxiliary/plotFeatures.m  (DriverBoundedFit call)
%   fig 80086  Physiology (input parameter) bounds, embedded by plotFeatures
%              -> plotPhysiologyDashboard.m
%   fig 80087  Physiology bounds, standalone -> plotPhysiologyDashboard.m
%              (DriverBoundedFit call)

DriverBoundedFit;     % runs the fit + draws all figures (clears workspace internally)

outDir = 'C:\home\git\ATP-depletion-and-heart-failure\Figures\fit_plots';
if ~exist(outDir, 'dir'); mkdir(outDir); end

figs = findall(0, 'Type', 'figure');
[~, ord] = sort([figs.Number]);
figs = figs(ord);

fprintf('\n=== Exporting %d figures to %s ===\n', numel(figs), outDir);
for i = 1:numel(figs)
    f = figs(i);
    set(f, 'Position', [60 60 1500 950]);   % enlarge so dense dashboards are legible
    drawnow;
    fname = sprintf('fig%05d.png', f.Number);
    try
        exportgraphics(f, fullfile(outDir, fname), 'Resolution', 150);
        fprintf('  %-14s (fig %d)\n', fname, f.Number);
    catch ME
        fprintf('  FAILED fig %d: %s\n', f.Number, ME.message);
    end
end
fprintf('Done — %d PNGs in %s\n', numel(figs), outDir);
