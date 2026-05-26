% Update features_captions figure: set YLim and save

close all;
fig = openfig('Figures/features_captions.fig');

% Set YLim for each axis
ax = findall(fig, 'type', 'axes');
fprintf('Processing %d axes:\n', length(ax));

for i = 1:length(ax)
    if i == 4  % steady state - valley2 (negative values)
        ax(i).YLim = [-inf 0];
        fprintf('  Axis %d: YLim = [-inf, 0]\n', i);
    else
        ax(i).YLim = [0 inf];
        fprintf('  Axis %d: YLim = [0, inf]\n', i);
    end
end

% Apply compact spacing
tcl = fig.Children(1);
tcl.Padding = 'compact';
tcl.TileSpacing = 'compact';
fprintf('\nApplied compact spacing\n');

% Save .fig
savefig(fig, 'Figures/features_captions_compact.fig');
fprintf('Saved: Figures/features_captions_compact.fig\n');

% Save PNG at 300 DPI
exportgraphics(fig, 'Figures/features_captions_compact.png', 'Resolution', 300);
fprintf('Saved: Figures/features_captions_compact.png\n');

% Show file info
fig_info = dir('Figures/features_captions_compact.fig');
png_info = dir('Figures/features_captions_compact.png');
fprintf('\nFile sizes:\n');
fprintf('  .fig: %.1f KB\n', fig_info.bytes / 1024);
fprintf('  .png: %.1f MB\n', png_info.bytes / (1024*1024));

fprintf('\nDone!\n');
