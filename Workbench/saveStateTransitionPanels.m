% Save each panel from plotStateTransitions_subpanels as separate high-res PNG
% Exports directly without copying (preserves all details and formatting)
%
% Usage:
%   plotStateTransitions_subpanels;
%   saveStateTransitionPanels(params, 'output_dir/');

function saveStateTransitionPanels(params, output_dir)

if nargin < 2
    output_dir = './Figures/TransitionPanels/';
end

% Create output directory if it doesn't exist
if ~isfolder(output_dir)
    mkdir(output_dir);
end

% Get the current figure with tiledlayout
fig = gcf();
tlo = fig.Children(1);

% Check that we have a tiledlayout with 8 tiles
if ~isa(tlo, 'matlab.graphics.layout.TiledChartLayout')
    error('Current figure does not have a TiledChartLayout. Run plotStateTransitions_subpanels first.');
end

% Panel names
panel_names = {
    'Panel_1_Hydrolysis'
    'Panel_2_Attachment_P1_Detachment'
    'Panel_3_Power_Stroke'
    'Panel_4_P2_Detachment'
    'Panel_5_SRX_ATP'
    'Panel_6_SRX_ADP'
    'Panel_7_Inter_SRX'
    'Panel_8_Summary'
};

% Get all axes from tiledlayout
axes_list = findobj(tlo, 'type', 'axes');
axes_list = flip(axes_list);

fprintf('Exporting %d panels directly at 600 DPI to %s\n', length(axes_list), output_dir);

% Export directly from each axis (no copying = perfect fidelity)
for i = 1:length(axes_list)
    ax = axes_list(i);
    output_file = fullfile(output_dir, [panel_names{i} '.png']);

    % Export directly from the axis at ultra-high resolution
    exportgraphics(ax, output_file, 'Resolution', 600, 'Colorspace', 'rgb');

    fprintf('  [%d/8] Saved: %s\n', i, output_file);
end

fprintf('Done! All panels exported at 600 DPI.\n');

end
