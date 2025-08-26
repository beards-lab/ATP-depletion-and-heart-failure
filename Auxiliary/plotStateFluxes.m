function plotStateFluxes(out, t, max_flux)

if t < out.t(1) || t > out.t(end)
    return;
end

i_pos = find(out.t >= t, 1, 'first');

% Define state positions for rectangular arrangement
x = [0, 1, 1, 0, 0, 1];
y = [1, 1, 0, 0, 2, 2];

% Define forward fluxes between states
ST2UT = out.RSR2PT(i_pos);
SD2UD = out.RSRD2PD(i_pos);
UT2UD = out.RTD(i_pos);  % Previously d2t, then s12s2
UD2A1 = out.RD1(i_pos);  % Previously t2a1, then s22a1
A12A2 = out.R12(i_pos);
A22UT = out.R2T(i_pos);  % Previously a22d, then a22s1

% Define backward fluxes between states
UT2ST = out.RPT2SR(i_pos);
UD2SD = out.RPD2SRD(i_pos);
UD2UT = 0;  % Previously t2d, then s22s1
A12UD = out.R1D(i_pos);  % Previously a12t, then a12s2
A22A1 = out.R21(i_pos);
UT2A2 = out.RT2(i_pos);

% Define flux between SD and ST
SD2ST = out.RSRD2SR(i_pos);

% Assemble forward and backward fluxes into arrays
forward_fluxes = [UT2UD, UD2A1, A12A2, A22UT, ST2UT, SD2UD];
backward_fluxes = [UD2UT, A12UD, A22A1, UT2ST, UD2SD, SD2ST, UT2A2];

% Define the maximal flux
if nargin < 3
    max_flux = 50;
end

% Set scaling factor for arrow thickness
% scaling_factor = 1;

% Plot states
hold on;
plot(x, y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10);

% Add state labels
labels = {'UT', 'UD', 'A1', 'A2', 'ST', 'SD'};
text(x + 0.1, y +0.1, labels, 'FontSize', 12, 'FontWeight', 'bold');

% Mapping forward fluxes to their corresponding start and end indices
forward_pairs = [1 2; 2 3; 3 4; 4 1; 5 1; 6 2];

% Mapping backward fluxes to their corresponding start and end indices
backward_pairs = [2 1; 3 2; 4 3; 1 5; 2 6; 6 5; 1 4];

% Offset for parallel arrows
vertical_offset = 0.02;
horizontal_offset = 0.02;

% Function to center arrows between states
centered_arrow = @(x1, y1, x2, y2, flux, offset_x, offset_y, color, thickness) quiver(...
        (x1 + x2) / 2 + offset_x, ...
        (y1 + y2) / 2 + offset_y, ...
        (x2 - x1) * (flux / max_flux) / 2, ...
        (y2 - y1) * (flux / max_flux) / 2, ...
        0, 'MaxHeadSize', 1.5, 'LineWidth', thickness, 'Color', color);

% Plotting forward fluxes
for i_pos = 1:size(forward_pairs, 1)
    si = forward_pairs(i_pos, 1);
    ei = forward_pairs(i_pos, 2);
    
    % Get starting and ending coordinates
    x_start = x(si);
    y_start = y(si);
    x_end = x(ei);
    y_end = y(ei);
    
    % Determine the correct offset
    if x_start == x_end
        offset_x = 0;
        offset_y = vertical_offset;
    else
        offset_x = horizontal_offset;
        offset_y = 0;
    end
    
    % Draw centered arrows for forward fluxes
    if forward_fluxes(i_pos) > max_flux        
        centered_arrow(x_start, y_start, x_end, y_end, max_flux, offset_x, offset_y, 'b', 1 + min(4, (forward_fluxes(i_pos) - max_flux)));
    else
        centered_arrow(x_start, y_start, x_end, y_end, forward_fluxes(i_pos), offset_x, offset_y, 'b', 1);
    end
end

% Plotting backward fluxes
for i_pos = 1:size(backward_pairs, 1)
    si = backward_pairs(i_pos, 1);
    ei = backward_pairs(i_pos, 2);

    % Get starting and ending coordinates
    x_start = x(si);
    y_start = y(si);
    x_end = x(ei);
    y_end = y(ei);

    % Determine the correct offset
    if x_start == x_end
        offset_x = 0;
        offset_y = -vertical_offset;
    else
        offset_x = -horizontal_offset;
        offset_y = 0;
    end
    
    % Draw centered arrows for forward fluxes
    if backward_fluxes(i_pos) > max_flux
        centered_arrow(x_start, y_start, x_end, y_end, max_flux, offset_x, offset_y, 'r', 1 + min(4, (backward_fluxes(i_pos) - max_flux)));
    else
        centered_arrow(x_start, y_start, x_end, y_end, backward_fluxes(i_pos), offset_x, offset_y, 'r', 1);
    end
    
end

% Set plot properties
% axis equal;
title('Fluxes between States with Forward and Backward Flows');
xlabel('X');
ylabel('Y');
grid off;
hold off;