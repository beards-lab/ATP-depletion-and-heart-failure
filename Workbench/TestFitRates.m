%% main.m  — Piecewise projection example
clear; clc; close all;

% ------------------------------------------------------------
% Mock data (can be replaced with real data later)
x_data = [1, 2, 3];
y_data = [1, 4, 9];
Y_top  = 100;

% Handle to model function
model_handle = @(x) MODEL(x);

% Compute top intersections
[t_vals, nodes] = find_t_intersections(x_data, y_data, Y_top, model_handle);

% Display results
disp('Top intersection X positions (t_vals):');
disp(t_vals);
disp('All nodes (X, Y):');
disp(nodes);

% Plot the final polyline
figure(1); hold on; grid on;
plot(x_data, y_data, 'bo-', 'LineWidth', 2);
plot(nodes(:,1), nodes(:,2), 'r--o', 'LineWidth', 1.5);
yline(Y_top, ':k');
xlabel('X'); ylabel('Y');
legend('Data points','Piecewise projection','Y_{top}');
title('Piecewise Projection Construction');

%% test_piecewise_smooth_variants.m
clear; clc; close all;

x = linspace(0, 10, 200);

% Example parameter sets
params1 = [0, 2, 5, 29, 200];           % Variant 1 (fixed x-grid)
params2 = [1.0, 2.0, 4.0, 7.0, 200];   % Variant 2 (variable grid)

% Evaluate
y1 = piecewise_smooth_variant(x, params1, true);
y2 = piecewise_smooth_variant(x, params2, false);

% Plot
figure(1); 
hold on; grid on;
plot(x, y1, 'r-|', 'LineWidth', 2);
title('Variant 1: Fixed X, 4 variable Y');
xlabel('x'); ylabel('y');

plot(x, y2, 'b-+', 'LineWidth', 2);
title('Variant 2: Variable Δx, 3 variable Y');
xlabel('x'); ylabel('y');

sgtitle('Smooth Monotonic Piecewise Function Variants');

function y = MODEL(x)
    % MODEL - Placeholder for a complex model
    % For now, returns quadratic values for testing
    % In a real use-case, this could be any numerical simulation.
    y = x.^3 + 2*x;
end

function [t_vals, nodes] = find_t_intersections(x_data, y_data, Y_top, model_handle)
% FIND_T_INTERSECTIONS
% Computes X positions where each segment's model output reaches Y_top
% using fzero, showing progress after each segment (line growth visualization).
%
% Inputs:
%   x_data, y_data : base datapoints
%   Y_top          : target Y value (e.g. 100)
%   model_handle   : function handle y = MODEL(x)
%
% Outputs:
%   t_vals         : X positions of intersections with Y_top
%   nodes          : all (X,Y) pairs for visualization

    N = numel(x_data);
    t_vals = nan(1, N-1);

    % Prepare figure for live updates
    figure(100); clf; hold on; grid on;
    xlabel('X'); ylabel('Y');
    title('Piecewise Line Growth');
    yline(Y_top, 'k--', 'Y_{top}');
    plot(x_data, y_data, 'bo-', 'LineWidth', 1.5);
    drawnow;

    fprintf('\n=== Starting piecewise intersection search ===\n');
    fprintf('Target Y_top = %.2f\n', Y_top);
 
    for i = 1:N-1
        fprintf('\n--- Segment %d/%d ---\n', i, N-1);

        x1 = x_data(i);
        x2 = x_data(i+1);
        y1 = model_handle(x1);
        y2 = model_handle(x2);

        fprintf('Base points: (%.3f, %.3f) → (%.3f, %.3f)\n', x1, y1, x2, y2);

        % Define root function
        f = @(x) model_handle(x) - Y_top;

        % Configure fzero options to limit model calls
        opts = optimset('Display', 'off', 'MaxFunEvals', 10);

        try
            % Try within [x1, x2], but if needed expand upward in X
            x_guess = [x1, x2];
            if model_handle(x2) < Y_top
                x_guess = [x1, x2 * 10]; % allow extension if model too low
            end

            [x_intersect, fval, exitflag, output] = fzero(f, x_guess, opts);

            if exitflag <= 0
                warning('Segment %d: fzero did not converge, using linear extrapolation.', i);
                slope = (y2 - y1) / (x2 - x1);
                x_intersect = x1 + (Y_top - y1) / slope;
            else
                fprintf('Converged: %.3f → fval = %.3e, evals = %d\n', ...
                        x_intersect, fval, output.funcCount);
            end

            t_vals(i) = x_intersect;

        catch ME
            warning('Segment %d failed (%s). Using extrapolation.', i, ME.message);
            slope = (y2 - y1) / (x2 - x1);
            t_vals(i) = x1 + (Y_top - y1) / slope;
        end

        % Build nodes incrementally
        top_nodes = [t_vals(1:i)', Y_top * ones(i, 1)];
        base_nodes = [x_data(:), y_data(:)];
        nodes = [base_nodes; top_nodes];
        nodes = sortrows(nodes, 1);

        % Update plot incrementally
        plot(nodes(:,1), nodes(:,2), 'r-', 'LineWidth', 1.2);
        drawnow;

        fprintf('Added top intersection X = %.3f, Y = %.3f\n', t_vals(i), Y_top);
    end

    fprintf('\n=== Intersection search complete ===\n');

    % Final node set
    top_nodes = [t_vals(:), Y_top * ones(length(t_vals),1)];
    base_nodes = [x_data(:), y_data(:)];
    nodes = [base_nodes; top_nodes];
    nodes = sortrows(nodes, 1, 'ascend');
    [~, idx] = unique(nodes(:,1), 'stable');
    nodes = nodes(idx,:);
end


function y = piecewise_smooth(x, params, fixedGrid)

    if fixedGrid
        % VARIANT 1: Fixed x-grid, variable y-values (monotonic, smooth)
        % params: [y1, y2, y3, y4]
    
        % Fixed X grid (can be adjusted)
        x_ctrl = linspace(0, 5, length(params)); 
        y_ctrl = params(:)';
    
    else
        % VARIANT 2: Variable equidistant x-grid, 3 variable y, smooth & monotonic
    % params: [Δx, y2, y3, y4] or we fix y1=0 for simplicity
    %
    % Equidistant x grid: [0, Δx, 2Δx, 3Δx]
    % Smooth monotonic connection with pchip
    
        dx = abs(params(1));          % enforce positive spacing
        x_ctrl = dx * (0:length(params)-1);          % equidistant grid
    
        % variable y values (first fixed to 0 for simplicity)
        y_ctrl = [0, params(2:end)];
    
    end

    % Ensure monotonic increasing y (clip or enforce)
    for i = 2:length(y_ctrl)
        if y_ctrl(i) <= y_ctrl(i-1)
            y_ctrl(i) = y_ctrl(i-1) + 1e-3;
        end
    end

    % Monotonic cubic interpolation
    y = pchip(x_ctrl, y_ctrl, x);   
end

