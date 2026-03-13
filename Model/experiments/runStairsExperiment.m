function [E_stairs, out_stairs] = runStairsExperiment(params0)
%RUNSTAIRSEXPERIMENT Run the ramp-up (stairs) protocol and return error + outputs.
%
%   [E_stairs, out_stairs] = runStairsExperiment(params0)
%
%   Loads Baker ramp-up data, simulates the prescribed velocity protocol, and
%   computes a mean-squared error between simulated and experimental force
%   traces (scaled by 10 to match the magnitude of other error terms).
%
%   Inputs:
%     params0  - Base parameter struct. Must include params0.modelFcn (string).
%
%   Outputs:
%     E_stairs   - Scalar cost for the ramp-up protocol (MSE * 10).
%     out_stairs - Output struct from evaluateModel.

    modelFcn = str2func(params0.modelFcn);
    params = params0;
    datastruct = load('data/bakers_rampup8.mat');
    datatable = datastruct.datatable;
    velocitytable = datastruct.velocitytable;
    velocitytable(1, 1) = -1; % enough time to get to steady state
    params.Velocity = velocitytable(1:end-1, 2);

    params.SL0 = 2.0;
    params.Slim_l = 1.9;
    params.Slim_r = 2.3;
    % update params with new N and Slims
    params = getParams(params, params.g, true);

    [~, out_stairs] = evaluateModel(modelFcn, velocitytable(:, 1), params);

    Fi = interp1(out_stairs.t, out_stairs.Force, datatable(:, 1));
    e = (datatable(:, 3) - Fi).^2;
    E_stairs = mean(e) * 10;

    if ~isempty(params.ghostSave)
        ghost = [out_stairs.t; out_stairs.Force]';
        save(['Ghost_' params.ghostSave '_rampup'], 'ghost');
    end

    if params.PlotEachSeparately
        if ~params.PlotFullscreen
            nexttile;
        end

        hold on;
        yyaxis right;
        plot(datatable(:, 1), datatable(:, 2), out_stairs.t, out_stairs.SL);
        yyaxis left;

        % manage GHOST
        if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_rampup.mat'],'file')
            ghost = load(['Ghost_' params.ghostLoad '_rampup']);
            ghost = ghost.ghost;
            gp = plot(ghost(:, 1), ghost(:, 2), '-', 'Linewidth', 3, 'Color', [0.5843    0.8157    0.9882]);
        else
            clear gp;
        end

        plot(datatable(:, 1), datatable(:, 3), 'k-', 'linewidth', 1);
        plot(out_stairs.t, out_stairs.Force, 'b-', 'linewidth', 1.5);
        xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
        ylabel('Force (rel.)','interpreter','latex','fontsize',16);
        set(gca,'fontsize',14, 'xlim', [-0.05 0.35]); box on;
        title('Ramp-up');

        if exist('gp', 'var') && isvalid(gp)
            legend(['Ghost ' params.ghostLoad], 'F data', 'F sim', 'SL data*', 'SL sim*', 'Location', 'Best');
        else
            legend('F data', 'F sim', 'SL data*', 'SL sim*', 'Location', 'Best');
        end
    end
end
