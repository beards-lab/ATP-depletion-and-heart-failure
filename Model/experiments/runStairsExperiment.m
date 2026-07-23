function [E_stairs, out_stairs] = runStairsExperiment(params0)
%RUNSTAIRSEXPERIMENT Run the ramp-up (stairs) protocol and return error + outputs.
%
%   [E_stairs, out_stairs] = runStairsExperiment(params0)
%
%   If params.stairs_velocitytableonfile is set (non-empty), loads a protocol
%   file (velocitytable + datatable) and drives the model with it, computing
%   MSE against normalized experimental force. Falls back to Baker ramp-up
%   data (bakers_rampup8.mat) when the field is empty or absent.
%
%   Inputs:
%     params0  - Base parameter struct. Must include params0.modelFcn (string).
%
%   Outputs:
%     E_stairs   - Scalar cost (MSE * 10).
%     out_stairs - Output struct from evaluateModel.

    modelFcn = str2func(params0.modelFcn);
    params = params0;

    if isfield(params, 'stairs_velocitytableonfile') && ~isempty(params.stairs_velocitytableonfile)
        % ── File-driven protocol ───────────────────────────────────────────────
        datastruct    = load(['data/' params.stairs_velocitytableonfile]);
        datatable     = datastruct.datatable;
        velocitytable = datastruct.velocitytable;

        params.SL0    = 2.0;
        params.Slim_l = 1.4;   % accommodate L=0.8 (SL=1.6 µm) and L=1.1 (SL=2.2 µm)
        params.Slim_r = 2.3;
        params = getParams(params, params.g, true);
        params.Velocity = velocitytable(:, 2);

        [~, out_stairs] = evaluateModel(modelFcn, velocitytable(:, 1), params);

        % Fss: model isometric force at the onset of the first length change.
        % Interpolated (not window-averaged) so it is robust to the sparse solver
        % output during the long quiet warm-start hold. Works for staircase
        % stretches and ktr-style shortening alike.
        i_move = find(velocitytable(:, 2) ~= 0, 1);
        t_move = velocitytable(i_move, 1);
        Fm_ss  = interp1(out_stairs.t, out_stairs.Force, t_move, 'linear', 'extrap');
        if Fm_ss == 0; Fm_ss = 1; end
        F_norm = out_stairs.Force / Fm_ss;

        Fi      = interp1(out_stairs.t, F_norm, datatable(:, 1), 'linear', 'extrap');
        E_stairs = mean((datatable(:, 3) - Fi).^2) * 10;

        if ~isempty(params.ghostSave)
            ghost = [out_stairs.t; F_norm']';
            save(['Ghost_' params.ghostSave '_rampup'], 'ghost');
        end

        if params.PlotEachSeparately
            if ~params.PlotFullscreen
                nexttile;
            end
            hold on;
            yyaxis right;
            plot(datatable(:, 1), datatable(:, 2), '-', out_stairs.t, out_stairs.SL, 'o-', 'Linewidth', 2, 'MarkerSize', 3);
            yyaxis left;

            if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_rampup.mat'],'file')
                ghost = load(['Ghost_' params.ghostLoad '_rampup']); ghost = ghost.ghost;
                gp = plot(ghost(:,1), ghost(:,2), '-', 'Linewidth', 3, 'Color', [0.5843 0.8157 0.9882]);
            else
                clear gp;
            end

            plot(datatable(:, 1), datatable(:, 3), 'k-', 'linewidth', 1);
            plot(out_stairs.t, F_norm, 'b-', 'linewidth', 1.5);
            xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
            ylabel('Force (norm.)','interpreter','latex','fontsize',16);
            title(sprintf('Stairs (file): E=%.3f', E_stairs));
            xlim([datatable(1,1), datatable(end,1)]);

            if exist('gp', 'var') && isvalid(gp)
                legend(['Ghost ' params.ghostLoad], 'F data', 'F sim', 'SL data*', 'SL sim*', 'Location', 'Best');
            else
                legend('F data', 'F sim', 'SL data*', 'SL sim*', 'Location', 'Best');
            end
        end

    else
        % ── Fallback: Baker ramp-up data ──────────────────────────────────────
        datastruct    = load('data/bakers_rampup8.mat');
        datatable     = datastruct.datatable;
        velocitytable = datastruct.velocitytable;
        velocitytable(1, 1) = -1;
        params.Velocity = velocitytable(1:end-1, 2);

        params.SL0    = 2.0;
        params.Slim_l = 1.9;
        params.Slim_r = 2.3;
        params = getParams(params, params.g, true);

        [~, out_stairs] = evaluateModel(modelFcn, velocitytable(:, 1), params);

        Fi       = interp1(out_stairs.t, out_stairs.Force, datatable(:, 1));
        e        = (datatable(:, 3) - Fi).^2;
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

            if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_rampup.mat'],'file')
                ghost = load(['Ghost_' params.ghostLoad '_rampup']); ghost = ghost.ghost;
                gp = plot(ghost(:,1), ghost(:,2), '-', 'Linewidth', 3, 'Color', [0.5843 0.8157 0.9882]);
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
end
