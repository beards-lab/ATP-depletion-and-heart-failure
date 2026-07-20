function [E_ktr, out_ktr] = runKtrExperiment(params0, Ktr_mean)
%RUNKTRTEXPERIMENT Run the ktr (rate of force redevelopment) protocol.
%
%   [E_ktr, out_ktr] = runKtrExperiment(params0, Ktr_mean)
%
%   Replicates the ktr experimental protocol: rapid shortening to 80% ML,
%   brief hold, then re-stretch to 105% ML, hold, return to 100% ML.
%   The rate of force redevelopment (Ktr) is extracted as 1/t(Frel=1-1/e).
%
%   Inputs:
%     params0   - Base parameter struct. Must include params0.modelFcn (string).
%     Ktr_mean  - Scalar (or vector) of target Ktr values (s^-1) from data.
%
%   Outputs:
%     E_ktr   - Scalar cost: squared deviation of simulated Ktr from Ktr_mean(1).
%     out_ktr - Output struct from evaluateModel (time-shifted so t=0 at
%               the start of force redevelopment).

    modelFcn = str2func(params0.modelFcn);
    params = params0;
    params.SL0 = 2.0;
    params = getParams(params, params.g, true);

    if isfield(params, 'ktr_velocitytableonfile') && ~isempty(params.ktr_velocitytableonfile)
        datastruct    = load(['data/' params.ktr_velocitytableonfile]);
        datatable     = datastruct.datatable;
        velocitytable = datastruct.velocitytable;
        params.Velocity = velocitytable(:, 2);

        [~, out_ktr] = evaluateModel(modelFcn, velocitytable(:, 1), params);

        % Normalize model force by pre-shortening isometric
        t_pre   = velocitytable(2, 1);   % first real protocol row (after warm-start)
        Fm_ss   = mean(out_ktr.Force(out_ktr.t < t_pre));
        if Fm_ss == 0; Fm_ss = 1; end
        F_norm  = out_ktr.Force / Fm_ss;

        % MSE against normalized experimental force trace
        Fi    = interp1(out_ktr.t, F_norm, datatable(:, 1), 'linear', 'extrap');
        E_ktr = mean((datatable(:, 3) - Fi).^2) * 10;

        Ktr = NaN;   % not extracted for file-driven protocol
    else
        % Fallback: hardcoded rapid ktr protocol, scalar Ktr cost
        v = 500; % ML/s
        pos_ML = [1   , 1, 0.8, 0.8, 1.05, 1.05, 1, 1];
        times = cumsum([0, 1.0004, 0.2/v, 0.01005 - 0.2/v, 0.25/v, 0.0045 - 0.25/v, 0.05/v, 1]);
        params.Velocity = diff(pos_ML) ./ diff(times);
        [~, out_ktr] = evaluateModel(modelFcn, times - times(end-1) + 1, params);
        out_ktr.t = out_ktr.t - 1;

        i_0 = find(out_ktr.t > 0 & out_ktr.FXB > 0, 1);
        Frel = out_ktr.FXB(i_0:end) ./ out_ktr.FXB(end);
        i_ktr = find(Frel >= 1 - exp(-1), 1);
        if ~isempty(i_ktr)
            Ktr   = 1 / out_ktr.t(i_ktr + i_0);
            E_ktr = abs(Ktr - Ktr_mean(1)).^2;
        else
            Ktr   = Inf;
            E_ktr = Inf;
        end
    end

    if ~isempty(params.ghostSave)
        ghost = [out_ktr.t; out_ktr.Force / out_ktr.Force(end)]';
        save(['Ghost_' params.ghostSave '_ktr'], 'ghost');
    end

    if params.PlotEachSeparately
        if params.PlotFullscreen
            clf;
        else
            nexttile;
        end
        hold on;

        if isfield(params, 'ktr_velocitytableonfile') && ~isempty(params.ktr_velocitytableonfile)
            % File-driven protocol: plot normalized force traces
            t_pre  = velocitytable(2, 1);
            Fm_ss  = mean(out_ktr.Force(out_ktr.t < t_pre));
            if Fm_ss == 0; Fm_ss = 1; end

            yyaxis right;
            plot(datatable(:, 1), datatable(:, 2), '-', out_ktr.t, out_ktr.SL, 'o-', 'Linewidth', 2, 'MarkerSize', 3);
            yyaxis left;

            if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_ktr.mat'],'file')
                ghost = load(['Ghost_' params.ghostLoad '_ktr']);
                ghost = ghost.ghost;
                gp = plot(ghost(:, 1), ghost(:, 2), '-', 'Linewidth', 3, 'Color', [0.5843    0.8157    0.9882]);
            else
                clear gp;
            end

            plot(datatable(:, 1), datatable(:, 3), 'k-', 'linewidth', 1);
            plot(out_ktr.t, out_ktr.Force / Fm_ss, 'b-', 'linewidth', 1.5);
            xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
            ylabel('Force (norm.)','interpreter','latex','fontsize',16);
            title(sprintf('Ktr protocol (file): E=%.3f', E_ktr));
            xlim([datatable(1,1), datatable(end,1)]);
        else
            % Fallback: bakers data, scalar Ktr display
            datastruct_plot = load('data/bakers_ktr_8.mat');
            datatable_plot  = datastruct_plot.datatable;
            yyaxis right;
            plot(datatable_plot(:, 1), datatable_plot(:, 2), '-', out_ktr.t, out_ktr.SL, 'o-', out_ktr.t, out_ktr.LXB, ':', 'Linewidth', 2, 'MarkerSize', 3);
            yyaxis left;

            if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_ktr.mat'],'file')
                ghost = load(['Ghost_' params.ghostLoad '_ktr']);
                ghost = ghost.ghost;
                gp = plot(ghost(:, 1), ghost(:, 2), '-', 'Linewidth', 3, 'Color', [0.5843    0.8157    0.9882]);
            else
                clear gp;
            end

            plot(datatable_plot(:, 1), datatable_plot(:, 3), 'k-', 'linewidth', 1);
            plot(out_ktr.t, out_ktr.Force, 'b-', 'linewidth', 1.5);
            xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
            ylabel('Force (rel.)','interpreter','latex','fontsize',16);
            title(sprintf('Speed of the transient: %1.1f s^{-1}', Ktr));
            xlim([-0.02, 0.1]);

            if exist('gp', 'var') && isvalid(gp)
                legend(['Ghost ' params.ghostLoad], 'F data', 'F sim','SL data*', 'SL sim*', 'LXB sim*', 'Location', 'southeast');
            else
                legend('F data', 'F sim','SL data*', 'SL sim*', 'LXB sim*', 'Location', 'southeast');
            end
        end
    end
end
