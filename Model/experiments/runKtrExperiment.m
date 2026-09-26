function [E_ktr, out_ktr, features_model, features_data] = runKtrExperiment(params0, Ktr_mean)
%RUNKTRTEXPERIMENT Run the ktr (rate of force redevelopment) protocol.
%
%   [E_ktr, out_ktr, features_model, features_data] = runKtrExperiment(params0, Ktr_mean)
%
%   Replicates the ktr experimental protocol: rapid shortening to 80% ML,
%   brief hold, then re-stretch to 105% ML, hold, return to 100% ML.
%   The rate of force redevelopment (Ktr) is extracted as 1/t(Frel=1-1/e).
%
%   Inputs:
%     params0   - Base parameter struct. Must include params0.modelFcn (string).
%     Ktr_mean  - Scalar (or vector) of target Ktr values (s^-1) from data.
%                 Optional: only the hardcoded-protocol fallback uses it, since
%                 the file-driven path scores against the data trace instead.
%
%   Outputs:
%     E_ktr   - Scalar cost: MSE vs the data trace (file-driven protocol), or
%               squared deviation of simulated Ktr from Ktr_mean(1) (fallback).
%     out_ktr - Output struct from evaluateModel (time-shifted so t=0 at
%               the start of force redevelopment).
%     features_model, features_data - redevelopment-rate features
%               (extractKtrAttributes). Populated only for the file-driven path
%               and only when the protocol actually ENDS on a release into an
%               isometric hold, i.e. is a ktr manoeuvre -- the FV2 ramp files
%               are driven through this same runner but end on a LENGTHENING
%               ramp, so they emit nothing. Empty structs otherwise.
%
%               These are what lets OptimizeOn='Feats' see this protocol at
%               all: without them the ktr block contributes only E_ktr, which
%               the feature objective never reads. See
%               Analyses/KtrProtocolVsSlackKtr/conclusions.md.

    if nargin < 2; Ktr_mean = []; end

    features_model = struct();
    features_data  = struct();

    modelFcn = str2func(params0.modelFcn);
    params = params0;
    params.SL0 = 2.0;

    if isfield(params, 'ktr_velocitytableonfile') && ~isempty(params.ktr_velocitytableonfile)
        datastruct    = load(['data/' params.ktr_velocitytableonfile]);
        datatable     = datastruct.datatable;
        velocitytable = datastruct.velocitytable;

        % Start at the length the protocol actually starts from, rather than
        % assuming 2.0 um. The ktr files begin at L=1.0 ML so this is a no-op
        % for them, but the FV2 ramp files begin at L=1.1 ML (SL 2.2) and the
        % fixed 2.0 offset the whole simulated trajectory by -0.2 um.
        % (col 4 is the length in ML; SL = 2 x ML.)
        params.SL0 = 2 * velocitytable(1, 4);
        params = getParams(params, params.g, true);
        params.Velocity = velocitytable(:, 2);

        [~, out_ktr] = evaluateModel(modelFcn, velocitytable(:, 1), params);

        % Normalize model force by pre-shortening isometric force, taken at the
        % onset of shortening (interpolated, not window-averaged, so the long
        % quiet warm-up's initial force-buildup transient does not bias it).
        t_pre   = velocitytable(2, 1);   % first real protocol row (after warm-start)
        Fm_ss   = interp1(out_ktr.t, out_ktr.Force, t_pre, 'linear', 'extrap');
        if Fm_ss == 0; Fm_ss = 1; end
        F_norm  = out_ktr.Force / Fm_ss;

        % MSE against normalized experimental force trace
        Fi    = interp1(out_ktr.t, F_norm, datatable(:, 1), 'linear', 'extrap');
        E_ktr = mean((datatable(:, 3) - Fi).^2) * 10;

        % Redevelopment-rate features, so OptimizeOn='Feats' can see this
        % protocol. Only for a real ktr manoeuvre: the last ramp must be a
        % RELEASE into an isometric hold. The FV2 ramp files share this runner
        % but end on a lengthening ramp, so they fall through and emit nothing.
        isKtrShaped = size(velocitytable, 1) >= 4 && velocitytable(end-2, 2) < -1 ...
                      && velocitytable(end-1, 2) == 0 ...
                      && (velocitytable(end, 1) - velocitytable(end-1, 1)) > 0.03;
        if isKtrShaped
            if isfield(params, 'FeatureReadGrid') && ischar(params.FeatureReadGrid) ...
                    && strcmpi(params.FeatureReadGrid, 'data')
                % Read the model at the DATA's own time stamps (opt-in): the raw
                % output is solver steps or a 0.5 ms grid, and extractKtrAttributes'
                % smoother is sized in samples. Later chunk wins where segments
                % overlap. Analyses/FeatureRobustness.
                km = makeMonotonous(out_ktr.t(:)); tk = out_ktr.t(:); tk = tk(km);
                Fk = F_norm(:); Fk = Fk(km);
                tq = datatable(datatable(:, 1) >= tk(1) & datatable(:, 1) <= tk(end), 1);
                features_model = extractKtrAttributes(tq, interp1(tk, Fk, tq, 'linear'), velocitytable);
            else
                features_model = extractKtrAttributes(out_ktr.t, F_norm, velocitytable);
            end
            features_data  = extractKtrAttributes(datatable(:, 1), datatable(:, 3), velocitytable);

            % POOLED TARGET. The protocol file is ONE preparation, and 03/27 is
            % both the strongest (Fss22 80.1 kPa) and the fastest (k63 49.4
            % against a 3-day mean of 43.2). Scoring against it while taking the
            % deadband from the 3-day spread mixes two populations and makes the
            % term roughly twice as permissive as intended. Under
            % UsePooledTarget, swap in the 3-preparation mean -- the same thing
            % pooledTargetFeatures does for the slack block, which deliberately
            % leaves this protocol alone. Measured by
            % Analyses/RefitStrategy/PlanRefit2.m from the three curated files:
            %   03/27 49.4 | 04/03 42.1 | 04/10 38.1  -> mean 43.19, rel SD 0.132
            % Now computed from the three curated files by pooledKtrTarget, so
            % every ktr-protocol feature (k63, the scored k80, F20) is pooled the
            % same way and cannot go stale after a re-curation.
            if isfield(params, 'UsePooledTarget') && ~isequal(params.UsePooledTarget, false)
                PK = pooledKtrTarget();
                for jk = 1:numel(PK.name)
                    features_data.(PK.name{jk}) = PK.val(jk);
                end
            end
        end

        Ktr = NaN;   % not extracted for file-driven protocol
    else
        params = getParams(params, params.g, true);
        % Fallback: hardcoded rapid ktr protocol, scalar Ktr cost
        if isempty(Ktr_mean)
            error('runKtrExperiment:missingKtrMean', ...
                ['The hardcoded ktr fallback needs Ktr_mean (run LoadData first), ' ...
                 'or set params0.ktr_velocitytableonfile to use a protocol file.']);
        end
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
            Fm_ss  = interp1(out_ktr.t, out_ktr.Force, t_pre, 'linear', 'extrap');
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
