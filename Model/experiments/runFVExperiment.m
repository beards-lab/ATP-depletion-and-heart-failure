function [E_fv, outs, features_model, features_data] = runFVExperiment(params0, ATP_c, Data_ATP)
%RUNFVEXPERIMENT Run the force-velocity protocol and return error + outputs.
%
%   [E_fv, outs, features_model] = runFVExperiment(params0, ATP_c, FV_velocities, Data_ATP)
%
%   Inputs:
%     params0       - Base parameter struct (from getParams). Must include
%                     params0.modelFcn (string name of the ODE function).
%     ATP_c         - Vector of ATP concentrations to evaluate (mM).
%     FV_velocities - Row vector of shortening velocities (ML/s); 0 = isometric.
%     Data_ATP      - Experimental force-velocity data matrix; column 1 is
%                     velocity (negative = shortening), columns 2+ are forces
%                     at each ATP concentration.
%
%   Outputs:
%     E_fv           - Scalar cost for the FV protocol.
%     outs           - Struct array (one element per velocity) of evaluateModel
%                      output structs.
%     features_model - Struct of simulated FV features (populated only when
%                      params0.EvalFeatures is true; otherwise empty struct).
%     features_data  - Struct of reference data FV features (populated only when
%                      params0.EvalFeatures is true; otherwise empty struct).
%                      Uses hardcoded Baker 8 mM values unless
%                      params0.recalculateDataFeats is true.

    modelFcn = str2func(params0.modelFcn);

    recalculateDataFeats = params0.recalculateDataFeats;
    features_model = struct();
    features_data  = struct();

    params = params0;
    t_ss = [0 1];
    t_sl0 = [0 0.1];

    params.UseForceOnsetShift = false;
    F_active = [];
    clear outs;

    if isfield(params, 'PU0')
        params = rmfield(params, 'PU0');
    end

    FV_velocities = params.FV_velocities;

    % this is disabled ATM
    % for a = params.EvalAtp
    % params.MgATP = ATP_c(i_atp);
    LoadData;
    if params0.MgATP == 8
        i_atp = 1;
    elseif params0.MgATP == 2
        i_atp = 3;
    end
    
    N = length(FV_velocities);

    if isempty(gcp('nocreate')) || params.justPlotStateTransitionsFlag
        parple = [];
    else
        parple = gcp('nocreate');
        futures = cell(1, N);
    end

    for j = 1:N
        if FV_velocities(j) == 0
            params.SL0 = 2.0;
            params.Velocity = 0;

            % optimize for speed, force is up to 80
            dsl = 80/params.kSE;
            params.Slim_r = 2.0 + params.A1AttachmentWidth + params.dS;
            params.Slim_l = 2.0 - dsl - params.A1AttachmentWidth - params.dS;

            if isempty(parple)
                [F_active(i_atp, j), out_fv] = evaluateModel(modelFcn, t_ss, params);
                outs(j) = out_fv;
            else
                futures{j} = parfeval(parple, @evaluateModel, 2, modelFcn, t_ss, params);
            end
        else
            % optimize for speed, force is up to 80
            dsl = 80/params.kSE;
            params.Slim_r = 2.2 + params.A1AttachmentWidth + params.dS;
            params.Slim_l = 2.0 - dsl - params.A1AttachmentWidth - params.dS;

            params.SL0 = 2.2;
            if ~isfield(params, 'PU0') && params.OptimizeFVInit
                %% speed things up by storing the initialization
                params.Velocity = 0;
                [~, out_fv] = evaluateModel(modelFcn, t_ss, params);
                params.PU0 = out_fv.PU(end, :);
                params.LXBpivot = out_fv.params.LXBpivot;
            end
            params.Velocity = FV_velocities(j);
            if isempty(parple)
                [F_active(i_atp, j), out_fv] = evaluateModel(modelFcn, t_sl0/abs(FV_velocities(j)), params);
                outs(j) = out_fv;
            else
                futures{j} = parfeval(parple, @evaluateModel, 2, modelFcn, t_sl0/abs(FV_velocities(j)), params);
            end
            if abs(FV_velocities(j)) >= 3
                breakpointIsHappening = 1; %#ok<NASGU> % only to place a bp
            end
        end
    end

    if ~isempty(parple)
        for j = 1:N
            if ~isempty(futures{j})
                [F_active(i_atp, j), out_fv] = fetchOutputs(futures{j});
                outs(j) = out_fv;
            end
        end
    end

    % we do eval for a 1 atp level only this time
    % atp_levels = params.EvalAtp;
    atp_levels = i_atp;


    % cost function
    [found, idx_ATP] = ismember(-FV_velocities, Data_ATP(:,1));
    valid_idx = idx_ATP(found);
    
    E_fv = sum((F_active(atp_levels,found)./Data_ATP(valid_idx,atp_levels+1)' - 1).^2, 'all');
    % normalize by number of data points
    E_fv = E_fv / size(valid_idx, 2) / length(atp_levels);

    if params.EvalFeatures
        if recalculateDataFeats
            features_data = extractForceVelocityAttributes(-Data_ATP(:,1)', Data_ATP(:,i_atp+1)', struct(), FV_velocities);
        else
            % Reference FV curve, taken from the LOADED data at THIS [ATP].
            % Previously these were hardcoded 8 mM numbers, which silently scored
            % the 2 mM condition against the 8 mM curve. The data for every
            % [ATP] is already in LoadData (Data_ATP columns are ATP_c = 8/4/2,
            % i_atp selects), so index it rather than transcribe it.
            % NB there are 8 forces, so only 8 velocities. The list previously
            % carried a 9th entry (-7) with no matching force, which would have
            % indexed past the end had FV_velocities ever asked for it.
            AllVelocities = -[0, 0.5, 1, 2, 3, 4, 5, 6];
            FV_SETS = { 'Baker2022',               Data_ATP
                        'IsovelocityForFilip2021', Data_ATP_AB2021 };
            iset = find(strcmp(FV_SETS(:,1), params0.FV_dataset), 1);
            if isempty(iset)
                error('runFVExperiment:unknownFVdataset', ...
                      'params0.FV_dataset = "%s" is not one of: %s', ...
                      params0.FV_dataset, strjoin(FV_SETS(:,1)', ', '));
            end
            AllForces = FV_SETS{iset,2}(1:numel(AllVelocities), i_atp+1)';
            if any(~isfinite(AllForces))
                error('runFVExperiment:noDataAtThisATP', ...
                      'dataset "%s" has no force-velocity data at %g mM MgATP.', ...
                      params0.FV_dataset, params0.MgATP);
            end

            vsel = find(ismember(AllVelocities, FV_velocities));
            features_data.FV_v     = -AllVelocities(vsel)';
            features_data.FV_f     = AllForces(vsel)';
            features_data.FV_fnorm = AllForces(vsel)'/AllForces(1);

            % ---- POWER, from BOTH datasets ------------------------------
            % The two FV curves come from different preparations/sources, so
            % their absolute forces are not comparable (isometric 56.4 vs 67.4
            % kPa) but their SHAPES are. Each is therefore normalised by its
            % OWN isometric force before power is formed, and the two are
            % averaged. The per-point variance between them is the honest
            % statement of how well-determined each velocity is, and is what
            % the 'w:FV_normpowerVar' weighting in the objective consumes.
            i0 = find(AllVelocities == 0, 1);
            if isempty(i0)
                error('runFVExperiment:noIsometric', ...
                      'FV reference set has no v=0 point to normalise by.');
            end
            % v = 0 excluded — power there is identically 0 in every dataset by
            % definition, so it can never contribute error, and including it
            % would only dilute the inverse-variance weights on the points that
            % do carry information. Must match extractForceVelocityAttributes.
            ipw  = features_data.FV_v > 0;
            vpos = features_data.FV_v(ipw);
            vsp  = vsel(ipw);
            features_data.FV_powerv = vpos;
            % Each source set is read at THIS [ATP] too. A set with no data at
            % this [ATP] (AB2021 has no 4 mM arm) drops out of the average
            % instead of poisoning it with NaN.
            P = nan(numel(vsp), size(FV_SETS,1));
            for q = 1:size(FV_SETS,1)
                Fq = FV_SETS{q,2}(1:numel(AllVelocities), i_atp+1)';
                if any(~isfinite(Fq)); continue; end
                features_data.(sprintf('FV_power%d', q))     = Fq(vsp)' .* vpos;
                features_data.(sprintf('FV_normpower%d', q)) = (Fq(vsp)'/Fq(i0)) .* vpos;
                P(:,q) = features_data.(sprintf('FV_normpower%d', q));
            end
            if all(all(~isfinite(P)))
                error('runFVExperiment:noPowerData', ...
                      'no force-velocity source set has data at %g mM MgATP.', params0.MgATP);
            end
            features_data.FV_normpowerAvg = mean(P, 2, 'omitnan');
            % omitnan so a dropped source set gives the remaining one's variance
            % (0 for a single set) rather than NaN. evalFeatureCost applies a
            % VAR_FLOOR_FRAC shrinkage floor, so a 0 here cannot become an
            % infinite inverse-variance weight.
            features_data.FV_normpowerVar = var(P, 0, 2, 'omitnan');
        end
        features_model = extractForceVelocityAttributes(FV_velocities, F_active(i_atp, :), features_model, FV_velocities);
    end

    better = false;
    if params.SaveBest
        e0 = Inf;
        if exist([params.ghostSave '_params.mat'], 'file')
            e0 = load([params.ghostSave '_params.mat']).E;
        end
        if E_fv < e0
            save([params.ghostSave '_params.mat'], 'params', 'E_fv');
            better = true;
        end
    end

    if params.PlotEachSeparately || better
        if ~params.PlotFullscreen
            nexttile;
        end
        hold on;

        if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_FV.mat'],'file')
            ghost = load(['Ghost_' params.ghostLoad '_FV']);
            ghost = ghost.ghost;
            gp = plot(ghost(:, 1), ghost(:, 2), 'x-', 'Linewidth', 4, 'Color', [0.5843    0.8157    0.9882]);
        else
            clear gp;
        end

        ls = []; ld = [];
        for ea = atp_levels
            set(gca,'ColorOrderIndex',ea);
            ld = [ld plot(Data_ATP(:,ea+1),Data_ATP(:,1),'o','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1])]; %#ok<AGROW>
        end
        for ea = atp_levels
            set(gca,'ColorOrderIndex',ea);
            ls = [ls plot(F_active(ea, :), -FV_velocities,'x-','linewidth',2, 'MarkerSize', 20)]; %#ok<AGROW>
        end
        ylabel('Velocity (ML/s)','interpreter','latex','fontsize',16);
        xlabel('Force (kPa)','interpreter','latex','fontsize',16);
        set(gca,'fontsize',14);
        axis([-10 65 0 6]);
        title('Force-velocity')
        set(gca,'fontsize',16); xlim([0 70])
        box on; grid on;

        if exist('gp', 'var') && isvalid(gp)
            legend(['Ghost ' params.ghostLoad], ['Sim'  params.SimTitle] , 'Data', 'interpreter','none');
        end
    end

    if ~isempty(params.ghostSave)
        if size(FV_velocities, 1) > size(FV_velocities, 2)
            ghost = [F_active(1, :), -FV_velocities(:)'];
        else
            ghost = [F_active(1, :)', -FV_velocities(:)];
        end
        save(['Ghost_' params.ghostSave '_FV'], 'ghost', 'params0');
    end
end
