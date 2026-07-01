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
    if params.MgATP == 8
        a = 1;
    elseif params.MgATP == 2
        a = 3;
    end
    params.MgATP = ATP_c(a);
    
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
                [F_active(a, j), out_fv] = evaluateModel(modelFcn, t_ss, params);
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
                [F_active(a, j), out_fv] = evaluateModel(modelFcn, t_sl0/abs(FV_velocities(j)), params);
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
                [F_active(a, j), out_fv] = fetchOutputs(futures{j});
                outs(j) = out_fv;
            end
        end
    end

    % cost function
    [found, idx_ATP] = ismember(-FV_velocities, Data_ATP(:,1));
    valid_idx = idx_ATP(found);

    E_fv = sum((F_active(params.EvalAtp,found)./Data_ATP(valid_idx,params.EvalAtp+1)' - 1).^2, 'all');
    % normalize by number of data points
    E_fv = E_fv / size(valid_idx, 2) / length(params.EvalAtp);

    if params.EvalFeatures
        if recalculateDataFeats
            features_data = extractForceVelocityAttributes(-Data_ATP(:,1)', Data_ATP(:,a+1)', struct(), FV_velocities);
        else
            % Hardcoded Baker lab 8 mM reference values (avoids re-extracting each iteration)
            AllVelocities = -[0, 0.5, 1, 2, 3, 4, 5, 6, 7];
            % if strcmp(params0.FV_dataset, 'Baker2022')
                % AllForces     = [56.40, 51.8120, 37.4459, 17.8025, 11.4430, 6.2643, 3.2759, 2.2120];
                % AllForces = Data_ATP(:, a + 1)';
            % elseif strcmp(params0.FV_dataset,'IsovelocityForFilip2021')
                % AllForces = [67.3942   60.9885   42.7059   18.8539 12.5956    6.4426    3.5136    1.8556];
                % AllForces = Data_ATP(:, a + 1)';
                
            % end
            AllForces = Data_ATP(:, a + 1)';
            vsel = find(ismember(AllVelocities, FV_velocities));
            features_data.FV_v     = -AllVelocities(vsel)';
            features_data.FV_f     = AllForces(vsel)';
            features_data.FV_fnorm = AllForces(vsel)'/AllForces(1);
        end
        features_model = extractForceVelocityAttributes(FV_velocities, F_active(a, :), features_model, FV_velocities);
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
        for a = params.EvalAtp
            set(gca,'ColorOrderIndex',a);
            ld = [ld plot(Data_ATP(:,a+1),Data_ATP(:,1),'o','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1])]; %#ok<AGROW>
        end
        for a = params.EvalAtp
            set(gca,'ColorOrderIndex',a);
            ls = [ls plot(F_active(a, :), -FV_velocities,'x-','linewidth',2, 'MarkerSize', 20)]; %#ok<AGROW>
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
