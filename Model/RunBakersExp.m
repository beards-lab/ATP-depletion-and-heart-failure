E = [];
features_model = struct();
features_data = struct();
% modelFcn = @dPUdTCaSimpleAlternative2State;
modelFcn = @dPUdT_CombinedTransitions;
recalculateDataFeats = false;

if isfield(params0, 'modelFcn')
    modelFcn = str2func(params0.modelFcn);
end

%% FORCE VELOCITY

if params0.RunForceVelocity
    params = params0;
    t_ss = [0 0.1];
    t_sl0 = [0 0.1];

    params = getParams(params, params.g, true, true);
    params.UseForceOnsetShift = false;
    F_active = [];
    clear outs;

    % params.UseTitinModel = false;
    % params.UseSerialStiffness = false;
    if isfield(params, 'PU0')
        params = rmfield(params, 'PU0');
    end

    if isfield(params, 'FV_velocities')
        FV_velocities = params.FV_velocities;
    end    

    % this is disabled ATM
    % for a = params.EvalAtp
    a = 1;
    params.MgATP = ATP_c(a);
    N = length(FV_velocities);
    % paralelize?
    if isempty(gcp('nocreate')) || params.justPlotStateTransitionsFlag
        parple = [];
    else
        parple = gcp('nocreate');
        % Launch asynchronous tasks
        futures = cell(1, N);
    end

    for j = 1:N
        % try
            if FV_velocities(j) == 0
                params.SL0 = 2.0;
                params.Velocity = 0;
                
                % optimize for speed, force is up to 80
                dsl = 80/params.kSE;
                params.Slim_r = 2.0 + params.A1AttachmentWidth + params.dS;
                params.Slim_l = 2.0 - dsl - params.A1AttachmentWidth - params.dS;
                
                % params = getParams(params);
                
                if isempty(parple)
                    [F_active(a, j) out] = evaluateModel(modelFcn, t_ss, params);
                    outs(j) = out;
                else
                     futures{j} = parfeval(parple, @evaluateModel, 2, modelFcn, t_ss, params);
                end
            else
                % optimize for speed, force is up to 80
                dsl = 80/params.kSE;
                params.Slim_r = 2.2 + params.A1AttachmentWidth + params.dS;
                params.Slim_l = 2.0 - dsl - params.A1AttachmentWidth - params.dS;
                
                % experimentally cut in half
                % params.Slim_l = params.Slim_r - (params.Slim_r-params.Slim_l)/2;
                % params.Slim_l = params.Slim_r - 2*params.dS*params.MaxStrainArraySize;

                params.SL0 = 2.2;
                % params = getParams(params);
                % true to start from 2.2um steady state isntead from scratch.
                % Neither is perfect though
                if ~isfield(params, 'PU0') && params.OptimizeFVInit
                    %% speed things up by storing the initialization
                    % tic
                    params.Velocity = 0;
                    [~, out] = evaluateModel(modelFcn, t_ss, params);
                    params.PU0 = out.PU(end, :);
                    % toc
                end
                params.Velocity = FV_velocities(j);
                if isempty(parple)
                    %%
                    % tic
                    % rmfield(params, 'PU0');
                    % params.Velocity = vel(j);
                    [F_active(a, j) out] = evaluateModel(modelFcn, t_sl0/abs(FV_velocities(j)), params);                        
                    % toc
                    % F_active
                    outs(j) = out;
                    % StatesInTime
                else
                     futures{j} = parfeval(parple, @evaluateModel, 2, modelFcn, t_sl0/abs(FV_velocities(j)), params);
                end
                if abs(FV_velocities(j)) >= 3
                    breakpointIsHappening = 1; % only to place a bp
                end
            end
        % catch e 
            % wraps additional error for the optimizer
            % handleAndRethrowCostException(e, t_ss*(length(vel) - j));
        % end
        
    end
    if ~isempty(parple)
        for j = 1:N 
            if ~isempty(futures{j})
                [F_active(a, j), out] = fetchOutputs(futures{j}); 
                outs(j) = out;
            end
        end
    end

    out.Slim_r = params.Slim_r;
    out.Slim_l = params.Slim_l;

    % end

    % cost function
    [found, idx_ATP] = ismember(-FV_velocities, Data_ATP(:,1));
    valid_idx = idx_ATP(found);

    % E(1) = sum((F_active(params.EvalAtp,:) - Data_ATP(:,params.EvalAtp+1)').^2, 'all');
    E(1) = sum((F_active(params.EvalAtp,found)./Data_ATP(valid_idx,params.EvalAtp+1)' - 1).^2, 'all');
    % normalize by number of data points
    E(1) = E(1)/size(valid_idx, 2)/length(params.EvalAtp);
%%    
    if params.EvalFeatures
        %%
        
        if recalculateDataFeats
            features_data = extractForceVelocityAttributes(-Data_ATP(:,1)', Data_ATP(:,a+1)', features_data, FV_velocities);
        else
            % features_data = extractForceVelocityAttributes(-Data_ATP(:,1)', Data_ATP(:,a+1)', features_data, -[0.5, 1, 2, 3, 4, 5, 6, 7])
            AllVelocities = -[0, 0.5, 1, 2, 3, 4, 5, 6, 7];
            AllForces = [56.40, 51.8120, 37.4459, 17.8025, 11.4430, 6.2643, 3.2759, 2.2120];
            vsel = find(ismember(AllVelocities, FV_velocities));
            features_data.FV_v = -AllVelocities(vsel)';
            features_data.FV_f = AllForces(vsel)';
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
        if sum(E) < e0
            ss.params = params; % save struct
            ss.E = sum(E);
            save([params.ghostSave '_params.mat'], 'params', 'E');
            better = true;
        end
    end

    %
    % figure(101); clf; axes('position',[0.1 0.6 0.35 0.35]); hold on;
    % figure(); clf; axes('position',[0.1 0.6 0.35 0.35]); hold on;

    if params.PlotEachSeparately || better
        %     figure(102);hold on;
        if ~params.PlotFullscreen
            % axes('position',[0.05 0.6 0.4 0.35]);
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

        %%
        ls = [];ld = []; % line legend for sim and for data
        for a = params.EvalAtp
            set(gca,'ColorOrderIndex',a);
            ld = [ld plot(Data_ATP(:,a+1),Data_ATP(:,1),'o','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1])];
        end
        for a = params.EvalAtp
            set(gca,'ColorOrderIndex',a);
            ls = [ls plot(F_active(a, :), -FV_velocities,'x-','linewidth',2, MarkerSize=20)];
        end
        % legend('8mM', '4mM', '2mM');
        ylabel('Velocity (ML/s)','interpreter','latex','fontsize',16);
        xlabel('Force (kPa)','interpreter','latex','fontsize',16);
        set(gca,'fontsize',14);
        % axis([0 65 0 6]);
        axis([-10 65 0 6]);
        title('Force-velocity')
        % presentation stuff
        set(gca,'fontsize',16);xlim([0 70])

        box on;grid on;


        % plot(Data_ATP(:,3),Data_ATP(:,1),'go','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
        % plot(Data_ATP(:,4),Data_ATP(:,1),'ro','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);

        if exist('gp', 'var') && isvalid(gp)
            legend(['Ghost ' params.ghostLoad], ['Sim'  params.SimTitle] , 'Data', 'interpreter','none');
        else
            %         legend(['Sim'  params.SimTitle] , 'Data');
            % legend([ld ls(1)], '8mM ATP', '4mM ATP', '2mM ATP', 'model');
        end
    end

    if ~isempty(params.ghostSave)
        if size(FV_velocities, 1) > size(FV_velocities, 2) 
            ghost = [F_active(1, :), -FV_velocities(:)'];
        else
            ghost = [F_active(1, :)', -FV_velocities(:)];
        end
        save(['Ghost_' params.ghostSave '_FV'], 'ghost', "params0");
    end
end
% return;
%% KTR EXPERIMENT
if params0.RunKtr
    params = params0;
    params.SL0 = 2.0;
    % params.dS = 0.008;
    % params.Slim_l = 1.8;
    % params.Slim_r = 2.0;
    % params.UseSlack = true;
    % params.PlotFullscreen = true;
    % params.LXBpivot = 2.0;
    params = getParams(params, params.g, true, true);

    % params.UseTitinModel = false;

    % replicate the ktr protocol
    v = 500; % ML/s
    %       times  = [-1e3, 0,  2, 20, 22.5, 25 , 25.5, 1e3]/1000 - 25.5e-3;
    pos_ML = [1   , 1,0.8,0.8, 1.05,1.05,     1, 1 ];

    % putting the numbers as a difference
    times = cumsum([0 , 1.0004,0.2/v,0.01005 - 0.2/v,0.25/v, 0.0045 - 0.25/v, 0.05/v, 1]);
    params.Velocity = diff(pos_ML)./diff(times);
    [t, out] = evaluateModel(modelFcn, times - times(end-1) + 1, params);
    out.t = out.t - 1;

    % calculate ktr
    i_0 = find(out.t > 0 & out.FXB > 0, 1);
    Frel = out.FXB(i_0:end)./out.FXB(end);
    % relative to max value from data to run optimizer faster
    % Frel = out.FXB(i_0:end)./55;
    i_ktr = find(Frel >= 1-exp(-1), 1);
    if ~isempty(i_ktr) 
        Ktr = 1/out.t(i_ktr+i_0);
        E(2) = abs(Ktr-Ktr_mean(1)).^2;
    else
        Ktr = Inf;
        E(2) = Inf;
    end

    if ~isempty(params.ghostSave)
        ghost = [out.t;out.Force/out.Force(end)]';
        save(['Ghost_' params.ghostSave '_ktr'], 'ghost');
    end
    if params.PlotEachSeparately
        %
        % clf;
        if params.PlotFullscreen
            clf;
        else
            % axes('position',[0.55 0.6 0.4 0.35]);
            nexttile;
        end
        hold on;
        datastruct = load('data/bakers_ktr_8.mat');
        datatable = datastruct.datatable;
        yyaxis right;
        plot(datatable(:, 1) ,datatable(:, 2), '-', out.t, out.SL, 'o-', out.t, out.LXB, ':', 'Linewidth', 2, 'MarkerSize', 3);
        yyaxis left;


        % manage GHOST
        if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_ktr.mat'],'file')
            ghost = load(['Ghost_' params.ghostLoad '_ktr']);
            ghost = ghost.ghost;
            gp = plot(ghost(:, 1), ghost(:, 2), '-', 'Linewidth', 3, 'Color', [0.5843    0.8157    0.9882]);
        else
            clear gp;
        end
        scaleData = 1/datatable(end, 3);
        scaleData = 1;

        plot(datatable(:, 1),datatable(:, 3)*scaleData,'k-','linewidth',1);

        scaleModel = 1/out.Force(end);
        scaleModel = 1;
        plot(out.t,out.Force*scaleModel,'b-','linewidth',1.5);
        % plot(Tspan,F_active,'linewidth',1.5);
        xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
        ylabel('Force (rel.)','interpreter','latex','fontsize',16);
        % set(gca,'fontsize',14,'ylim',[0 1.1], 'xlim', [-0.05 0.45]);  box on;
        title(sprintf('Speed of the transient: %1.1f s^{-1}', Ktr));
        xlim([-0.02, 0.1]);
        % xlim([1, 1.0008])

        if exist('gp', 'var') && isvalid(gp)
            legend(['Ghost ' params.ghostLoad], 'F data', 'F sim','SL data*', 'SL sim*', 'LXB sim*', 'Location', 'southeast');
        else
            legend('F data', 'F sim','SL data*', 'SL sim*', 'LXB sim*', 'Location', 'southeast');
        end
    end
end
%% RAMP UP
if params0.RunStairs
    params = params0;
    datastruct = load('data/bakers_rampup8.mat');
    datatable = datastruct.datatable;
    velocitytable = datastruct.velocitytable;
    velocitytable(1, 1) = -1; % enough time to get to steady state
    params.Velocity = velocitytable(1:end-1, 2);

    params.SL0 = 2.0;
    % params.LXBpivot = 2.0;
    params.Slim_l = 1.9;
    params.Slim_r = 2.3;
    % params.dS
    % params.WindowsOverflowStepCount
    % params.N = 30;
    % update params with new N and Slims
    params = getParams(params, params.g, true, true);

    [F out] = evaluateModel(modelFcn, velocitytable(:, 1), params);

    Fi = interp1(out.t, out.Force, datatable(:, 1));
    e = (datatable(:, 3) - Fi).^2;
    E(3) = mean(e)*10;

    if ~isempty(params.ghostSave)
        ghost = [out.t; out.Force]';
        save(['Ghost_' params.ghostSave '_rampup'], 'ghost');
    end
    if params.PlotEachSeparately
        if ~params.PlotFullscreen
            % axes('position',[0.05 0.1 0.4 0.35]);
            nexttile;
        end

        hold on;
        yyaxis right;
        plot(datatable(:, 1),datatable(:, 2), out.t, out.SL);
        yyaxis left;

        % manage GHOST
        if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_rampup.mat'],'file')
            ghost = load(['Ghost_' params.ghostLoad '_rampup']);
            ghost = ghost.ghost;
            gp = plot(ghost(:, 1), ghost(:, 2), '-', 'Linewidth', 3, 'Color', [0.5843    0.8157    0.9882]);
        else
            clear gp;
        end

        plot(datatable(:, 1),datatable(:, 3),'k-','linewidth',1);
        plot(out.t,out.Force,'b-','linewidth',1.5);
        % plot(Tspan,F_active,'linewidth',1.5);
        xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
        ylabel('Force (rel.)','interpreter','latex','fontsize',16);
        set(gca,'fontsize',14, 'xlim', [-0.05 0.35]);  box on;
        title('Ramp-up');


        if exist('gp', 'var') && isvalid(gp)
            legend(['Ghost ' params.ghostLoad],'F data', 'F sim','SL data*', 'SL sim*', 'Location', 'Best');
        else
            legend('F data', 'F sim','SL data*', 'SL sim*', 'Location', 'Best');
        end
    end
end
%% SLACK
if params0.RunSlack
    params = params0;

    if isfield(params, 'velocitytableonfile')
    
    % datastruct = load("data\bakers_slack8mM_all_fix.mat");
    % datastruct = load("data\bakers_slack8mM_all_fudge.mat");
        datastruct = load(['data/' params.velocitytableonfile ]);
    else
        datastruct = load('data/bakers_slack8mM_all.mat');
    end
    datatable = datastruct.datatable;
    validZone = datatable(:, 1) > 1;
    
    
    datastruct.velocitytable(1, 1) = -20;
    PU0 = [];
    par_velocitytable = [];
    parple = [];

    % show the indexes
    % [(1:length(datastruct.velocitytable))' datastruct.velocitytable]

    switch params.RunSlackSegments
        % first slack
        case 'First'
    velocitytable = datastruct.velocitytable(1:7, :);
    % validZone = datatable(:, 1) > 1;
    validZone = datatable(:, 1) > datastruct.velocitytable(2, 1) & datatable(:, 1) < datastruct.velocitytable(5, 1);
        case 'FirstTwo'
    % two slacks
    velocitytable = datastruct.velocitytable(1:11, :);
    validZone = datatable(:, 1) > 1;
        case 'Fourth-rampuponly'
            params.SL0 = 1.9194;
            % only the last slack
            velocitytable = [datastruct.velocitytable(1, :); datastruct.velocitytable(16:19, :)];
            % ramp0up and peak
            validZone = datatable(:, 1) > datastruct.velocitytable(17, 1) & datatable(:, 1) < datastruct.velocitytable(19, 1) - 0.1;
            % % only ramp-up at first - 
            % validZone = datatable(:, 1) > datastruct.velocitytable(16, 1) & datatable(:, 1) < datastruct.velocitytable(17, 1);            
        case 'Fourth'
            params.SL0 = 2.2;
            % only the pre-last slack
            velocitytable = [datastruct.velocitytable(1, :); datastruct.velocitytable(15:19, :)];
            % ramp0up and peak
            validZone = datatable(:, 1) > datastruct.velocitytable(15, 1) - 0.1 & datatable(:, 1) < datastruct.velocitytable(19, 1) - 0.1;

        case 'FirstAndLast'
            velocitytable = datastruct.velocitytable([1:6 19:23], :);
            % validZone = (datatable(:, 1) > 1 & datatable(:, 1) < 1.45) | (datatable(:, 1) > 2.7);
            validZone = datatable(:, 1) > datastruct.velocitytable(2, 1) & datatable(:, 1) < datastruct.velocitytable(5, 1) | datatable(:, 1) > datastruct.velocitytable(19, 1)-.1 & datatable(:, 1) < datastruct.velocitytable(21, 1);
        case 'FirstAndLastShort'
            velocitytable = datastruct.velocitytable([1:6 19:23], :);
            velocitytable([5:6 9:10], 1) = velocitytable([5:6, 9:10], 1) + 1;
            % validZone = (datatable(:, 1) > 1 & datatable(:, 1) < 1.45) | (datatable(:, 1) > 2.7);
            validZone = datatable(:, 1) > datastruct.velocitytable(3, 1) - 0.1 & datatable(:, 1) < datastruct.velocitytable(3, 1)+0.04 ...
                | datatable(:, 1) > datastruct.velocitytable(19, 1)-.1 & datatable(:, 1) < datastruct.velocitytable(19, 1) + 0.04;            
            
        case 'FirstAndLastExtended'
            velocitytable = datastruct.velocitytable([1:6 19:23], :);
            velocitytable([5:6 9:10], 1) = velocitytable([5:6, 9:10], 1) + 1;
            % validZone = (datatable(:, 1) > 1 & datatable(:, 1) < 1.45) | (datatable(:, 1) > 2.7);
            validZone = datatable(:, 1) > datastruct.velocitytable(2, 1) & datatable(:, 1) < datastruct.velocitytable(5, 1) ...
                | datatable(:, 1) > datastruct.velocitytable(18, 1) & datatable(:, 1) < datastruct.velocitytable(21, 1);            
        case 'AllButLast'
            % all but the last
            velocitytable = datastruct.velocitytable(1:19, :);
            validZone = datatable(:, 1) > 1;
        case 'Last'
            % only the last slack
            velocitytable = datastruct.velocitytable(18:end, :);
            validZone = datatable(:, 1) > datastruct.velocitytable(19, 1)-.1 & datatable(:, 1) < datastruct.velocitytable(23, 1);
            velocitytable(1, 1) = -2;
            % validZone = datatable(:, 1) > 2;
        case 'All'
            % all
            velocitytable = datastruct.velocitytable(1:end, :);
            validZone = datatable(:, 1) > 1;
            validZone = datatable(:, 1) > datastruct.velocitytable(2, 1);
        case 'AllPar'
            % all, but run all in parallel
            velocitytable = datastruct.velocitytable(1:end, :);
            validZone = datatable(:, 1) > 1;
            validZone = datatable(:, 1) > datastruct.velocitytable(2, 1);
            
            % 1. precalculate the init to be reused
            if isfield(params, 'PU0')
                params = rmfield(params, 'PU0');
            end        
            % reset the PU0
            % params = getParams(params, params.g, true);
            params.Velocity = 0;
            [~, out] = evaluateModel(modelFcn, [-10 velocitytable(3, 1)], params);
            PU0 = out.PU(end, :);

            wins = 3:4:(size(velocitytable, 1) -1);
            % velocitytable(wins, :)
            
            % Preallocate cell array
            par_velocitytable = cell(length(wins), 1);
            
            % Fill the cells
            for k = 1:length(wins)
                par_velocitytable{k} = velocitytable(wins(k):wins(k)+4, :);
            end

            if isempty(gcp('nocreate')) || params.justPlotStateTransitionsFlag
                parple = [];
            else
                parple = gcp('nocreate');
                % Launch asynchronous tasks
                futures = cell(1, N);
            end
            
        case 'AllNoBump'
            % all, but evaluate only force onset, not the restretch
            velocitytable = datastruct.velocitytable(1:end, :);
            validZone = datatable(:, 1) > 1;
            validZone = datatable(:, 1) > datastruct.velocitytable(2, 1) & datatable(:, 1) < datastruct.velocitytable(5, 1) |...
            datatable(:, 1) > datastruct.velocitytable(7, 1) - 0.05 & datatable(:, 1) < datastruct.velocitytable(9, 1) |...
            datatable(:, 1) > datastruct.velocitytable(11, 1) - 0.05 & datatable(:, 1) < datastruct.velocitytable(13, 1) |...
            datatable(:, 1) > datastruct.velocitytable(15, 1) - 0.05 & datatable(:, 1) < datastruct.velocitytable(17, 1) |...
            datatable(:, 1) > datastruct.velocitytable(19, 1)-.1 & datatable(:, 1) < datastruct.velocitytable(21, 1);
        case 'ramp-up'
            velocitytable = [-5, 0;0, 1;.15,0;1,0];
            params.SL0 = 1.9;
            validZone = datatable(:, 1) > 1;
        case 'ramp-down'
            params.SL0 = 2.2
            velocitytable = [-2, 0;0, -.01;15,0;20,0];
        case 'stairs-down'
            %%
            rampvel = -0.05;
            rampup = 0.5;
            ramphold = 1;
            params.SL0 = 2.2
            velocitytable = [-2, 0;0,0];
            for iv = 2:5
                velocitytable = [velocitytable;velocitytable(end, 1) + ramphold, rampvel/rampup; velocitytable(end, 1) + rampup+ramphold, 0];
            end
            velocitytable = [velocitytable;velocitytable(end, 1) + ramphold*2, 0];
        case 'stairs-up'
            %%
            rampvel = 0.05;
            rampup = 0.5;
            ramphold = 10;
            params.SL0 = 1.9;
            velocitytable = [-2, 0;0,0];
            for iv = 2:5
                velocitytable = [velocitytable;velocitytable(end, 1) + ramphold, rampvel/rampup; velocitytable(end, 1) + rampup+ramphold, 0];
            end
            velocitytable = [velocitytable;velocitytable(end, 1) + ramphold*2, 0];

    end     
%%
    if isempty(par_velocitytable)
        par_velocitytable = {velocitytable};
    end
    N = length(par_velocitytable);
    outs = cell(1, N);    
    params = getParams(params, params.g, true, true);    

    for i_par_chunk = 1:N

        velocitytable_chunk = par_velocitytable{i_par_chunk};
        params.Velocity = velocitytable_chunk(:, 2);
        % reset the PU0
        params.PU0 = PU0;
    
        
    
        % [F out] = evaluateModel(modelFcn, velocitytable(:, 1), params);
        % [F out] = evaluateModel(modelFcn, velocitytable(:, 1), params);
        if isempty(parple)
            [F, out] = evaluateModel(modelFcn, velocitytable_chunk(:, 1), params);
            outs{i_par_chunk} = out;
        else
            % [F, out] = evaluateModel(modelFcn, velocitytable(:, 1), params);
            % futures{i_par_chunk} = [F, out];
            futures{i_par_chunk} = parfeval(parple, @evaluateModel, 2, modelFcn, velocitytable_chunk(:, 1), params);
        end
    end
%%
    if ~isempty(parple)
        for j = 1:N 
            if ~isempty(futures{j})
                % [~, out] = futures{j}; 
                [~, out] = fetchOutputs(futures{j});                 
                outs{j} = out;
            end
        end        
    end
    out = mergeOutStructs([outs{:}]);

    % i_0 = find(datatable(:, 1) > 2.77, 1);
    % i_e = length(datatable(:, 1));
    % i_0 = find(datatable(:, 1) > 2.762, 1); % start a bit earlier
    % i_e = find(datatable(:, 1) > 2.9, 1); % not all the way in
    % i_e = length(datatable(:, 1));
    
    % limiting the first slack only
    % validZone = datatable(:, 1) > velocitytable(2, 1) - 0.1 & datatable(:, 1) < velocitytable(2, 1) ... % pre-slack steady state
    %     | datatable(:, 1) > velocitytable(3, 1) + 0.002 ... & datatable(:, 1) < velocitytable(4, 1) + 0.01 ... Redevelopment zone
    %     | datatable(:, 1) > velocitytable(5, 1) + 0.01; % after the titin transient
    % validZone = true([1 length(datatable(:, 1))]);
    
%%
    nonrepeating = makeMonotonous(out.t);
    Fi = interp1(out.t(nonrepeating), out.Force(nonrepeating), datatable(validZone, 1));
    % plot(datatable(validZone, 1), datatable(validZone, 3));hold on;
    % plot(datatable(validZone, 1), Fi, '|');

    e = (datatable(validZone, 3) - max(0, Fi)).^2;
    % e(isnan(e)) = 10;
    E(4) = mean(e(~isnan(e)))*20;


    % tet = [tet; params.N, et]

    if params.PlotEachSeparately
        % axes('position',[0.55 0.1 0.4 0.35]); hold on;
        nexttile;
        yyaxis right;
        pl1 = plot(datatable(:, 1),datatable(:, 2), '-', out.t, out.SL, '--', out.t, out.LXB, ':', 'Linewidth', 2, 'MarkerSize', 3);
        % pl2 = plot(out.t, out.SL, '-', out.t, out.LXB, ':', 'Linewidth', 2, 'MarkerSize', 3);
        
        yyaxis left;hold on;

        % manage GHOST
        if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_slack.mat'],'file')
            ghost = load(['Ghost_' params.ghostLoad '_slack']);
            ghost = ghost.ghost;
            gp = plot(ghost(:, 1), ghost(:, 2), '-', 'Linewidth', 3, 'Color', [0.5843    0.8157    0.9882]);
        else
            clear gp;
        end

        % params0.ghostSave = 'NiceFit_Compliant';

        plot(datatable(:, 1),datatable(:, 3),'k-','linewidth',0.5);
        pl5 = plot(datatable(validZone, 1),datatable(validZone, 3),'k-','linewidth',1.5);        

        pl3 = plot(out.t,out.Force,'b-','linewidth',1.5);
        pl4 = plot(out.t,out.FXBPassive,'b--','linewidth',1.5);
        
        % plot(Tspan,F_active,'linewidth',1.5);
        xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
        ylabel('Force (rel.)','interpreter','latex','fontsize',16);
        set(gca,'fontsize',14);  box on;
        title('Slack');
        xlim([par_velocitytable{1}(2, 1)-0.02 velocitytable(end, 1)])
        yl = ylim;
        % plot(out.t, out.XB_TORs, '-')
        ylim(yl)
        xl = xlim();


        if exist('gp', 'var') && isvalid(gp)
            legend([gp; pl1(1:3); pl5; pl3;pl4], ['Ghost ' params.ghostLoad],'MLx2.0', 'SL*', 'LXB*','Force (data)', 'Force (sim)', 'Passive (sim)', 'Location', 'best');
        else
            legend([pl1(1:3); pl5; pl3;pl4],'MLx2.0', 'SL*', 'LXB*','Force (data)', 'Force (sim)', 'Passive (sim)', 'Location', 'best');
        end

        if params.ShowStatePlots
            nexttile;
            if params.NumberOfStates == 2
                plot(out.t, out.p1_0, '-', out.t, out.p2_0, '-', ...
                    out.t, out.PuATP, '-',out.t, out.PuR, '-', out.t, out.SR, out.t, out.SRD, LineWidth=1.5, LineStyle='-')
                legend('A1','A2',...
                    'UT','UD', 'SR', 'SRD')
            elseif params.NumberOfStates == 3
                plot(out.t, out.p1_0, '-', out.t, out.p2_0, '-', out.t, out.p3_0, '-', ...
                    out.t, out.PuATP, '-',out.t, out.PuR, '-', out.t, out.SR, out.t, out.SRD, LineWidth=1.5, LineStyle='-')
                legend('A1','A2','A3', ...
                    'UT','UD', 'SR', 'SRD')
            end
            xlim(xl);
        end
    end
    
    if params.ShowResidualPlots
        % nexttile;
        figure(1001);
        
        % if there is a ghost, then subtrack that
        if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_slack.mat'],'file')
            ghost = load(['Ghost_' params.ghostLoad '_slack']);
            ghost = ghost.ghost;
            Gi = interp1(ghost(:, 1), ghost(:, 2), datatable(validZone, 1));
            % Fi = interp1(out.t(nonrepeating), out.Force(nonrepeating), datatable(validZone, 1));
            % gp = plot(ghost(:, 1), ghost(:, 2), '-', 'Linewidth', 3, 'Color', [0.5843    0.8157    0.9882]);        
            plot(datatable(validZone, 1), (Fi - Gi), 'linewidth', 2); hold on;
            gp = plot(datatable(validZone, 1), (Gi - datatable(validZone, 3)), 'k');
        else        
            plot(datatable(validZone, 1), (Fi - datatable(validZone, 1)));
        end
    end

    if ~isempty(params.ghostSave)
        ghost = [out.t; out.Force]';
        
        % do not save ghostSave command!
        gs = params0.ghostSave; params0.ghostSave = '';
        save(['Ghost_' params.ghostSave '_slack'], 'ghost', "params0");
        params0.ghostSave = gs;
    end

    %% SAVE FIG
    if params.PlotEachSeparately
        fig = gcf;
        % saveas(fig, ['XBBakersDataFit.png']);
        if ~isempty(params.ghostSave)
            saveas(fig, ['XBBakersDataFit_' params.ghostSave '.png']);
        end
    end


    if params.EvalFitSlackOnset    
        try
            % figure(8989);clf;
            [e_dt e_ktr] = fitSlackForceOnset(datatable, velocitytable, out.t, out.SL, out.Force, params.PlotEachSeparately & params.drawForceOnset);
            E(4) = E(4);
            E(5) = 0*5e6*e_dt; E(6) = 0*e_ktr;
        catch e
            E(5) = 1e3;
            warning('Eval slack onset failed');
        end
    end

    if params.EvalFeatures
        %%
        % figure(42);clf;hold on;        
        if recalculateDataFeats
            features_data = extractSlackAttributes(datatable(:, 1), datatable(:, 3), datatable(:, 2), velocitytable, features_data, [], false);
            %% print out to store later for speed-up
            fieldNames = fieldnames(features_data);
            for k = 1:numel(fieldNames)
                fname = fieldNames{k};
                val   = features_data.(fname);
                if ~isnumeric(val)
                    continue;
                end
                s = mat2str(round(val, 4, 'significant'));
                fprintf('features_data.%s = %s;\n', fname, s);
            end
        else
            features_data.ktr = [39.76 34.18 31.63 28.36 27.6];
            features_data.A = [68.4 65.47 58.92 52.2 43.51];
            features_data.t0 = [0.002165 0.004643 0.007818 0.01117 0.01589];
            features_data.Am = [58.18 55.29 52.1 47.33 41.98];
            features_data.SLslack = [2.04 2 1.96 1.92 1.88];
            features_data.SLdiff = [0.162 0.202 0.242 0.282 0.322];
            features_data.restretchSlopeStart = [1131 1011 1006 907 953.4];
            features_data.restretchSlopeEnd = [111.3 122.2 108.7 119.7 135.3];
            features_data.v_restretch = [4 5 6 7 3];            
            features_data.peak1_y = [77.59 76.02 74.2 71.11 56.89];
            features_data.peak1_t = [0.0035 0.003 0.0025 0.002 0.005];
            features_data.peak1_SL = [2.068 2.034 1.99 1.954 1.912];
            features_data.peak1_dSL = [0.028 0.034 0.03 0.034 0.032];
            features_data.vall_y = [68.97 66.95 65.33 62.49 54.37];
            features_data.vall_t = [0.0085 0.006 0.006 0.005 0.0095];
            features_data.peak2 = [77.22 78.47 80.01 82.42 61.85];
            features_data.steady = [76.91 76.48 75.98 75.53 56.75];
            features_data.vall2_dy = [-13.56 -13.84 -13.6 -13.08 -8.826];
            features_data.vall2_t = [0.0085 0.0055 0.008 0.0065 0.0095];
            features_data.ovrsht_dy = [0.431 1.392 1.244 1.552 0.3205];
            features_data.ovrsht_t = [0.2116 0.225 0.247 0.239 0.2676];
            features_data.XTOR = [10 10 10 10 10];
        end
        % figure(43);clf;hold on;
        features_model = extractSlackAttributes(out.t, out.Force, out.SL, velocitytable, features_model, out, params0.PlotFeatureFitting);
    end
    

end

%% MINIMAL STRETCH

if isfield(params0, 'RunMinStretch') && params0.RunMinStretch

fname = ['../data/PassiveCaSrc2/' '20241121' '/02_03_Fmax_stiff_ktr.txt'];
fmaxdata = readtable(fname, ReadVariableNames=false);

datatable = table2array([1 2.0 1].*fmaxdata(:, 1:3));    
%
    h = [2.5, 5, 7.5]*1e-3;
    dt = 1e-3;
    s = 5e-4;
    
    velocitytable = [
    -2, 0, 0, 2.0;
    s+ 0.1, h(1)/dt, 2*h(1)*dt, 2.0;
    s+ 0.1 + dt, 0, 0, 2.0 + h(1)*2;
    s+ 0.1+0.2, -h(1)/dt, -2*h(1)*dt, 2.0 + h(1)*2;
    s+ 0.1+0.2+dt, 0, 0, 2.0;

    s+ 0.5, h(2)/dt, 2*h(2)*dt, 2.0;
    s+ 0.5 + dt, 0, 0, 2.0 + h(2)*2;
    s+ 0.5+0.2, -h(2)/dt, -2*h(2)*dt, 2.0 + h(2)*2;
    s+ 0.5+0.2+dt, 0, 0, 2.0;

    s+ 0.9, h(3)/dt, 2*h(3)*dt, 3.0;
    s+ 0.9 + dt, 0, 0, 2.0 + h(3)*2;
    s+ 0.9+0.2, -h(3)/dt, -2*h(3)*dt, 2.0 + h(3)*2;
    s+ 0.9+0.2+dt, 0, 0, 2.0;
        
    1.5, 0, 0, 2.0;
    ];
    
    v = 250;
    vt_ktr =[
   -0.0146 ,-v
   -0.0146 + 0.2/v         0
   -0.0046   v
   -0.0046 + 0.25/v         0
   -0.0001, -v
   -0.0001 + 0.05/v         0
    1.0000         0    ];

    velocitytable = [velocitytable ; [vt_ktr + [1.6,0],vt_ktr(:, 2)/2, vt_ktr(:, 1)*0]];

    params = params0;

    params.Velocity = velocitytable(:, 2);

    params.SL0 = 2.0;
    params.Slim_l = 1.9;
    params.Slim_r = 2.0;   
    params = getParams(params, params.g, true);
    params.UseSuperRelaxed = false;
    params.UseSuperRelaxedADP = false;
    params.UseA2MechanicalRecocking = false;    

    [F out] = evaluateModel(modelFcn, velocitytable(:, 1), params);
%
% clf;    hold on;
%     plot(datatable(:, 1), datatable(:, 2),datatable(:, 1), datatable(:, 3)/datatable(1, 3), 'k');
%     plot(out.t, out.SL, out.t, out.Force/out.Force(end), LineWidth=2);

mask = out.t >= 0;

% 2. Create the new table with the exact column names needed
% new_table = table(out.t(mask)', out.SL(mask)', out.Force(mask)', ...
%     'VariableNames', {'x_s_', 'x_Lo_', 'x_kPa_'});

% 3. Append it to the fmaxdata cell array
%% This adds it as the last element
fmaxdataCurrent = {};
fmaxdataCurrent{1} = datatable;
fmaxdataCurrent{end+1} = [out.t(mask)', out.SL(mask)', out.Force(mask)'];
ProcessFmaxData;
% StatesInTime
end

%% FORCE VELOCITY RESIMULATION
if params0.RunForceVelocityTime


    isovelocity = load('../data/bakers_isovelocity.mat', 'datatable', 'velocitytable');
    datatable = isovelocity.datatable;
    velocitytable = isovelocity.velocitytable;

    params = params0;
    % here the isovelocity was in ML isntead of SL
    % datatable(:, 2) = 2*datatable(:, 2);
    % here the isovelocity was in ms isntead of s
    % datatable(:, 1) = datatable(:, 1)/1000;
    params.datatable = datatable;
    params.Velocity = velocitytable(:, 2);

    params.SL0 = 2.2;
    params.Slim_l = 1.5;
    params.Slim_r = 2.25;
    if isfield(params, 'PU0')
        params = rmfield(params, 'PU0');
    end

    % reset the PU0
    % params = getParams(params, params.g, true);
    % velocitytable(1, 1) = 2;
    % velocitytable = velocitytable(1:4, 1);
    [F out] = evaluateModel(modelFcn, velocitytable(:, 1), params);
    %%  postprocess the isovelocity velocities
    
    constantSpeedSegment = find(velocitytable(:, 2) > 0) - 3;
    speeds = -velocitytable(constantSpeedSegment, 2);
    
    % x = velocitytable(constantSpeedSegment(1), 1);
    smashed2_0_data = arrayfun(@(x) find(datatable(:, 2) < 2.0 & datatable(:, 1) > x, 1, 'first'), velocitytable(constantSpeedSegment, 1), UniformOutput=true);
    % we are not sure it won
    smashed2_0_model = arrayfun(@(x) find(out.SL <= 2.0 & out.t > x, 1, 'first'), velocitytable(constantSpeedSegment, 1), UniformOutput=true);
    % add a final one at speed 0
    speeds = [speeds;0];
    smashed2_0_data = [smashed2_0_data; length(datatable(:, 1))];
    smashed2_0_model = [smashed2_0_model; length(out.t)];
    %%
    clf;
    nexttile;hold on;
    plot(datatable(:, 1), datatable(:, 2)/2*100,'k-', datatable(:, 1), datatable(:, 3), 'k--'); 
    plot(datatable(smashed2_0_data, 1), datatable(smashed2_0_data, 3), 'o', LineWidth=4, MarkerSize=12);
    
    [x, i] = ismember([0.5, 1], speeds)

    
    % isValid = cellfun(@(x) ~isempty(x),  smashed2_0_model);
    % isValid = smashed2_0_model)
    plot(out.t, out.Force, out.t, out.SL/2*100);
    plot(out.t(smashed2_0_model), out.Force(smashed2_0_model), 'x', LineWidth=4, MarkerSize=12);
    plot(out.t(smashed2_0_model), out.Force(smashed2_0_model), 'x', LineWidth=4, MarkerSize=12);
    
    %%
    nexttile;hold on;
    cla;hold on;
    modelForce = out.Force(smashed2_0_model);    
    
    [~, speed_sorted] = sort(speeds)
    a = 1; % ATP 8 for now
    plot(modelForce(speed_sorted), speeds(speed_sorted), 'x-', LineWidth=2, MarkerSize=12);
    
    plot(Data_ATP(:,a+1),Data_ATP(:,1),'o','linewidth',2.5,'Markersize',12,'markerfacecolor',[1 1 1])
    plot(datatable(smashed2_0_data(speed_sorted), 3), speeds(speed_sorted), 's', LineWidth=2.5, MarkerSize=12);

    if exist('F_active', 'var')
        % we have simulated the force-constvelocity, so lets compare
        plot(F_active, -FV_velocities(1:end), '+--', LineWidth=2, MarkerSize=12);
        legend('SIM - isovelocity timecourse', 'Data - isovelocity', 'Data - isovelocity timecourse (Baker 2021)', 'SIM isovelocity steady');
    else
        legend('data', 'Force-velocity timecourse');
    end
end

%% FORCE LENGTH ESTIMATED relation
if params0.RunForceLengthEstim
    params = params0;
    params.PlotEachSeparately = true;
    gp = ones(5, 1);
    Es = simulateForceLengthEstim(ones(5, 1), params, modelFcn, params.PlotEachSeparately);
    E(end+1) = Es*1e0;
end


%% READJUST ERROR FUNCTION
if isfield(params0, 'ErrorMultiplier')
    erl = 1:min(length(E), length(params0.ErrorMultiplier));
    E(erl) = E(erl).*params0.ErrorMultiplier(erl);
end


function mono = makeMonotonous(a)
%% generated by GPT
    % Reverse the array
    arev = flip(a);
    
    % Find the longest monotonically increasing subsequence
    max_value = arev(1);
    keep_rev = false(size(arev));
    keep_rev(1) = true;
    
    for i = 2:length(arev)
        if arev(i) < max_value
            keep_rev(i) = true;
            max_value = arev(i);
        end
    end
    
    % Reverse the logical array to match the original indexing
    mono = flip(keep_rev);
    
    % % Filter the array based on the 'keep' logical indices
    % result = a(keep);
end


function merged = mergeOutStructs(outs, dim)
%MERGEOUTSTRUCTS Concatenate arrays in a struct array field-by-field.
%   merged = mergeOutStructs(outs, dim)
%   outs: 1xN struct array, all elements share the same fields
%   dim : concatenation dimension (default = 1)
    dbg = 1;

    if nargin < 2 || isempty(dim)
        dim = 1; % default: stack along rows
    end
    if ~isstruct(outs)
        error('Input must be a struct array.');
    end
    if isempty(outs)
        merged = struct(); 
        return;
    end

    % Ensure all structs share the same fields
    baseFields = fieldnames(outs(1));
    for i = 2:numel(outs)
        if ~isequal(fieldnames(outs(i)), baseFields)
            error('All structs must have identical field names.');
        end
    end

    merged = struct();
    for f = 1:numel(baseFields)
        name = baseFields{f};
        vals = {outs.(name)};  % collect values from each struct for this field

        % Decide concatenation based on class of the first element
        firstVal = vals{1};

        if iscell(firstVal)
            % For cell arrays: use bracket concatenation in the chosen dimension
            % dim==1 → vertical [vals{:}]'; dim==2 → horizontal [vals{:}]
            % fprintf('Cell: %s \n', name);
            if dim == 1
                % Ensure column: wrap each cell array in a column if needed
                merged.(name) = vertcat(vals{:});
            else
                merged.(name) = horzcat(vals{:});
            end

        elseif isnumeric(firstVal) || islogical(firstVal)
            % Numeric/logical arrays: use cat(dim, ...)
            try
                % fprintf('Numeric: %s \n', name);
                merged.(name) = cat(dim, vals{:});
            catch ME
                % Fallback: try vertical then horizontal
                try
                    merged.(name) = vertcat(vals{:});
                catch
                    merged.(name) = horzcat(vals{:});
                end
            end

        elseif isstring(firstVal) || ischar(firstVal)
            % Strings/char arrays
            if isstring(firstVal)
                % fprintf('String: %s \n', name);
                merged.(name) = vertcat(vals{:});   % strings usually stack by rows
            else
                % for char, concatenate along dim if sizes permit
                merged.(name) = cat(dim, vals{:});
            end

        elseif isstruct(firstVal)
            % Generic fallback: put pieces in a cell and flatten
            % merged.(name) = [vals{:}];
            merged.(name) = firstVal;
            % fprintf('Skipping struct %s \n', name);
            
        else
            warning('Field %s unknown while merging parallel chunks \n', name)
        end
    end
end
