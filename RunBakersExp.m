E = [];
modelFcn = @dPUdTCaSimpleAlternative2State;
%% FORCE VELOCITY

if params0.RunForceVelocity
    params = params0;
    params.Slim_l = 1.9;
    t_ss = [0 1];
    t_sl0 = [0 0.1];

    params = getParams(params);
    F_active = [];

    params.UseTitinModel = false;
    % params.UseSerialStiffness = false;
    if isfield(params, 'PU0')
        params = rmfield(params, 'PU0');
    end
    for a = params.EvalAtp
        params.MgATP = ATP_c(a);
        for j = 1:length(vel)
            if vel(j) == 0
                params.SL0 = 2.0;
                params.Velocity = 0;
                [F_active(a, j) out] = evaluateModel(modelFcn, t_ss, params);
            else
                params.SL0 = 2.2;
                % true to start from 2.2um steady state isntead from scratch.
                % Neither is perfect though
                if ~isfield(params, 'PU0')
                    % speed things up by storing the initialization
                    params.Velocity = 0;
                    [~, out] = evaluateModel(modelFcn, t_ss, params);
                    params.PU0 = out.PU(end, :);
                end
                params.Velocity = vel(j);
                [F_active(a, j) out] = evaluateModel(modelFcn, t_sl0/abs(vel(j)), params);
                if abs(vel(j)) >= 3
                    breakpointIsHappening = 1; % only to place a bp
                end
            end
        end
    end
    % cost function
    E(1) = sum((F_active(params.EvalAtp,:) - Data_ATP(:,params.EvalAtp+1)').^2, 'all');
    % normalize by number of data points
    E(1) = E(1)/size(Data_ATP, 1)/length(params.EvalAtp);

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
            gp = plot(ghost(:, 1), ghost(:, 2), '-', 'Linewidth', 3, 'Color', [0.5843    0.8157    0.9882]);
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
            ls = [ls plot(F_active(a, :), -vel,'-','linewidth',1)];
        end
        legend('8mM', '4mM', '2mM');
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
            legend([ld ls(1)], '8mM ATP', '4mM ATP', '2mM ATP', 'model');
        end
    end

    if ~isempty(params.ghostSave)
        ghost = [F_active(1, :), -vel(:)];
        save(['Ghost_' params.ghostSave '_FV'], 'ghost');
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
    params = getParams(params, params.g, true);

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
    Ktr = 1/out.t(i_ktr+i_0);
    E(2) = abs(Ktr-Ktr_mean(1)).^2;

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
    params = getParams(params, params.g, true);

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
    datastruct = load('data/bakers_slack8mM_all.mat');
    datatable = datastruct.datatable;
    
    % first slack
    % velocitytable = velocitytable(1:4, 1);

    % two slacks
    velocitytable = datastruct.velocitytable(1:11, :);

    % all but the last
    % velocitytable = datastruct.velocitytable(1:19, :);
    
    % only the last slack
    % velocitytable = datastruct.velocitytable(18:end, :);
    
    % all
    velocitytable = datastruct.velocitytable(1:end, :);
    
    velocitytable(1, 1) = -2;
    
    
    params.Velocity = velocitytable(:, 2);
    params.datatable = datatable;

    % params.SL0 = 2.2;
    params.Slim_l = 1.85;
    % params.Slim_r = 2.2;
    % params.LXBpivot = 2.2;
    % params.dS = 0.0025;
    
    if isfield(params, 'PU0')
        params = rmfield(params, 'PU0');
    end

    % reset the PU0
    params = getParams(params, params.g, true);
    
    % [F out] = evaluateModel(modelFcn, velocitytable(:, 1), params);
    [F out] = evaluateModel(modelFcn, velocitytable(:, 1), params);

    % i_0 = find(datatable(:, 1) > 2.77, 1);
    % i_e = length(datatable(:, 1));
    % i_0 = find(datatable(:, 1) > 2.762, 1); % start a bit earlier
    % i_e = find(datatable(:, 1) > 2.9, 1); % not all the way in
    % i_e = length(datatable(:, 1));
    validZone = datatable(:, 1) > velocitytable(2, 1) - 0.1 & datatable(:, 1) < velocitytable(2, 1) ... % pre-slack steady state
        | datatable(:, 1) > velocitytable(3, 1) + 0.002 ... & datatable(:, 1) < velocitytable(4, 1) + 0.01 ... Redevelopment zone
        | datatable(:, 1) > velocitytable(5, 1) + 0.01; % after the titin transient
    nonrepeating = diff(out.t) ~= 0;
    Fi = interp1(out.t(nonrepeating), out.Force(nonrepeating), datatable(validZone, 1));
    % plot(datatable(validZone, 1), datatable(validZone, 3));hold on;
    % plot(datatable(validZone, 1), Fi, '|');

    e = (datatable(validZone, 3) - Fi).^2;
    % e(isnan(e)) = 10;
    E(4) = mean(e(~isnan(e)))*20;


    % tet = [tet; params.N, et]

    if ~isempty(params.ghostSave)
        ghost = [out.t; out.Force]';
        save(['Ghost_' params.ghostSave '_slack'], 'ghost');
    end
    if params.PlotEachSeparately
        % axes('position',[0.55 0.1 0.4 0.35]); hold on;
        nexttile;
        yyaxis right;
        plot(datatable(:, 1),datatable(:, 2), '-', out.t, out.SL, 'o-', out.t, out.LXB, ':', 'Linewidth', 2, 'MarkerSize', 3);
        yyaxis left;hold on;

        % manage GHOST
        if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_slack.mat'],'file')
            ghost = load(['Ghost_' params.ghostLoad '_slack']);
            ghost = ghost.ghost;
            gp = plot(ghost(:, 1), ghost(:, 2), '-', 'Linewidth', 3, 'Color', [0.5843    0.8157    0.9882]);
        else
            clear gp;
        end


        plot(datatable(:, 1),datatable(:, 3),'k-','linewidth',0.5);
        plot(out.t,out.Force,'b-','linewidth',1.5);
        % plot(Tspan,F_active,'linewidth',1.5);
        xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
        ylabel('Force (rel.)','interpreter','latex','fontsize',16);
        set(gca,'fontsize',14);  box on;
        title('Slack');
        xlim([velocitytable(2, 1) velocitytable(end, 1)])
        yl = ylim;
        plot(out.t, out.XB_TORs, '-')
        ylim(yl)
        % xl = xlim();


        if exist('gp', 'var') && isvalid(gp)
            legend(['Ghost ' params.ghostLoad],'F data', 'F sim','SL data*', 'SL sim*', 'Location', 'southwest');
        else
            legend('F data', 'F sim','SL data*', 'SL sim*', 'Location', 'southwest');
        end
        nexttile;
        if params.NumberOfStates == 2
            plot(out.t, out.p1_0, '-', out.t, out.p2_0, '-', out.t, out.PuATP, '-',out.t, out.PuR, '-', out.t, out.SR, LineWidth=1.5, LineStyle='-')
            legend('P1','P2','PuATP','PuR', 'SR')
        elseif params.NumberOfStates == 3
            plot(out.t, out.p1_0, '-', out.t, out.p2_0, '-', out.t, out.p3_0, '-',out.t, out.PuATP, '-',out.t, out.PuR, '-', out.t, 1 - out.SR, LineWidth=1.5, LineStyle='-')
            legend('P1','P2','P3','PuATP','PuR', 'SR')
        end
        % xlim(xl);

    end
    %% SAVE FIG
    if params.PlotEachSeparately
        fig = gcf;
        % saveas(fig, ['XBBakersDataFit.png']);
        if ~isempty(params.ghostSave)
            saveas(fig, ['XBBakersDataFit_' params.ghostSave '.png']);
        end
    end
end
%% FORCE VELOCITY RESIMULATION
if params0.RunForceVelocityTime
    params = params0;

    % datafile = "data/2021 06 15 isovelocity fit Filip.xlsx";
    % datatable = readtable(datafile, ...
    %     "filetype", 'spreadsheet', ...
    %     'VariableNamingRule', 'modify', ...
    %     'Sheet', '8 mM', ...
    %     'Range', 'A5:C86004');
    % 
    % datatable.Properties.VariableNames = {'Time', 'L', 'F'};
    % datatable.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};
    % 
    % params.datatable = table2array(datatable);
    % datatable = table2array(datatable);
    % save('data/isovelocity.mat', 'datatable')
    datatable = load('data/isovelocity.mat').datatable;
    datatable(:, 2) = 2*datatable(:, 2);
    datatable(:, 1) = datatable(:, 1)/1000;
    params.datatable = datatable;

    params.UseSLInput = true;

    params.SL0 = 2.2;
    params.Slim_l = 1.5;
    params.Slim_r = 2.25;
    % params.LXBpivot = 2.2;
    % params.dS = 0.0025;
    if isfield(params, 'PU0')
        params = rmfield(params, 'PU0');
    end

    % reset the PU0
    params = getParams(params, params.g, true);
    % velocitytable(1, 1) = 2;
    % velocitytable = velocitytable(1:4, 1);
    % [F out] = evaluateModel(modelFcn, velocitytable(:, 1), params);
    [F out] = evaluateModel(modelFcn, [0 1], params);
end