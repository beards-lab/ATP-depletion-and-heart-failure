% driver code
g = ones(30, 1);
LoadData;

params = getParams();
params.g = g;
params.SL0 = 2.2;
% params.Slim = 0.18;
params.Slim = 0.3;
params.N = 60;
params.MgATP = 8;
t_ss = [0 1];
t_sl0 = [0 0.1];

figure(1);clf;
ghostSave = '';
% ghostSave = 'beardsOrig_passive';
% ghostSave = 'ShiftingStrain40_Slim0_3';
% ghostSave = 'ShiftingStrain80_Slim0_3';
ghostSave = 'operatorSplittingPU020';
ghostSave = 'operatorSplittingPU0';
ghostSave = 'operatorSplittingPU080';
ghostSave = 'ShiftingStrain160_Slim0_3';
ghostSave = 'beardsOrig_all60';
ghostSave = '';

ghostLoad = '';
% ghostLoad = 'beardsOrig';
% ghostLoad = 'ShiftingStrain40';
% ghostLoad = 'operatorSplitting'
ghostLoad = 'operatorSplittingPU020';
ghostLoad = 'operatorSplittingPU0';% N = 40
ghostLoad = 'operatorSplittingPU080';
ghostLoad = 'ShiftingStrain40_Slim0_3';
ghostLoad = 'ShiftingStrain80_Slim0_3';
ghostLoad = 'ShiftingStrain160_Slim0_3';
ghostLoad = 'ShiftingStrainTest20';
ghostLoad = 'ShiftingStrainTest40';
ghostLoad = 'ShiftingStrainTest80';
ghostLoad = 'ShiftingStrainTest160';
ghostLoad = 'ShiftingStrainTest240';
ghostLoad = 'beardsOrig_all40';

% ghostLoad = '';

% testing setup
% params.UseOverlap = true;
% params.UsePassive = true;

params.UseTORNegShift = false;
params.UseMutualPairingAttachment = false;
% set as a default to be modified
% params0 = params;
params.PlotEachSeparately = true;
tic
%% FORCE VELOCITY

params = getParams(params);
F_active = [];
params.UseSerialStiffness = true;
if isfield(params, 'PU0')
    params = rmfield(params, 'PU0');
end
for j = 1:length(vel)
    if vel(j) == 0 
        params.SL0 = 2.0;
        params.Velocity = 0;
        [F_active(j) out] = evaluateModel(@dPUdTCa, t_ss, params);        
    else
        params.SL0 = 2.2;
        if ~isfield(params, 'PU0')
            % speed things up by storing the initialization
            params.Velocity = 0;
            [~, out] = evaluateModel(@dPUdTCa, t_ss, params);
            params.PU0 = out.PU(end, :);
        end
        params.Velocity = vel(j);
        [F_active(j) out] = evaluateModel(@dPUdTCa, t_sl0/abs(vel(j)), params);
    end        
end

% cost function
E(1) = sum(abs(F_active(:,1)-Data_ATP(:,1+1)).^2);
% normalize by number of data points
E(1) = E(1)/size(Data_ATP, 1);
  

%
% figure(101); clf; axes('position',[0.1 0.6 0.35 0.35]); hold on;
% figure(); clf; axes('position',[0.1 0.6 0.35 0.35]); hold on;

if ~isempty(ghostSave)
    ghost = [F_active(:), -vel(:)];
    save(['Ghost_' ghostSave '_FV'], 'ghost');
end

if params.PlotEachSeparately
%     figure(102);hold on;
    axes('position',[0.05 0.6 0.4 0.35]); 
    hold on;

    if ~isempty(ghostLoad) && exist(['Ghost_' ghostLoad '_FV.mat'],'file')
        ghost = load(['Ghost_' ghostLoad '_FV']);
        ghost = ghost.ghost;
        gp = plot(ghost(:, 1), ghost(:, 2), '-', 'Linewidth', 3, 'Color', [0.5843    0.8157    0.9882]);
    else
        clear gp;
    end

    %
    plot(F_active(:), -vel,'b-^','linewidth',1);
    ylabel('Velocity (ML/s)','interpreter','latex','fontsize',16);
    xlabel('Force (kPa)','interpreter','latex','fontsize',16);
    set(gca,'fontsize',14); 
    % axis([0 65 0 6]);
    axis([-10 65 0 6]);
    title('Force-velocity')
    box on;grid on;
    plot(Data_ATP(:,2),Data_ATP(:,1),'bo','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
    % plot(Data_ATP(:,3),Data_ATP(:,1),'go','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
    % plot(Data_ATP(:,4),Data_ATP(:,1),'ro','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);

    if exist('gp', 'var') && isvalid(gp)
        legend(['Ghost ' ghostLoad], 'Sim', 'Data');
    else
        legend('Sim', 'Data');
    end
end                            

%% KTR EXPERIMENT
params.SL0 = 2.0;
params.UseSlack = true;
% takes like 20mins
% params.UseSerialStiffness = true;
% params.UseSerialStiffness = false;
params.N = 40;
params.LXBpivot = 1.8;
params = getParams(params, params.g, true);

% replicate the ktr protocol
v = 150; % ML/s
%       times  = [-1e3, 0,  2, 20, 22.5, 25 , 25.5, 1e3]/1000 - 25.5e-3;
pos_ML = [1   , 1,0.8,0.8, 1.05,1.05,     1, 1 ];

% putting the numbers as a difference
times = cumsum([0 , 1,0.2/v,0.01 - 0.2/v,0.25/v, 0.005 - 0.25/v, 0.05/v, 1]);
params.Velocity = diff(pos_ML)./diff(times);
[t, out] = evaluateModel(@dPUdTCa, times - times(end-1), params);

% calculate ktr
i_0 = find(out.t > 0 & out.FXB > 0, 1);
Frel = out.FXB(i_0:end)./out.FXB(end);
% relative to max value from data to run optimizer faster
% Frel = out.FXB(i_0:end)./55;
i = find(Frel >= 1-exp(-1), 1);
Ktr = 1/out.t(i+i_0);
E(2) = abs(Ktr-Ktr_mean(1)).^2;

if ~isempty(ghostSave)
    ghost = [out.t;out.Force/out.Force(end)]';
    save(['Ghost_' ghostSave '_ktr'], 'ghost');
end
if params.PlotEachSeparately
    
    axes('position',[0.55 0.6 0.4 0.35]); hold on;
    datastruct = load('data/bakers_ktr_8.mat');
    datatable = datastruct.datatable;
    yyaxis right;
    plot(datatable(:, 1),datatable(:, 2), out.t, out.SL, out.t, out.LXB);
    yyaxis left;

    % manage GHOST
    if ~isempty(ghostLoad) && exist(['Ghost_' ghostLoad '_ktr.mat'],'file')
        ghost = load(['Ghost_' ghostLoad '_ktr']);
        ghost = ghost.ghost;
        gp = plot(ghost(:, 1), ghost(:, 2), '-', 'Linewidth', 3, 'Color', [0.5843    0.8157    0.9882]);
    else
        clear gp;
    end

    plot(datatable(:, 1),datatable(:, 3)/datatable(end, 3),'k-','linewidth',1);

    plot(out.t,out.Force/out.Force(end),'b-','linewidth',1.5);
    % plot(Tspan,F_active,'linewidth',1.5);
    xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
    ylabel('Force (rel.)','interpreter','latex','fontsize',16);
    set(gca,'fontsize',14,'ylim',[0 1.1], 'xlim', [-0.05 0.45]);  box on;
    title('Speed of the transient');

    if exist('gp', 'var') && isvalid(gp)
        legend(['Ghost ' ghostLoad], 'F data', 'F sim','SL data*', 'SL sim*', 'LXB sim*', 'Location', 'southeast');
    else
        legend('F data', 'F sim','SL data*', 'SL sim*', 'LXB sim*', 'Location', 'southeast');
    end
end
%% RAMP UP
datastruct = load('data/bakers_rampup2_8.mat');
datatable = datastruct.datatable;    
velocitytable = datastruct.velocitytable;
velocitytable(1, 1) = -1; % enough time to get to steady state
params.Slim = 0.25;
params.Velocity = velocitytable(1:end-1, 2);
params.SL0 = 2.0;
params.LXBpivot = 2.1;
% params.N = 30;
% update params with new N and Slims
params = getParams(params, g, true);

[F out] = evaluateModel(@dPUdTCa, velocitytable(:, 1), params);

Fi = interp1(out.t, out.Force, datatable(:, 1));
e = (datatable(:, 3) - Fi).^2;
E(3) = mean(e)*10;

if ~isempty(ghostSave)
    ghost = [out.t; out.Force]';
    save(['Ghost_' ghostSave '_rampup'], 'ghost');
end
if params.PlotEachSeparately

    axes('position',[0.05 0.1 0.4 0.35]); hold on;
    yyaxis right;
    plot(datatable(:, 1),datatable(:, 2), out.t, out.SL);
    yyaxis left;

    % manage GHOST
    if ~isempty(ghostLoad) && exist(['Ghost_' ghostLoad '_rampup.mat'],'file')
        ghost = load(['Ghost_' ghostLoad '_rampup']);
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
        legend(['Ghost ' ghostLoad],'F data', 'F sim','SL data*', 'SL sim*', 'Location', 'Best');
    else
        legend('F data', 'F sim','SL data*', 'SL sim*', 'Location', 'Best');
    end
end
%% SLACK

datastruct = load('data/bakers_slack8mM.mat');
datatable = datastruct.datatable;    
velocitytable = datastruct.velocitytable(4:end, :);
params.Velocity = velocitytable(:, 2);
params.datatable = datatable;
params.UseSLInput = false;
params.UseSlack = true;
params.SL0 = 2.2;
params.Slim = 0.3;
params.UseSlack = true;
params.LXBpivot = 2.0;
if isfield(params, 'PU0')
    params = rmfield(params, 'PU0');
end

% reset the PU0
params = getParams(params, g, true);

[F out] = evaluateModel(@dPUdTCa, velocitytable(:, 1), params);

i_0 = find(datatable(:, 1) > 2.77, 1);
nonrepeating = diff(out.t) ~= 0;
Fi = interp1(out.t(nonrepeating), out.Force(nonrepeating), datatable(i_0:end, 1));
e = (datatable(i_0:end, 3) - Fi).^2;
e(isnan(e)) = 10;
E(4) = mean(e)*20;


% tet = [tet; params.N, et]

if ~isempty(ghostSave)
    ghost = [out.t; out.Force]';
    save(['Ghost_' ghostSave '_slack'], 'ghost');
end
if params.PlotEachSeparately
    axes('position',[0.55 0.1 0.4 0.35]); hold on;
    yyaxis right;
    plot(datatable(:, 1),datatable(:, 2), '-', out.t, out.SL,'|--', 'linewidth',1);
    yyaxis left;hold on;

    % manage GHOST
    if ~isempty(ghostLoad) && exist(['Ghost_' ghostLoad '_slack.mat'],'file')
        ghost = load(['Ghost_' ghostLoad '_slack']);
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
    set(gca,'fontsize',14, 'xlim', [2.4 3.1]);  box on;
    title('Slack');

    if exist('gp', 'var') && isvalid(gp)
        legend(['Ghost ' ghostLoad],'F data', 'F sim','SL data*', 'SL sim*', 'Location', 'southwest');
    else
        legend('F data', 'F sim','SL data*', 'SL sim*', 'Location', 'southwest');
    end
end
%% SAVE FIG
toc
if params.PlotEachSeparately
    fig = gcf;
    % saveas(fig, ['XBBakersDataFit.png']);
    if ~isempty(ghostSave)
     saveas(fig, ['XBBakersDataFit_' ghostSave '.png']);
    end
end

