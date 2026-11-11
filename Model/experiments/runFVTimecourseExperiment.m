function out = runFVTimecourseExperiment(params0, Data_ATP, FV_velocities, F_active)
%RUNFVTIMECOURSEEXPERIMENT Re-simulate and plot Baker isovelocity timecourse.
%
%   out = runFVTimecourseExperiment(params0, Data_ATP, FV_velocities, F_active)
%
%   Loads bakers_isovelocity.mat, runs a single full-timecourse ODE simulation
%   across all shortening velocities, extracts the force value at the moment
%   SL crosses 2.0 um for each velocity, and plots the result against the
%   experimental timecourse and the steady-state FV curve.
%
%   Inputs:
%     params0       - Base parameter struct. Must include params0.modelFcn (string).
%     Data_ATP      - Experimental FV data matrix (from LoadBakersExp); column 1
%                     is velocity, columns 2+ are force at each ATP concentration.
%     FV_velocities - Velocity vector used in the steady-state FV run (for overlay).
%     F_active      - Force matrix from runFVExperiment (ATP x velocity); pass []
%                     if the FV experiment was not run in this session.
%
%   Outputs:
%     out - Output struct from evaluateModel (fields: t, SL, Force, etc.)
%
%   See also: runFVExperiment, evaluateModel

    modelFcn = str2func(params0.modelFcn);

    isovelocity = load('bakers_isovelocity.mat', 'datatable', 'velocitytable');
    datatable    = isovelocity.datatable;
    velocitytable = isovelocity.velocitytable;

    params = params0;
    params.datatable = datatable;
    params.Velocity  = velocitytable(:, 2);

    params.SL0    = 2.2;
    params.Slim_l = 1.5;
    params.Slim_r = 2.25;
    if isfield(params, 'PU0')
        params = rmfield(params, 'PU0');
    end

    [~, out] = evaluateModel(modelFcn, velocitytable(:, 1), params);

    constantSpeedSegment = find(velocitytable(:, 2) > 0) - 3;
    speeds = -velocitytable(constantSpeedSegment, 2);

    smashed2_0_data  = arrayfun(@(x) find(datatable(:, 2) < 2.0 & datatable(:, 1) > x, 1, 'first'), velocitytable(constantSpeedSegment, 1), 'UniformOutput', true);
    smashed2_0_model = arrayfun(@(x) find(out.SL <= 2.0 & out.t > x, 1, 'first'), velocitytable(constantSpeedSegment, 1), 'UniformOutput', true);
    % add a final one at speed 0
    speeds           = [speeds; 0];
    smashed2_0_data  = [smashed2_0_data;  length(datatable(:, 1))];
    smashed2_0_model = [smashed2_0_model; length(out.t)];

    clf;
    nexttile; hold on;
    plot(datatable(:, 1), datatable(:, 2)/2*100, 'k-', datatable(:, 1), datatable(:, 3), 'k--');
    plot(datatable(smashed2_0_data, 1), datatable(smashed2_0_data, 3), 'o', 'LineWidth', 4, 'MarkerSize', 12);

    [~, ~] = ismember([0.5, 1], speeds); %#ok<ASGLU>

    plot(out.t, out.Force, out.t, out.SL/2*100);
    plot(out.t(smashed2_0_model), out.Force(smashed2_0_model), 'x', 'LineWidth', 4, 'MarkerSize', 12);

    nexttile; hold on;
    cla; hold on;
    modelForce = out.Force(smashed2_0_model);

    [~, speed_sorted] = sort(speeds);
    a = 1; % ATP 8 mM
    plot(modelForce(speed_sorted), speeds(speed_sorted), 'x-', 'LineWidth', 2, 'MarkerSize', 12);

    plot(Data_ATP(:,a+1), Data_ATP(:,1), 'o', 'linewidth', 2.5, 'Markersize', 12, 'markerfacecolor', [1 1 1])
    plot(datatable(smashed2_0_data(speed_sorted), 3), speeds(speed_sorted), 's', 'LineWidth', 2.5, 'MarkerSize', 12);

    if ~isempty(F_active)
        % overlay the steady-state FV curve for comparison
        plot(F_active, -FV_velocities(1:end), '+--', 'LineWidth', 2, 'MarkerSize', 12);
        legend('SIM - isovelocity timecourse', 'Data - isovelocity', 'Data - isovelocity timecourse (Baker 2021)', 'SIM isovelocity steady');
    else
        legend('data', 'Force-velocity timecourse');
    end
end
