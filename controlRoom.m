%% Control room for evaluating the simulations
% init
load gopt;
LoadData;

% Set up environment
MgADP = 0; 
Pi    = 0;
F_active_0 = zeros(length(vel), length(MgATP));
t_sl0 = [0 0.11]; % time at which the SL = 2.2 - time to stop the experiment
t_ss = [0 0.3]; %% steady state time

fcn = @dPUdT;
opts = struct('N', 40, 'Slim', 0.05, 'PlotProbsOnFig', 0, 'ValuesInTime', 1);
%%

[F outs] = evaluateModel(fcn, -1, t_ss, MgATP(3),Pi,MgADP,g, opts);

figure(4);clf;
subplot(211);
plot(outs.t, outs.F);
title('Force');
xlabel('time')

subplot(223);
plot(outs.t, outs.p1_0,outs.t, outs.p2_0, outs.t, outs.p3_0);
title('Zero moment to time')

subplot(224);
plot(outs.t, outs.p1_1,outs.t, outs.p2_1, outs.t, outs.p3_1);
title('First moment to time');



if isfield(outs, 'ps0_t')

    figure(5);clf;
    plot(outs.t, outs.ps0_t);
    title('S(1) maximal value (checking the boundary)')
end

