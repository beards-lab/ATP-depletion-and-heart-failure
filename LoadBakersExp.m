%% length-velocity protocol
% load Anthony Baker's experiments
datafile = "data/2021 06 15 isovelocity fit Filip.xlsx";
load gopt;

%% load length-force data for 2 mM
data_table = readtable(datafile, ...
    "filetype", 'spreadsheet', ...
    'VariableNamingRule', 'modify', ...
    'Sheet', '2 mM', ...
    'Range', 'A5:C86004');

data_table.Properties.VariableNames = {'Time', 'L', 'F'};
data_table.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};
datalabel = "2 mM";
%% load length-force data for 8 mM
data_table = readtable(datafile, ...
    "filetype", 'spreadsheet', ...
    'VariableNamingRule', 'modify', ...
    'Sheet', '8 mM', ...
    'Range', 'A5:C86004');

data_table.Properties.VariableNames = {'Time', 'L', 'F'};
data_table.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};
datalabel = "8 mM";
ts = [490, 500.25, 500.9, 509.25, 510, 519.5,539.8, 550];
SL0 = 2.2*1.1;

%% load step-up data for 2 mM
datafile = "data/06 21 21 Ramps 2 mM ATP.xlsx";
data_table = readtable(datafile, ...
    "filetype", 'spreadsheet', ...
    'VariableNamingRule', 'modify', ...
    'Sheet', 'Sheet1', ...
    'Range', 'A5:C10005');

data_table.Properties.VariableNames = {'Time', 'L', 'F'};
data_table.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};
datalabel = "step-up 2 mM";
ts = [-200, 20.6, 40.7, 62.1, 80.1, 101.3, 121.4, 141.7, 161.5, 181.7, 201.8, 221.9, 241.9, 261.9, 281.9, 300.5, 321.9, 500];
% ts = ts(1:4);
stopTime = ts(end);    
SL0 = 2.2;

%% plot
figure(1); clf; 
subplot(211);hold on;
plot(data_table.Time, data_table.L, 'x-');
% xlim(xl);
subplot(212);hold on;
plot(data_table.Time, data_table.F);
% xlim(xl);

%% Plot 8 mM
subplot(211);
plot(data_table.Time, data_table.L);
subplot(212);
plot(data_table.Time, data_table.F);

%% Relabel and downsample 
dsf = 5;
maxvel = 50;
imax = find(data_table.Time >= stopTime, 1);
t = data_table.Time(1:imax);
td = downsample(t, dsf);
l = data_table.L(1:imax);
ld = downsample(l, dsf);
f = data_table.F(1:imax);
fd = downsample(f, dsf);
dL = [diff(data_table.L);0];
ddL = [diff(max(min(1000*diff(movmean(ld,[5 5])), maxvel), -maxvel));0;0];
dT = [diff(data_table.Time);0];
%% Extract the speed
% dL = [0;diff(ld)];
% dT = [0;diff(td)];
sp8 = dL./dT*1000;
sp8a = movmean(max(min(sp8, maxvel), -maxvel),[5 5]);
sp8ad = downsample(movmean(max(min(sp8, maxvel), -maxvel),[5 5]), dsf);
% sp8a = smooth(sp8, 1000, 'lowess');
xl = [500 540];
% yl = [-10 10];

figure(1);clf;
% subplot(311);
plot(data_table.Time, data_table.L, td, ld)
xlim(xl);

% subplot(312);plot(t, sp8a,t, sp8a, td, sp8ad, 'x-')
% xlim(xl)
% 
% subplot(313);plot(td, ddL,'x-')
% xlim(xl)
% 
% ddl_th = 1e-4;
% ylim([-ddl_th ddl_th])
%% get points of interests
% ramping down
% pif = ddl > 1e-4 && circshift(ddl, -1) < 1e-4; %inflection points
% TODO

%% Split it into segments with const velocities

vs = [];
for it = 1:(length(ts)-1)
    t1 = find(t >= ts(it), 1);
    t2 = find(t >= ts(it + 1), 1);
    vs(it) = round((l(t1) - l(t2))/(t(t1) - t(t2))*1000, 1);
end
[ts;[vs 0]]

%% reconstruct from vel

%% Throw it into evaluateModel
fcn = @dPUdT;
simulateTimes = ts/1000;
velocities = vs;
params = struct('Pi', 0, 'MgATP', 8, 'MgADP', 0, 'Ca', 1000,...
    'Velocity', velocities, ...
    'UseCa', false,'UseOverlap', false);

% TODO breaking and slacking velocities
opts = struct('N', 50, 'Slim', 0.08, 'PlotProbsOnFig', 0, 'ValuesInTime', 1, ...
    'BreakingVelocity', -10, 'SlackVelocity', 10, 'SL0', SL0, 'MatchTimeSegments', 1, 'ML', 2.2);


[F outs] = evaluateModel(fcn, simulateTimes, params, g, opts);


% vel = 1:20;
% dt = 1e-4;
% N = opts.N; % space (strain) discretization--number of grid points in half domain
% Slim = opts.Slim; 
% % dS = Slim/N;
% dS = dt*vel;
% Slim = dS*N
% dt = dS./abs(vel)
% tend = T(end)./abs(vel); % ending time of simulation
% Nstep = round(tend/dt);% = Tspan(end)/dS
%
% figure;
%
clf;
subplot(211);hold on;
plot(t, l, '-');
plot(outs.t*1000, outs.SL/2.2, '|-', 'MarkerSize', 2)
plot([simulateTimes;simulateTimes]*1000, repmat([min(ld);max(ld)], [1 size(simulateTimes, 2)]))
legend('Muscle length (-), Data', 'Muscle length (-), Simulation')
% xlim([0.45 0.55])
xlim([0 inf])

subplot(212);hold on;
plot(td, fd);
plot(outs.t*1000, outs.F, '|-', 'MarkerSize', 2)
legend('Force (kPa?), Data', 'Force (mmHg), Simulation')
xlim([0 inf])
ylim([-50, Inf])

