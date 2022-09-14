%% length-velocity protocol
% load Anthony Baker's experiments
datafile = "data/2021 06 15 isovelocity fit Filip.xlsx";

%% load length-force data for 2 mM
data_lf2 = readtable(datafile, ...
    "filetype", 'spreadsheet', ...
    'VariableNamingRule', 'modify', ...
    'Sheet', '2 mM', ...
    'Range', 'A5:C86004');

data_lf2.Properties.VariableNames = {'Time', 'L', 'F'};
data_lf2.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};

%% load length-force data for 8 mM
data_lf8 = readtable(datafile, ...
    "filetype", 'spreadsheet', ...
    'VariableNamingRule', 'modify', ...
    'Sheet', '8 mM', ...
    'Range', 'A5:C86004');

data_lf8.Properties.VariableNames = {'Time', 'L', 'F'};
data_lf8.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};

%% plot 2 mM
figure(1); clf; 
subplot(211);hold on;
plot(data_lf2.Time, data_lf2.L);
% xlim(xl);
subplot(212);hold on;
plot(data_lf2.Time, data_lf2.F);
% xlim(xl);

%% Plot 8 mM
subplot(211);
plot(data_lf8.Time, data_lf8.L);
subplot(212);
plot(data_lf8.Time, data_lf8.F);

%% Extract the speed
dsf = 5;
maxvel = 50;
t = data_lf8.Time;
td = downsample(t, dsf);
l = data_lf8.L;
ld = downsample(l, dsf);
fd = downsample(data_lf8.F, dsf);
dL = [diff(data_lf8.L);0];
ddL = [diff(max(min(1000*diff(movmean(ld,[5 5])), maxvel), -maxvel));0;0];
dT = [diff(data_lf2.Time);0];
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
plot(data_lf8.Time, data_lf8.L, td, ld)
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
ts = [490, 500.25, 500.9, 509.25, 510, 519.5,539.8, 550];
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
    'BreakingVelocity', -10, 'SlackVelocity', 10, 'SL0', 2.2*1.1, 'MatchTimeSegments', 1, 'ML', 2.2);


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
plot(t/1000, l, '|-');
plot([simulateTimes;simulateTimes], repmat([min(ld);max(ld)], [1 size(simulateTimes, 2)]))
plot(outs.t, outs.SL/2.2, '*-')
xlim([0.45 0.55])
subplot(212);hold on;
plot(td/1000, fd);
plot(outs.t, outs.F*7.5, '*-')
xlim([0.45 0.55])
ylim([-50, Inf])

