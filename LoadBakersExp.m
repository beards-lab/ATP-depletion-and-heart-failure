%% length-velocity protocol
% load Anthony Baker's experiments
load gopt;

%% load length-force data for 2 mM
datafile = "data/2021 06 15 isovelocity fit Filip.xlsx";
dt2 = readtable(datafile, ...
    "filetype", 'spreadsheet', ...
    'VariableNamingRule', 'modify', ...
    'Sheet', '2 mM', ...
    'Range', 'A5:C86004');

dt2.Properties.VariableNames = {'Time', 'L', 'F'};
dt2.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};
datalabel = "2 mM";
SL0 = 2.2*1.1;
%% load length-force data for 8 mM
dt8 = readtable(datafile, ...
    "filetype", 'spreadsheet', ...
    'VariableNamingRule', 'modify', ...
    'Sheet', '8 mM', ...
    'Range', 'A5:C86004');

dt8.Properties.VariableNames = {'Time', 'L', 'F'};
dt8.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};
datalabel = "8 mM";
ts = [490, 500.25, 500.9, 509.25, 510, 519.5,539.8, 550];
SL0 = 2.2*1.1;
stopTime = ts(end);
%% plot
figure(1); clf; 
subplot(211);hold on;
plot(dt2.Time, dt2.L, '-');
plot(dt8.Time, dt8.L, '-');
yyaxis right
plot(dt8.Time, [0;diff(dt8.L)./diff(dt8.Time)], '-');
title('Length')
xlabel('time (ms)');
ylabel('Length (ML)');
legend('ATP 2mM', 'ATP 8mM');
% 0.85 0.95

% xlim(xl);
subplot(212);hold on;
plot(dt2.Time, dt2.F, '-');
plot(dt8.Time, dt8.F, '-');
title('force')
xlabel('time')
ylabel('force (kPa)')
% xlim(xl);

%% load step-up data for 2 mM
% datafile = "data/06 21 21 Ramps 2 mM ATP.xlsx";
datafile = "data/06 21 21 Ramps 8mM ATP.xlsx";
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
% ts(end) = 200;
stopTime = ts(end);    
SL0 = 2.0;

%% Proof that the velocities are in ML/s and that the ML = SL0
poi = [506,	509, 1075	1085, 1519	1523, 2007.5, 2010.5, 2670	2690 3110.5	3115.5 3627	3635 4100, 4200];
for i = 1:length(poi)/2
i_pois(i) = find(dt2.Time > poi(i*2-1), 1);
i_poie(i) = find(dt2.Time > poi(i*2), 1);
inds = i_pois(i):i_poie(i);
v(i) = mean([diff(dt8.L(inds))./diff(dt8.Time(inds))*1e3]);
end
v
% clf;plot(dt8.Time, v)

%% Plot 8 mM
clf;
subplot(211);hold on;
plot(data_table.Time, data_table.L);
plot([ts;ts], repmat([min(data_table.L);max(data_table.L)], 1, length(ts)))
subplot(212);hold on;
plot(data_table.Time, data_table.F);
plot([ts;ts], repmat([min(data_table.F);max(data_table.F)], 1, length(ts)))
legend('2 mM ATP', '8 mM ATP')

%% Relabel and downsample 
dsf = 1;
% maxvel = 50;

imax = find(data_table.Time >= stopTime, 1);
t = data_table.Time(1:imax);
td = downsample(t, dsf);
l = data_table.L(1:imax);
lf = movmean(l,[dsf dsf]); % l filtered
ld = downsample(lf, dsf);

f = data_table.F(1:imax);
ff = movmean(f,[5 5]); % force filtered
fd = downsample(ff, dsf);

figure(1);clf;
subplot(211);plot(t, l, t, lf, td, ld, 'x-');
subplot(212);plot(t, f, t, ff, td, fd, 'x-', 'Linewidth', 2, 'MarkerSize', 10);

% fff = fft(f);
% clf;plot(1:length(fff), fff)

% dL = [diff(data_table.L);0];
% ddL = [diff(max(min(1000*diff(movmean(ld,[5 5])), maxvel), -maxvel));0;0];
% dT = [diff(data_table.Time);0];

%% Split it into segments with const velocities

vs = [];
for it = 1:(length(ts)-1)
    t1 = find(t >= ts(it), 1);
    t2 = find(t >= ts(it + 1), 1);
    vs(it) = round((l(t1) - l(t2))/(t(t1) - t(t2))*1000, 1);
end
vsum = vs*SL0; % 
% time (s), velocity ML/s, velocity um/s
velocitytable = [ts/1000;[vs 0];[vsum 0]]'; 

% ts(end) = 150;
%% Export the data into modelica-readable format and for identificatoin
datatable = [td/1000, ld*SL0, fd];
save data/bakers_rampup8.mat datatable velocitytable;
%% Extract the speed
% dL = [0;diff(ld)];
% dT = [0;diff(td)];
% sp8 = dL./dT*1000;
% sp8a = movmean(max(min(sp8, maxvel), -maxvel),[5 5]);
% sp8ad = downsample(movmean(max(min(sp8, maxvel), -maxvel),[5 5]), dsf);
% % sp8a = smooth(sp8, 1000, 'lowess');
% % xl = [500 540];
% % yl = [-10 10];
% 
% figure(1);clf;
% % subplot(311);
% plot(data_table.Time, data_table.L, td, ld)
% xlim(xl);

% subplot(312);plot(t, sp8a,t, sp8a, td, sp8ad, 'x-')
% xlim(xl)
% 
% subplot(313);plot(td, ddL,'x-')
% xlim(xl)
% 
% ddl_th = 1e-4;
% ylim([-ddl_th ddl_th])

%% overlap the segments on top of each other
segle = 800;
segment = zeros(800, 3);
s = 0;
    segment(:, 1) = t(1:segle);
for it = 2:2:(length(ts)-1)
    ind = find(t >= ts(it), 1);
    segment(:, 2) = segment(:, 2) + l(ind:ind+segle-1);
    force = f(ind:ind+segle-1);
    fmi = min(force);   fma = max(force); norm = 1/fma;
    segment(:, 3) = segment(:, 3) + force/fma;
    s = s+1;
end
segment(:, 2:3) = segment(:, 2:3)./s;
clf;subplot(211)
plot(segment(:, 1), segment(:, 2));
subplot(212)
plot(segment(:, 1), segment(:, 3));
%% get points of interests
% ramping down
% pif = ddl > 1e-4 && circshift(ddl, -1) < 1e-4; %inflection points
% TODO

%% test velocities
params.Velocity = [0 1];
opts.Slim = 0.06;
opts.N = 30;
params.UseOverlap = false;
[F out] = evaluateModelUpwind(fcn, [0 1 2], params, g, opts);

%% Throw it into evaluateModel
fcn = @dPUdTCa;
simulateTimes = ts/1000;
velocities = vs;

params = struct('Pi', 0, 'MgATP', 8, 'MgADP', 0, 'Ca', 1000,...
    'Velocity', velocities, ...
    'UseCa', false,'UseOverlap', false,'kSE', 1000, 'mu', 10);

% TODO breaking and slacking velocities
opts = struct('N', 40, 'Slim', 0.025, 'PlotProbsOnFig', 0, 'ValuesInTime', 1, ...
    'BreakingVelocity', -10, 'SlackVelocity', 10, 'SL0', SL0, ...
    'MatchTimeSegments', 1, 'ML', 2.0, 'PlotProbsOnStep', false, 'ReduceSpace', false,...
    'OutputAtSL', 2.1);

% figure(1);clf;

[F out] = evaluateModelUpwind(fcn, simulateTimes, params, g, opts);


%% Plot the lengths and forces
clf;
subplot(211);hold on;
% plot(t, l, '-');
plot(out.t*1000, out.SL/params.ML, '|-', 'MarkerSize', 2)
plot(out.t*1000, out.LXB/params.ML, '|-', 'MarkerSize', 2)
% plot([simulateTimes;simulateTimes]*1000, repmat([min(ld);max(ld)], [1 size(simulateTimes, 2)]))
legend('Muscle length (-), Data', 'Muscle length (-), Simulation', 'Sarcomere length (-), Simulation')
% xlim([0.45 0.55])
xlim([0 inf])

subplot(212);hold on;
% plot(td, fd);
plot(out.t*1000, out.Force, '|-', 'MarkerSize', 2)
plot(out.t*1000, out.FXB, '|-', 'MarkerSize', 2)
% plot([simulateTimes;simulateTimes]*1000, repmat([min([fd;out.F']);max([fd;out.F'])], [1 size(simulateTimes, 2)]))
legend('Force (kPa?), Data', 'Force (mmHg), Simulation', 'XB Force (mmHg), Simulation')

xlim([0 inf])
ylim([-50, Inf])

%% plot states
figure;clf; hold on;
Pus = 1 - out.p1_0 - out.p2_0 - out.p3_0;% PU substitute

leg = plot(out.t, Pus, out.t, out.p1_0, out.t, out.p2_0, out.t, out.p3_0, out.t, out.NR);
% plot([simulateTimes;simulateTimes], repmat([0; 1], [1 size(simulateTimes, 2)]))
legend(leg)
legend('Pu', 'P1', 'P2', 'P3', 'NR')



%% Animate the probs
% out = outs;


% times2 = linspace(out.t(1), out.t(end), 1000);
times2 = out.t;
ut = find(out.t ~= circshift(out.t, 1));
PU = interp1(out.t(ut), out.PU(ut, :), times2);


figure(1);clf;
% N = opts.N;

% params.dS = opts.Slim/opts.N;
%     params.Slim = opts.Slim;
%     params.s = (0:1:opts.N)*params.dS; % strain 
%     params.s = (-para.N:opts.N)*params.dS; % strain space

for i = 1:size(PU, 1)
    
    ss = length(params.s);
    p1 = PU(i, 1:ss);
    p2 = PU(i, ss +1:2*ss);
    p3 = PU(i, 2*ss+1:3*ss);
    NR = PU(i, 3*ss + 1);
    NP = PU(i, 3*ss +2);
    SL = PU(i, 3*ss +3);    
    LSe = PU(i, 3*ss +4);    
    subplot(211);cla;hold on;
    plot(params.s, p1, params.s, p2, params.s, p3);
    plot([out.Vxb(i) out.Vxb(i)]/100, [0 max([p1, p2, p3])])
%     ylim([0 500]);
    subplot(212);cla;hold on;
    plot(out.t, out.SL)
    plot(out.t, out.LXB)
    plot([times2(i) times2(i)], [2 2.2])
    xlim([out.t(1), out.t(end)])
    
%     drawnow;
    pause(1/30);
end
