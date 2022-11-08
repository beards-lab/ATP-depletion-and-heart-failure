%% load Anthony Baker's experiments
load gopt;
% decimation sampling (each x)
dsf = 10;
ML = 2.0;
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
%%
% time segment data
% ts_d = [490, 500.25, 500.9, 502, 509, 509.25, 510, 519.5,539.8, 550];
ts_d = [950.0    1000.3    1001.0    1080.9 1081.75 1091.2    1111.6    1121.8   1200];
% time segment simulation
ts_s = [500 ts_d(2:end)];
clf;
[datatable, velocitytable] = DownSampleAndSplit(dt2, ts_d, ts_s, ML, dsf, 1, 'ForceLength2mM');
velocitytable
%% load length-force data for 8 mM
datafile = "data/2021 06 15 isovelocity fit Filip.xlsx";
dt8 = readtable(datafile, ...
    "filetype", 'spreadsheet', ...
    'VariableNamingRule', 'modify', ...
    'Sheet', '8 mM', ...
    'Range', 'A5:C86004');

dt8.Properties.VariableNames = {'Time', 'L', 'F'};
dt8.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};
datalabel = "8 mM";
%%
% ts_d = [450 500.25, 500.9, 509.25, 510, 519.5,539.8, 590];
% ts_s = [400 ts_d(2:end)]
ts_d = [950.0    1000.3    1001.0    1091.2   1092 1101.6    1121.8   1150];
ts_s = [600 ts_d(2:end)]

clf;
[datatable, velocitytable] = DownSampleAndSplit(dt8, ts_d, ts_s, ML, 5, 1, 'ForceLength8mM');
subplot(211);title('Length (ML)');xlabel('Time (ms)');ylabel('ML');
subplot(212);title('Force (kPa)');xlabel('Time (ms)');ylabel('kPa');
velocitytable
%% load step-up data for 8 mM
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
ts_d = [10, 20.6, 40.7, 62.1, 80.1, 101.3, 121.4, 141.7, 161.5, 181.7, 201.8, 221.9, 241.9, 261.9, 281.9, 300.5, 321.9, 500];
ts_s = [-200 ts_d(2:end)]
% clf;
% [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, dsf*10, 93/70, 'bakers_rampup8');
[datatable, velocitytable] = DownSampleAndSplit(data_table, [0 0], ts_s, ML, dsf*10, 1, 'bakers_rampup8');
velocitytable

%% Proof that the velocities are in ML/s and that the ML = SL0
% poi = [506,	509, 1075	1085, 1519	1523, 2007.5, 2010.5, 2670	2690 3110.5	3115.5 3627	3635 4100, 4200];
% for i = 1:length(poi)/2
% i_pois(i) = find(dt8.Time > poi(i*2-1), 1);
% i_poie(i) = find(dt8.Time > poi(i*2), 1);
% inds = i_pois(i):i_poie(i);
% v(i) = mean([diff(dt8.L(inds))./diff(dt8.Time(inds))*1e3]);
% end
% v
% % clf;plot(dt8.Time, v)
% ts(end) = 150;

%% overlap the segments on top of each other and average that
% segle = 800;
% segment = zeros(800, 3);
% s = 0;
%     segment(:, 1) = t(1:segle);
% for it = 2:2:(length(ts)-1)
%     ind = find(t >= ts(it), 1);
%     segment(:, 2) = segment(:, 2) + l(ind:ind+segle-1);
%     force = f(ind:ind+segle-1);
%     fmi = min(force);   fma = max(force); norm = 1/fma;
%     segment(:, 3) = segment(:, 3) + force/fma;
%     s = s+1;
% end
% segment(:, 2:3) = segment(:, 2:3)./s;
% clf;subplot(211)
% plot(segment(:, 1), segment(:, 2));
% subplot(212)
% plot(segment(:, 1), segment(:, 3));

%% Simulate the step-up experiment
% TODO update this

fcn = @dPUdTCa;
simulateTimes = velocitytable(:, 1);
velocities = velocitytable(1:end-1, 2);

params = struct('Pi', 0, 'MgATP', 8, 'MgADP', 0, 'Ca', 1000,...
    'Velocity', velocities, ...
    'UseCa', false,'UseOverlap', false,'kSE', 1000, 'mu', 10,...
    'N', 40, 'Slim', 0.025, 'ValuesInTime', 1, ...
    'BreakingVelocity', -10, 'SlackVelocity', 10, 'SL0', 2.0, ...
    'MatchTimeSegments', 1, 'ML', ML, 'PlotProbsOnStep', false, ...
    'ReduceSpace', false...
    );

params = getParams(params, g)
% figure(1);clf;

[F out] = evaluateModel(fcn, simulateTimes, params);


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

%% Animate states
AnimateStateProbabilities(out);


%% function definition
function [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, dsf, scaleF, saveAs)

    if isempty(ts_d)
        ts_d = [0 data_table.Time(end)];
    end

    % Relabel and downsample 
    imin_d = find(data_table.Time >= ts_d(1), 1);
    imax_d = find(data_table.Time >= ts_d(end), 1);

    imin_s = find(data_table.Time >= ts_s(1), 1);
    imax_s = find(data_table.Time >= ts_s(end), 1);    

    t = data_table.Time(imin_s:imax_s);
    tf = data_table.Time(imin_d:imax_d);
    td = downsample(tf, dsf);
    l = data_table.L(imin_s:imax_s);
    lf = movmean(data_table.L(imin_d:imax_d),[dsf/2 dsf/2]); % l filtered
    ld = downsample(lf, dsf);

    f = data_table.F(imin_s:imax_s)*scaleF;
    ff = movmean(data_table.F(imin_d:imax_d)*scaleF,[dsf/2 dsf/2]); % force filtered
    fd = downsample(ff, dsf);
    
    datatable = [td/1000, ld*ML, fd];

    % Split it into segments with const velocities
    vs = [];
    for it = 1:(length(ts_s)-1)
        t1 = find(t >= ts_s(it), 1);
        t2 = find(t >= ts_s(it + 1), 1);
        vs(it) = round((l(t1) - l(t2))/(t(t1) - t(t2))*1000, 1);
    end
    vsum = vs*ML; % 
    % time (s), velocity ML/s, velocity um/s
    velocitytable = [ts_s/1000;[vs 0];[vsum 0]]'; 
    

    if ~isempty(saveAs)
        % Export the data into modelica-readable format and for identificatoin
        fn = ['data/' saveAs '.mat'];
        save(fn,  'datatable', 'velocitytable');
        disp(['Saved as ' fn])
    end

%     figure(1);clf;
    subplot(211);hold on;
    plot(t, l, tf, lf, td, ld, '|-');
    plot([ts_d;ts_d], repmat([min(data_table.L);max(data_table.L)], 1, length(ts_d)))
    subplot(212);hold on;
    plot(t, f, tf, ff, td, fd, '|-', 'Linewidth', 2, 'MarkerSize', 10);
    plot([ts_d;ts_d], repmat([min(data_table.F);max(data_table.F)], 1, length(ts_d)))
end