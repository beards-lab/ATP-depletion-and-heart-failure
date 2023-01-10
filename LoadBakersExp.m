%% load Anthony Baker's experiments
load gopt;
% decimation sampling (each x)
dsf = 10;
ML = 2.0;
% normalized force multiplier
nf = 56;
clf;close all;
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
ts_d = [450 500.25, 500.9, 509.25, 510, 519.5,539.8, 650];
% ts_s = [400 ts_d(2:end)]
% ts_d = [ts_d 950.0    1000.3    1001.0    1091.2   1092 1101.6    1121.8   1150];
% ts_d = [1500.3 1500.7  1524.1 1524.9]
% ts_s = [-50 ts_d(2:end) dt8.Time(end-1)]
ts_s = [-50 ts_d(2:end) 2200]

% clf
% figure(202);clf
% [datatable, velocitytable] = DownSampleAndSplit(dt8, ts_d, ts_s, ML, 5, nf/67, 'ForceLength8mM');
[datatable, velocitytable] = DownSampleAndSplit(dt8, [0 4250], [-5000, 500, 2000, 4250], ML, dsf/10, nf/67, 'ForceLength8mM_all');

% subplot(211);title('Length (ML)');xlabel('Time (ms)');ylabel('ML');
% subplot(212);title('Force (kPa)');xlabel('Time (ms)');ylabel('kPa');
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
datalabel = "step-up 8 mM";
ts_d = [10, 20.6, 40.7, 62.1, 80.1, 101.3, 121.4, 141.7, 161.5, 181.7, 201.8, 221.9, 241.9, 261.9, 281.9, 300.5, 321.9, 500];
ts_s = [-500 ts_d(2:end)]
% clf;
% figure
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, dsf*10, nf/55, 'bakers_rampup8');
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, dsf*10, 1, '');
% velocitytable

%% Load stretch step-up data
% clf;
data_table = readtable('data/8 mM stretch.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [ts_d(1:end-1) 339.95], [ts_s(1:end-1) 339.95], ML, dsf*5, nf/54, 'bakers_rampup2_8');

%% slack 8 mM
% clf;
data_table = readtable('data/8 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
o = 1150 - 100 + 9.4;
ts_s = [0 1070 1159 2259 2759 3058]; % to prevent skipping events with large integrator step
ts_s = [2500, 2759.6, 2760.4, 2910.4, 2930, 3050]
% [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_s([1, end])-o, ts_s -o, ML, dsf, nf/54, 'bakers_slack8mM', o);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
% subplot(211)
% title('Slack experiment for different ATP concentrations')
% legend('8 mM', '2 mM', '0.2 mM')
%% get the ktr of the zones
zones = [1162, 1209;1464 1519;1816 1889;2269 2359.5;2774 2900];
clf;    
fitRecovery(datatable, zones);

%% 8 mM long scope data
% figure(101);clf;
% tss_d = [118555, 126800]
% tss_s = [118555, 121890, 121900, 121910,121920, ts_d(1) + (122070+710)-10, ts_d(end-1) + (122070+710), 123910, 123930, 123960, 124000, ...
%     124210, 124230, 124270, 124310, 124560, 124580, 124640, 124680, 125010, 125030, 125110, 125150, 125510, 125530, 125660, 125700, 126800]
data_table = readtable('data/8 mM ATP scope.txt', 'filetype', 'text', 'NumHeaderLines',4);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, 1, '', -(122070+710)+20);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, tss_d, tss_s, ML, 1, nf/54, 'bakers_rampup2_8_long', 0);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, 'bakers_rampup2_8_long', -(122070+710)+20);

legend('ForceLength8mM_all', 'bakers_rampup8', 'bakers_rampup2_8', 'bakers_rampup2_8_long', 'nevim', 'Interpreter', 'None')
% xlim([1800, 2000])
%% 8mM ATP
clf;
data_table = readtable('data/8mM ATP 2ktr.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, '');
fitRecovery(datatable, [100, 700;]); 
%%
data_table = readtable('data/2 mM ATP ktr.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, '');
fitRecovery(datatable, [100, 700;]); 
%%
data_table = readtable('data/0.2 mM ATP ktr.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, '');
fitRecovery(datatable, [100, 700;]); 
%% load length-force data for 2 mM
figure('2 mM');
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
% clf;
[datatable, velocitytable] = DownSampleAndSplit(dt2, ts_d, ts_s, ML, dsf, nf/65, 'ForceLength2mM');
% [datatable, velocitytable] = DownSampleAndSplit(dt2, [], [], ML, dsf, nf/65, '');
velocitytable

%% same preparation as 8mM, using the same scale
% clf;
data_table = readtable('data/2 mM stretch.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, dsf*10, nf/54, 'bakers_rampup2_2');
%% same preparation as 8mM, using the same scale
data_table = readtable('data/02 mM ATP scope.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, 'bakers_rampup2_2_long', -(122070+710)-2700);
%%
% clf;
data_table = readtable('data/0.2 mM stretch.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, dsf*10, nf/54, 'bakers_rampup2_02');
%%
data_table = readtable('data/0.2 mM ATP scope.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, 'bakers_rampup2_02_long', -(122070+710)-2700 + 1280 + 20 - 5);
%%
data_table = readtable('data/relaxed stretch.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, dsf*10, nf/54, 'bakers_rampup2_rel');

return;
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
% 
% fcn = @dPUdTCa;
% simulateTimes = velocitytable(:, 1);
% velocities = velocitytable(1:end-1, 2);
% 
% params = struct('Pi', 0, 'MgATP', 8, 'MgADP', 0, 'Ca', 1000,...
%     'Velocity', velocities, ...
%     'UseCa', false,'UseOverlap', false,'kSE', 1000, 'mu', 10,...
%     'N', 40, 'Slim', 0.025, 'ValuesInTime', 1, ...
%     'BreakingVelocity', -10, 'SlackVelocity', 10, 'SL0', 2.0, ...
%     'MatchTimeSegments', 1, 'ML', ML, 'PlotProbsOnStep', false, ...
%     'ReduceSpace', false...
%     );
% 
% params = getParams(params, g)
% % figure(1);clf;
% 
% [F out] = evaluateModel(fcn, simulateTimes, params);
% 

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
figure(9);clf; hold on;
Pus = 1 - out.p1_0 - out.p2_0 - out.p3_0;% PU substitute

leg = plot(out.t, Pus, out.t, out.p1_0, out.t, out.p2_0, out.t, out.p3_0, out.t, out.NR);
% plot([simulateTimes;simulateTimes], repmat([0; 1], [1 size(simulateTimes, 2)]))
legend(leg)
legend('Pu', 'P1', 'P2', 'P3', 'NR')

%% Animate states
animateStateProbabilities(out, params);


%% function definition
function [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, dsf, scaleF, saveAs, offset)
% ts_d - time segment data for cost function
% ts_s - time segment simulation - broke by constant velocity segments
% dsf - downsample factor

% offset in ms
if nargin < 8
    offset = 0;
end

    data_table.Properties.VariableNames = {'Time', 'L', 'F'};
    data_table.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};

    if isempty(ts_d)
        ts_d = [data_table.Time(1) data_table.Time(end)];
    end
    
    if isempty(ts_s)
        ts_s = [data_table.Time(1) data_table.Time(end)];
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
    % round to limit the oscillations
    lf = round(lf, 3);
    
    ld = downsample(lf, dsf);


    f = data_table.F(imin_s:imax_s)*scaleF;
    ff = movmean(data_table.F(imin_d:imax_d)*scaleF,[dsf/2 dsf/2]); % force filtered
    fd = downsample(ff, dsf);
    
    datatable = [td/1000  + offset/1000, ld*ML, fd];

    % Split it into segments with const velocities
    vs = [];
    pos = [data_table.L(find(t >= ts_s(1), 1))]*ML;
    for it = 1:(length(ts_s)-1)
        t1 = find(t >= ts_s(it), 1);
        t2 = find(t >= ts_s(it + 1), 1);
        vs(it) = round((l(t1) - l(t2))/(t(t1) - t(t2))*1000, 1);
        pos(it+1) = pos(it) + vs(it)*ML*(t(t2)-t(t1))/1000;
    end
    vsum = vs*ML; % 
    % time (s), velocity ML/s, velocity um/s
    velocitytable = [(ts_s + offset)/1000;[vs 0];[vsum 0];pos(1:end)]'; 
    

    if ~isempty(saveAs)
        % Export the data into modelica-readable format and for identificatoin
        fn = ['data/' saveAs '.mat'];
        save(fn,  'datatable', 'velocitytable');
        disp(['Saved as ' fn])
    end

%     figure();clf;
    subplot(211);hold on;title(saveAs, 'Interpreter', 'none');
    plot(td + offset, ld, '-');
    
    
%     plot(t + offset, l, tf + offset, lf, td + offset, ld, '|-');   
%     plot(velocitytable(:, 1)*1000, velocitytable(:, 4)/ML, 'x-', 'Linewidth', 1, 'MarkerSize', 10)
%     plot([ts_d;ts_d], repmat([min(data_table.L);max(data_table.L)], 1, length(ts_d)))
    xlabel('time (ms)');
    ylabel('Length (ML)')

    
    subplot(212);hold on;
    plot(td + offset, fd, '-', 'Linewidth', 2, 'MarkerSize', 10);
    
%     plot(t + offset, f, tf + offset, ff, td + offset, fd, '|-', 'Linewidth', 2, 'MarkerSize', 10);
%     plot([ts_d;ts_d], repmat([min(data_table.F);max(data_table.F)], 1, length(ts_d)))
    xlabel('time (ms)');
    ylabel('Force (kPa)')
    
end

