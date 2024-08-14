%% load Anthony Baker's experiments
% load gopt;
% decimation sampling (each x)
dsf = 10;
ML = 2.0;
% normalized force multiplier
nf = 56;
ts_s = []; ts_d = [];
% clear;
% close all;
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
dt2 = readtable(datafile, ...
    "filetype", 'spreadsheet', ...
    'VariableNamingRule', 'modify', ...
    'Sheet', '2 mM', ...
    'Range', 'A5:C86004');

dt2.Properties.VariableNames = {'Time', 'L', 'F'};
dt2.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};
% datalabel = "8 mM";
%%
% clf;
plot(dt8.Time, dt8.L, '--b');hold on;
plot(dt2.Time, dt2.L, '--r');
yyaxis right;
plot(dt8.Time, dt8.F, '-b');hold on;
plot(dt2.Time, dt2.F, '-r');

%%
ts_d = [450 500.25, 500.9, 509.25, 510, 519.5,539.8, 650];
% ts_s = [400 ts_d(2:end)]
% ts_d = [ts_d 950.0    1000.3    1001.0    1091.2   1092 1101.6    1121.8   1150];
% ts_d = [1500.3 1500.7  1524.1 1524.9]
% ts_s = [-50 ts_d(2:end) dt8.Time(end-1)]
ts_s = [-50 ts_d(2:end) 2200]

% clf
% figure(202);clf

% [datatable, velocitytable] = DownSampleAndSplit(dt8, [], [], ML, 10, nf/67, '');
% % gwet the velocities
% dSL = [0;diff(datatable(:, 2))]; dT = [-datatable(2, 1);diff(datatable(:, 1))]; dSLdT = dSL./dT;
% subplot(211);hold on;yyaxis right;plot(datatable(:, 1)*1000, dSLdT/2);ylim([-10, 10]);
% velocities = [-6, -1, -3, -5, -0.5, -4, -2]

% [datatable, velocitytable] = DownSampleAndSplit(dt8, ts_d, ts_s, ML, 5, nf/67, 'ForceLength8mM');
[datatable, velocitytable] = DownSampleAndSplit(dt8, [0 4250], [-5000, 500, 2000, 4250], ML, dsf/10, nf/67, 'ForceLength8mM_all');

% subplot(211);title('Length (ML)');xlabel('Time (ms)');ylabel('ML');
% subplot(212);title('Force (kPa)');xlabel('Time (ms)');ylabel('kPa');
velocitytable
%% load step-up data for 8 mM
datafile = "data/06 21 21 Ramps 2 mM ATP.xlsx";
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
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, dsf/5, nf/55, '');
% [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, dsf*10, nf/55, 'bakers_rampup8');
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, dsf*10, 1, '');
% velocitytable

%% Load stretch step-up data
clf;
data_table = readtable('data/0.2 mM stretch.txt', 'filetype', 'text', 'NumHeaderLines',4);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [ts_d(1:end-1) 339.95], [ts_s(1:end-1) 339.95], ML, dsf*5, nf/54, 'bakers_rampup2_8');
[datatable, velocitytable] = DownSampleAndSplit(data_table, [ts_d(1:end-1) 339.95], [ts_s(1:end-1) 339.95], ML, dsf, nf/54, '');


%% slack 8 mM
figure(1);
 clf;hold on;

data_table = readtable('data/8 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
o = 1150 - 100 + 9.4;
% ts_s = [0 1070 1159 2259 2759.6 2760.4 2910 2930 3058]; % to prevent skipping events with large integrator step
ts_s = [0 1070 1159 2259 2759 3058]; % to prevent skipping events with large integrator step
ts_s = [0 1070 1159]
ts_s = [0 1070 1159.4 1160.3 1209.9 1229.8 1459.4 1460.4 1519.95 1540.1 1809.4 1810.4 1889.9 1909.9 2259.4 2260.4 2360 2380.0, 2759.4, 2760.3, 2910.4, 2930.35, 3050]
% [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_s([1, end])-o, ts_s -o, ML, dsf, nf/54, 'bakers_slack8mM', o);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, 1, nf/54, 'bakers_slack8mM_all', o);
% subplot(211)
% title('Slack experiment for different ATP concentrations')
% legend('8 mM', '2 mM', '0.2 mM')
%%
figure(1);
data_table = readtable('data/2 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, 10, nf/54, 'bakers_slack2mM', o);
%%
data_table = readtable('data/0.2 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, 1, nf/54, 'bakers_slack02mM', o);


yyaxis right;
plot(datatable(:, 1)*1000, datatable(:, 2));

title('Slack experiment')
xlabel('Time (ms)')
ylabel('Sarcomere length (um)')
xlim([2250, 2500]);
set(gca, 'fontsize', 22);
xt = xticks;
% xticks(xt);
xticklabels(xt-xt(1));
legend('8mM', '2mM', '0.2mM', 'Musc. L*');
set(gca,'YColor',[0.49 0.18 0.56])
%% get the ktr of the zones
figure(4);
zones = [1162, 1209;1464 1519;1816 1889;2269 2359.5;2774 2900];
clf;    
[dSLpc, ktr, df, del, E, SL, x0lin]  = fitRecovery(datatable, zones, 0);
plot([1 3], [0 0], 'k-')
% times of start of the SL drop
dropstart = velocitytable([3, 7, 11, 15, 19], 1);
dt = x0lin' - dropstart; 
dL = 2.2 - SL;

v = dL'./dt;
%%
figure(5);clf;
nexttile;hold on;
plot(datatable(:, 1)-dropstart', datatable(:, 2));
plot([3e-4 dt'], [2.2 SL], '*-', LineWidth=2)
% plot([3e-4 dt'], [2.2 SL], '*--', LineWidth=1)
title('Time of force rebound, plotted over muscle length');
ylabel('ML (normalized to 2\mum)');xlabel('Time after slack');

xlim([-0.05, 0.25])
nexttile;hold on;
plot(datatable(:, 1)-dropstart', datatable(:, 3));
hold on; set(gca,'ColorOrderIndex',1);
% plot(datatable(:, 1)-dropstart', datatable(:, 3), '--', linewidth=2);
xlim([-0.05, 0.25])

slack_x = [3e-4 dt'] - 3e-4;
slack_y = [2.2 SL];
%% first peak and SL length dependence - is there a connection?
x_peak1_8 = [53.95 64 83.45 103 155.45]/1000;
x_peak1_2 = [58, 66.5, 85.7, 105, 164]/1000;
x_peak1 = x_peak1_2; 
hold on;
plot([x_peak1;x_peak1], x_peak1*0 +  [0;100], '--');

for i_sp = 1:5
   i_dt = find(datatable(:, 1) > x_peak1(i_sp) + dropstart(i_sp), 1, 'first');
   SL_peak1(i_sp) = datatable(i_dt, 2);
end
title('Time at first force peak')
%% switch to panel 1
plot([0;0.25] + SL_peak1*0, [SL_peak1;SL_peak1], '--');
plot(x_peak1, SL_peak1, 'x',LineWidth=2);
title('Sarcomere length at first peak')
%%
figure(3);clf;nexttile;
plot(SL, SL_peak1 - SL, '*-',SL(end-1:end), SL_peak1(end-1:end) - SL(end-1:end), 'o-', LineWidth=2);
legend('Velocities increasing from 4', 'Velocity 3')
title('Base sarcomere length to force peak1')
xlabel('base SL')
ylabel('SL at 1st force peak');
nexttile;
velocities = velocitytable([5, 9, 13, 17, 21], 2);
i_sorted = [5, 1 2 3 4];
i_sorted = 1:4;
plot(-velocities(i_sorted), SL_peak1(i_sorted) - SL(i_sorted), '*-', -velocities([1,5]), SL_peak1([1,5]) - SL([1,5]), 'o-', LineWidth=2)
legend('SL increasing from 1.92', 'SL 1.88')
xlabel('Negative velocity (SL extension) (\mum)', Interpreter='tex')
title('Ramp-up velocity to force peak1')

%%

% %% Summary graph    
% 8mM	
% X dT (ms),	Y dML (% ML)
bp8 = [2.735152222	8;
4.380311224	10;
6.466291958	12;
9.056710864	14;
12.05655714	16];

% 2mM	
% X dT (ms),	Y dML (% ML)
bp2 = [1.191034575	7.954971857;
3.025150961	9.981238274;
5.552426411	11.96998124;
7.953651442	13.9587242;
11.48933025	15.94746717];

% 02mM	
% X dT (ms),	Y dML (% ML)
bp02 = [-9.491620288	7.99249531
-6.018966686	9.981238274
-3.302497359	12.00750469
-1.090347959	13.99624765
0.2394485	15.98499062];

% figure;
clf;hold on;
plot(bp8(:, 2), bp8(:, 1),'o', 'linewidth', 2);
plot(bp2(:, 2), bp2(:, 1),'s', 'linewidth', 2);
plot(bp02(:, 2), bp02(:, 1),'d', 'linewidth', 2);

bp = bp8;


y_line = @(k, x0, x)k.*(x-x0);
% [ae be] = fit(bp(:, 2), bp(:, 1), y_line, 'StartPoint', [1, 0]);
% 
timebase = (0:0.1:20);
% plot(timebase, y_line(ae.k, ae.x0, timebase), '--', 'Linewidth', 2);
        
% y_exp = @(df, ktr, s, x)df*(1-exp(-(x-s)*ktr));        
% y_line = @(k, x0, x)k.*(x-x0);
    y_exp = @(a, b, c, x)a.*x.^2 +b.*x + c;        
    [ae be] = fit(bp(:, 2), bp(:, 1), y_exp, 'StartPoint', [1, 1, 0]);        
    plot(timebase, y_exp(ae.a, ae.b, ae.c, timebase), '--', 'Linewidth', 2);
    disp('8 mM exp');disp(ae);
    [ael bel] = fit(bp(:, 2), bp(:, 1), y_line, 'StartPoint', [1, 0]);
    plot(timebase, y_line(ael.k, ael.x0, timebase), '--', 'Linewidth', 2);

    
    bp = bp2;
    [ae be] = fit(bp(:, 2), bp(:, 1), y_exp, 'StartPoint', [1, 1, 0]);        
    plot(timebase, y_exp(ae.a, ae.b, ae.c, timebase), '--', 'Linewidth', 2);
    disp('2 mM');disp(ae);
    [ael bel] = fit(bp(:, 2), bp(:, 1), y_line, 'StartPoint', [1, 0]);
    plot(timebase, y_line(ael.k, ael.x0, timebase), '--', 'Linewidth', 2);
    
    bp = bp02;
%     y_exp = @(a, b, c, x)(a.*x-b).^2 + c;
    [ae be] = fit(bp(:, 2), bp(:, 1), y_exp, 'StartPoint', [-1, -20, 0]);        
    plot(timebase, y_exp(ae.a, ae.b, ae.c, timebase), '--', 'Linewidth', 2);
    disp('0.2 mM');disp(ae);
    [ael bel] = fit(bp(:, 2), bp(:, 1), y_line, 'StartPoint', [1, 0]);
    plot(timebase, y_line(ael.k, ael.x0, timebase), '--', 'Linewidth', 2);
    
    xlabel('Step down size (% ML)')
    ylabel('Time delay (ms)')
    
    %zero crossing with linear at stepdown of ML%
    sd = [6, 7.3, 15];
    
%% 8 mM long scope data
% figure(101);clf;
% tss_d = [118555, 126800]
% tss_s = [118555, 121890, 121900, 121910,121920, ts_d(1) + (122070+710)-10, ts_d(end-1) + (122070+710), 123910, 123930, 123960, 124000, ...
%     124210, 124230, 124270, 124310, 124560, 124580, 124640, 124680, 125010, 125030, 125110, 125150, 125510, 125530, 125660, 125700, 126800]
data_table = readtable('data/8 mM ATP scope.txt', 'filetype', 'text', 'NumHeaderLines',4);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, 1, '', -(122070+710)+20);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, tss_d, tss_s, ML, 1, nf/54, 'bakers_rampup2_8_long', 0);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, 'bakers_rampup2_8_long', -(122070+710)+20);
datatable8s = datatable;
% legend('ForceLength8mM_all', 'bakers_rampup8', 'bakers_rampup2_8', 'bakers_rampup2_8_long', 'nevim', 'Interpreter', 'None')
% xlim([1800, 2000])
%% 8mM ATP
% clf;
data_table = readtable('data/8mM ATP 2ktr.txt', 'filetype', 'text', 'NumHeaderLines',4);
% match to the scope
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, 'bakers_ktr_8', -870-85 + 4.5);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, 'bakers_ktr_8', -100);

% fitRecovery(datatable, [100, 700;], 0); 
%%
data_table = readtable('data/2 mM ATP ktr.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, '', -100);
% fitRecovery(datatable, [100, 700;],0); 
%%
data_table = readtable('data/0.2 mM ATP ktr.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, '', -100);
% fitRecovery(datatable, [100, 700;],0); 
%% load length-force data for 2 mM
% figure(2);
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
% [datatable, velocitytable] = DownSampleAndSplit(dt2, ts_d, ts_s, ML, dsf, nf/65, 'ForceLength2mM');
[datatable, velocitytable] = DownSampleAndSplit(dt2, [], [], ML, dsf, nf/67, '');
velocitytable

%% same preparation as 8mM, using the same scale
clf;
data_table = readtable('data/8 mM stretch.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, 'bakers_rampup2_2');
%% same preparation as 8mM, using the same scale
data_table = readtable('data/02 mM ATP scope.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, 'bakers_rampup2_2_long', -(122070+710)-2700);
datatable2s = datatable;
%%
% clf;
data_table = readtable('data/0.2 mM stretch.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, dsf*10, nf/54, 'bakers_rampup2_02');
%%
legend('ATP 8mM', 'ATP 2mM', 'ATP 0.2mM')
set(gca,'fontsize',16);
%%
title('Active force ramp-up')
ylabel('Muscle length(um)')
figure(1);
yyaxis right;
plot(datatable(:, 1)*1000, datatable(:, 2), 'r--');

%%
data_table = readtable('data/0.2 mM ATP scope.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, 'bakers_rampup2_02_long', -(122070+710)-2700 + 1280 + 20 - 5);
datatable02s = datatable;
%% Show timepoints for ??? reasons
timepoints = [-1020, 0, 3300:4100];
% indx = [];
for i = 1:length(timepoints)
    indx_1(i) = find(datatable(:, 1) >= timepoints(i)/1000, 1);
end
plot(datatable(indx_1, 1)*1000, datatable(indx_1, 3), '*', 'Linewidth', 2)

timepoints = [800, 1000, 1305, 1664, 2105, 2613];
% indx = [];
for i = 1:length(timepoints)
    indx_1_1(i) = find(datatable(:, 1) >= timepoints(i)/1000, 1);
end
plot(datatable(indx_1_1, 1)*1000, datatable(indx_1_1, 3), '^', 'Linewidth', 2)

%% passive stretch

data_table = readtable('data/relaxed stretch.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, dsf, nf/54, 'bakers_rampup2_rel');

return;
%% plot active scope
% Atp 0.2 is shifted. Extend before 0 and after 1000ms by x
% datatable8s = datatable;
i0 = find(datatable8s(:, 1) > 0, 1);
datatable8s(1:i0, 1) = datatable8s(1:i0, 1) - 0.14;
datatable2s(1:i0, 1) = datatable2s(1:i0, 1) - 0.14;
% datatable02s; stays
% Extend after 800ms
datatable02s = datatable;
i0 = find(datatable02s(:, 1) > 0.8, 1);
% 8 and 2 stays
datatable02s(i0:end, 1) = datatable02s(i0:end, 1) + 0.13;

%%
clf;subplot(211);
plot(datatable8s(1:end, 1), datatable8s(1:end, 2)/ML, 'Linewidth', 2);
xlim([-1.2, 4]);
xlabel('Time'); ylabel('Muscle length (-)')
set(gca,'fontsize',16);
subplot(212);hold on;
plot(datatable8s(1:end, 1), datatable8s(1:end, 3),    datatable2s(:, 1), datatable2s(:, 3), datatable02s(:, 1), datatable02s(:, 3), 'Linewidth', 2);
xlim([-1.2, 4]);
xlabel('Time'); ylabel('Tension (kPa)');
set(gca,'fontsize',16);
legend('8 mM ATP', '2 mM ATP','0.2 mM ATP');

%% New stretch experiments


dsf = 10;
ML = 2.0;
% normalized force multiplier
nf = 56;
ts_s = []; ts_d = [];

figure(12);clf;
subplot(121);hold on;xlabel('ML');ylabel('Tension');title('Tension on Muscle length')
subplot(122);hold on;xlabel('ML (um)');ylabel('SL (um)');title('SL on Muscle length')
figure(13);clf;hold on;title('Semilog axis');ylabel('log (Tension)');xlabel('time');
figure(14);clf;hold on;
figure(11);clf;
%%
el = 20; % experiment length in ms
data_table = readtable('data/20ms_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):2000+el*2 + el*40];
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):200000];
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, 10, 1, 'bakers_passiveStretch_20ms');

figure(12); 
subplot(121); plot(datatable(:, 2), datatable(:, 3), 'o-');
subplot(122); plot(datatable(:, 2), datatable(:, 4), 'o-');
figure(13);semilogy(datatable(:, 1), datatable(:, 3), 'LInewidth', 2);
figure(14);plot(datatable(1:end-1, 1), diff(abs(log10(datatable(:, 3)))), 'LInewidth', 2);
figure(11);

%%
el = 100; % experiment length in ms
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):2000+el*2 + el*40];
ts_d = [-5000, 0:10:2000, 2000:ceil(el/200):2000 + el*2, 2000 + el*2:ceil(el/20):200000];
data_table = readtable('data/100ms_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s, ML, 1, 1, '');
% [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, 10, 1, 'bakers_passiveStretch_100ms');
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, 10, 1, '');

figure(12); 
subplot(121); plot(datatable(:, 1), datatable(:, 3), 'o-');
subplot(122); plot(datatable(:, 2), datatable(:, 4), 'o-');
figure(13);semilogy(datatable(:, 1), datatable(:, 3), 'LInewidth', 2);
figure(14);plot(datatable(1:end-1, 1), diff(abs(log10(datatable(:, 3)))), 'LInewidth', 2);
figure(11);
%%
ML = 2.0;
el = 1000; % experiment length in ms
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):2000+el*2 + el*40];
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):200000];
ts_s = [0, 200000];
% ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el, 2000 + el:ceil(el/100):2000 + el*2, 2000 + el*2:ceil(el/2):200000];
data_table = readtable('data/1s_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, 5, 1, 'bakers_passiveStretch_1000ms');
% explanatory figure
figure(10);clf;hold on;
plot(datatable(:, 1), datatable(:, 2), '-', 'Linewidth', 2); title({"Stretching ramp", "protocol"});
xlabel('Time (s)');ylabel('Sarcomere length (um)');
plot([2 2], [1.5 2.5], 'r--');
plot([3 3], [1.5 2.5], 'r--');
plot([4.95,4.75],[2.44 2.36], 'k', 'Linewidth', 4)
plot([5.2, 5.00],[2.444 2.36], 'k', 'Linewidth', 4)
xlim([0, 6]);
set(gca,'fontsize',16);



figure(12); 
subplot(121); plot(datatable(:, 2), datatable(:, 3), 'o-');
subplot(122); plot(datatable(:, 2), datatable(:, 4), 'o-');
figure(13);semilogy(datatable(:, 1), datatable(:, 3), 'LInewidth', 2);
figure(14);plot(datatable(1:end-1, 1), diff(abs(log10(datatable(:, 3)))), 'LInewidth', 2);
figure(11);
%%
el = 10000; % experiment length in ms
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):2000+el*2 + el*15];
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):200000];
data_table = readtable('data/10s_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, 1, 1, 'bakers_passiveStretch_10000ms');
% ts_d = [-5000, 0:10:2000, 2000:ceil(el/200):2000 + el*2, 2000 + el*2:ceil(el/20):200000];
% [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, 1, 1, '');


figure(12); 
subplot(121); plot(datatable(:, 1), datatable(:, 3), 'o-');
subplot(122); plot(datatable(:, 2), datatable(:, 4), 'o-');
figure(13);semilogy(datatable(:, 1), datatable(:, 3), 'LInewidth', 2);
figure(14);plot(datatable(1:end-1, 1), diff(abs(log10(datatable(:, 3)))), 'LInewidth', 2);
figure(11);
%%
el = 100000; % experiment length in ms
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el, 2000+el:ceil(el/100):2000 + el*1.1,2000 + el*1.1:ceil(el/10):2000+el*2 + el*0.9];
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):200000];
data_table = readtable('data/100s_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, [], ML, 40, 1, 'bakers_passiveStretch_100000ms');
slowest = datatable(:, 3);

figure(12); 
subplot(121); plot(datatable(:, 1), datatable(:, 3), 'o-');
subplot(122); plot(datatable(:, 2), datatable(:, 4), 'o-');
figure(13);semilogy(datatable(:, 1), datatable(:, 3), 'LInewidth', 2);
figure(14);plot(datatable(1:end-1, 1), diff(abs(log10(datatable(:, 3)))), 'LInewidth', 2);
figure(11);

figure(12); 
subplot(121); legend('20 ms', '100 ms', '1000 ms', '10 s', '100 s', 'location', 'Northwest')

%% Fit decay of the stretch data - loop for any param
% start with the slowest one


% ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el):200000];
rd = 1; % experiment length in s
% ts_d = [-5000, 0:500:2000, 2000:ceil(el/100):2000 + el*2, 2000 + el*2:ceil(el/10):200000];
% data_table = readtable('data/20ms_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
% % [datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s, ML, 1, 1, '');
% [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, [], ML, 5, 1, '');

% rds = [0.02, 0.1, 1, 10, 100];
% new ramp-ups with Ca
rds = [0.1, 1, 10];
colors = colormap(lines(length(rds)));
ss_rmse = [];
sa = [0:0.1:2.1];
% sa = [0.1:0.005:0.16];
for ss = sa
%% Loop for all ramps only
figure(1);clf;
rmse = [];
    for rd_i = 1:length(rds)
    rd = rds(rd_i);
% datatable = load(['data/bakers_passiveStretch_' num2str(rd*1000) 'ms.mat']).datatable;
data_table = readtable(['Data/PassiveCa_1/bakers_passiveStretch_pCa11_' num2str(1000*rds(rd_i)) 'ms.csv']);
datatable = table2array(data_table);

[~, i_peak] = max(datatable(:, 3));
dataplotpoints = 1:1:length(datatable(:, 1));
zone = i_peak:length(datatable(:, 1));
timebase = datatable(zone, 1) - datatable(zone(1), 1);
extrap_time = [timebase; timebase(end) + (1:10:(60*30))'];
fbase = datatable(zone, 3);
% this fits the fastest one pretty well
% y_dec = @(a, b, c, d, e, x)a*x.^(-b) + c + 0*a*b*c*d*e;
% exploring the steady state param
y_dec = @(a, b, c, d, e, x)a*x.^(-b) +ss+ 0*c*d*e;
% exploring the exponent param with constsnt offset
% y_dec = @(a, b, c, d, e, x)a*x.^(-b) +2.04 + 0*d*exp(-x*e)+ 0*a*b*c*d*e;


[ae be] = fit(timebase(2:end), fbase(2:end), y_dec, 'StartPoint', [1, 1, 3.750, 0, 1]);
rmse(rd_i) = be.rmse;
fitparam.a(rd_i) = ae.a;fitparam.b(rd_i) = ae.b;fitparam.c(rd_i) = ae.c;fitparam.d(rd_i) = ae.d;fitparam.e(rd_i) = ae.e;
%%
% t = datatable(dataplotpoints, 1);
% f = datatable(dataplotpoints, 3);
% plot data
% clf;
% hold on;
% plot(t, f);
% loglog(t, f-2.04);
% semilog(t, f);
%%
% plot(datatable(dataplotpoints, 1), datatable(dataplotpoints, 3), 'Color', colors(rd_i, :));hold on;
semilogx(datatable(dataplotpoints, 1), datatable(dataplotpoints, 3), 'Color', colors(rd_i, :));hold on;
% plot fit
semilogx(extrap_time+datatable(zone(1), 1), y_dec(ae.a, ae.b,ae.c, ae.d,ae.e, extrap_time), '--','Linewidth', 4, 'Color', max(colors(rd_i, :) - [0.2 0.2 0.2], [0 0 0]));

% clf;loglog(timebase, fbase-2.04);hold on;
%% redo as differential equation
a = fitparam.a(rd_i);
b = fitparam.b(rd_i);
c = 2.04;

dxdt = @(t, x) -b*(x-c)/t;
dxd =  @(~, y) -b/(a^(1/b))*(y-c)^(1 + 1/b); 
% dxdt = -b*x/t;
% clf;
% semilogx(datatable(dataplotpoints, 1), datatable(dataplotpoints, 3), 'Color', colors(rd_i, :));hold on;


% for m = 2.8:0.1:3.2
    % m = 3;
x0 = a*extrap_time(2)^(-b) + c;
[t, x] = ode15s(dxdt, [extrap_time(2), extrap_time(end)], x0);
[t2, x2] = ode15s(dxd, [extrap_time(2), extrap_time(end)], x0);

% semilogx(t + 2 + rds(rd_i), x, 'k:')
% semilogx(t2 + 2 + rds(rd_i), x2, 'k--', LineWidth=2);

% end


end
% legend(...
%     '0.02s data', sprintf('0.02s fit, rmse %0.2f', rmse(1)),...
%     '0.1s data', sprintf('0.1s fit, rmse %0.2f', rmse(2)),...
%     '1s data', sprintf('1s fit, rmse %0.2f', rmse(3)),...
%     '10s data', sprintf('10s fit, rmse %0.2f', rmse(4)),...
%     '100s data', sprintf('100s fit, rmse %0.2f', rmse(5))...
%     );
legend(...
    '0.1s data', sprintf('0.1s fit, rmse %0.2f', rmse(1)),...
    '1s data', sprintf('1s fit, rmse %0.2f', rmse(2)),...
    '10s data', sprintf('10s fit, rmse %0.2f', rmse(3))...
    );

title(sprintf('OVerall fit: %0.2f', sum(rmse)));
figure(2);clf;
% semilogx(rds, fitparam.a, 'o-',rds, fitparam.b, 'o-',rds, fitparam.c, 'o-');
semilogx(rds, fitparam.a, 'o-',rds, fitparam.b, 'o-',rds, fitparam.d, 'o-',rds, fitparam.e, 'o-');
legend('Optimized a', 'Optimized b', 'Optimized c', 'd', 'e')
xlabel('Ramp duration');
title('Optimized parameter values for different ramp durations');
ss_rmse = [ss_rmse sum(rmse)];
% ss
end
figure(3);clf;
plot(sa(1:length(ss_rmse)), ss_rmse, 'x-', 'LineWidth',2);
% fix the 'a' param
% [ae be] = fit(timebase, fbase, @(b, x)y_dec(0.95, b, x), 'StartPoint', [0.05]);
% plot(datatable(zone, 1)*1000, y_dec(0.95, ae.b, timebase) + datatable(end, 3), 'Linewidth', 2)
%% Fit the pCa4 passive experiments
datatable = readtable('data/pCa4_008.txt', 'filetype', 'text', 'NumHeaderLines',4);
% datatable = readtable('data/Relax_018.txt', 'filetype', 'text', 'NumHeaderLines',4);
datatable.Properties.VariableNames = {'t', 'L','F'};

t_rd = [0.1, 1, 10, 100]; % ramp duration time

% ramp start time
% for pCa4_008
t_rs = fliplr([186.32, 430.73, 585.82, 730.81]);
% for Relax_018
% t_rs = fliplr([26, 270.45, 425.4, 570.42]);
t_dd = ones(1, 4)*30; % decay duration + ramp-down duration
d_rd = 10; % duration ramp-down
d_zh = 10; % duration zero-hold
t_re = t_rs + t_rd; % ramp end time
t_rec = t_re + 30; % end of recover

% zero-offset for pCa4_008
t0 = [0, 45, fliplr(t_rs) + fliplr(t_rd) + t_dd + d_rd];
% zero-offset for Relax_018
% t0 = fliplr(t_rs) + fliplr(t_rd) + t_dd + d_rd;
i_fzero = any(datatable.t > t0 & datatable.t < t0 + d_zh, 2);
t_fzero = datatable.t(i_fzero);
f_fzero = datatable.F(i_fzero);

% fit a func
y_f0 = @(a, b, x)a*x + b;
    
[ae be] = fit(t_fzero, f_fzero, y_f0, 'StartPoint', [.001, 2]);

% zero level drift
f_0 = y_f0(ae.a, ae.b, datatable.t);
datatable.F = datatable.F - f_0;

figure(1);clf;
subplot(211);plot(datatable.t, datatable.L);hold on;
plot(datatable.t(i_fzero), datatable.L(i_fzero), '.');

legend('Length (L0)');
subplot(212);hold on;
plot(datatable.t, datatable.F+f_0);
plot(datatable.t(i_fzero), datatable.F(i_fzero)+f_0(i_fzero));
plot(datatable.t, f_0, 'Linewidth', 2);
plot(datatable.t, datatable.F);legend('Tension raw', 'zeros', 'zero drift approx', 'Tension clear') ;

xlabel('Time (s)')
colors = colormap(lines(length(t_rs)));

pca = '4'; % pCa level
figure(2);clf;
figure(3);clf;
for i_ramp = 1:length(t_rs)
   %% 
    saveZone = find(datatable.t >= t_rs(i_ramp) - 2 & datatable.t < t_rec(i_ramp)); % time span to save wioth 2s steady state prior to ramp-up
    el = t_rd(i_ramp)*1000;% experiment length in ms
    ts_d = [0:200:2000, 2000:ceil(el/20):2000 + el*2, 2000 + (el*2:500:el+30000)];
    % saving as current ramp's datatable
    
    datatable_cur = datatable(saveZone, :);
    datatable_cur.t = datatable_cur.t - datatable.t(saveZone(1));
    % datatable_cur.F = datatable_cur.F*1.5; % adjust the peak heights
    filename = ['bakers_passiveStretch_pCa' num2str(pca) '_' num2str(t_rd(i_ramp)*1000) 'ms'];
    % writetable(datatable_cur, ['data/' filename '.csv']);

    DownSampleAndSplit(datatable_cur, ts_d/1000, [], 2.0, 1, 1.5, filename, 0, 1);

    rampZone = find(datatable.t >= t_rs(i_ramp) & datatable.t < t_rec(i_ramp));
    peakZone = find(datatable.t >= t_re(i_ramp) & datatable.t < t_rec(i_ramp));
    timebase = datatable.t(rampZone) - datatable.t(rampZone(1));

    peakTimebase = datatable.t(peakZone) - datatable.t(peakZone(1));
    % extrap_time = [timebase; timebase(end) + (1:10:(60*30))'];
    fbase = datatable.F(rampZone);
    peakFbase = datatable.F(peakZone);
    
    % fit
    % same function as without PNB, does not fit well
    % y_dec = @(a, b, c, d, e, x)a*x.^(-b) + 3.3 + 2.04 + 0*a*b*c*d*e;
    % not a great fit using the params identified by the resting experiment
    % y_dec = @(a, b, c, d, e, x)fitparam.a(i_ramp+1)*x.^(-fitparam.b(i_ramp+1)) + 2.04 + 3.3 + d*exp(-x*e) + 0*a*b*c*d*e;
    % y_dec = @(a, b, c, d, e, x)fitparam.a(i_ramp+1)*x.^(-fitparam.b(i_ramp+1)) + 2.04 + 3.3 + d*x.^(-e) + 0*a*b*c*d*e;
    % best fit, with clamped offset
    offset =  2.04;
    y_dec = @(a, b, c, d, e, x)a*x.^(-b) + offset + d*exp(-x*e) + 0*a*b*c*d*e;
    
    [ae be] = fit(peakTimebase(2:end), peakFbase(2:end), y_dec, 'StartPoint', [1, 1, 0, 0.1, 1], 'Lower', [0 0 1, 0, 0], 'Upper', [10, 1, 5, 10, 1]);
    rmse(i_ramp) = be.rmse;
    fitparam.pnba(i_ramp) = ae.a;fitparam.pnbb(i_ramp) = ae.b;fitparam.pnbc(i_ramp) = ae.c;fitparam.pnbd(i_ramp) = ae.d;fitparam.pnbe(i_ramp) = ae.e;
    
    % plot data  - decay only
    figure(2);   
    semilogx(peakTimebase, peakFbase, 'Color', colors(i_ramp, :));hold on;
    % semilogx(timebase, fbase, 'Color', colors(i_ramp, :));hold on;

    % plot fit
    % semilogx(peakTimebase + datatable.t(peakZone(1)) - datatable.t(rampZone(1)), y_dec(ae.a, ae.b,ae.c, ae.d,ae.e, peakTimebase), '--','Linewidth', 4, 'Color', max(0, colors(i_ramp, :)-0.2));
    semilogx(peakTimebase, y_dec(ae.a, ae.b,ae.c, ae.d,ae.e, peakTimebase), '--','Linewidth', 4, 'Color', max(0, colors(i_ramp, :)-0.2));

    figure(3);
    semilogx(datatable_cur.t, datatable_cur.F);hold on;
    semilogx(peakTimebase + 2 + t_rd(i_ramp), y_dec(ae.a, ae.b,ae.c, ae.d,ae.e, peakTimebase), '--', 'Color', colors(i_ramp, :), 'linewidth', 2);hold on;

end
legend('0.1s data', sprintf('0.1s fit, rmse %0.2f', rmse(1)),...
    '1s data', sprintf('1s fit, rmse %0.2f', rmse(2)),...
    '10s data', sprintf('10s fit, rmse %0.2f', rmse(3)),...
    '100s data', sprintf('100s fit, rmse %0.2f', rmse(4))...
    );
title(sprintf('PNB overall fit: %0.2f', sum(rmse)));
xlabel('time (s)');ylabel('Tension (kPa)');
sum(rmse)

%% plot param fit
figure(3);clf;
loglog(rds, fitparam.a, rds, fitparam.b, Marker="s", LineWidth=2);hold on;
set(gca, "colororderindex", 1);
loglog(t_rd, fitparam.pnba, t_rd, fitparam.pnbb, t_rd, fitparam.pnbc, t_rd, fitparam.pnbd, t_rd, fitparam.pnbe, Marker="o", LineStyle="--", LineWidth=2);hold on;
legend('Resting a', 'Resting b','PNB a','PNB b','PNB c','PNB d','PNB e');
xlabel('Ramp duration');title('Values of identified parameters');

figure(4);clf;
for i_rd =1:length(t_rd)
    rd = t_rd(i_rd);
    t = rd:rd:min(200, rd*10);
    f_pow = y_dec(fitparam.pnba(i_rd), fitparam.pnbb(i_rd), 0, 0, 0, max(1e-2, t-rd)) - offset;
    f_exp = y_dec(0, 0, 0, fitparam.pnbd(i_rd), fitparam.pnbe(i_rd), max(1e-1, t-rd)) - offset;
    % area(repmat(t+rd, [3, 1]), [f_exp; f_pow; repmat(offset, [1, length(t)])]);
    
    
    area(t, offset + f_pow + f_exp, FaceColor=colors(1, :));
    hold on;
    plot(t, offset + f_pow + f_exp, 'kv--');
    area(t, offset + f_pow, FaceColor=colors(2, :));
    area(t, repmat(offset, [length(t), 1]), FaceColor=colors(3, :));
    plot([rd, rd], [0,offset + f_pow(1) + f_exp(1)], 'k-', LineWidth=2);
end
set (gca, 'Xscale', 'log');
title('Composition of decay response');
legend('d.exp(-e.x)', 'a.x^{-b}', 'const')
%% load the pCa experiments and save the peaks
pca = [Inf, 4];
t_rd = [0.1, 1, 10, 100];
peaktable = table();
i = 1;
figure(45);clf;hold on;
% figure(44);clf;hold on;
for i_pca = 1:length(pca)
    for i_ramp = 1:length(t_rd)
        filename = ['data/bakers_passiveStretch_pCa' num2str(pca(i_pca)) '_' num2str(t_rd(i_ramp)*1000) 'ms.csv'];
        datatable = readtable(filename);
        [peak, i_peak] = max(datatable.F);
        mm = movmean(datatable.F, [32, 4]);
        ss = mm(datatable.t > datatable.t(i_peak) + 29);
             % RampDuration,Peak,SteadyState
        peaktable(i, :) = {pca(i_pca), t_rd(i_ramp), peak, ss(1)};
        i = i+1;
        plot(datatable.t, datatable.F);
        plot(datatable.t, mm, 'Linewidth', 2);
        plot(datatable.t(i_peak), datatable.F(i_peak), '*', 'Linewidth', 2);
        plot(datatable.t(i_peak) + 29, ss(1), 'd', 'Linewidth', 2);
    end
end
figure(44);clf;
peaktable.Properties.VariableNames  = {'pCa', 'RampDuration','Peak','SteadyState'};
semilogx(peaktable.RampDuration(1:4), peaktable.Peak(1:4), 'o-', 'LineWidth',2);
hold on;
semilogx(peaktable.RampDuration(5:8), peaktable.Peak(5:8), 'd-', 'LineWidth',2);
semilogx(peaktable.RampDuration(1:4), peaktable.SteadyState(1:4), 'o:', 'LineWidth',2);
semilogx(peaktable.RampDuration(5:8), peaktable.SteadyState(5:8), 'd:', 'LineWidth',2);

% compare with the higher step protocol
% peakTable = readtable("data/bakers_passiveStretch_Peaks.csv");
% semilogx(peakTable.RampDuration, peakTable.Peak, 'x--', peakTable.RampDuration, peakTable.SteadyState, 'x--', 'LineWidth',2)

title('Maximal peaks by ramp duration (pCa protocol)')
legend('Resting (no Ca)', 'pCa4', '30s after peak (no Ca)','30s after peak (pCa4)')

xlabel('ramp duration (s)');ylabel('Tension (kPa)');
%% plot semilog

ts_d = [-5000, 0:500:2000, 2000:ceil(el/100):2000 + el*2, 2000 + el*2:ceil(el/10):200000];
data_table = readtable('data/100s_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s, ML, 1, 1, '');
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, [], ML, 5, 1, '');
datatable_slow = datatable;

ts_d = [-5000, 0:500:2000, 2000:ceil(el/1):2000 + el*3, 2000 + el*3:ceil(el/10):200000];
data_table = readtable('data/100ms_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s, ML, 1, 1, '');
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, [], ML, 5, 1, '');
datatable_fast = datatable;
%%
figure(2); clf;hold on;
plot(datatable_slow(:, 1), log10(datatable_slow(:, 3)));
plot(datatable_fast(:, 1), log10(datatable_fast(:, 3)));
semilogy(datatable_slow(:, 1), (datatable_slow(:, 3)));
semilogy(datatable_fast(:, 1), (datatable_fast(:, 3)));

%% Fit the fast one
clf;
el = 100; % experiment length in ms
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*3, 2000 + el*3:ceil(el):200000];
data_table = readtable('data/100ms_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s, ML, 1, 1, '');
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, [], ML, 5, 1, '');
[~, i_peak] = max(datatable(:, 3));
zone = i_peak:length(datatable(:, 1));
timebase = datatable(zone, 1) - datatable(zone(1), 1);
fbase = datatable(zone, 3) - datatable(end, 3);

y_dec = @(a, b, a2, b2, x)a*exp(-max(b, 0).*x) + (a2)*exp(-max(b2, 0).^x);


% all params free
% [ae be] = fit(timebase, fbase, @(a, b, a2, b2,x)y_dec(a, b, a2, b2, x), 'StartPoint', [3, 0.05, 5, 0.5]);
% plot(datatable(zone, 1)*1000, y_dec(ae.a, ae.b, ae.a2, ae.b2, timebase) + datatable(end, 3),'o-', 'Linewidth', 2)

% fix the slow decay
% got from the slow ramp
% a = ae.a;b = ae.b;
[ae be] = fit(timebase, fbase, @(a2, b2,x)y_dec(a, b, a2, b2, x), 'StartPoint', [0.05, 0.5]);
% together
plot(datatable(zone, 1)*1000, y_dec(a, b, ae.a2, ae.b2, timebase) + datatable(end, 3), 'o-', 'Linewidth', 2);
% slow
plot(datatable(zone, 1)*1000, y_dec(a, b, 0, 0, timebase) + datatable(end, 3), '--', 'Linewidth', 1);
% fast
plot(datatable(zone, 1)*1000, y_dec(0, 0, ae.a2, ae.b2, timebase) + datatable(end, 3), '--', 'Linewidth', 1);


ae

%% Animate states
% animateStateProbabilities(out, params);