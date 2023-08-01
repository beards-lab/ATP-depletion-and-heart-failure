%% load Anthony Baker's experiments
% load gopt;
% decimation sampling (each x)
dsf = 10;
ML = 2.0;
% normalized force multiplier
nf = 56;
ts_s = []; ts_d = [];
clear;
close all;
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
 % clf;hold on;

data_table = readtable('data/0.2 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
o = 1150 - 100 + 9.4;
ts_s = [0 1070 1159 2259 2759.6 2760.4 2910 2930 3058]; % to prevent skipping events with large integrator step
% ts_s = [2500, 2759.6, 2760.4, 2910.4, 2930, 3050]
% [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_s([1, end])-o, ts_s -o, ML, dsf, nf/54, 'bakers_slack8mM', o);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, 10, nf/54, 'bakers_slack8mM', o);
% subplot(211)
% title('Slack experiment for different ATP concentrations')
% legend('8 mM', '2 mM', '0.2 mM')

data_table = readtable('data/2 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, 10, nf/54, 'bakers_slack2mM', o);

data_table = readtable('data/0.2 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, 10, nf/54, 'bakers_slack02mM', o);


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
zones = [1162, 1209;1464 1519;1816 1889;2269 2359.5;2774 2900];
clf;    
fitRecovery(datatable, zones, 0);

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
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, '', -870-85 + 4.5);
% fitRecovery(datatable, [100, 700;],0); 
%%
data_table = readtable('data/0.2 mM ATP ktr.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, '', -870-85 + 4.5);
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
% clf;
data_table = readtable('data/2 mM stretch.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, dsf*10, nf/54, 'bakers_rampup2_2');
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

el = 1000; % experiment length in ms
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):2000+el*2 + el*40];
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):200000];
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
rd = 0.02; % experiment length in s
% ts_d = [-5000, 0:500:2000, 2000:ceil(el/100):2000 + el*2, 2000 + el*2:ceil(el/10):200000];
% data_table = readtable('data/20ms_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
% % [datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s, ML, 1, 1, '');
% [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, [], ML, 5, 1, '');

rds = [0.02, 0.1, 1, 10, 100];
colors = colormap(lines(length(rds)));
ss_rmse = [];
% sa = [2:0.01:2.1];
sa = [0.1:0.005:0.16];
for ss = sa
%% Loop for all ramps only
figure(1);clf;
rmse = [];
    for rd_i = 1:length(rds)
    rd = rds(rd_i);
datatable = load(['data/bakers_passiveStretch_' num2str(rd*1000) 'ms.mat']).datatable;

[~, i_peak] = max(datatable(:, 3));
dataplotpoints = 1:1:length(datatable(:, 1));
zone = i_peak:length(datatable(:, 1));
timebase = datatable(zone, 1) - datatable(zone(1), 1);
extrap_time = [timebase; timebase(end) + (1:10:(60*30))'];
fbase = datatable(zone, 3);
% this fits the fastest one pretty well
% y_dec = @(a, b, c, d, e, x)a*x.^(-b) + c + 0*a*b*c*d*e;
% exploring the steady state param
% y_dec = @(a, b, c, d, e, x)a*x.^(-b) +ss+ 0*c*d*e;
% exploring the exponent param with constsnt offset
y_dec = @(a, b, c, d, e, x)a*x.^(-b) +2.04 + 0*d*exp(-x*e)+ 0*a*b*c*d*e;


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
semilogx(extrap_time+datatable(zone(1), 1), y_dec(ae.a, ae.b,ae.c, ae.d,ae.e, extrap_time), '--','Linewidth', 4, 'Color', colors(rd_i, :));

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
legend('0.02s data', sprintf('0.02s fit, rmse %0.2f', rmse(1)),...
    '0.1s data', sprintf('0.1s fit, rmse %0.2f', rmse(2)),...
    '1s data', sprintf('1s fit, rmse %0.2f', rmse(3)),...
    '10s data', sprintf('10s fit, rmse %0.2f', rmse(4)),...
    '100s data', sprintf('100s fit, rmse %0.2f', rmse(5))...
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
plot(sa, ss_rmse, 'x-', 'LineWidth',2);
% fix the 'a' param
% [ae be] = fit(timebase, fbase, @(b, x)y_dec(0.95, b, x), 'StartPoint', [0.05]);
% plot(datatable(zone, 1)*1000, y_dec(0.95, ae.b, timebase) + datatable(end, 3), 'Linewidth', 2)
%% Fit the pCa4 passive experiments
datatable = readtable('data/pCa4_008.txt', 'filetype', 'text', 'NumHeaderLines',4);
datatable.Properties.VariableNames = {'t', 'L','F'};
t_rs = fliplr([186.32, 430.73, 585.82, 730.81]); % ramp start time
t_rd = [0.1, 1, 10, 100]; % ramp duration time
t_re = t_rs + t_rd; % ramp end time
t_rec = t_re + 30; % end of recover

figure(1);clf;
plot(datatable.t, datatable.F, datatable.t, datatable.L*10)
colors = colormap(lines(length(t_rs)));
%%
figure(2);clf;
for i_ramp = 1:length(t_rs)
    rampZone = find(datatable.t >= t_rs(i_ramp) & datatable.t < t_rec(i_ramp));
    peakZone = find(datatable.t >= t_re(i_ramp) & datatable.t < t_rec(i_ramp));
    timebase = datatable.t(rampZone) - datatable.t(rampZone(1));
    peakTimebase = datatable.t(peakZone) - datatable.t(peakZone(1));
    % extrap_time = [timebase; timebase(end) + (1:10:(60*30))'];
    fbase = datatable.F(rampZone);
    peakFbase = datatable.F(peakZone);
    
    % plot data
    semilogx(peakTimebase, peakFbase, 'Color', colors(i_ramp, :));hold on;
    % semilogx(timebase, fbase, 'Color', colors(i_ramp, :));hold on;

    % fit
    % same function as without PNB, does not fit well
    % y_dec = @(a, b, c, d, e, x)a*x.^(-b) + 3.3 + 2.04 + 0*a*b*c*d*e;
    % not a great fit using the params identified by the resting experiment
    % y_dec = @(a, b, c, d, e, x)fitparam.a(i_ramp+1)*x.^(-fitparam.b(i_ramp+1)) + 2.04 + 3.3 + d*exp(-x*e) + 0*a*b*c*d*e;
    % y_dec = @(a, b, c, d, e, x)fitparam.a(i_ramp+1)*x.^(-fitparam.b(i_ramp+1)) + 2.04 + 3.3 + d*x.^(-e) + 0*a*b*c*d*e;
    % best fit, with clamped offset
    offset =  3.3 + 2.04 ;
    y_dec = @(a, b, c, d, e, x)a*x.^(-b) + offset+ d*exp(-x*e) + 0*a*b*c*d*e;
    
    [ae be] = fit(peakTimebase(2:end), peakFbase(2:end), y_dec, 'StartPoint', [1, 1, 0, 0.1, 1], 'Lower', [0 0 1, 0, 0], 'Upper', [10, 1, 5, 10, 1]);
    rmse(i_ramp) = be.rmse;
    fitparam.pnba(i_ramp) = ae.a;fitparam.pnbb(i_ramp) = ae.b;fitparam.pnbc(i_ramp) = ae.c;fitparam.pnbd(i_ramp) = ae.d;fitparam.pnbe(i_ramp) = ae.e;
    
    % plot fit
    % semilogx(peakTimebase + datatable.t(peakZone(1)) - datatable.t(rampZone(1)), y_dec(ae.a, ae.b,ae.c, ae.d,ae.e, peakTimebase), '--','Linewidth', 4, 'Color', max(0, colors(i_ramp, :)-0.2));
    semilogx(peakTimebase, y_dec(ae.a, ae.b,ae.c, ae.d,ae.e, peakTimebase), '--','Linewidth', 4, 'Color', max(0, colors(i_ramp, :)-0.2));
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


%% function definition
function [datatable, velocitytable] = DownSampleAndSplit(data_table, dwnsmpl, ts_s, ML, dsf, scaleF, saveAs, offset)
% df / ts_d - empty: no downsampling. 
%             Scalar: downsampling by said scalar
%             vector: time segment data points for cost function, filtered by dsf.
%
% ts_s - time segment simulation - broke by constant velocity segments
% dsf - data smoothing factor (if ts_d is a two-point vector) factor

% offset in ms
if nargin < 8
    offset = 0;
end

    if length(data_table.Properties.VariableNames) > 3
        % passive recording with SL
        data_table.Properties.VariableNames = {'Time', 'L', 'F', 'SL'};
%         data_table.Properties.VariableUnits = {'s', 'Lo', 'kPa', 'um'};
        % fix the unit
        data_table.Time = data_table.Time * 1000;
        data_table.Properties.VariableUnits = {'ms', 'Lo', 'kPa', 'um'};
        SL = true;
    else
        % active recording without SL
        data_table.Properties.VariableNames = {'Time', 'L', 'F'};
        data_table.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};
        SL = false;
    end
    
    % Downsample 
    if isempty(dwnsmpl)
        % no downsampling
        dwnsmpl = [data_table.Time(1) data_table.Time(end)];
        imin_d = find(data_table.Time >= dwnsmpl(1), 1);
        imax_d = find(data_table.Time >= dwnsmpl(end), 1);
        i_data = imin_d:imax_d;
        % Oh, Hi, Mark! We should not use marker unless they provide datapoitns
        hiMark = '-';
    elseif length(dwnsmpl) == 1
        % downsample by a factor
        dwnsmpl = 1:dwnsmpl:data_table.Time(end);
        i_data = interp1(data_table.Time, 1:length(data_table.Time), dwnsmpl, 'nearest');
        % Oh, Hi, Mark!
        hiMark = '|-';
    else
        % get the nearest neighbors to prevent missing the data
        i_data = interp1(data_table.Time, 1:length(data_table.Time), dwnsmpl, 'nearest', 'extrap');
        
        % make sure the bounds are where neeeded 
        % TODO FIX ME OR KILL ME
        data_table.Time(1) = dwnsmpl(1);data_table.Time(end) = dwnsmpl(end);
        % Oh, Hi, Mark!
        hiMark = '|-';
    end
    
    if isempty(ts_s)
        ts_s = [data_table.Time(1) data_table.Time(end)];
    end
    
    imin_s = find(data_table.Time >= ts_s(1), 1);
    imax_s = find(data_table.Time >= ts_s(end), 1);    

    t = data_table.Time(imin_s:imax_s);
    tf = data_table.Time; % we do not filter time    
    td = tf(i_data);
    
    l = data_table.L(imin_s:imax_s);
    lf = movmean(data_table.L,[dsf/2 dsf/2]); % l filtered
    % round to limit the oscillations
    lf = round(lf, 3);
    ld = lf(i_data);

    f = data_table.F(imin_s:imax_s)*scaleF;
    ff = movmean(data_table.F*scaleF,[dsf/2 dsf/2]); % force filtered
    fd = ff(i_data);
    
    if SL
        SL = data_table.SL(imin_s:imax_s);
        SLf = movmean(data_table.SL,[dsf/2 dsf/2]); % l filtered
        SLd = SLf(i_data);
        datatable = [td/1000  + offset/1000, ld*ML, fd, SLd];
    else
        datatable = [td/1000  + offset/1000, ld*ML, fd];
    end
    
    

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
        T = array2table(datatable);
        T.Properties.VariableNames(1:3) = {'Time','ML','Force'};
        writetable(T, ['data/' saveAs '.csv']);
        disp(['Saved as ' fn ' and csv'])
    end

%     figure();clf;
    subplot(211);hold on;title(saveAs, 'Interpreter', 'none');
    plot(td + offset, ld, hiMark);
    
    if SL
        plot(td + offset, SLd, '-');
        legend('ML', 'SL','AutoUpdate','off');
    end
    
    
%     plot(t + offset, l, tf + offset, lf, td + offset, ld, '|-');   
%     plot(velocitytable(:, 1)*1000, velocitytable(:, 4)/ML, 'x-', 'Linewidth', 1, 'MarkerSize', 10)
    plot([ts_s;ts_s]+offset, repmat([min(data_table.L);max(data_table.L)], 1, length(ts_s)))
    
    xlabel('time (ms)');
    ylabel('Length (ML)')

    
    subplot(212);hold on;
    plot(td + offset, fd, hiMark, 'Linewidth', 2, 'MarkerSize', 10);
    
%     plot(t + offset, f, tf + offset, ff, td + offset, fd, '|-', 'Linewidth', 2, 'MarkerSize', 10);
%     plot([ts_d;ts_d], repmat([min(data_table.F);max(data_table.F)], 1, length(ts_d)))
    xlabel('time (ms)');
    ylabel('Force (kPa)')
    
end

