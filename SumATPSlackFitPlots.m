% prepares the plot for the physiome abstract

% decimation sampling (each x)
dsf = 10;
ML = 2.0;
% normalized force multiplier
nf = 56;
%% fitting figure
clf
data_table = readtable('data/8 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
o = -100;
ts_s = [];
% ts_s = [2500, 2759.6, 2760.4, 2910.4, 2930, 3050]
% [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_s([1, end])-o, ts_s -o, ML, dsf, nf/54, 'bakers_slack8mM', o);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [80 200], ts_s -o, ML, dsf, nf/54, '', o);
zones = [10, 50;];

figure(1);clf;    
[dSLpc, ktr, df, del] = fitRecovery(datatable, zones, 0);

%% different atp levels figure
ts_s = [];
o = 1150 - 100 + 9.4;
zones = [1162, 1209;1464 1519;1816 1889;2269 2359.5;2774 2900];

% 8 mM
data_table = readtable('data/8 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
figure(2);clf; 
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
figure(2);clf;
[dSLpc_8, ktr_8, df_8, del_8] = fitRecovery(datatable, zones, 0);

% 2 mM
data_table = readtable('data/2 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
figure(2);clf; 
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
figure(2);clf;
[dSLpc_2, ktr_2, df_2, del_2] = fitRecovery(datatable, zones, 0);

% 0.2 mM
data_table = readtable('data/0.2 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
figure(2);clf; 
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
figure(2);clf;
[dSLpc_02, ktr_02, df_02, del_02] = fitRecovery(datatable, zones, 0);

%% plotting
figure(2); clf;
subplot(121);hold on;
plot(-dSLpc_8, df_8, '*-',-dSLpc_2, df_2, 's-',-dSLpc_02, df_02, 'o-');
title({'Predicted maximal force', 'at step-down'});
xlabel('% of ML decrease');
ylabel('Force (kPa)')


subplot(122);hold on;
plot(-dSLpc_8, ktr_8, '*-',-dSLpc_2, ktr_2, 's-',-dSLpc_02, ktr_02, 'o-');
title({'Fitted exponential', 'factor k_{tr}'});
xlabel('% of ML decrease');
ylabel('k_{tr} (s_{-1})');



% subplot(133);hold on;
% plot(-dSLpc_8, del_8*1000, '*-',-dSLpc_2, del_2*1000, 's-',-dSLpc_02, del_02*1000, 'o-');
% title({'Fitted exponential', 'factor k_{tr}'});
% xlabel('% of ML decrease');
% ylabel('k_{tr} (s_{-1})')