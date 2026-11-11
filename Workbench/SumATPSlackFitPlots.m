% prepares the plot for the physiome abstract

% decimation sampling (each x)
dsf = 1;
ML = 2.0; 
% normalized force multiplier
nf = 56;
%% Stitch the experiments together

figure(101);clf;
% tss_d = [118555, 126800]
% tss_s = [118555, 121890, 121900, 121910,121920, ts_d(1) + (122070+710)-10, ts_d(end-1) + (122070+710), 123910, 123930, 123960, 124000, ...
%     124210, 124230, 124270, 124310, 124560, 124580, 124640, 124680, 125010, 125030, 125110, 125150, 125510, 125530, 125660, 125700, 126800]
data_table = readtable('data/8 mM ATP scope.txt', 'filetype', 'text', 'NumHeaderLines',4);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, 1, '', -(122070+710)+20);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, tss_d, tss_s, ML, 1, nf/54, 'bakers_rampup2_8_long', 0);
[datatable8s, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, '', -(122070+710)+20);

data_table = readtable('data/02 mM ATP scope.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable2s, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, '', -(122070+710)-2700);

data_table = readtable('data/0.2 mM ATP scope.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable02s, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, '', -(122070+710)-2700 + 1280 + 20-6);

% Atp 0.2 is shifted. Extend before 0 and after 1000ms by x
% datatable8s = datatable;
i0 = find(datatable8s(:, 1) > 0, 1);
datatable8s(1:i0, 1) = datatable8s(1:i0, 1) - 0.14;
datatable2s(1:i0, 1) = datatable2s(1:i0, 1) - 0.14;
% datatable02s; stays
% Extend after 800ms
% datatable02s = datatable;
i0 = find(datatable8s(:, 1) > 1, 1);
% 8 and 2 stays
% datatable02s(i0:end, 1) = datatable02s(i0:end, 1) + 0.13;
datatable8s(i0:end, 1) = datatable8s(i0:end, 1) - 0.13;

clf;
% subplot(211);
% plot(datatable8s(:, 1), datatable8s(:, 2)/ML,datatable2s(:, 1), datatable2s(:, 2)/ML,datatable02s(:, 1), datatable02s(:, 2)/ML, 'Linewidth', 2);
% xlim([-1.2, 4]);
% xlabel('Time'); ylabel('Muscle length (-)')
% set(gca,'fontsize',16);
% subplot(212);

yyaxis right; 
c = lines;
p1 = plot(datatable8s(:, 1), datatable8s(:, 2)/ML, '-', 'Linewidth', 2, Color=c(4, :));
ylim([0.8 1.5])
ax = gca;
ax.YAxis(2).Color = c(4, :);
ylabel('Muscle length (ML*, L/L_0)')
yyaxis left;
hold on;

p2 = plot(datatable8s(:, 1), datatable8s(:, 3),'-',  'Linewidth', 2, Color=c(1,  :));
p3 = plot(datatable2s(:, 1), datatable2s(:, 3),'-', 'Linewidth', 2, Color=c(2,  :));
p4 = plot(datatable02s(:, 1), datatable02s(:, 3),'-', 'Linewidth', 2, Color=c(3, :));
xlim([-1.2, 3]);
xlabel('Time (s)'); ylabel('Tension (kPa)');
set(gca,'fontsize',16);
leg = legend([p1 p2 p3 p4], 'ML*', '8 mM', '2 mM','0.2 mM', 'Location','northwest');
leg.ItemTokenSize = [20, 10]
ylim([-10 95])
fontsize(14, 'points')
c = [0.1 0.1 0.6] ;
rectangle(Position=[-1.1 -5 0.8 65], EdgeColor=c, FaceColor="none", LineStyle=":", LineWidth=4);
text(-1.1 + 0.3, 20, 'k_{tr}', 'FontSize',18, FontWeight='bold', FontName='Arial', Color=c, HorizontalAlignment='left', VerticalAlignment='top')


rectangle(Position=[-0.05 48 0.55 45], EdgeColor=c, FaceColor="none", LineStyle=":", LineWidth=4);
text( -0.2, 48, 'Ramp-up', 'FontSize',18, FontWeight='bold', FontName='Arial', Color=c, HorizontalAlignment='left', VerticalAlignment='top')

rectangle(Position=[0.95 0 1.95 85], EdgeColor=c, FaceColor="none", LineStyle=":", LineWidth=4);
text(0.95, 85, 'Slack', 'FontSize',18, FontWeight='bold', FontName='Arial', Color=c, HorizontalAlignment='left', VerticalAlignment='bottom')

set(gcf, "Position", [488.2000  195.4000  612.0000  466.4000]);

%% fitting figure
clf
data_table = readtable('data/8 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
o = -100;
ts_s = [];
% ts_s = [2500, 2759.6, 2760.4, 2910.4, 2930, 3050]
% [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_s([1, end])-o, ts_s -o, ML, dsf, nf/54, 'bakers_slack8mM', o);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [80 200], ts_s -o, ML, dsf, nf/54, '', o);
Fm_8 = mean(datatable(datatable(1:100, 1) < 0, 3));
zones = [5, 50;];

figure(1);clf;    
[dSLpc, ktr_s, df, del, e] = fitRecovery(datatable, zones, 0);
%% Ktr protocol
fign = 10;
clf;
data_table = readtable('data/8mM ATP 2ktr.txt', 'filetype', 'text', 'NumHeaderLines',4);
figure(fign+10);clf; 
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
zones = [1170, 1500];
figure(fign);
[dSLpc, ktr_ktr(1), df(1), del, e] = fitRecovery(datatable, zones, 0);
%
data_table = readtable('data/2 mM ATP ktr.txt', 'filetype', 'text', 'NumHeaderLines',4);
figure(fign+10);clf; 
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
zones = [1170, 1500];
figure(fign);
[dSLpc, ktr_ktr(2), df(2), del, e] = fitRecovery(datatable, zones, 0, []);
%
data_table = readtable('data/0.2 mM ATP ktr.txt', 'filetype', 'text', 'NumHeaderLines',4);
figure(fign+10);clf; 
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
zones = [1170, 1500];
figure(fign);
[dSLpc, ktr_ktr(3), df(3), del, e] = fitRecovery(datatable, zones, 0, []);
%% different atp levels figure
fign = 3;
ts_s = [];
o = 1150 - 100 + 9.4;
zones = [1162, 1209;1464 1519;1816 1889;2269 2359.5;2774 2900];

% 8 mM
data_table = readtable('data/8 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
figure(fign);clf; 
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
figure(fign);clf;
[dSLpc_8, ktr_sla8, df_8, del_8, e_8, SL] = fitRecovery(datatable, zones, 0);
% df_8 = [];
% 2 mM
data_table = readtable('data/2 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
figure(fign+100);clf; 
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
figure(fign);
[dSLpc_2, ktr_sla2, df_2, del_2, e_2] = fitRecovery(datatable, zones, 0, df_8);

% 0.2 mM
data_table = readtable('data/0.2 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
figure(fign+100);clf; 
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
figure(fign);
[dSLpc_02, ktr_sla02, df_02, del_02, e_02] = fitRecovery(datatable, zones, 0, df_8);

plot([1 3], [0 0], 'k', LineWidth=0.5)
%% cutting out unnecessary lines manually, then jsut setting the legend
% figure(4); clf;
zones = [2269 2359.5];
fign = 1010;
% 8 mM
data_table = readtable('data/8 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
figure(fign);clf; 
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
figure(fign);clf;
[dSLpc_8, ktr_sla8, df_8, del_8, e_8, SL] = fitRecovery(datatable, zones, 0);
% df_8 = [];
% 2 mM
data_table = readtable('data/2 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
figure(fign+100);clf; 
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
figure(fign);
[dSLpc_2, ktr_sla2, df_2, del_2, e_2] = fitRecovery(datatable, zones, 0, df_8);

% 0.2 mM
data_table = readtable('data/0.2 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
figure(fign+100);clf; 
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, dsf, nf/54, '', o);
figure(fign);
[dSLpc_02, ktr_sla02, df_02, del_02, e_02] = fitRecovery(datatable, zones, 0, df_8);
%%
legend('ATP 8 mM','Fitted 8 mM','ATP 0.2 mM','Fitted 0.2 mM')
fontsize(16, 'points');
% xticklabels = {'a', 'b'}
currentTicks = xticks;

% Calculate the new tick labels (starting from 0)
newTickLabels = currentTicks - currentTicks(1);

% Apply new tick labels
xticklabels(newTickLabels);
ylim([-10 100])
%% results:
% Free df:
% e_2: 0.3241    0.1784    0.2072    0.1893    0.1835
% e_02:0.2093    0.2007    0.1802    0.1782    0.1975
% sum(e_2):    1.0825
% sum(e_02):   0.9658
% 2.0484
% Fixed df:
% e_2:  0.6620    0.4679    0.4892    0.5106    0.6342
% e_02: 0.3572    0.2785    0.2975    0.2304    0.2145
% sum(e_2):    2.7639
% sum(e_02):    1.3782
% 4.1421
%%
%SLm = max(data_table(:, 2).*ML);
SLm = 2.2;
dSL = SLm - SL;
figure(4);clf;hold on;
SLs = [SLm SL];
Fs = [Fm_8 df_8];
plot(SLs, Fs, 'x-', LineWidth=2); xlabel('SL (um)');ylabel('Fmax (est) (kPa)');
%% plotting against drop
figure(2); clf;
subplot(121);hold on;
% plot(-[0 dSLpc_8], [Fm_8 df_8], '*-',-dSLpc_2, df_2, 's-',-dSLpc_02, df_02, 'o-', LineWidth=2);
plot(-[0 dSLpc_8], [Fm_8 df_8], 'x-', LineWidth=2);
title({'Predicted maximal force', 'at step-down'});
xlabel('% of ML decrease');
ylabel('Force (kPa)')


subplot(122);hold on;
plot(-dSLpc_8, ktr_sla8, '*-',-dSLpc_2, ktr_sla2, 's-',-dSLpc_02, ktr_sla02, 'o-', LineWidth=2);
title({'Fitted exponential', 'factor k_{tr}'});
xlabel('% of ML decrease');
ylabel('k_{tr} (s_{-1})');

%% plotting against SL
figure(2); clf;
subplot(121);hold on;
% plot(-[0 dSLpc_8], [Fm_8 df_8], '*-',-dSLpc_2, df_2, 's-',-dSLpc_02, df_02, 'o-', LineWidth=2);
plot(SLs, [Fm_8 df_8], 'x-', LineWidth=2);
set(gca,'ColorOrderIndex',1)
plot(SLs(1), Fs(1), 's', LineWidth=2); 
set(gca,'ColorOrderIndex',1)
plot(2.0, df, 'o', LineWidth=2); 
title({'Predicted maximal force', 'at step-down'});
xlabel('SL');
ylabel('Force (kPa)')
legend('Slack protocol', '2.2 steady state', ...
    '2.0 steady state - 8mM', '2.0 steady state - 2mM', '2.0 steady state - 0.2mM',Location='best')

subplot(122);hold on;
set(gca,'ColorOrderIndex',1)
plot(SL, ktr_sla8, '*-',SL, ktr_sla2, 's-',SL, ktr_sla02, 'o-', LineWidth=2);
set(gca,'ColorOrderIndex',1)
plot(2.0, ktr_ktr, 'o', LineWidth=2);
title({'Fitted exponential', 'factor k_{tr}'});
xlabel('SL');
ylabel('k_{tr} (s_{-1})');
legend('Slack protocol - 8mM', 'Slack protocol - 2mM','Slack protocol - 0.2mM',...
    'Ktr protocol - 8mM','Ktr protocol - 2mM','Ktr protocol - 0.2mM',Location='best')



% subplot(133);hold on;
% plot(-dSLpc_8, del_8*1000, '*-',-dSLpc_2, del_2*1000, 's-',-dSLpc_02, del_02*1000, 'o-');
% title({'Fitted exponential', 'factor k_{tr}'});
% xlabel('% of ML decrease');
% ylabel('k_{tr} (s_{-1})')

%% This is weird

data_table_ktr = readtable('data/8mM ATP 2ktr.txt', 'filetype', 'text', 'NumHeaderLines',4);
data_table_slack = readtable('data/8 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
%%
figure(101);clf;
o= 301.5;
subplot(211);hold on;
plot(data_table_slack.Var1, data_table_slack.Var2, LineWidth=2);
plot(data_table_ktr.Var1 + o, data_table_ktr.Var2, LineWidth=2);
xlim([90+o-10 o+100+100]);
legend('slack protocl', 'ktr protocol', Location='best')
subplot(212);hold on;
plot(data_table_slack.Var1, data_table_slack.Var3, LineWidth=2);
plot(data_table_ktr.Var1 + o, data_table_ktr.Var3, LineWidth=2);
legend('slack protocl', 'ktr protocol', Location='best')
xlim([90+o-10 o+100+100]);