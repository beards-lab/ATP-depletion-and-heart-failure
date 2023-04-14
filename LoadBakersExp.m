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
 clf;
data_table = readtable('data/8 mM ATP slack.txt', 'filetype', 'text', 'NumHeaderLines',4);
o = 1150 - 100 + 9.4;
ts_s = [0 1070 1159 2259 2759.6 2760.4 2910 2930 3058]; % to prevent skipping events with large integrator step
% ts_s = [2500, 2759.6, 2760.4, 2910.4, 2930, 3050]
% [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_s([1, end])-o, ts_s -o, ML, dsf, nf/54, 'bakers_slack8mM', o);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s -o, ML, 1, nf/54, 'bakers_slack8mM', o);
% subplot(211)
% title('Slack experiment for different ATP concentrations')
% legend('8 mM', '2 mM', '0.2 mM')
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

legend('ForceLength8mM_all', 'bakers_rampup8', 'bakers_rampup2_8', 'bakers_rampup2_8_long', 'nevim', 'Interpreter', 'None')
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
%%
% clf;
data_table = readtable('data/0.2 mM stretch.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, dsf*10, nf/54, 'bakers_rampup2_02');
%%
data_table = readtable('data/0.2 mM ATP scope.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, 1, nf/54, 'bakers_rampup2_02_long', -(122070+710)-2700 + 1280 + 20 - 5);
%%
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
%%
data_table = readtable('data/relaxed stretch.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, [], [], ML, dsf, nf/54, 'bakers_rampup2_rel');

return;


%% New stretch experiments
figure(12);clf;
subplot(121);hold on;xlabel('ML');ylabel('Tension');title('Tension on Muscle length')
subplot(122);hold on;xlabel('ML (um)');ylabel('SL (um)');title('SL on Muscle length')

figure(11);clf;

el = 20; % experiment length in ms
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):2000+el*2 + el*20];
data_table = readtable('data/20ms_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s, ML, 1, 1, '');
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, 10, 1, 'bakers_passiveStretch_20ms');

figure(12); 
subplot(121); plot(datatable(:, 2), datatable(:, 3), 'o-');
subplot(122); plot(datatable(:, 2), datatable(:, 4), 'o-');
figure(11);

el = 100; % experiment length in ms
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):2000+el*2 + el*20];
data_table = readtable('data/100ms_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
% [datatable, velocitytable] = DownSampleAndSplit(data_table, [], ts_s, ML, 1, 1, '');
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, 10, 1, 'bakers_passiveStretch_100ms');

figure(12); 
subplot(121); plot(datatable(:, 2), datatable(:, 3), 'o-');
subplot(122); plot(datatable(:, 2), datatable(:, 4), 'o-');
figure(11);

el = 1000; % experiment length in ms
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):2000+el*2 + el*20];
data_table = readtable('data/1s_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, 1, 1, 'bakers_passiveStretch_1000ms');

figure(12); 
subplot(121); plot(datatable(:, 2), datatable(:, 3), 'o-');
subplot(122); plot(datatable(:, 2), datatable(:, 4), 'o-');
figure(11);

el = 10000; % experiment length in ms
ts_d = [-5000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):2000+el*2 + el*15];
data_table = readtable('data/10s_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, 1, 1, 'bakers_passiveStretch_10000ms');

figure(12); 
subplot(121); plot(datatable(:, 2), datatable(:, 3), 'o-');
subplot(122); plot(datatable(:, 2), datatable(:, 4), 'o-');
figure(11);
%%
el = 100000; % experiment length in ms
ts_d = [-500000, 0:500:2000, 2000:ceil(el/20):2000 + el*2, 2000 + el*2:ceil(el/2):2000+el*2 + el];
data_table = readtable('data/100s_4.txt', 'filetype', 'text', 'NumHeaderLines',4);
[datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, [], ML, 10, 1, 'bakers_passiveStretch_100000ms');
slowest = datatable(:, 3)

figure(12); 
subplot(121); plot(datatable(:, 2), datatable(:, 3), 'o-');
subplot(122); plot(datatable(:, 2), datatable(:, 4), 'o-');
figure(11);

figure(12); 
subplot(121); legend('20 ms', '100 ms', '1000 ms', '10 s', '100 s', 'location', 'Northwest')

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

