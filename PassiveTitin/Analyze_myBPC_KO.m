% load Kerr McDonalds data to check influence of myBP-C

%% load
clear;
data = readtable('..\data\myBP-C KO\WT RAMP STRETCH IN 45 M&B 1 Kerry S. McDonald.TXT');
datatable = data(:, [3 2 4 6]);
datatable.Properties.VariableNames = {'Time', 'Motor', 'SL', 'Force'};
datatableWT = datatable;

data = readtable('..\data\myBP-C KO\KO RAMP STRETCH IN 45 MAVA BLEBB 1 Kerry S. McDonald.TXT');
datatable = data(:, [3 2 4 6]);
datatable.Properties.VariableNames = {'Time', 'Motor', 'SL', 'Force'};
datatableKO = datatable;

%% get the peaks to zero
[~, imax] = max(datatableWT.Force);
datatableWT.Time = datatableWT.Time - datatableWT.Time(imax) + 0.01;

[~, imax] = max(datatableKO.Force);
datatableKO.Time = datatableKO.Time - datatableKO.Time(imax) + 0.01;

%% normalize force to tension
% from personal communication
% WT: 0.74 kPa prior to ramp
% KO: 0.70 kPa
b_avg = (datatableWT.Time < -0.02 & datatableWT.Time > -2.02);
avg0WT = mean(datatableWT.Force(b_avg));
datatableWT.Force = datatableWT.Force*0.74/avg0WT;
% limit to 25s
datatableWT = datatableWT(datatableWT.Time < 25, :);

b_avg = (datatableKO.Time < -0.02 & datatableKO.Time > -2.02);
avg0KO = mean(datatableKO.Force(b_avg))
datatableKO.Force = datatableKO.Force*0.7/avg0KO;
% limit to 25s
datatableKO = datatableKO(datatableKO.Time < 25, :);

%% plot & compare
clf;nexttile;title('Data from Kerry McDonald');hold on;
plot(datatableKO.Time, datatableKO.Force, datatableWT.Time, datatableWT.Force, 'LineWidth',2);

set(gca,'ColorOrderIndex',1)
[~, imax] = max(datatableKO.Force);
plot(datatableKO.Time(imax), datatableKO.Force(imax), '*', LineWidth=4)

[~, imax] = max(datatableWT.Force);
plot(datatableWT.Time(imax), datatableWT.Force(imax), '*', LineWidth=4)

Fp_KO = datatableKO.Force(imax) - 0.7;
Fp_WT = datatableWT.Force(imax) - 0.74;

legend('KO', 'WT', sprintf('Peak KO %.1d', Fp_KO),sprintf('Peak WT %.1d', Fp_WT));
xlim([-5 25]);
xlabel('Time (s)');ylabel('Tension (kPa?)');

nexttile;title('Data from Anthony Baker');hold on;
load pca11data.mat
plot(Tarr{4}, Farr{4},Tarr{3}, Farr{3},Tarr{2}, Farr{2});
legend('Ramp 0.1s (averaged)','Ramp 1s (averaged)','Ramp 10s (averaged)')
xlim([-5 25])
xlabel('Time (s)');ylabel('Tension (kPa)');


%% power law ident
Farr = cell([4, 1]);
Tarr = cell([4, 1]);

% KO init
x = [0.0769    0.2394    0.6799]
Farr{4} = datatableKO.Force;
Tarr{4} = datatableKO.Time;

% WT init
% x = [0.0133    0.1746    0.7382]
% Farr{5} = datatableWT.Force;
% Tarr{5} = datatableWT.Time;

% Farr{2} = []; Tarr{2} = [];
% Farr{3} = []; Tarr{3} = [];
% Farr{1} = []; Tarr{1} = [];

options = optimset('Display','iter', 'TolFun', 1e-4, 'Algorithm','sqp', 'UseParallel', true, ...
    'TolX', 0.0001, 'PlotFcns', @optimplotfval, 'MaxIter', 150);
fitfunOpt = @(init) evalPowerFit(init, Farr, Tarr, false);
x = fminsearch(fitfunOpt, x, options)
[c rampShift] = evalPowerFit(x, Farr, Tarr, true, [], false);

%% plot KO
clf;
x = [0.0769    0.2394    0.6799]
Farr{4} = datatableKO.Force; Tarr{4} = datatableKO.Time;
[c rampShift] = evalPowerFit(x, Farr, Tarr, true, [], true);
ylim('auto');
yticks('auto')
yl = ylim();
ylim(yl)
xlim([1e-2, 35])
%% plot WT
x = [0.0133    0.1746    0.7382]
Farr{4} = datatableWT.Force; Tarr{4} = datatableWT.Time;
[c rampShift] = evalPowerFit(x, Farr, Tarr, true, [], true);