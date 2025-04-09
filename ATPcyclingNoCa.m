% Based on ATP decay data from Toepfer 2020, 
% aka "Myosin Sequestration Regulates Sarcomere Function, Cardiomyocyte Energetics, and Metabolism, Informing the Pathogenesis of Hypertrophic Cardiomyopathy"

atpd = [0.252525, 1.00126
10.8586, 0.715006
20.4545, 0.491803
30.3030, 0.388398
40.6566, 0.327869
49.7475, 0.283733
60.3535, 0.248424
70.4545, 0.228247
81.0606, 0.201765
90.1515, 0.184111
100.758, 0.172762
109.848, 0.155107
120.960, 0.142497
130.556, 0.136192
140.404, 0.126103
150.505, 0.119798
160.606, 0.110971
170.455, 0.107188
180.808, 0.0996217
191.162, 0.0933165
200.758, 0.0870113
210.859, 0.0870113
220.707, 0.0794451
230.808, 0.0756620
240.404, 0.0693569
250.505, 0.0668348
260.606, 0.0643127
270.455, 0.0630517
281.061, 0.0580076
290.909, 0.0580076
300.758, 0.0517024
];

tab = readtable("data/ATP_control.csv");
% tab = readtable("data/ATP_Mava.csv");

% save as mat
save('data/tables.mat', 'tab', '-mat')
% save("data/ATP_control.mat", "tab")

tab_i = 1;

switch tab_i
    % case 3
    % atpd = table2array(tab(tab.t>40, :), [1, 1+ tab_i]);
    otherwise
atpd = table2array(tab(:, [1, 1+ tab_i]));
end

% Separate the data into x and y
x = [0 ; atpd(:,1)];
y = [1 ; atpd(:,2)];
x = atpd(:,1);
y = atpd(:,2);

% plot just the raw data

figure(2);clf;
% hold on;
% plot(atpd(:, 1), atpd(:, 2), '-*')



% Define the model as a sum of two exponential decays
model = fittype('a*exp(-b*(x - 0*x0)) + (y_n -a)*exp(-d*(x-0*x0)) + y0', ...
                'independent', 'x', ...
                'coefficients', {'a', 'b', 'd', 'x0', 'y_n', 'y0'});

% Set initial guesses for the parameters
initialGuess = [0.5, 0.001, 0.01, 0, 1, 0]; % Adjust these based on your data

opts = fitoptions('StartPoint', initialGuess, 'Method', 'NonlinearLeastSquares','Lower',[0 0, 0, 0, 0.5, -1],'Upper',[1, 1, 1, 40, 1.5, 1]);

% Perform the fit
[fitResult, gof] = fit(x, y, model, opts);
fitResult

% sqrErr = (fitResult(x) - y).^2;
% plot(x, cumsum(sqrErr))
fprintf('SRX %0.1f, rate %0.3e \n', fitResult.a*100, fitResult.b)

% Display the fit results and goodness of fit
% disp('Fit Coefficients:');
% disp(fitResult);
% disp('Goodness of Fit:');
% disp(gof);

% Plot the data and the fitted curve
% figure;
hold on;
scatter(x, y, 'o', 'DisplayName', 'Data'); % Scatter plot of the data
xFit = linspace(fitResult.x0, max(x), 100); % Generate fine-grained x values for the fit
yFit = feval(fitResult, xFit); % Evaluate the fit at these x values
plot(xFit, yFit, 'k-', 'DisplayName', 'Fit', LineWidth=2); % Plot the fit
xlabel('x');
ylabel('y');
legend;
title('Two Exponential Decay Fit');
grid on;

% Evaluate the two exponential components separately
% fitResult_1 = fitResult; fitResult_1.c = 0;
yFit1 = fitResult.a * exp(-fitResult.b * xFit); % First exponential component
yFit2 =  (1 - fitResult.a)* exp(-fitResult.d * xFit); % Second exponential component

% % Evaluate the two exponential components separately
% yFit1_s = 0.70 * exp(-0.052 * xFit); % First exponential component
% yFit2_s = 0.3 * exp(-0.0061* xFit); % Second exponential component


% % Plot the individual components
plot(xFit, yFit1, 'g--', 'DisplayName', 'First Exponential Component');
plot(xFit, yFit2, 'b--', 'DisplayName', 'Second Exponential Component');


% plot(xFit, yFit1_s + yFit2_s, 'k--', 'DisplayName', 'Second Exponential Component', 'LineWidth',2);


% Define the single exponential decay model
singleExpModel = fittype('a*exp(-b*(x-0*x0))', ...
                         'independent', 'x', ...
                         'coefficients', {'a', 'b', 'x0'});

% Set initial guesses for the parameters
initialGuessSingle = [0.8, 0.01, 0.1]; % Adjust these based on your data

% Perform the fit for the single exponential
[fitResultSingle, gofSingle] = fit(x, y, singleExpModel, 'StartPoint', initialGuessSingle);

% Display the fit results and goodness of fit
% disp('Single Exponential Fit Coefficients:');
% disp(fitResultSingle);
% disp('Single Exponential Goodness of Fit:');
% disp(gofSingle);

% Evaluate the single exponential fit
yFitSingle = feval(fitResultSingle, xFit);

% Plot the single exponential fit
plot(xFit, yFitSingle, 'm-.', 'DisplayName', 'Single Exponential Fit', 'LineWidth',2); % Magenta dashed-dotted line

%% Rates from Alison's model fit
% mava
names = ["k_{SR^-}", "k_{SR^+}", "{TOR_0}"];

f = figure(1);clf; 
aspect = 1;
f.Position = [300 200 7.2*96 7.2*96/aspect];

tiledlayout(2, 2, TileSpacing='tight');
nexttile([1 2]);
imshow(imread('mantatpschema.png'))

nexttile;hold on;
ctrl = readtable('ctrl1.csv');mava = readtable('Mava1.csv');
plot(ctrl.Time, ctrl.integratedSquaredDeviation_u1, 'b-');plot(ctrl.Time, ctrl.integratedSquaredDeviation_u2, 'b--', ctrl.Time, ctrl.SRX_solute, 'b:', LineWidth=2)
plot(mava.Time, mava.integratedSquaredDeviation_u1, 'r-');plot(mava.Time, mava.integratedSquaredDeviation_u2, 'r--', mava.Time, mava.SRX_solute, 'r:', LineWidth=2)
xlim([-10 inf]);
legend('Control: data', 'Control: model', 'Control: predicted SRX', 'Mava: data', 'Mava: simulation', 'Mava: predicted SRX')
ylim([0, 1.1]);
ylabel('Normalized (-)');xlabel('Time (s)');
fontsize(12, "points");
title('Fluorescence');

nexttile;
mava = [0.2586141, 16.004807, 1.5708864];
ctrl = [0.25814757, 6.1799517, 1.1841985];

bar([ctrl;mava]')
% set(gca, 'YScale', 'log')
ylabel('Transition rate (s^{-1})'); xticklabels(names);
% figure(2);
legend('Control', 'Mavacamten');
fontsize(12, "points");

%% Rates from Hoojiman & Cooke 2011
% Define the data
Incubation = ["mantATP"; "ATP"; "mantADP"; "ADP*"; "mantATP + ATP"; "mantATP"; "mantATP"];
Chase = ["ATP"; "mantATP"; "ADP"; "mantATP"; "ATP"; "ADP"; "ATP/Ca"];
N = [64; 9; 15; 17; 8; 17; 22];
P1 = [0.66; 0.63; 0.83; 0.92; 0.51; 0.74; 0.64];
% Standard error
SE_P1 = [0.01; 0.03; 0.01; 0.02; 0.06; 0.02; 0.04];
T1 = [14; 12; 10; 9; 17; 12; 26];
SE_T1 = [1; 1; 1; 1; 4; 1; 2];
P2 = [0.27; 0.29; 0.09; 0.04; 0.08; 0.19; 0.26];
SE_P2 = [0.01; 0.02; 0.01; 0.02; 0.01; 0.02; 0.02];
T2 = [144; 138; 261; 75; 123; 76; 228];
SE_T2 = [10; 21; 49; 16; 12; 10; 22];

% Create the table
T = table(Incubation, Chase, N, P1, SE_P1, T1, SE_T1, P2, SE_P2, T2, SE_T2);

% Display the table
% disp(T);

figure(1);clf;
t = 0:0.1:600;
for i = [1:4 6:length(N)]
    nexttile;hold on;
    
    combinedModelHoojiman = @(t) 1 - P1(i)*(1 - exp(-t/T1(i))) - P2(i)*(1 - exp(-t/T2(i)));
    plot(t, combinedModelHoojiman(t), ...
    'DisplayName', sprintf("Hoojiman's fit: 1 - %0.2g*(1 - exp(-%0.2ft)) - %0.2g*(1 - exp(-%0.3ft))", ...
    P1(i), 1/T1(i), P2(i), 1/T2(i)), LineWidth=3);

    % singleExpModel_slow = @(t) P1(i) - P1(i)*(1 - exp(-t/T1(i)));
    % singleExpModel_fast = @(t) P2(i) - P2(i)*(1 - exp(-t/T2(i)));
    % plot(t, singleExpModel_slow(t), '--','DisplayName','Hoojiman slow',LineWidth=2)
    % plot(t, singleExpModel_fast(t), '--', 'DisplayName','Hoojiman fast',LineWidth=2);
    
    
    % combinedModelExpForm = @(t) P1(i)*(exp(-t/T1(i))) + P2(i)*(exp(1-t/T2(i)));
    % plot(t, combinedModelExpForm(t), '--', 'DisplayName', sprintf("Exp form: 1 - %0.2g*(exp(-%0.2ft)) - %0.2g*(exp(-%0.3ft))", ...
    % P1(i), 1/T1(i), P2(i), 1/T2(i)), LineWidth=2);
    
    model = fittype('a*exp(-b*(x)) + (1 -a)*exp(-d*(x))', ...
                'independent', 'x', ...
                'coefficients', {'a', 'b', 'd'});

    % Set initial guesses for the parameters
    initialGuess = [0.5, 0.01, 0.001]; % Adjust these based on your data
    
    opts = fitoptions('StartPoint', initialGuess, 'Method', 'NonlinearLeastSquares','Lower',[0 0, 0],'Upper',[1, 1, 1]);
    
    % Perform the fit
    [fitResult, gof] = fit(t', combinedModelHoojiman(t)', model, opts);
    fr1 = fitResult;fr1.d = 1e9;
    fr2 = fitResult;fr2.b = 1e9;
    plot(t, fitResult(t), ':', 'DisplayName', sprintf("Exp refit: 1 - %0.2g*(exp(-%0.2ft)) - %0.2g*(exp(-%0.3ft))", ...
    fitResult.a, fitResult.b, 1-fitResult.a, fitResult.d), LineWidth=3);
    t_offset = t(2:end);
    plot(t_offset, fr1(t_offset), ':', t_offset, fr2(t_offset), ':', LineWidth=2)
    title(sprintf("%s -> %s, SRX %0.2f Hoojiman, %0.2f alternatively", Incubation(i), Chase(i), P1(i)/(P1(i)+P2(i)), fitResult.a))

legend show;
end


%% Load newest data from Julia
% Specify the filename and the sheet name  
filename = '../data/EHT_Mant_Data_Filip.xlsx'; % Replace with your Excel file's path  
sheetname = 'EHT_Control'; % Name of the sheet to load  

% Load the data into a table  
dataTable = readtable(filename, 'Sheet', sheetname);  

% Display the table (optional)  
% disp(dataTable);


figure(10);clf; hold on;
for i_Pos = [8     7     6 5     4     3     2     1]
nexttile; hold on;
% cut from 60s
x = dataTable.Time(2:end) - dataTable.Time(2);
y = dataTable(2:end, 2);

% raw
x = dataTable.Time;
y = dataTable(:, i_Pos + 1).Variables;


y = (y- min(y))./(max(y) - min(y));

% Define the model as a sum of two exponential decays
model = fittype('a*exp(-b*(x - 0*x0)) + (y_n -a)*exp(-d*(x-0*x0)) + y0', ...
                'independent', 'x', ...
                'coefficients', {'a', 'b', 'd', 'x0', 'y_n', 'y0'});

modelHooj = fittype('1 - a*(1 - exp(-b*(x - 0*x0))) - c*(1 - exp(-d*(x-0*x0))) + y0', ...
                'independent', 'x', ...
                'coefficients', {'a', 'b', 'c', 'd', 'x0', 'y0'});

% Set initial guesses for the parameters
initialGuess = [0.5, 0.001, 0.01, 0, 1, 0]; % Adjust these based on your data

opts = fitoptions('StartPoint', initialGuess, 'Method', 'NonlinearLeastSquares','Lower',[0 0, 0, 0, 0.5, -.1],'Upper',[1, 1, 1, 40, 1.5, .1]);

% Perform the fit
[fitResult, gof] = fit(x, y, model, opts);
fitResult

[fitResultHooj, gof] = fit(x, y, modelHooj, opts);


params(i_Pos, 1) = fitResult.a;
params(i_Pos, 2) = fitResult.b;
params(i_Pos, 3) = fitResult.d;

t = 0:.1:600;
plot(x, y, '|-', 'DisplayName', 'Data', LineWidth=2);
plot(t, fitResult(t), '--', 'DisplayName', sprintf("Exp refit: %0.2g*(exp(-%0.2ft)) - %0.2g*(exp(-%0.3ft))", ...
    fitResult.a, fitResult.b, 1-fitResult.a, fitResult.d), LineWidth=3);
plot(t, fitResultHooj(t), '--', 'DisplayName', sprintf("Hoojiman fit: %0.2g*(exp(-%0.2ft)) - %0.2g*(exp(-%0.3ft))", ...
    fitResultHooj.a, fitResultHooj.b, fitResultHooj.c, fitResultHooj.d), LineWidth=3);

    fr1 = fitResult;fr1.d = 1e9;
    fr2 = fitResult;fr2.b = 1e9;
    
    
    t_offset = t(2:end);
    plot(t_offset, fr1(t_offset), ':', 'DisplayName', sprintf("Slow mode: %0.2g*(exp(-%0.2ft))", ...
    fitResult.a, fitResult.b), LineWidth=2);  
    plot(t_offset, fr2(t_offset), ':', 'DisplayName', sprintf("Fast mode: %0.2g*(exp(-%0.2ft))", ...
    fitResult.y_n - fitResult.a, fitResult.d), LineWidth=2)
    legend show;
    ylim([0 1]);xlim([0 600]);
    
    % title(sprintf("%s SRX %0.2f Hoojiman, %0.2f alternatively", dataTable.Properties.VariableNames{i_Pos+1}, P1(i_Pos)/(P1(i_Pos)+P2(i_Pos))))
    fontsize(12, 'points')
end