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

figure(2);clf;
plot(atpd(:, 1), atpd(:, 2), '-*')

% Separate the data into x and y
x = atpd(:,1);
y = atpd(:,2);

% Define the model as a sum of two exponential decays
model = fittype('a*exp(-b*x) + c*exp(-d*x)', ...
                'independent', 'x', ...
                'coefficients', {'a', 'b', 'c', 'd'});

% Set initial guesses for the parameters
initialGuess = [1, 0.1, 0.5, 0.01]; % Adjust these based on your data

% Perform the fit
[fitResult, gof] = fit(x, y, model, 'StartPoint', initialGuess);

% Display the fit results and goodness of fit
disp('Fit Coefficients:');
disp(fitResult);
disp('Goodness of Fit:');
disp(gof);

% Plot the data and the fitted curve
figure;
hold on;
scatter(x, y, 'o', 'DisplayName', 'Data'); % Scatter plot of the data
xFit = linspace(min(x), max(x), 100); % Generate fine-grained x values for the fit
yFit = feval(fitResult, xFit); % Evaluate the fit at these x values
plot(xFit, yFit, 'r-', 'DisplayName', 'Fit'); % Plot the fit
xlabel('x');
ylabel('y');
legend;
title('Two Exponential Decay Fit');
grid on;

% Evaluate the two exponential components separately
yFit1 = fitResult.a * exp(-fitResult.b * xFit); % First exponential component
yFit2 = fitResult.c * exp(-fitResult.d * xFit); % Second exponential component

% Evaluate the two exponential components separately
yFit1_s = 0.70 * exp(-0.052 * xFit); % First exponential component
yFit2_s = 0.3 * exp(-0.0061* xFit); % Second exponential component


% Plot the individual components
plot(xFit, yFit1, 'g--', 'DisplayName', 'First Exponential Component');
plot(xFit, yFit2, 'b--', 'DisplayName', 'Second Exponential Component');
plot(xFit, yFit1_s + yFit2_s, 'k--', 'DisplayName', 'Second Exponential Component', 'LineWidth',2);


% Define the single exponential decay model
singleExpModel = fittype('a*exp(-b*x)', ...
                         'independent', 'x', ...
                         'coefficients', {'a', 'b'});

% Set initial guesses for the parameters
initialGuessSingle = [1, 0.1]; % Adjust these based on your data

% Perform the fit for the single exponential
[fitResultSingle, gofSingle] = fit(x, y, singleExpModel, 'StartPoint', initialGuessSingle);

% Display the fit results and goodness of fit
disp('Single Exponential Fit Coefficients:');
disp(fitResultSingle);
disp('Single Exponential Goodness of Fit:');
disp(gofSingle);

% Evaluate the single exponential fit
yFitSingle = feval(fitResultSingle, xFit);

% Plot the single exponential fit
plot(xFit, yFitSingle, 'm-.', 'DisplayName', 'Single Exponential Fit'); % Magenta dashed-dotted line
