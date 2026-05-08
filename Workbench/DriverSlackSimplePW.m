% DriverSlackSimplePW.m
% Compare full vs simplified piecewise strain-dependent rates (R1D 7→5, R2 10→5).
% Uses justPlotStateTransitionsFlag to visualise the curves side-by-side,
% then runs the slack experiment with the simplified parameterisation.

clear; close all; clc;

%% Load latest optimised params
params0 = getParams();
ModelOptParams_SRXTD_iter_9

params0.RunForceVelocity    = false;
params0.RunKtr              = false;
params0.RunStairs           = false;
params0.RunSlack            = true;
params0.RunForceLengthEstim = false;
params0.RunSlackSegments    = 'All';
params0.EvalFeatures        = true;
params0.PlotEachSeparately  = true;
params0.BreakOnODEUnstable  = true;
params0.MaxRunTime          = 600;

%% Simplified piecewise definitions
%
% R1D (kd multiplier, p1→PD): 7→5 pts
%   dropped: X=-0.015 (P=4000, redundant with far-left clamp)
%            X=+0.00417 (P=8.6, smooth PCHIP covers 0→0.01 adequately)
R1D_X_full  = [-0.05  -0.015  -0.00366  0       0.00417  0.01048  0.1];
R1D_P_full  = [4000    4000    43.891   4.187   8.578   135.308  2000];
R1D_X_simp  = [-0.05          -0.00366  0                0.01048  0.1];
R1D_P_simp  = [4000             43.891  4.187           135.308  2000];

% R2 (k2 multiplier, p2→PD): 10→5 pts
%   collapsed left plateau (X=-0.01,-0.0105,-0.0085, all P≈50) into single far-left clamp
%   kept left corner, center, right corner, far-right clamp
R2_X_full   = [-0.01  -0.0105  -0.0085  -0.00555  -0.00427  0       0.00618  0.01  0.0233  0.1];
R2_P_full   = [50      50       50.608   1.133     0.968     0.467   1        20    50      50];
R2_X_simp   = [-0.01                    -0.00555            0       0.00618  0.01          0.1];
R2_P_simp   = [50                        1.133              0.467   1        20            50];

%% Figure: compare full vs simplified curves
s_eval = linspace(-0.06, 0.12, 500);

figure(200); clf;
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(s_eval, ppval(pchip(R1D_X_full, R1D_P_full), s_eval), 'b-',  'LineWidth', 2, 'DisplayName', 'full (7 pts)');
plot(s_eval, ppval(pchip(R1D_X_simp, R1D_P_simp), s_eval), 'r--', 'LineWidth', 2, 'DisplayName', 'simplified (5 pts)');
plot(R1D_X_full, R1D_P_full, 'bo', 'MarkerFaceColor', 'b');
plot(R1D_X_simp, R1D_P_simp, 'rs', 'MarkerFaceColor', 'r');
set(gca, 'YScale', 'log'); xlim([-0.06 0.12]); grid on; box on;
xlabel('Strain'); ylabel('R1D multiplier'); title('R1D  (kd, p1→PD)');
legend('Location', 'best');

nexttile; hold on;
plot(s_eval, ppval(pchip(R2_X_full, R2_P_full), s_eval), 'b-',  'LineWidth', 2, 'DisplayName', 'full (10 pts)');
plot(s_eval, ppval(pchip(R2_X_simp, R2_P_simp), s_eval), 'r--', 'LineWidth', 2, 'DisplayName', 'simplified (5 pts)');
plot(R2_X_full, R2_P_full, 'bo', 'MarkerFaceColor', 'b');
plot(R2_X_simp, R2_P_simp, 'rs', 'MarkerFaceColor', 'r');
set(gca, 'YScale', 'log'); xlim([-0.06 0.12]); grid on; box on;
xlabel('Strain'); ylabel('R2 multiplier'); title('R2  (k2, p2→PD)');
legend('Location', 'best');

sgtitle('Full vs simplified piecewise — shape comparison');

%% Apply simplified curves
params0.PieceWiseStrainDepR1DX      = R1D_X_simp;
params0.PieceWiseStrainDepR1DParams = R1D_P_simp;
params0.PieceWiseStrainDep2X        = R2_X_simp;
params0.PieceWiseStrainDep2Params   = R2_P_simp;

%% Run slack experiment
figure(201); clf;
tic
RunBakersExp;
toc

if params0.EvalFeatures
    plotFeatures(features_data, features_model, [], params0.fn);
end
