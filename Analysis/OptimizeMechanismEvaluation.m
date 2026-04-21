% OptimizeMechanismEvaluation.m
% Optimizes the top K parameters identified by Sensitivity Analysis
% against the experimental data.

clear; close all; clc;

%% 1. Initialize Baseline Model
params0 = getParams();
ModelParamsSlackKtrOpt;

% Apply Mechanism Corrections (Actions 1-3)
params0.UseTitinIdentifiedPassive = 0;
params0.k_pas = 77*2*2*2;
params0.gamma = 7.9;
params0.UseOverlap = true;
params0.UseOverlapFactor = true;

params0.L_thick = 1.65;
params0.L_thin = 1.15;

params0.UseNegativeKstiff = 0;
params0.kstiff2 = 15000;
params0.kstiff1 = "=kstiff2";
params0.MaxSlackNegativeForce = -1;
params0.kSE = 2000;

params0.ka = 67*2*2;
params0.kd = 9.55*2;
params0.k1 = 270*2;
params0.k2 = 133*0.85;
params0.kah = 30*1.5;
params0.A1AttachmentWidth = 0.004;
params0.PieceWiseStrainDepR1DX = [-0.05, -0.015, -0.004, 0, 0.004, 0.01, 0.1];
params0.PieceWiseStrainDepR1DParams = [1000, 1000, 10.147, 1, 2, 30, 500]*4; 
params0.dr2 = "=dr";
params0.mu = 0.633*0.5;
params0.mu_neg = "=mu";
params0.PieceWiseStrainDepX = [-1, -0.0083, -0.0023, 0.0046, 0.02];
params0.PieceWiseStrainDepParams = [10.0000, 10.8001, 2.8824, 1.1554, 0.015000];
params0.PieceWiseStrainDep2X = [-0.01, -0.0105, -0.0085, -0.0055, -0.004, -0.000, 0.006, 0.01, 0.0233, 0.1];
params0.PieceWiseStrainDep2Params = [50.0000, 50.0000, 50.6847, 1, 1, 0.5, 1, 20, 50, 50.0000];

% Use freshly constructed velocity table
params0.velocitytableonfile = 'protocol_03_27_2026_velocitytable.mat';

% Configuration
params0.PlotEachSeparately = 0;
params0.RunSlackSegments = 'All';
params0.RunSlack = true;
params0.RunForceVelocity = false;
params0.RunForceLengthEstim = false;
params0.BreakOnODEUnstable = true;

% Features to target
params0.fn = {'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};
params0.LoadPassiveMat = 'protocol_03_27_2026_PassiveCa_slack.mat';
% params0.LoadPassiveMat = '';
tic;RunBakersExp;
toc
%% 2. Define Parameters for Optimization
% Top identifiable parameters from previous sensitivity analysis:
% 1: k2, 2: kSE, 3: kstiff2, 4: k_pas, 5: dr, 6: k1, 7: sigma1, 8: kd
params0.mods = {'k2', 'kSE', 'kstiff2', 'k_pas', 'dr', 'k1', 'sigma1', 'kd'};
params0.g = ones(1, length(params0.mods)); 

%% 3. Run Initial Evaluation
disp('Evaluating Pre-Optimization Baseline...');
[InitResiduals, InitJac] = ResidualAndJacobian(params0.g, params0, true);
initCost = sum(InitResiduals.^2);
fprintf('Initial Cost: %.4f\n', initCost);

%% 4. Optimization
disp('Starting fminsearch Optimization...');
options = optimset('Display', 'iter', 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter', 500);

% Define objective function using internal Jacobian handling to prevent wasteful finite differences
costFun = @(g) sum(ResidualAndJacobian(g, params0, true).^2);

[g_opt, fval] = fminsearch(costFun, params0.g, options);

fprintf('\n=== Optimization Complete ===\n');
fprintf('Final Cost: %.4f\n', fval);
for i=1:length(params0.mods)
    fprintf('  %s multiplier: %.4f\n', params0.mods{i}, g_opt(i));
end

%% 5. Evaluate Post-Optimization & Plot
disp('Generating final plots...');
params0.PlotEachSeparately = 1;
params0.g = g_opt;
ResidualAndJacobian(g_opt, params0, true);
disp('Done!');
