% Mechanism Evaluation Runner
% Implements Actions 1-3 from Docs/mechanism-evaluation.md
%
% Generates params0, applies the listed mechanism changes, and evaluates the model.

% 1. Load base parameters
addpath('ModelOptParams');
addpath('Model');
addpath('Auxiliary');
ModelParamsSlackKtrOpt;

% === Action 1: Activate correct passive force exponent ===
params0.UseTitinIdentifiedPassive = 1;

% === Action 2: Fix filament geometry ===
% Set to published values for cardiac muscle
params0.L_thick = 1.65;
params0.L_thin = 1.15;

% === Action 3: Replace negative stiffness asymmetry with fast negative-strain detachment ===
params0.UseNegativeKstiff = 0; % Disable unphysical asymmetry
% Increase baseline kd
params0.kd = params0.kd * 2; 
% Make PCHIP for R1D asymmetric with steep rise at s < -dr
params0.PieceWiseStrainDepR1DX = [-0.05, -0.015, -0.005, 0.006, 0.05];
params0.PieceWiseStrainDepR1DParams = [100, 50, 1.147, 1.595, 50]; 

% Let RunBakersExp know we want to see the plots
params0.PlotEachSeparately = 1;
params0.PlotFullscreen = 0;
params0.ghostSave = '';
params0.ghostLoad = '';

% Load requisite data variables (e.g. ATP_c)
LoadData;

% Run the experiment scripts
RunBakersExp;

params0.fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};
% params0.fn = {'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};

% params0.fn = {'restretchSlopeStart'};
plotFeatures(features_data, features_model, [], params0.fn);

totalCost = sum(evalFeatureCost(features_data, features_model, params0.fn, 1));
