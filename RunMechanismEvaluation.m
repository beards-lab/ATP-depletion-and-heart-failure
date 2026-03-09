% Mechanism Evaluation Runner
% Implements Actions 1-3 from Docs/mechanism-evaluation.md
%
% Generates params0, applies the listed mechanism changes, and evaluates the model.

% p = parpool('threads', 5);
figure(101);clf;

% 1. Load base parameters
addpath('ModelOptParams');
addpath('Model');
addpath('Auxiliary');
ModelParamsSlackKtrOpt;

% === Action 1: Activate correct passive force exponent ===
params0.UseTitinIdentifiedPassive = 0;
params0.k_pas = 77*2*2;
params0.gamma = 7.9;

params0.UseOverlapFactor = false;

% === Action 2: Fix filament geometry ===
% Set to published values for cardiac muscle
params0.L_thick = 1.65;
params0.L_thin = 1.15;

% === Action 3: Replace negative stiffness asymmetry with fast negative-strain detachment ===
params0.UseNegativeKstiff = 0; % Disable unphysical asymmetry

params0.kstiff2 = 50000;
params0.kstiff1 = params0.kstiff2;
params0.MaxSlackNegativeForce =-1;
params0.kSE = 1000;
% Increase baseline kd
% params0.kd = params0.kd * 2; 


% Make PCHIP for R1D asymmetric with steep rise at s < -dr
% old 
% params0.PieceWiseStrainDepR1DX = [-0.0500   -0.0050  0.0063    0.0500];
% params0.PieceWiseStrainDepR1DParams = [50.0000    10.1474   1.5953   50.0000];


params0.ka = 67*2;
params0.kd = 9.55*2;
params0.k1 = 270*2;
params0.k2 = 133;
params0.kah = 30;

params0.A1AttachmentWidth = 0.002;

% R1D
params0.PieceWiseStrainDepR1DX = [-0.05, -0.015, -0.002, 0, 0.006, 0.05] + 0e-3;
params0.PieceWiseStrainDepR1DParams = [5000, 5000, 10.147, 1, 1.595, 50]*4; 
params0.dr2 = 0.01;
params0.mu = 0.633*0.5;
params0.mu_neg = "=mu";

% R1
params0.PieceWiseStrainDepX = [0.1000    0.0079         0   -0.0050   -1.0000] - params0.dr/3;
params0.PieceWiseStrainDepParams = 1*[0.5000    1.1554    2.8824   10.8001   10.0000];

% R2
params0.PieceWiseStrainDep2X = [-0.0100   -0.0080   -0.0079         0    0.0127    0.0258    0.1000]-params0.dr/4;
params0.PieceWiseStrainDep2Params = [50.0000   50.0000   50.6847    1.0383    0.9402    2.9058   50.0000];

% R21
% params0.PieceWiseStrainDepR21X = [-0.0580   -0.0049    0.0062    0.0500];
% params0.PieceWiseStrainDepR21Params = [50.0000    0.9940    0.9299   50.0000];
% disable for now
params0.PieceWiseStrainDepR21X = [-1 1];
params0.PieceWiseStrainDepR21Params = [1 1];

% Let RunBakersExp know we want to see the plots
params0.PlotEachSeparately = 1;
params0.PlotFullscreen = 0;
params0.RunSlackSegments = 'AllPar';
params0.ghostSave = '';
params0.ghostLoad = '';
params0.xrate = 1;


% params0.justPlotStateTransitionsFlag = true;
params0.RunForceVelocity = false;
params0.UseSuperRelaxed = false;
params0.UseSuperRelaxedADP = false;

params0.dS = 0.001;

% Load requisite data variables (e.g. ATP_c)
LoadData;

% Run the experiment scripts
RunBakersExp;
ylim([-20, 100]);
return
%%
params0.fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};
params0.RunForceVelocity = false;
params0.fn = {'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};

% params0.fn = {'restretchSlopeStart'};
plotFeatures(features_data, features_model, [], params0.fn);

totalCost = sum(evalFeatureCost(features_data, features_model, params0.fn, 1));
