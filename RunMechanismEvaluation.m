% Mechanism Evaluation Runner
% Implements Actions 1-3 from Docs/mechanism-evaluation.md
%
% Generates params0, applies the listed mechanism changes, and evaluates the model.

%% Test lattice spacing
SL = 1.8:0.01:2.4;
d10 = params.d10_ref .* sqrt(params.SL_ref_lattice ./ SL);
d_surface = d10 - params.R_thick - params.R_thin;
excess = max(0, d_surface - params.d_optimal); % only penalize when lattice is wider than optimal
f_lattice = exp(-(excess.^2) ./ (2 * params.sigma_lattice^2));
plot(SL, f_lattice)

%%

if isempty(gcp('nocreate'))
    p = parpool('threads', 5);
end
figure(101);clf;

% 1. Load base parameters
addpath('ModelOptParams');
addpath('Model');
addpath('Auxiliary');
ModelParamsSlackKtrOpt;



% === Action 1: Activate correct passive force exponent ===
params0.UseTitinIdentifiedPassive = 0;
params0.k_pas = 77*2*2*2;
params0.gamma = 7.9;

params0.UseOverlap = true;

params0.UseOverlapFactor = true;

% === Action 2: Fix filament geometry ===
% Set to published values for cardiac muscle
params0.L_thick = 1.65;
params0.L_thin = 1.15;

% === Action 3: Replace negative stiffness asymmetry with fast negative-strain detachment ===
params0.UseNegativeKstiff = 0; % Disable unphysical asymmetry

params0.kstiff2 = 40000;
params0.kstiff1 = "=kstiff2";
params0.MaxSlackNegativeForce =-1;
params0.kSE = 2000;
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
params0.kah = 30*2;

params0.A1AttachmentWidth = 0.004;

% R1D
params0.PieceWiseStrainDepR1DX = [-0.05, -0.015, -0.001, 0, 0.006, 0.05] + 0e-3;
params0.PieceWiseStrainDepR1DParams = [5000, 5000, 10.147, 1, 1.595, 50]*4; 
params0.dr2 = "=dr";
params0.mu = 0.633*0.5;
params0.mu_neg = "=mu";

% R1
params0.PieceWiseStrainDepX =        [-1   -0.0083   -0.0023    0.0046    0.0967];
params0.PieceWiseStrainDepParams = 1*[10.0000   10.8001    2.8824    1.1554    0.5000];

% R2
params0.PieceWiseStrainDep2X = [-0.01   -0.0105   -0.0085   -0.0055 -0.004  -0.000    0.0233    0.1];
params0.PieceWiseStrainDep2Params = [50.0000   50.0000   50.6847    1 1         1    2.9058   50.0000];

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
params0.RunSlackSegments = 'FirstAndLastExtended';
params0.RunSlack = true;
% params0.RunSlackSegments = 'stairs-down';
params0.ghostSave = 'oj';
params0.ghostLoad = 'oj';
params0.xrate = 1;
% params0.justPlotStateTransitionsFlag = true;
params0.FV_velocities = -[0, 0.5, 1, 3, 4];
params0.RunForceVelocity = true;
params0.UseSuperRelaxed = true;
params0.UseSuperRelaxedADP = true;

params0.UseA2AttachmentShift = false;
params0.BreakOnODEUnstable = false;

params0.dS = 0.001;
params0.MaxStrainArraySize = 40;
params0.UseLatticeSpacing = true;

params0.sigma_lattice = 0.0006;

% Load requisite data variables (e.g. ATP_c)
LoadData;

% Run the experiment scripts
tic
RunBakersExp;
toc
% ylim([-20, 100]);
% return
%%
params0.fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};
params0.RunForceVelocity = false;
params0.fn = {'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};
params0.fn = {'FV_f|FV_v'};
% params0.fn = {'restretchSlopeStart'};
plotFeatures(features_data, features_model, [], params0.fn);

totalCost = sum(evalFeatureCost(features_data, features_model, params0.fn, 1))

%% frank -starling investigation??
steps = [0 2 3.5, 5, 7];
% SL_Steps = interp1(out.t, out.SL, steps, 'nearest');
% ForceSteps_SrxOFLat = interp1(out.t, out.Force, steps, 'nearest');
% ForceSteps_SrxOF = interp1(out.t, out.Force, steps, 'nearest');
% ForceSteps_OF = interp1(out.t, out.Force, steps, 'nearest');
% ForceSteps_OV = interp1(out.t, out.Force, steps, 'nearest');
% ForceSteps_nonOV = interp1(out.t, out.Force, steps, 'nearest');
% ForceSteps_passive = interp1(out.t, out.FXBPassive, steps, 'nearest');

clf; hold on;

plot(SL_Steps, ForceSteps_SrxOFLat, 'x-', LineWidth=2);
plot(SL_Steps, ForceSteps_SrxOF, 'x-', LineWidth=2);
plot(SL_Steps, ForceSteps_OF, 'x-', LineWidth=2);
plot(SL_Steps, ForceSteps_OV, 'x-', LineWidth=2);
plot(SL_Steps, ForceSteps_nonOV, 'x-', LineWidth=2);
plot(SL_Steps, ForceSteps_nonOV(1) - ForceSteps_passive(1) + ForceSteps_passive, 'x-', LineWidth=2);
legend('ForceSteps_SrxOFLat','ForceSteps_SrxOF','ForceSteps_OF','ForceSteps_OV','ForceSteps_nonOV', interpreter='none');
% ForceSteps_SrxOFLat = interp1(out.t, out.Force, steps, 'nearest');
% ForceSteps_SrxOF = interp1(out.t, out.Force, steps, 'nearest');
% ForceSteps_OF = interp1(out.t, out.Force, steps, 'nearest');
% ForceSteps_OV = interp1(out.t, out.Force, steps, 'nearest');
% ForceSteps_nonOV = interp1(out.t, out.Force, steps, 'nearest');