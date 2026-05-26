% OptimizeMechanismEvaluation.m
% Optimizes the top K parameters identified by Sensitivity Analysis
% against the experimental data.

clear; close all; clc;

if isempty(gcp('nocreate'))
    p = parpool('threads', 5);
end
figure(102);clf;
%%

% Get the parent directory of your current folder
parentDir = fileparts(pwd); 

% Generate the path string for the parent and all its subfolders
allPaths = genpath(parentDir); 

% Add these paths to the MATLAB search path
addpath(allPaths);

% Optional: Save these changes for future sessions
% savepath; 
%% 1. Initialize Baseline Model
params0 = getParams();
ModelParamsWPassive_PreOpt


%% Configuration
params0.PlotEachSeparately = 1;
params0.RunSlackSegments = 'AllPar';
params0.RunSlack = true;
params0.RunForceVelocity = true;
params0.FV_velocities = -[0 0.5, 1, 2, 4];
params0.RunForceLengthEstim = false;
params0.BreakOnODEUnstable = true;
params0.UseFastPPval = true;

% Features to target
params0.fn = {'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};
params0.MaxRunTime = 10;
params0.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
%% Define the output metrics
params0 = getParams();
ModelParams_tuesdayLunch;
% ModelOptParams_TL_iter_4
ModelOptParams_TL3_iter_17
ModelOptimParam_TL5_iter_81
% ModelOptParams_TL3_iter_17_SRXstart_v5

params0.UseA2AttachmentShift = 1;
params0.slope = 100;

params0.UseStretchActivation = 1;
params0.k_SA = .1;

params0.UseSuperRelaxed = false;
params0.UseSuperRelaxedADP = false;
params0.SRXD_0 = 0.08;
params0.SRXT_0 = 0.08;
params0.RunForceVelocity = true;

params0.fn = {
    'FV_fnorm|FV_v|10', 'ktr|2', 'A|50', ...
    'ktr_rmse|0-0.2|.1', ...
    'XTOR[1]|2-8|10,XTOR_vmax[1]|4-30,SRX_ss[1]|0.1-0.8,attached_ss[1]|0.3-0.6|10', ...
    't0_crossing|SLdiff|2', ...
    'restretchSlopeStart|0.1',  ...
    'peak1_y|10', 'peak1_dSL|0.2', ...
    'vall_y|10', 'vall_t|0.2', 'peak2|5', ...
    'steady|50', 'vall2_dy|0.1', 'ovrsht_dy|0.1', ...
    'ovrsht_dy|0-2|1','ovrsht_t|0.4-100|100' ...
    'AssertParams|1'...
};  
%
params0.MaxRunTime = 100;
params0.justPlotStateTransitionsFlag = false;
figure(1);clf;
tic;RunBakersExp;
toc
%
sum(plotFeatures(features_data, features_model, [], params0.fn))
features_ghost = features_model;

%% old ghost
            % features_ghost.ktr                 = [39.76 34.18 31.63 28.36 27.6];
            % features_ghost.A                   = [68.4 65.47 58.92 52.2 43.51];
            % features_ghost.t0                  = [0.002165 0.004643 0.007818 0.01117 0.01589];
            % features_ghost.Am                  = [58.18 55.29 52.1 47.33 41.98];
            % features_ghost.SLslack             = [2.04 2 1.96 1.92 1.88];
            % features_ghost.SLdiff              = [0.162 0.202 0.242 0.282 0.322];
            % features_ghost.restretchSlopeStart = [1131 1011 1006 907 953.4];
            % features_ghost.restretchSlopeEnd   = [111.3 122.2 108.7 119.7 135.3];
            % features_ghost.v_restretch         = [4 5 6 7 3];
            % features_ghost.peak1_y             = [77.59 76.02 74.2 71.11 56.89];
            % features_ghost.peak1_t             = [0.0035 0.003 0.0025 0.002 0.005];
            % features_ghost.peak1_SL            = [2.068 2.034 1.99 1.954 1.912];
            % features_ghost.peak1_dSL           = [0.028 0.034 0.03 0.034 0.032];
            % features_ghost.vall_y              = [68.97 66.95 65.33 62.49 54.37];
            % features_ghost.vall_t              = [0.0085 0.006 0.006 0.005 0.0095];
            % features_ghost.peak2               = [77.22 78.47 80.01 82.42 61.85];
            % features_ghost.steady              = [76.91 76.48 75.98 75.53 56.75];
            % features_ghost.vall2_dy            = [-13.56 -13.84 -13.6 -13.08 -8.826];
            % features_ghost.vall2_t             = [0.0085 0.0055 0.008 0.0065 0.0095];
            % features_ghost.ovrsht_dy           = [0.431 1.392 1.244 1.552 0.3205];
            % features_ghost.ovrsht_t            = [0.2116 0.225 0.247 0.239 0.2676];
            % features_ghost.XTOR               = [10 10 10 10 10];

%% 2. Define Parameters for Optimization

% XB scalar rates
xb_rates = {'ka', 'kd', 'k1', 'k_1', 'k2', 'kah', 'kamh'};

% Force generation
force = {'kstiff1', 'kstiff2', 'dr', 'estiff'};

% Serial elastic element
se = {'kSE', 'ekSE', 'mu', 'mu_neg'};
% SE is fixed
se = {};

% Passive / titin viscoelastic
titin = {'k_pas', 'gamma', 'kSE_M', 'eta_M'};
% titin is fixed
titin = {};

% A2 hopping
if params0.UseA2AttachmentShift
    a2 = {'slope', 's_threshold_R'};
else 
    a2 = {};
end

% Filament overlap geometry
overlap = {'L_thick', 'L_hbare', 'L_thin'};
overlap = {};

% PWSD shape — p1→p2 (R12, multiplies k1)
% values at X = [-1  -0.0083  -0.0023  0.0046  0.02]
pw_k1_vals = {'PieceWiseStrainDepParams__2', 'PieceWiseStrainDepParams__3', 'PieceWiseStrainDepParams__4'};
pw_k1_x    = {'PieceWiseStrainDepX__4', ...
              'PieceWiseStrainDepX__3'};  % __1=-1 and __5=0.02 are far endpoints
pw_k1_x = {};


% PWSD shape — p1→PD (R1D, multiplies kd)
% values at X = [-0.05  -0.015  -0.004  0  0.004  0.01  0.1]
% __1=4000, __2=4000, __7=2000 are "clamp to rapid detachment" endpoints
pw_kd_vals = {'PieceWiseStrainDepR1DParams__3', 'PieceWiseStrainDepR1DParams__4', ...
              'PieceWiseStrainDepR1DParams__5', 'PieceWiseStrainDepR1DParams__6'};
pw_kd_x    = {'PieceWiseStrainDepR1DX__3', ...
              'PieceWiseStrainDepR1DX__5', ...
              'PieceWiseStrainDepR1DX__6'};  % __1=-0.05 and __7=0.1 are far endpoints
pw_kd_x = {};
pw_kd_vals = {'PieceWiseStrainDepR1DParams__3', 'PieceWiseStrainDepR1DParams__2'};


% PWSD shape — p2→PD (R2, multiplies k2)
% values at X = [-0.01 -0.0105 -0.0085 -0.0055 -0.004  0  0.006  0.01  0.0233  0.1]
% __1,__2,__9,__10 ≈ 50 are high-rate endpoints
pw_k2_vals = {'PieceWiseStrainDep2Params__2', 'PieceWiseStrainDep2Params__3', 'PieceWiseStrainDep2Params__4'};
pw_k2_x    = {'PieceWiseStrainDep2X__4', 'PieceWiseStrainDep2X__5', ...
              'PieceWiseStrainDep2X__7'};
pw_k2_x = {};

% % PWSD shape — p2→p1 (R21, multiplies k_1) — low impact, include if needed
% pw_k_1_vals = {'PieceWiseStrainDepR21Params__2', 'PieceWiseStrainDepR21Params__3'};
% pw_k_1_x    = {'PieceWiseStrainDepR21X__2', 'PieceWiseStrainDepR21X__3'};
        
offsets  = {'PieceWiseStrainDep2X_logOffset', 'PieceWiseStrainDepR1DX_logOffset', 'PieceWiseStrainDepR21X_logOffset', 'PieceWiseStrainDepX_logOffset'};

if params0.UseStretchActivation

    stretchAct = {'k_SA'};
else
    stretchAct  = {};
end

if params0.UseLatticeSpacing 
    lattice = {'d10_ref', 'R_thick', 'R_thin', 'd_optimal'};
else
    lattice = {};
end

% Super-relaxed (SRX): NR↔SR transition rates and strain sensitivity
if params0.UseSuperRelaxed
    sr = {'ksr0', 'kmsr', 'sigma1', 'sigma2'};
else
    sr = {};
end

% Super-relaxed ADP (SRD): SRD state rates
if params0.UseSuperRelaxedADP
    srd = {'kmsrd', 'sigma_srd1', 'sigma_srd2'};
else
    srd = {};
end


if params0.UseNegativeKstiff
    nonlinKstiff = {'kstiff2_n', 'kstiff1_n'};
else
    nonlinKstiff  = {};
end

othrs = {'MaxSlackNegativeForce', 'kah', 'kamh', 'xrate'};
%% Groups for random parameter draw
% Add/remove groups here to control the pool
active_groups = [xb_rates, force, se, titin, a2, overlap, ...
                 pw_k1_vals, pw_k1_x, pw_kd_vals, pw_kd_x, ...
                 pw_k2_vals, pw_k2_x, offsets, ...
                 sr, srd, othrs, nonlinKstiff, lattice, stretchAct];
active_groups = unique(active_groups, 'stable');
compulsory_params = {'ksr0', 'kmsr', 'sigma1', 'sigma2', 'ksrd', 'kmsrd', 'sigma_srd1', 'sigma_srd2'};
% compulsory_params = {'kstiff2_n', 'kstiff1_n', 'kstiff2', 'kstiff1', 'estiff'};
compulsory_params = {'k_SA', 'slope'};

% initial setup
% params0.kstiff1_n = params0.kstiff1;
% params0.kstiff2_n = params0.kstiff2;
% params0.kstiff1 = params0.kstiff2;
% params0.ksr2srd = params0.kah;
% params0.ksrd2sr = params0.kamh;
% params0.UseLatticeSpacing = true;
% 
% params0.UseNegativeKstiff = false;
% params0.UseSuperRelaxed = true;
% params0.UseSuperRelaxedADP = true;

% Top identifiable parameters from previous sensitivity analysis:
% 1: k2, 2: kSE, 3: kstiff2, 4: k_pas, 5: dr, 6: k1, 7: sigma1, 8: kd
params0.mods = {'ka', 'kd', 'k2', 'kstiff2', 'slope', 's_threshold_R', 'PieceWiseStrainDepR21X_logOffset', 'PieceWiseStrainDepR1DX_logOffset', 'PieceWiseStrainDep2X_logOffset'};
params0.g = ones(1, length(params0.mods));
optimTag = 'kSA-slope';


%% 3. Run Initial Evaluation
disp('Evaluating Pre-Optimization Baseline...');
params0.OptimizeOn = 'Feats';
params0.EvalFeatures = 1;
params0.BreakOnODEUnstable = 1;
params0.MaxRunTime = 30;
params0.RunSlackSegments = 'AllPar';
params0.UseFastPPval = true;

tic
[InitResiduals] = ResidualAndJacobian(params0.g, params0, true);
toc
initCost = sum(InitResiduals);
fprintf('Initial Cost: %.4f\n', initCost);

%% 4. Optimization
disp('Starting fminsearch Optimization...');
options = optimset('Display', 'iter', 'TolFun', 1e-3, 'TolX', 1e-2, 'MaxIter', 15, 'PlotFcns',@optimplotfval, 'MaxFunEvals', 60);
options = optimset('Display', 'iter', 'TolFun', 1e-3, 'TolX', 1e-2, 'MaxIter', 30, 'PlotFcns',@optimplotfval, 'MaxFunEvals', 90);
g_opt = params0.g; 
iter_num = 0;
fval = initCost;
% optimTag = 'TL6';
% compulsory_params = {'ksr0', 'kmsr', 'sigma1', 'sigma2', 'ksrd', 'kmsrd', 'sigma_srd1', 'sigma_srd2'};
% compulsory_params = {'ksr0', 'kmsr'};

for iter_num = 1:100    
    try
        params0 = draw10(params0, g_opt, active_groups, 10, compulsory_params);
        % Define objective function using internal Jacobian handling to prevent wasteful finite differences
        costFun = @(g) sum(ResidualAndJacobian(g, params0, true));
        fval_pre = fval;
        [g_opt, fval, a, optim_outputs] = fminsearch(costFun, params0.g, options);
        save("envOptIter4.mat", 'params0', 'g_opt');
        params0.g = g_opt;
        writeParamsToMFile(sprintf('../params/ModelOptimParam_%s_iter_%g.m', optimTag, iter_num), params0, [], ...
        sprintf("Iteration cost: %0.1f (from %0.1f) \r\n%% params0.mods = {'%s'};\r\n%% params0.g = [%s]", fval,fval_pre, strjoin(params0.mods, "', '"), num2str(g_opt)));
        captureOptimIter(params0, iter_num, fval, fval_pre, optimTag);
    catch e
        disp(e)
        disp("Going on...")
    end

end

fprintf('\n=== Optimization Complete ===\n');
fprintf('Final Cost: %.4f\n', fval);
for i=1:length(params0.mods)
    fprintf('  %s multiplier: %.4f\n', params0.mods{i}, g_opt(i));
end

%% 5. Evaluate Post-Optimization & Plot
disp('Generating final plots...');
parback = params0;
params0.g = g_opt;
params0 = getParams(params0, params0.g, false, true);
features_ghost = features_model;
%%
xxx = @(x, y) [x/x(1);y/y(1)];

xxx([1 2 3], [4 5 6])
StatesInTime
%% Test it 
figure(1); clf;
params0 = getParams();
% ModelOptParams_iter_27
% ModelOptParams_NonLinKstiff_iter_4
% ModelOptParams_SRXTD_iter_3

params0.PlotEachSeparately = 1;
% params0.MaxSlackNegativeForce = 0;
% % params0.UseMaxwellDashpot = 1;
% % params0.UseDynamicPassive = 1;
% % params0.kSE_M = 50;
% % params0.eta_M  = 0.05;
% % params0.mu  = 0.22;
% % params0.mu_neg = 0.22;
% params0.xrate = 3.05;
% params0.kstiff2 = 21000;
% params0.UseFastPPval = true;
params0.justPlotStateTransitionsFlag = false;
% params0.PieceWiseStrainDepR1DX_logOffset = 1;
% params0.PieceWiseStrainDep2X_logOffset = 1;
params0.PlotFeatureFitting = true;
% params0.fn = {'ktr_rmse|0-0.2', 'AssertParams|1'};
ModelOptParams_TL3_iter_17
ModelOptParams_TL3_iter_17_SRXstart_v5

params0.UseSuperRelaxed = 1;
params0.UseSuperRelaxedADP = 1;

tic
RunBakersExp;
toc
% ResidualAndJacobian(g_opt, params0, true);

% features_ghost = features_model;
sum(plotFeatures(features_data, features_model, features_ghost, params0.fn))
disp('Done!');
%% Experiment
params0 = getParams();
% features_ghost = features_model;
ModelOptParams_iter_27
params0.PlotEachSeparately = 1;
params0.RunForceVelocity = false;

params0.justPlotStateTransitionsFlag = false;
% params0.PieceWiseStrainDepR1DX_logOffset = 1;
% params0.PieceWiseStrainDep2X_logOffset = 1;

tic
RunBakersExp;
toc
%% ResidualAndJacobian(g_opt, params0, true);


sum(plotFeatures(features_data, features_model, features_ghost, params0.fn))
disp('Done!');


function params = draw10(params0, g, active_groups, n_draw, compulsory_groups)
    if nargin < 5
        compulsory_groups = {};
    end

    % Remove compulsory from the pool to avoid duplicates, then draw remainder
    pool = active_groups(~ismember(active_groups, compulsory_groups));
    n_rand = max(0, n_draw - numel(compulsory_groups));
    idx = randperm(numel(pool), min(n_rand, numel(pool)));
    drawn = [compulsory_groups, pool(idx)];

    params0.g = g;
    params = getParams(params0, params0.g, false, true);
    params.mods = drawn;
    params.g = ones(1, numel(params.mods));
    fprintf('Drew %d params: %s\n', numel(params.mods), strjoin(params.mods, ', '));
end
