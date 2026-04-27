% OptimizeMechanismEvaluation.m
% Optimizes the top K parameters identified by Sensitivity Analysis
% against the experimental data.

clear; close all; clc;

if isempty(gcp('nocreate'))
    p = parpool('threads', 5);
end
figure(102);clf;

%% 1. Initialize Baseline Model
params0 = getParams();
ModelParamsWPassive_PreOpt


%% Configuration
params0.PlotEachSeparately = 1;
params0.RunSlackSegments = 'AllPar';
params0.RunSlack = true;
params0.RunForceVelocity = true;
params0.FV_velocities = -[0.5, 1, 2, 4];
params0.RunForceLengthEstim = false;
params0.BreakOnODEUnstable = true;
params0.UseFastPPval = true;

% Features to target
params0.fn = {'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};
params0.MaxRunTime = 900;
params0.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
tic;RunBakersExp;
toc
%% old ghost
            features_ghost.ktr                 = [39.76 34.18 31.63 28.36 27.6];
            features_ghost.A                   = [68.4 65.47 58.92 52.2 43.51];
            features_ghost.t0                  = [0.002165 0.004643 0.007818 0.01117 0.01589];
            features_ghost.Am                  = [58.18 55.29 52.1 47.33 41.98];
            features_ghost.SLslack             = [2.04 2 1.96 1.92 1.88];
            features_ghost.SLdiff              = [0.162 0.202 0.242 0.282 0.322];
            features_ghost.restretchSlopeStart = [1131 1011 1006 907 953.4];
            features_ghost.restretchSlopeEnd   = [111.3 122.2 108.7 119.7 135.3];
            features_ghost.v_restretch         = [4 5 6 7 3];
            features_ghost.peak1_y             = [77.59 76.02 74.2 71.11 56.89];
            features_ghost.peak1_t             = [0.0035 0.003 0.0025 0.002 0.005];
            features_ghost.peak1_SL            = [2.068 2.034 1.99 1.954 1.912];
            features_ghost.peak1_dSL           = [0.028 0.034 0.03 0.034 0.032];
            features_ghost.vall_y              = [68.97 66.95 65.33 62.49 54.37];
            features_ghost.vall_t              = [0.0085 0.006 0.006 0.005 0.0095];
            features_ghost.peak2               = [77.22 78.47 80.01 82.42 61.85];
            features_ghost.steady              = [76.91 76.48 75.98 75.53 56.75];
            features_ghost.vall2_dy            = [-13.56 -13.84 -13.6 -13.08 -8.826];
            features_ghost.vall2_t             = [0.0085 0.0055 0.008 0.0065 0.0095];
            features_ghost.ovrsht_dy           = [0.431 1.392 1.244 1.552 0.3205];
            features_ghost.ovrsht_t            = [0.2116 0.225 0.247 0.239 0.2676];
            features_ghost.XTOR               = [10 10 10 10 10];

%% 2. Define Parameters for Optimization

% XB scalar rates
xb_rates = {'ka', 'kd', 'k1', 'k_1', 'k2', 'kah', 'kamh'};

% Force generation
force = {'kstiff1', 'kstiff2', 'dr', 'estiff'};

% Serial elastic element
se = {'kSE', 'ekSE', 'mu', 'mu_neg'};

% Passive / titin viscoelastic
titin = {'k_pas', 'gamma', 'kSE_M', 'eta_M'};

% A2 hopping
a2 = {'slope', 's_threshold_R'};

% Filament overlap geometry
overlap = {'L_thick', 'L_hbare', 'L_thin'};

% PWSD shape — p1→p2 (R12, multiplies k1)
% values at X = [-1  -0.0083  -0.0023  0.0046  0.02]
pw_k1_vals = {'PieceWiseStrainDepParams__3', 'PieceWiseStrainDepParams__4'};
pw_k1_x    = {'PieceWiseStrainDepX__4', ...
              'PieceWiseStrainDepX__3'};  % __1=-1 and __5=0.02 are far endpoints

% PWSD shape — p1→PD (R1D, multiplies kd)
% values at X = [-0.05  -0.015  -0.004  0  0.004  0.01  0.1]
% __1=4000, __2=4000, __7=2000 are "clamp to rapid detachment" endpoints
pw_kd_vals = {'PieceWiseStrainDepR1DParams__3', 'PieceWiseStrainDepR1DParams__4', ...
              'PieceWiseStrainDepR1DParams__5', 'PieceWiseStrainDepR1DParams__6'};
pw_kd_x    = {'PieceWiseStrainDepR1DX__3', ...
              'PieceWiseStrainDepR1DX__5', ...
              'PieceWiseStrainDepR1DX__6'};  % __1=-0.05 and __7=0.1 are far endpoints

% PWSD shape — p2→PD (R2, multiplies k2)
% values at X = [-0.01 -0.0105 -0.0085 -0.0055 -0.004  0  0.006  0.01  0.0233  0.1]
% __1,__2,__9,__10 ≈ 50 are high-rate endpoints
pw_k2_vals = {'PieceWiseStrainDep2Params__3', 'PieceWiseStrainDep2Params__4', ...
              'PieceWiseStrainDep2Params__5', 'PieceWiseStrainDep2Params__6', ...
              'PieceWiseStrainDep2Params__7', 'PieceWiseStrainDep2Params__8'};
pw_k2_x    = {'PieceWiseStrainDep2X__2', 'PieceWiseStrainDep2X__3', ...
              'PieceWiseStrainDep2X__4', 'PieceWiseStrainDep2X__5', ...
              'PieceWiseStrainDep2X__6', 'PieceWiseStrainDep2X__7', ...
              'PieceWiseStrainDep2X__8', 'PieceWiseStrainDep2X__9'};

% PWSD shape — p2→p1 (R21, multiplies k_1) — low impact, include if needed
pw_k_1_vals = {'PieceWiseStrainDepR21Params__2', 'PieceWiseStrainDepR21Params__3'};
pw_k_1_x    = {'PieceWiseStrainDepR21X__2', 'PieceWiseStrainDepR21X__3'};

offsets  = {'PieceWiseStrainDep2X_logOffset', 'PieceWiseStrainDepR1DX_logOffset', 'PieceWiseStrainDepR21X_logOffset', 'PieceWiseStrainDepX_logOffset'};

%% Random parameter draw
% Add/remove groups here to control the pool
active_groups = [xb_rates, force, se, titin, a2, overlap, ...
                 pw_k1_vals, pw_k1_x, pw_kd_vals, pw_kd_x, ...
                 pw_k2_vals, pw_k2_x, pw_k_1_vals, pw_k_1_x, offsets];
active_groups = unique(active_groups, 'stable');

n_draw = 10;
idx = randperm(numel(active_groups), min(n_draw, numel(active_groups)));
params0.mods = active_groups(idx);
fprintf('Drew %d params: %s\n', numel(params0.mods), strjoin(params0.mods, ', '));

% connections
params0.kstiff1 = "=kstiff2";
params0.UseLatticeSpacing = true;

lattice = {'d10_ref', 'R_thick', 'R_thin', 'd_optimal'};


% Top identifiable parameters from previous sensitivity analysis:
% 1: k2, 2: kSE, 3: kstiff2, 4: k_pas, 5: dr, 6: k1, 7: sigma1, 8: kd
params0.mods = {'ka', 'kd', 'k2', 'kstiff2', 'slope', 's_threshold_R', 'PieceWiseStrainDepR21X_logOffset', 'PieceWiseStrainDepR1DX_logOffset', 'PieceWiseStrainDep2X_logOffset'};
params0.g = ones(1, length(params0.mods)); 


%% Define the output metrics
params0.fn = {'FV_f|FV_v|@(x,y)sum((x/65 - y/56.4).^2)', ...
'ktr|2' , ...               
'A|50' , ...                 
't0|5' , ...                
'ktr_rmse|0' , ...          
'restretchSlopeStart|0.1', ...
'restretchSlopeEnd|0.1' , ... 
'peak1_y|10' , ...   
'peak1_dSL|0.2' , ...                    
'vall_y|10' , ...            
'vall_t|0.2' , ...            
'peak2|5' , ...             
'steady|50' , ...            
'vall2_dy|0.1' , ...          
'ovrsht_dy|0.1' , ...              
'XTOR|@(X, Y_data) 0.01*sum((max(0, 2 - X) + max(0, X - 8)).^2)'}          

% c = @(X, Y_data) 0.1*sum((max(0, 2 - X) + max(0, X - 8)).^2);
plotFeatures(features_data, features_model, features_ghost, params0.fn)
%% 3. Run Initial Evaluation
disp('Evaluating Pre-Optimization Baseline...');
params0.OptimizeOn = 'Feats';
params0.EvalFeatures = 1;
params0.BreakOnODEUnstable = 1;
params0.MaxRunTime = 10;
params0.RunSlackSegments = 'AllPar';
params0.UseFastPPval = true;

tic
[InitResiduals] = ResidualAndJacobian(params0.g, params0, true);
toc
initCost = sum(InitResiduals);
fprintf('Initial Cost: %.4f\n', initCost);

%% 4. Optimization
disp('Starting fminsearch Optimization...');
options = optimset('Display', 'iter', 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter', 100, 'PlotFcns',@optimplotfval, 'MaxFunEvals', 100);

% Define objective function using internal Jacobian handling to prevent wasteful finite differences
costFun = @(g) sum(ResidualAndJacobian(g, params0, true));

[g_opt, fval] = fminsearch(costFun, params0.g, options);
params0.g = g_opt;
save envOpt3;
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
figure(1); clf;
params0.PlotEachSeparately = 1;
% params0.UseMaxwellDashpot = 1;
% params0.UseDynamicPassive = 1;
% params0.kSE_M = 50;
% params0.eta_M  = 0.05;
% params0.mu  = 0.22;
% params0.mu_neg = 0.22;
params0.xrate = 3.05;
params0.kstiff2 = 21000;
params0.UseFastPPval = true;
params0.justPlotStateTransitionsFlag = true;
% params0.PieceWiseStrainDep2X_logOffset = 0.98;

tic
RunBakersExp;
toc
% ResidualAndJacobian(g_opt, params0, true);


sum(plotFeatures(features_data, features_model, features_ghost, params0.fn))
disp('Done!');
