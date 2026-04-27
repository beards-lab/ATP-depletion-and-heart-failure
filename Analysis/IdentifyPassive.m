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
params0.kstiff2 = 18000;
params0.kstiff1 = "=kstiff2";
params0.MaxSlackNegativeForce = -1;
params0.kSE = 1000;

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
% ModelParamsSlackKtrOpt;
params0.kSE = 3.2028e+03;
params0.ekSE = 0.9247;

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
params0.LoadPassiveMat = '03 27 2026 M/06_Merged_8mM_Active_PNB_Mava_loess_slack_sig.mat';
params0.LoadPassiveMat = '';
params0.MaxRunTime = 900;
params0.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
% params0.velocitytableonfile = 'protocol_03_27_2026_ActivePNBMava_slack.mat';
params0.PlotEachSeparately = 1;
% params0.ka = 100;
% params0.MaxSlackNegativeForce = -5;

ModelParams_Passive

% params.k_pas*max(SL-LSE - Lsc0, 0)^params.gamma; 
params0.UsePieceWiseStrainDep = 1;
params0.DryRun = 0;

% figure(215);clf;
params0.UseSuperRelaxed = 0;
params0.UseSuperRelaxedADP = 0;
clf;
tic;RunBakersExp;
toc
%%
features_data = extractSlackAttributes(out.datatable(:, 1), out.datatable(:, 3), out.datatable(:, 2)*2, velocitytable, struct(), out, true);
features_model = extractSlackAttributes(out.t, out.Force, out.SL, velocitytable, struct(), out, true);

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
%% old passive ghost
features_ghost = [];
features_ghost = features_model;
%%

plotFeatures(features_data, features_model, features_ghost, params0.fn)
%% 2. Define Parameters for Optimization
% Top identifiable parameters from previous sensitivity analysis:
% 1: k2, 2: kSE, 3: kstiff2, 4: k_pas, 5: dr, 6: k1, 7: sigma1, 8: kd
% params.k_pas*max(SL-LSE - Lsc0, 0)^params.gamma; 
params0.mods = {'k_pas', 'gamma', 'mu', 'kSE_M', 'eta_M', 'Lsc0'};
params0.g = ones(1, length(params0.mods)); 
% features_oldPassive = features_model
disp('Evaluating Pre-Optimization Baseline...');
params0.OptimizeOn = 'Time';
params0.EvalFeatures = 0;
params0.BreakOnODEUnstable = 1;
params0.MaxRunTime = 10;
params0.UseTitinIdentifiedPassive = 0;
params0.EvalPeaks = true;

[InitResiduals] = ResidualAndJacobian(params0.g, params0, true);
initCost = sum(InitResiduals);
fprintf('Initial Cost: %.4f\n', initCost);

%% 4. Optimization
disp('Starting fminsearch Optimization...');
options = optimset('Display', 'iter', 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter', 100, 'PlotFcns',@optimplotfval, 'MaxFunEvals', 100);

% Define objective function using internal Jacobian handling to prevent wasteful finite differences
costFun = @(g) sum(ResidualAndJacobian(g, params0, true));

[g_opt, fval] = fminsearch(costFun, params0.g, options);
params0.g = g_opt;
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

%%
clf;
params0.PlotEachSeparately = 1;
params0.UseMaxwellDashpot = 1;
params0.UseDynamicPassive = 1;
params0.UseTitinIdentifiedPassive = 0;

% params0.kSE_M = 73;
% params0.eta_M  = 0.027;
% params0.mu  = 0.0808*0.5;
% params0.mu_neg = "=mu";
% params0.gamma = 7.95;
% params0.k_pas = 15000;
% params0.Lsc0 = 1.8;
tic;RunBakersExp;toc;
E
% ResidualAndJacobian(g_opt, params0, true);
disp('Done!');

%
nexttile;
% params = params0;
LXB = 1.8:0.01:2.2;
F_passive = params.k_pas*max(LXB -  params.Lsc0, 0).^params.gamma; 
plot(LXB, F_passive, LXB, params0.k_pas*max(LXB -  params0.Lsc0, 0).^params0.gamma)

% writeParamsToMFile('../params/ModelParamsWPassive_PreOpt', params0, [])


