figure(1);clf;
params0 = getParams();
% features_ghost = features_model;
% ModelOptParams_iter_27
ModelParamsFeats_FVSlackUpdateValley2

params0.PlotEachSeparately = 1;
params0.RunForceVelocity = true;

params0.justPlotStateTransitionsFlag = false;
% params0.PieceWiseStrainDepR1DX_logOffset = 1;
% params0.PieceWiseStrainDep2X_logOffset = 1;

params0.BreakOnODEUnstable = false;
% params0.UseVernierVelocity = true;
% params0.alpha_vernier = 1.5;
% params0.v_ref_vernier = 0.5;
% params0.ghostSave = 'baseline';
params0.ghostSave = '';
params0.ghostLoad = 'baseline';
params0.RunSlackSegments = 'AllPar';
params0.FV_velocities = [0  -0.5000   -1.0000   -2.0000   -4.0000];
% params0.PieceWiseStrainDepR1DX_logOffset = 1;
params0.MaxRunTime = 360;
% params0.xrate = 0.25;
params0.fn = {'FV_f|FV_v', 'ktr|2', 'A|50', 't0|5', 'ktr_rmse|0', 'restretchSlopeStart|0.1', 'restretchSlopeEnd|0.1', 'peak1_y|10', 'peak1_dSL|0.2', 'vall_y|10', 'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|0.1', 'ovrsht_dy|0.1', 'XTOR[3]|2-8', 'SRX_ss[3]|0.25-0.75,attached_ss[3]|0.02-0.45,PT_ss[3]|0.02-0.70'};
params0.OptimizeOn = 'Feats';    
params0.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
params0.RunSlackSegments = 'All';

params0.kstiff2 = 101185*1.3;
params0.UseVernierVelocity = false;
params0.alpha_vernier = 0.5;
params0.v_ref_vernier =0.5;

% params0.UsePieceWiseStrainDep = true;

% params0.PieceWiseStrainDepR1DX = [-0.05 -0.015 -0.00366440431163266 0 0.00416212250313296 0.0104827633189078 0.1];
% params0.PieceWiseStrainDepR1DParams = [150 150 43.8911085627699 4.08490160249196 8.53555029444009 150 150];


tic
RunBakersExp;
toc
% ResidualAndJacobian(g_opt, params0, true);

% features_ghost = features_model;
sum(plotFeatures(features_data, features_model, features_ghost, params0.fn))
disp('Done!');

return
%%
% v_hs = -params0.FV_velocities;
v_hs = [0:0.5:10];

params0.alpha_vernier = 0.9;
params0.v_ref_vernier =1;


f_sat = 1 + params0.alpha_vernier .* v_hs ./ (v_hs + params0.v_ref_vernier);

f_sat = params0.alpha_vernier + (1 - params0.alpha_vernier)./(1+params0.v_ref_vernier./v_hs);

figure(546);clf;
plot(f_sat*100, v_hs)
xlabel('Augmentation (%)');ylabel('v (um/s)');



% matchStructFields
%%
