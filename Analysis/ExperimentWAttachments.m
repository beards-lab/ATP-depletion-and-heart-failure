figure(1);clf;
params0 = getParams();
% features_ghost = features_model;
ModelOptParams_iter_27

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
params0.FV_velocities = [0  -0.5000   -1.0000   -2.0000   -4.0000];
% params0.PieceWiseStrainDepR1DX_logOffset = 1;
params0.MaxRunTime = 60;
% params0.xrate = 0.25;

params0.alpha_vernier = 0.5;
params0.v_ref_vernier =0.5;

params0.UsePieceWiseStrainDep = true;

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
