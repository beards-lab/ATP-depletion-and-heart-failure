% driver code
g = ones(30, 1);

params = getParams();
params.g = g;
params.SL0 = 2.2;
% params.Slim = 0.18;
params.Slim = 0.3;
params.dS = 10e-3;
params.N = 30;
params.MgATP = 8;

figure(12);clf;
ghostSave = '';
% ghostSave = 'beardsOrig_passive';
% ghostSave = 'ShiftingStrain40_Slim0_3';
% ghostSave = 'ShiftingStrain80_Slim0_3';
ghostSave = 'operatorSplittingPU020';
ghostSave = 'operatorSplittingPU0';
ghostSave = 'operatorSplittingPU080';
ghostSave = 'ShiftingStrain160_Slim0_3';
ghostSave = 'beardsOrig_all60';
ghostSave = 'beardsOrig_OV_Pas_SS';
ghostSave = '';

ghostLoad = '';
% ghostLoad = 'beardsOrig';
% ghostLoad = 'ShiftingStrain40';
% ghostLoad = 'operatorSplitting'
ghostLoad = 'operatorSplittingPU020';
ghostLoad = 'operatorSplittingPU0';% N = 40
ghostLoad = 'operatorSplittingPU080';
ghostLoad = 'ShiftingStrain40_Slim0_3';
ghostLoad = 'ShiftingStrain80_Slim0_3';
ghostLoad = 'ShiftingStrain160_Slim0_3';
ghostLoad = 'ShiftingStrainTest20';
ghostLoad = 'ShiftingStrainTest40';
% ghostLoad = 'ShiftingStrainTest80';
% ghostLoad = 'ShiftingStrainTest160';
% ghostLoad = 'ShiftingStrainTest240';
% ghostLoad = 'beardsOrig_all30';
ghostLoad = 'beardsOrig_all60';
ghostLoad = 'operatorSplittingPU0';% N = 40
% ghostLoad = '';

% testing setup
% params.UseOverlap = true;
% params.UsePassive = true;

params.UseTORNegShift = false;
params.UseMutualPairingAttachment = false;
params.UseOverlap = false;
params.UsePassive = false;
params.UseSerialStiffness = false;
% set as a default to be modified
% params0 = params;
params.PlotEachSeparately = true;
params.ghostLoad = ghostLoad;
params.ghostSave = ghostSave;
%% init
LoadData;
t_ss = [0 1];
t_sl0 = [0 0.1];
tic

RunBakersExp;
toc
