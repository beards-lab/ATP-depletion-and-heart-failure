% driver code
g = ones(30, 1);

params = getParams();
params.g = g;
params.SL0 = 2.2;
% params.Slim = 0.18;
params.Slim = 0.3;
params.dS = 10e-3;
% params.N = 30;
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

% ghostSave = 'beardsOrig50'; % strain shifting dS 50 nm
% ghostSave = 'beardsOrig5'; % strain shifting dS 5 nm
% ghostSave = 'beardsOrig5_BW';%bells and whistles :)

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
ghostLoad = 'beardsOrig5';% N = 40
ghostLoad = 'beardsOrig5_BW';% N = 40
% ghostLoad = '';

% testing setup
% params.UseOverlap = true;
% params.UsePassive = true;

params.UseTORNegShift = false;
params.UseMutualPairingAttachment = false;
params.UseSlack = true;
params.UseOverlap = true;
params.UsePassive = true;
params.UseSerialStiffness = true;

% params.alpha1 = 0;
% params.alpha2 = 0;

% need a ksttiff1 and kstiff2 parameter retune
params.F_act_UseP31 = true;
params.UseP31Shift = true;


params.PlotEachSeparately = true;
params.ghostLoad = ghostLoad;
params.ghostSave = ghostSave;

% save as default
params0 = getParams(params);
params0.N

%% init
LoadData;
t_ss = [0 1];
t_sl0 = [0 0.1];
tic

RunBakersExp;
toc

%% OPTIM
% return
params0.mods = {"kstiff1", "kstiff2"};

options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.1, 'PlotFcns', @optimplotfval, 'MaxIter', 1500);
g = [1, 1];
optimfun = @(g)evaluateBakersExp(g, params0);
x = fminsearch(optimfun, g, options)
g = x;
