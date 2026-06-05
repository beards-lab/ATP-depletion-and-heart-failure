params0 = getParams();
ModelParams_tuesdayLunch;
% ModelOptParams_TL_iter_4
ModelOptParams_TL3_iter_17
% ModelOptimParam_TL5_iter_81.m
% ModelOptParams_TL3_iter_17_SRXstart_v5
SRX_FastGentle;
% addpath(genpath('../'));

params0.ksr2srd = "=kah";
params0.ksrd2sr = "=kah/10";
params0.kah = 100;
params0.kamh = "=kah/10";


% --- Parking (into SRX): sets the pool size via ratio to mobilization ---
params0.ksr0  = 500;    % PT  -> SR
params0.kmsr  = 0;   % SR  -> PT

params0.kmsrd = 40;   % SRD -> PD
params0.ksrd  = 0;    % PD  -> SRD

params0.ksr2srd = "=kah";   % SR→SRD mixing = kah (100)
params0.ksr0    = 190;      % PT→SR parking rate  — sets pool to ~0.20 at max force
params0.kmsrd   = 60;       % SRD→PD outflow      — shapes SR/SRD split

% params0.kstiff2 = 6.9794e+04*0.95;

params0.RunForceVelocity = false;

params0.UseSuperRelaxed = true;
params0.UseSuperRelaxedADP = true;
   
% params0.SRXD_0 = 0.08;
% params0.SRXT_0 = 0.08;

params0.MaxRunTime = 100;
params0.justPlotStateTransitionsFlag = false;
tic
RunBakersExp;
toc
%%
if ~exist("features_ghost", 'var')
    features_ghost = features_model;
end
sum(plotFeatures(features_data, features_model, features_ghost, params0.fn))

%%
StatesInTime