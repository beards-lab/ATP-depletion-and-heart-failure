params0 = getParams();
ModelParams_tuesdayLunch;
% ModelOptParams_TL_iter_4
ModelOptParams_TL3_iter_17
% ModelOptParams_TL3_iter_17_SRXstart_v5
addpath(genpath('../'));

params0.UseSuperRelaxed = false;
params0.UseSuperRelaxedADP = false;
% params0.SRXD_0 = 0.08;
% params0.SRXT_0 = 0.08;

params0.MaxRunTime = 100;
params0.justPlotStateTransitionsFlag = false;
RunBakersExp;
