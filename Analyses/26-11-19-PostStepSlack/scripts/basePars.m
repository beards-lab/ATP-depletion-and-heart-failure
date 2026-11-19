function pP = basePars(core)
% basePars  The incumbent as RunFVBeard2022 sets it up: robust1_seed + the fit-A k2 gate.
if nargin < 1, core = 'robust1_seed'; end
pP = getParams(loadParams(core), [], true, false);
pP.mods = {};  pP.g = ones(1, 0);
pP.EvalFeatures = false;  pP.BreakOnODEUnstable = false;  pP.MaxRunTime = 600;
pP.PlotEachSeparately = 0;  pP.PlotFeatureFitting = 0;  pP.PlotFullscreen = 0;
pP.ghostSave = '';  pP.ghostLoad = '';  pP.SaveBest = 0;
pP.UseAtpPowerGates = true;  pP.MgATP_ref = 8;  pP.UseAtpDetachR2 = 0;  pP.UseAtpKmsrd = 0;
pP.pAtp_k2 = log(0.54)/log(1/4);
end
