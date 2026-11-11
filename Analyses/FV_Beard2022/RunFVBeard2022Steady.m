% RunFVBeard2022Steady.m -- the model's own step-free F-V curve at 8, 4 and 2 mM.
%
% Same model as RunFVBeard2022.m (robust1_seed + fit-A k2 gate), run through the
% standard runFVExperiment protocol (1 s isometric at SL 2.2, then a pure ramp
% to SL 2.0, force read at the end of the ramp; v = 0 read at SL 2.0). It has
% NO pre-ramp step, so it separates the cross-bridge F-V from the series-
% compliance mismatch the measured step-then-ramp protocol exposes.
%
%   matlab -batch "run('Analyses/FV_Beard2022/RunFVBeard2022Steady.m')"
%
% Output -> results/steady.mat  (vFine, Fs = 3 x numel(vFine), atpS = [8 4 2])

rootP = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(rootP);
addpath(genpath(rootP));
wt = strsplit(genpath(fullfile(rootP, '.claude')), pathsep); wt = wt(~cellfun(@isempty, wt));
if ~isempty(wt), rmpath(wt{:}); end
outP = fullfile(rootP, 'Analyses', 'FV_Beard2022', 'results');
if ~exist('CORE', 'var'); CORE = 'robust1_seed'; end

vFine = [0 0.1 0.25 0.5:0.5:7];
LoadData;                                               % Data_ATP (only passed through)
pP = getParams(loadParams(CORE), [], true, false);
pP.mods = {}; pP.g = ones(1,0);
pP.EvalFeatures = true; pP.BreakOnODEUnstable = false; pP.MaxRunTime = 900;
pP.PlotEachSeparately = 0; pP.PlotFeatureFitting = 0; pP.PlotFullscreen = 0;
pP.ghostSave = ''; pP.ghostLoad = ''; pP.SaveBest = 0;
pP.UseAtpPowerGates = true; pP.MgATP_ref = 8; pP.UseAtpDetachR2 = 0; pP.UseAtpKmsrd = 0;
pP.pAtp_k2 = log(0.54)/log(1/4);
q = pP; q.FV_velocities = -vFine; q.OptimizeFVInit = true;
q.RunSlack = false; q.RunStairs = false; q.RunKtr = false; q.RunForceVelocity = true;

atpS = [8 4 2]; Fs = nan(numel(atpS), numel(vFine));
for i = 1:numel(atpS)
    qi = q; qi.MgATP = atpS(i);
    % runFVExperiment takes [ATP] from the FIRST entry of its 2nd argument
    t0 = tic; [~, o] = runFVExperiment(qi, [atpS(i) setdiff([8 4 2], atpS(i))], Data_ATP);
    Fs(i, :) = arrayfun(@(x) x.Force(end), o);
    fprintf('steady FV %d mM: %.0f s, F(0) %.2f kPa\n', atpS(i), toc(t0), Fs(i, 1));
end
save(fullfile(outP, 'steady.mat'), 'vFine', 'Fs', 'atpS', 'CORE');
fprintf('wrote %s\n', fullfile(outP, 'steady.mat'));
