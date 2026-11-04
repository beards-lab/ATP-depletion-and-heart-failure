% RunPooledFigureSims.m -- the two simulations the poster-style figures need
% that are not already saved by Analyses/LowATP_Refit9K2Sweep:
%   1. a fine force-velocity sweep (0 : 0.25 : 7 ML/s) at 8 and 2 mM, so the
%      F-V and power panels show a curve, not the fit's five points;
%   2. the passive (PNB + Mava) slack protocol on the current core.
% Slack, ktr-protocol and feature values are reused from the saved runs of fit A
% (bestfit_k2only.mat, levers_ref.mat, k2sweep.mat).
%
% Model = the current incumbent robust1_seed (= refit9_opt + robust feature
% reader) at 8 mM; fit A at 2 mM = the k2 power gate, R2 x0.54. The parameter
% construction is copied from RunBestFitFigures.m so both legs are the same model.
%
%   matlab -batch "run('Analyses/PooledFitFigures/RunPooledFigureSims.m')"
%
% Output -> results/sims.mat

rootP = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(rootP); addpath(genpath(rootP));
outP = fullfile(rootP, 'Analyses', 'PooledFitFigures', 'results');
if ~isfolder(outP); mkdir(outP); end
if ~exist('CORE', 'var'); CORE = 'robust1_seed'; end
if isempty(gcp('nocreate')); parpool('local', 5); end

vFine = 0:0.25:7;
LoadData;                                               % Data_ATP, Data_ATP_AB2021

pP = getParams(loadParams(CORE), [], true, false);
pP.mods = {}; pP.g = ones(1,0);
pP.EvalFeatures = true; pP.BreakOnODEUnstable = false; pP.MaxRunTime = 900;
pP.PlotEachSeparately = 0; pP.PlotFeatureFitting = 0; pP.PlotFullscreen = 0;
pP.ghostSave = ''; pP.ghostLoad = ''; pP.SaveBest = 0;
pP.UseAtpPowerGates = true; pP.MgATP_ref = 8; pP.UseAtpDetachR2 = 0; pP.UseAtpKmsrd = 0;

%% ---- 1. fine F-V -------------------------------------------------------------
q = pP; q.FV_velocities = -vFine; q.OptimizeFVInit = true;
q.RunSlack = false; q.RunStairs = false; q.RunKtr = false; q.RunForceVelocity = true;
q8 = q; q8.MgATP = 8;
q2 = q; q2.MgATP = 2; q2.pAtp_k2 = log(0.54)/log(1/4);  % fit A: R2 x0.54 at 2 mM
t0 = tic; [~, o8] = runFVExperiment(q8, [8 4 2], Data_ATP); F8 = arrayfun(@(o) o.Force(end), o8);
fprintf('FV 8 mM: %.0f s, F_iso %.2f kPa\n', toc(t0), F8(1));
% runFVExperiment takes [ATP] from the FIRST entry of its 2nd argument, not from MgATP
t0 = tic; [~, o2] = runFVExperiment(q2, [2 4 8], Data_ATP); F2 = arrayfun(@(o) o.Force(end), o2);
fprintf('FV 2 mM: %.0f s, F_iso %.2f kPa (2/8 = %.3f)\n', toc(t0), F2(1), F2(1)/F8(1));

%% ---- 2. passive slack (scored on 03/27, the only PNB+Mava recording) ---------
pp = pP; pp.MgATP = 8;
[~, ~, op] = runPassiveExperiment(pp);
PAS = struct('t', op.t(:), 'F', op.Force(:), 'SL', op.SL(:), 'datatable', op.datatable);

%% ---- 3. staircase PREDICTION, each preparation's own protocol -------------------
% Never scored. 03/27 is an 8-step (0.1 s) staircase, 04/03 and 04/10 are 4-step
% (0.2 s), so the three are not poolable as one trace. No 2 mM staircase exists:
% the 2 mM curve is a prediction at an [ATP] never measured.
STD = {'03_27', '04_03', '04_10'}; STR = struct('day', STD, 't', [], 'F8', [], 'F2', [], 't0', [], 'dt', []);
for i = 1:3
    qs = pP; qs.stairs_velocitytableonfile = sprintf('protocol_%s_2026_velocitytable_stairs.mat', STD{i});
    S = load(fullfile('data', qs.stairs_velocitytableonfile)); vt = S.velocitytable;
    STR(i).t0 = vt(find(vt(:,2) ~= 0, 1), 1); STR(i).dt = S.datatable;
    s8 = qs; s8.MgATP = 8; [~, o] = runStairsExperiment(s8);
    STR(i).t = o.t(:); STR(i).F8 = o.Force(:) / interp1(o.t, o.Force, STR(i).t0, 'linear', 'extrap');
    s2 = qs; s2.MgATP = 2; s2.pAtp_k2 = log(0.54)/log(1/4); [~, o] = runStairsExperiment(s2);
    STR(i).t2 = o.t(:); STR(i).F2 = o.Force(:) / interp1(o.t, o.Force, STR(i).t0, 'linear', 'extrap');
    fprintf('stairs %s done\n', STD{i});
end

save(fullfile(outP, 'sims.mat'), 'vFine', 'F8', 'F2', 'PAS', 'STR', 'CORE');
fprintf('wrote %s\n', fullfile(outP, 'sims.mat'));
