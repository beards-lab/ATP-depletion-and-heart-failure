% Export2mMSlidePanels.m — 2 mM counterpart of Export8mMPanels + makeSlidePanels.
% Runs the low-ATP slack experiment once, caches it, and writes the two
% slide-legible panels into Analyses/LowATP_AbsoluteFit/results/panels_2mM/.

root = 'C:\home\git\ATP-wt-lowatp';
cd(root); addpath(genpath(root)); LoadData;
outdir = fullfile(root,'Analyses','LowATP_AbsoluteFit','results','panels_2mM');
if ~exist(outdir,'dir'); mkdir(outdir); end

BASE = 'lowATP_powergate2_opt';
p0 = getParams(loadParams(BASE), [], true, false);
p0.MgATP = 2; p0.MgADP = 0; p0.Pi = 0.8;
p0.EvalFeatures = true; p0.BreakOnODEUnstable = false; p0.OptimizeOn = 'Feats';
p0.RunSlackSegments = 'All'; p0.MaxRunTime = 400;
p0.PlotEachSeparately = 0; p0.PlotFeatureFitting = 0; p0.PlotFullscreen = 0;
p0.RunForceVelocity = false; p0.RunForceVelocityTime = false;
p0.RunStairs = false; p0.RunKtr = false; p0.RunSlackPassive = false;
p0.RunSlack = true;
p0.velocitytableonfile       = 'protocol_03_27_2026_2mM_slack.mat';
p0.ktr_velocitytableonfile   = 'protocol_03_27_2026_velocitytable_ktr_2mM.mat';
p0.stairs_velocitytableonfile= 'protocol_03_27_2026_velocitytable_stairs.mat';

fprintf('\n===== slack experiment, 2 mM =====\n');
[E_slack, out_slack, fm, fd] = runSlackExperiment(p0);
fprintf('E_slack = %s\n', mat2str(round(E_slack,4)));
save(fullfile(outdir,'slack2_cache.mat'), 'out_slack','fm','fd','E_slack','-v7.3');

close all force;
makeSlidePanels(fullfile(outdir,'slack2_cache.mat'), outdir, '2mM', '2 mM ATP');
disp('EXPORT 2mM DONE');
