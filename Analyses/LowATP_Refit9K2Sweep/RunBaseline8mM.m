% RunBaseline8mM.m -- the incumbent's 8 mM cost breakdown, recorded BEFORE and
% AFTER the floored-campaign merge, so the merge can be shown not to move it.
%
%   matlab -batch "LABEL='premerge'; run('Analyses/LowATP_Refit9K2Sweep/RunBaseline8mM.m')"
if ~exist('LABEL','var'); LABEL = 'premerge'; end
if ~exist('SNAPS','var'); SNAPS = {'refit9_opt','robust1_seed'}; end
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','LowATP_Refit9K2Sweep','results');
diary(fullfile(resdir, ['baseline8_' LABEL '.log'])); diary on;

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_b8_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 5); end

B = struct();
for i = 1:numel(SNAPS)
    t0 = tic;
    p = getParams(loadParams(SNAPS{i}), [], true, false);
    p.RunForceVelocity=true; p.RunSlack=true; p.RunStairs=false; p.RunKtr=true;
    p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
    p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
    p.RunSlackSegments='AllPar'; p.MaxRunTime=300; p.OptimizeOn='Feats';
    p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0; p.mods={}; p.g=ones(1,0);
    [cost, ~, featCost, featNames] = evaluateBakersExp(1, p, true, false);
    B(i).snap = SNAPS{i}; B(i).cost = cost; B(i).featCost = featCost(:)';
    B(i).featNames = featNames(:)'; B(i).sec = toc(t0);
    fprintf('%s [%s] = %.6f  (%.0f s)\n', SNAPS{i}, LABEL, cost, B(i).sec);
    save(fullfile(resdir, ['baseline8_' LABEL '.mat']), 'B');
end
fprintf('DONE\n'); diary off;
