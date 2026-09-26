% SaveSeedBreakdown.m -- the seed's full per-feature cost, once, so every
% check-in can report "seed -> now" instead of "first improvement -> latest".
%
% optimizeFeatures only records featCost for IMPROVED iterates, so the seed
% itself never appears in state.optIter -- and a checker that diffs the first
% improved snapshot against the latest reports a flat zero after round 1,
% which is exactly the false "*** not moving ***" alarm of 2026-09-23 13:13.
%
% Runs 'AllPar' on a small pool: the per-segment initial conditions are
% precalculated, so the result does not depend on pool size and reproduces the
% campaign's own seed number.
%
%   matlab -batch "SEEDSNAP='refit6_seed'; run('Analyses/KtrSegmentSlope/SaveSeedBreakdown.m')"
if ~exist('SEEDSNAP','var'); SEEDSNAP = 'refit6_seed'; end
if ~exist('TAG','var');      TAG      = 'refit6';      end
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_ssb_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 2); end

p = getParams(loadParams(SEEDSNAP), [], true, false);
p.RunForceVelocity=true; p.RunSlack=true; p.RunStairs=false; p.RunKtr=true;
p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
p.RunSlackSegments='AllPar'; p.MaxRunTime=300; p.OptimizeOn='Feats';
p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0; p.mods={}; p.g=ones(1,0);

[cost, ~, featCost, featNames] = evaluateBakersExp(1, p, true, false);
seed = struct('snap', SEEDSNAP, 'cost', cost, 'featCost', featCost(:)', ...
              'featNames', {featNames(:)'}, 'when', datestr(now,'yyyy-mm-dd HH:MM'));
out = fullfile(root,'Analyses','KtrSegmentSlope','results',[TAG '_seed_breakdown.mat']);
save(out, 'seed');
fprintf('%s = %.4f  -> %s\n', SEEDSNAP, cost, out);
fprintf('DONE\n');
