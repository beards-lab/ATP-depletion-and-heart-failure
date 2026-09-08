% CollectCampaign.m — re-score every campaign snapshot into one table.
%
% The campaign runs as several concurrent MATLAB processes, so no single
% process sees all the results and a shared .mat would be clobbered. This
% script is the authoritative collector: it re-scores params/ec<NNN>_seed.m and
% params/ec<NNN>_opt.m from disk, on the objective and on the held-out
% staircase, and writes results/campaign_table.mat.
%
% Re-scoring rather than trusting the shard logs also recovers exponents whose
% shard died: e = 1.40 tripped optimizeFeatures' final self-consistency guard
% (best_params re-evaluated to 1.3980 against a recorded 1.4321 — a 2.4 % drift
% against a 2 % tolerance, and in the BETTER direction), but its improving round
% had already written the snapshot, so the point is recoverable here.
%
% Run: matlab -batch "run('Analyses/StairsPrediction/CollectCampaign.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
if ~exist('NWORKERS','var'); NWORKERS = 5; end
if isempty(gcp('nocreate')); parpool('local', NWORKERS); end

% The optimiser rewrites these snapshots from another process and MATLAB caches
% parsed script files, so force a re-read before loading any of them.
rehash;

d  = load('data/protocol_03_27_2026_velocitytable_stairs.mat');
fD = extractStairsAttributes(d.datatable(:,1), d.datatable(:,3), ...
                             d.velocitytable, 'st_', struct(), true);

snaps = dir(fullfile(pwd,'params','ec*_opt.m'));
r = struct('tag',{},'estiff',{},'seedEt',{},'optEt',{}, ...
           'seedStair',{},'optStair',{},'kstiff2',{},'kSE',{});

for i = 1:numel(snaps)
    tag = erase(snaps(i).name, '_opt.m');
    if ~isfile(fullfile(pwd,'params',[tag '_seed.m']))
        fprintf('  %-6s no seed snapshot, skipping\n', tag); continue;
    end
    pO = getParams(loadParams([tag '_opt']),  [], true, false);
    pS = getParams(loadParams([tag '_seed']), [], true, false);

    r(end+1) = struct('tag', tag, 'estiff', pO.estiff, ...
        'seedEt',    evaluateBakersExp(pS.g, pS, true, false), ...
        'optEt',     evaluateBakersExp(pO.g, pO, true, false), ...
        'seedStair', iStair(pS, d, fD), ...
        'optStair',  iStair(pO, d, fD), ...
        'kstiff2',   pO.kstiff2, 'kSE', pO.kSE);                  %#ok<SAGROW>

    fprintf('%-6s e=%.3f | Et %8.4f -> %8.4f | stair %8.3f -> %8.3f | kstiff2 %.5g  kSE %.5g\n', ...
        tag, r(end).estiff, r(end).seedEt, r(end).optEt, ...
        r(end).seedStair, r(end).optStair, r(end).kstiff2, r(end).kSE);
end

[~,ix] = sort([r.estiff]); r = r(ix);
save(fullfile(here,'results','campaign_table.mat'), 'r');
fprintf('\nwrote results/campaign_table.mat (%d exponents)\n', numel(r));

function Cf = iStair(p, d, fD)
%ISTAIR  Held-out staircase feature cost, normalised to the force at the
%        instant the first ramp starts — the instrument that scored the data.
ps = p; ps.RunStairs = 1; ps.PlotEachSeparately = false;
ps.ghostSave = ''; ps.ghostLoad = '';
[~, os] = runStairsExperiment(ps);
tMov = d.velocitytable(find(d.velocitytable(:,2)~=0,1),1);
Fm   = interp1(os.t, os.Force, tMov, 'linear', 'extrap');
Cf   = stairsCost(os.t, os.Force/Fm, d.velocitytable, fD);
end
