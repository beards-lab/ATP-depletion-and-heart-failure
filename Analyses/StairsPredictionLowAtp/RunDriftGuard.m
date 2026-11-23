% RunDriftGuard.m — does THIS worktree reproduce the main tree's 8 mM model?
%
% RunLowAtpLeg's first attempt at this check ran evaluateBakersExp under iCfg,
% which switches OFF force-velocity and the passive protocol while the fn list
% still demands FV_* and PS_* features. Missing features are penalised, so both
% snapshots scored Et ~12531 and the guard reported "drift" that was entirely an
% artefact of the reduced experiment set.
%
% This is the check done properly: the SAME experiment configuration the
% snapshots were fitted under (FV on, passive on, slack on, stairs off), so the
% number is comparable with the main tree's.
%
%   main tree:  spreadw_opt   Et = 0.9786
%               estiff125_opt Et = 0.9559
%
% It also dumps the 8 mM slack FEATURE values, which is the check that actually
% matters for the ATP ratios: those ratios come from runSlackExperiment, not
% from the objective, so what has to match between trees is the features.
%
% Run: matlab -batch "run('Analyses/StairsPredictionLowAtp/RunDriftGuard.m')"

clear; clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if ~exist('NWORKERS','var'); NWORKERS = 3; end
if isempty(gcp('nocreate')); parpool('local', NWORKERS); end

SNAPS = { 'spreadw_opt', 0.9786; 'estiff125_opt', 0.9559 };
FEAT  = {'steady','A','peak1_y','vall_y','peak2','vall2_dy','ktr','rsK','rsT0','peak1_dSL','vall_t'};

out = struct();
for k = 1:size(SNAPS,1)
    nm = SNAPS{k,1};
    p  = getParams(loadParams(nm), [], true, false);

    % the configuration the snapshot was FITTED under (cf. RunEstiffCampaign)
    p.RunForceVelocity = true;  p.RunKtr = false;
    p.RunSlack = true;          p.RunStairs = false;
    p.RunForceVelocityTime = false; p.RunSlackPassive = true;
    p.RunSlackSegments = 'AllPar';
    p.EvalFeatures = true;      p.BreakOnODEUnstable = false;
    p.OptimizeOn = 'Feats';     p.MaxRunTime = 600;
    p.PlotEachSeparately = 0;   p.PlotFeatureFitting = 0; p.PlotFullscreen = 0;
    p.MgATP = 8; p.MgADP = 0; p.Pi = 0;
    p.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
    p.mods = {}; p.g = ones(1,0);

    [Et, ~, ct, cn] = evaluateBakersExp(1, p, true, false);
    fprintf('\n=== %s (estiff %.4f) ===\n', nm, p.estiff);
    fprintf('  Et here %.4f   main tree %.4f   drift %+.1f %%\n', Et, SNAPS{k,2}, 100*(Et/SNAPS{k,2}-1));
    [~,ord] = sort(ct,'descend');
    for i = ord(1:min(6,numel(ord)))
        fprintf('     %-24s %8.4f\n', cn{i}, ct(i));
    end

    ps = p; ps.RunForceVelocity = false; ps.RunSlackPassive = false;
    [~,~,fm,~] = runSlackExperiment(ps);
    fprintf('  8 mM slack features:\n');
    for j = 1:numel(FEAT)
        if isfield(fm,FEAT{j})
            v = mean(double(fm.(FEAT{j})),'omitnan');
            fprintf('     %-16s %10.5g\n', FEAT{j}, v);
            out.(nm).(FEAT{j}) = v;
        end
    end
    out.(nm).Et = Et;
end
save(fullfile(here,'results','drift_guard.mat'), 'out', 'SNAPS');
