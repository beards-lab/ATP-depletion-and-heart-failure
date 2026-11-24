% AuditBaselines8mM.m -- which 8 mM snapshot is actually the best?
%
% The pathway ladder was seeded from joint82_opt, whose 8 mM leg is 4.0094 --
% WORSE than realign4_opt's 3.8513, because the joint refit traded 8 mM away for
% 2 mM. That makes the ladder's absolute numbers wrong (its ranking survives,
% since every rung shared the same core, but the core was the wrong one).
%
% Scores every candidate on ONE objective, in ONE session, with the current code
% (i.e. including the force-velocity data-source fix), so the answer is
% like-for-like. The winner becomes the frozen core for the hypothesis tests.
%
% Run: matlab -batch "run('Analyses/LowATP_Realign4Transfer/AuditBaselines8mM.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if isempty(gcp('nocreate')); parpool('local',5); end

CAND = {'realign4_opt','realign3_opt','realign2_opt','realign1_opt', ...
        'joint82_opt','rskR2_w025_opt','lowATP_H2_chem_opt'};

% one objective: realign4's feature list, full 8 mM battery
ref = getParams(loadParams('realign4_opt'), [], true, false);
FN  = ref.fn;

fprintf('\n%-24s %10s %10s %10s %10s %8s\n', ...
        'snapshot','8 mM Et','FV power','rsK term','slack rest','rsK s^-1');
R = struct('name',{},'Et',{},'fv',{},'rsk',{},'rskAbs',{});
for i = 1:numel(CAND)
    try
        p = getParams(loadParams(CAND{i}), [], true, false);
        p.fn = FN; p.mods = {}; p.g = [];
        p.UseAtpJointFit = false;
        p.MgATP = 8; p.MgADP = 0; p.Pi = 0;
        p.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
        p.RunForceVelocity = true; p.RunSlack = true; p.RunSlackPassive = true;
        p.RunKtr = false; p.RunStairs = false; p.RunForceVelocityTime = false;
        p.EvalFeatures = true; p.BreakOnODEUnstable = false; p.OptimizeOn = 'Feats';
        p.RunSlackSegments = 'AllPar'; p.MaxRunTime = 400;
        p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0;

        [Et,~,fc,nm] = evaluateBakersExp(1, p, true, false);
        ifv = find(contains(nm,'FV_normpowerAvg'),1);
        irk = find(strncmp(nm,'rsK',3),1);
        fv  = 0; if ~isempty(ifv); fv = fc(ifv); end
        rk  = 0; if ~isempty(irk); rk = fc(irk); end

        % absolute rsK, for the "was rsK tackled" question
        q = p; q.RunForceVelocity=false; q.RunSlackPassive=false;
        q.RunSlackSegments='All';
        [~,~,fm,~] = runSlackExperiment(q);
        rka = mean(fm.rsK(:),'omitnan');

        R(end+1) = struct('name',CAND{i},'Et',Et,'fv',fv,'rsk',rk,'rskAbs',rka); %#ok<SAGROW>
        fprintf('%-24s %10.4f %10.4f %10.4f %10.4f %8.1f\n', ...
                CAND{i}, Et, fv, rk, Et-fv-rk, rka);
    catch ME
        fprintf('%-24s   FAILED: %s\n', CAND{i}, ME.message);
    end
end

[~,ib] = min([R.Et]);
fprintf('\nBEST 8 mM: %s  (Et = %.4f)\n', R(ib).name, R(ib).Et);
fprintf('data rsK at 8 mM is ~43.7 s^-1 for reference.\n');

save(fullfile(here,'results','baseline_audit.mat'), 'R', 'FN');
fprintf('\nAUDIT DONE\n');
