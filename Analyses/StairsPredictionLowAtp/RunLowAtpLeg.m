% RunLowAtpLeg.m — does the estiff flip survive low ATP?
%
% The staircase work (main tree, Analyses/StairsPrediction) found that freezing
% estiff at 1.25 and refitting everything else gives a BETTER scored objective
% and a much better held-out staircase than the incumbent at estiff = 0.965.
% But that whole campaign was scored at MgATP = 8 only: the main tree has no ATP
% gate code at all, so nothing in it could have tested the 2 mM leg.
%
% This script does. The ATP gates here are pure power laws
%       rate_eff = rate * (MgATP/MgATP_ref)^p ,   MgATP_ref = 8
% so they are EXACTLY 1 at 8 mM for any p. That has two consequences that make
% this a clean experiment:
%
%   1. Switching the gates on cannot change either snapshot's 8 mM behaviour,
%      so the 8 mM fit each snapshot was bought with is untouched.
%   2. Transplanting ONE gate vector onto BOTH snapshots therefore compares
%      nothing except how each parameter set RESPONDS to lowering ATP.
%
% Gate vector: the H4 hypothesis (params/atphyp_H4_opt.m), the best physiological
% low-ATP account in this worktree — pAtp_k2 0.44 (slower strong-bridge
% detachment), pAtp_kah 0.30 (slower re-priming), plus the small Maxwell/series
% terms. 2 mM additionally carries Pi = 0.8 mM, as in every other head-to-head
% here (iCfg in Analyses/AtpEffectHypotheses/RunAtpHeadToHead.m).
%
% BRANCH-DRIFT GUARD: this worktree is a different branch from the one the
% snapshots were fitted on. The script first re-scores each snapshot at 8 mM and
% reports it against the value recorded in the main tree. If those disagree the
% ratios below are contaminated by model drift, not by estiff, and the run says so.
%
% Run: matlab -batch "run('Analyses/StairsPredictionLowAtp/RunLowAtpLeg.m')"

clc;   % NB: no `clear` — an externally supplied SNAPS must survive
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if ~exist('NWORKERS','var'); NWORKERS = 5; end
if isempty(gcp('nocreate')); parpool('local', NWORKERS); end
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

% name, estiff, Et recorded in the main tree. Override by defining SNAPS first.
if ~exist('SNAPS','var')
    SNAPS = { 'spreadw_opt',   0.965322, 0.9786
              'estiff125_opt', 1.25,     0.9559 };
end

% H4 power-gate vector, transplanted verbatim onto both snapshots.
GATE = struct('pAtp_k2',0.44, 'pAtp_kah',0.30, 'pAtp_eta_M',-0.15, ...
              'pAtp_kSE',0.07, 'pAtp_kSE_M',-0.15);

DS = { '03/27 M','protocol_03_27_2026_%s_slack.mat'
       '04/03 F','protocol_04_03_2026_%s_slack.mat'
       '04/10 M','protocol_04_10_2026_%s_slack.mat' };

% Force-class features are normalised by each preparation's own steady force
% before ratioing, so the 2/8 ratio is an ATP effect and not a prep scale.
FEAT = { 'steady','F'; 'A','F'; 'peak1_y','F'; 'vall_y','F'; 'peak2','F'
         'ktr','R';    'rsK','R' };

%% ---------- 1. the measured 2 mM / 8 mM ratio, three preparations ---------
DR = struct();
for a = [8 2]
    tag = sprintf('%dmM', a);
    for i = 1:size(DS,1)
        S  = load(['data/' sprintf(DS{i,2}, tag)]);
        ft = extractSlackAttributes(S.datatable(:,1), S.datatable(:,3), ...
                                    S.datatable(:,2), S.velocitytable, [], [], false, true);
        sc = mean(double(ft.steady),'omitnan');
        for j = 1:size(FEAT,1)
            nm = FEAT{j,1};
            if ~isfield(ft,nm); continue; end
            v = mean(double(ft.(nm)),'omitnan');
            if FEAT{j,2}=='F' && ~strcmp(nm,'steady'); v = v/sc; end
            DR.(sprintf('a%d',a)).(nm)(i) = v; %#ok<AGROW>
        end
    end
end
fprintf('\n=== MEASURED 2 mM / 8 mM (raw, uncorrected for rundown) ===\n');
dataRatio = struct();
for j = 1:size(FEAT,1)
    nm = FEAT{j,1};
    if ~isfield(DR.a8,nm) || ~isfield(DR.a2,nm); continue; end
    r = DR.a2.(nm)(:)./DR.a8.(nm)(:);
    dataRatio.(nm) = r;
    fprintf('  %-14s %.3f %.3f %.3f   mean %.3f  CV %.0f %%\n', nm, r, mean(r), 100*std(r)/mean(r));
end
fprintf('  [rundown-corrected targets, ATPEffectReconciliation: force x1.36, ktr x0.547]\n');

%% ---------- 2. each snapshot, both conditions -----------------------------
R = struct();
for k = 1:size(SNAPS,1)
    nm = SNAPS{k,1};
    p0 = getParams(loadParams(nm), [], true, false);
    for f = fieldnames(GATE)'; p0.(f{1}) = GATE.(f{1}); end
    p0.UseAtpPowerGates = true;

    fprintf('\n===== %s (estiff %.4f) =====\n', nm, p0.estiff);
    assert(abs(p0.estiff - SNAPS{k,2}) < 1e-3, 'estiff mismatch in %s', nm);

    % --- branch-drift guard: does this tree reproduce the main tree's 8 mM Et?
    % It must run under the configuration the snapshot was FITTED with. The
    % first version of this guard used iCfg, which switches FV and the passive
    % protocol OFF while the fn list still demands FV_* and PS_* features; the
    % missing terms are penalised and BOTH snapshots scored Et ~12531, a
    % spurious "drift" that was purely an artefact of the reduced experiment
    % set. With the full set it reproduces to 0.1 % (see RunDriftGuard.m).
    pchk = p0;
    pchk.RunForceVelocity = true;  pchk.RunKtr = false;
    pchk.RunSlack = true;          pchk.RunStairs = false;
    pchk.RunForceVelocityTime = false; pchk.RunSlackPassive = true;
    pchk.RunSlackSegments = 'AllPar';
    pchk.EvalFeatures = true;      pchk.BreakOnODEUnstable = false;
    pchk.OptimizeOn = 'Feats';     pchk.MaxRunTime = 600;
    pchk.PlotEachSeparately = 0;   pchk.PlotFeatureFitting = 0; pchk.PlotFullscreen = 0;
    pchk.MgATP = 8; pchk.MgADP = 0; pchk.Pi = 0;
    pchk.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
    pchk.mods = {}; pchk.g = ones(1,0);
    Et8 = evaluateBakersExp(1, pchk, true, false);
    fprintf('  8 mM Et here %.4f   (main tree recorded %.4f)  drift %+.1f %%\n', ...
            Et8, SNAPS{k,3}, 100*(Et8/SNAPS{k,3}-1));

    for a = [8 2]
        ps = iCfg(p0, a);
        [~,~,fm,~] = runSlackExperiment(ps);
        R.(nm).(sprintf('a%d',a)) = fm;
    end

    f8 = R.(nm).a8; f2 = R.(nm).a2;
    sc8 = mean(double(f8.steady),'omitnan'); sc2 = mean(double(f2.steady),'omitnan');
    fprintf('  %-14s %8s %8s %8s %10s\n','feature','8 mM','2 mM','ratio','data mean');
    for j = 1:size(FEAT,1)
        fn = FEAT{j,1};
        if ~isfield(f8,fn) || ~isfield(f2,fn); continue; end
        v8 = mean(double(f8.(fn)),'omitnan'); v2 = mean(double(f2.(fn)),'omitnan');
        if FEAT{j,2}=='F' && ~strcmp(fn,'steady'); v8 = v8/sc8; v2 = v2/sc2; end
        dm = NaN; if isfield(dataRatio,fn); dm = mean(dataRatio.(fn)); end
        fprintf('  %-14s %8.4g %8.4g %8.3f %10.3f\n', fn, v8, v2, v2/v8, dm);
        R.(nm).ratio.(fn) = v2/v8;
    end
    R.(nm).Et8 = Et8; R.(nm).estiff = p0.estiff;
end

%% ---------- 3. side by side ----------------------------------------------
fprintf('\n=== ATP RESPONSE: every snapshot, one gate vector ===\n');
fprintf('%-14s', 'feature');
for k = 1:size(SNAPS,1); fprintf('%12s', sprintf('e=%.3f', SNAPS{k,2})); end
fprintf('%12s\n', 'data');
for j = 1:size(FEAT,1)
    fn = FEAT{j,1};
    have = all(cellfun(@(n) isfield(R,n) && isfield(R.(n).ratio,fn), SNAPS(:,1)'));
    if ~have; continue; end
    fprintf('%-14s', fn);
    for k = 1:size(SNAPS,1); fprintf('%12.3f', R.(SNAPS{k,1}).ratio.(fn)); end
    dm = NaN; if isfield(dataRatio,fn); dm = mean(dataRatio.(fn)); end
    fprintf('%12.3f\n', dm);
end

save(fullfile(resDir,'lowatp_leg.mat'), 'R','dataRatio','GATE','SNAPS');

function p = iCfg(p, atp)
p.UseAtpJointFit=false; p.EvalFeatures=true; p.BreakOnODEUnstable=false;
p.OptimizeOn='Feats'; p.MaxRunTime=600;
p.PlotEachSeparately=0; p.PlotFeatureFitting=0; p.PlotFullscreen=0;
p.RunKtr=false; p.RunStairs=false; p.RunForceVelocityTime=false;
p.RunSlackPassive=false; p.RunForceVelocity=false; p.RunSlack=true;
p.RunSlackSegments='AllPar';
p.MgATP=atp; p.MgADP=0;
if atp==8
    p.Pi=0; p.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
else
    p.Pi=0.8; p.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
end
end
