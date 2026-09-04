% RunBothAtpSlideRuns.m -- simulate the 2 mM slack and staircase protocols for
% EACH rung of the ATP hypothesis ladder, and export the traces for the slide
% figures that live in the OTHER worktree.
%
% WHY IT LIVES HERE. The ATP power-gate machinery (params0.pAtp_*, MgATP_ref) is
% implemented in THIS worktree's Model/getParams.m and Model/evaluateModel.m
% only. The Further-optim worktree cannot evaluate a gated parametrisation at
% all -- loadParams would set the pAtp_* fields and nothing would read them, so
% the "2 mM" run would silently come back as an 8 mM run. Hence: simulate here,
% draw there, and pass a plain .mat per rung between the two.
%
% THE RUNGS (RunAtpHypothesisOptim.m; the gate SET is what differs, and most
% gates in the wider sets fit to exactly zero):
%   H1   3 gates   chemistry only          k2 kah kamh
%   H2   7 gates   + SRX nucleotide states + kmsrd ksr0 ksr2srd ksrd2sr
%   H3  10 gates   + mechanical            + eta_M kSE mu_neg
%   H4  21 gates   any rate                + kstiff2 kSE_M dr dr2 P_bound_max
%                                            ka kd k1 k_1 tau_reg A0
% Every gate is exactly 1 at MgATP_ref = 8, so ALL FOUR rungs are the same
% 8 mM model -- the core, spreadw_iter_018. Only the 2 mM leg differs.
%
% THE STAIRCASE IS UNMATCHED. All three protocol days recorded the staircase at
% 8 mM only (DataCuration/CreateStairsProtocolVelocityTable.m drives all three
% from *_Merged_8mM_Active.txt). Every 2 mM staircase below is therefore a
% prediction with NO data to test it against -- the figures must say so.
%
% Run:  matlab -batch "run('Analyses/AtpEffectHypotheses/RunBothAtpSlideRuns.m')"
%       matlab -batch "RUNGS={'H2'}; run('Analyses/.../RunBothAtpSlideRuns.m')"

clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

if ~exist('RUNGS','var'); RUNGS = {'H1','H2','H3','H4'}; end
if ~exist('NWORK','var'); NWORK = 4; end
CORE = 'spreadw_iter_018';

% Where the drawing worktree expects to find these.
FIGROOT = fullfile(root,'..','ATP-depletion-and-heart-failure','Analyses', ...
                   'OptimResultReport','results','bothatp');

LoadData;
if isempty(gcp('nocreate')); parpool('local', NWORK); end

for r = 1:numel(RUNGS)
    rung = RUNGS{r};
    SNAP = sprintf('atphyp_%s_opt', rung);
    fprintf('\n================ %s (%s) ================\n', rung, SNAP);
    if ~isfile(fullfile('params',[SNAP '.m']))
        fprintf(2,'  MISSING params/%s.m -- skipped\n', SNAP); continue;
    end

    p = getParams(loadParams(SNAP), [], true, false);
    p.mods = {}; p.g = ones(1,0);
    p.EvalFeatures = true; p.OptimizeOn = 'Feats';
    p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0; p.PlotFullscreen = 0;
    p.BreakOnODEUnstable = false; p.MaxRunTime = 900;
    p.UseAtpJointFit = false; p.UseAtpPowerGates = true; p.MgATP_ref = 8;

    % 2 mM chemistry, exactly as the ladder scored it (iCfg in
    % RunLowAtpFitReport.m and RunAtpHypothesisOptim.m): the Pi = 0.8 mM is the
    % measured accompaniment of the low-ATP solution, not a free parameter.
    ATPMM = 2;
    p.MgATP = ATPMM; p.MgADP = 0; p.Pi = 0.8;
    p.velocitytableonfile        = 'protocol_03_27_2026_2mM_slack.mat';
    p.stairs_velocitytableonfile = 'protocol_03_27_2026_velocitytable_stairs.mat';

    nz = fieldnames(p); nz = nz(startsWith(nz,'pAtp_'));
    nz = nz(cellfun(@(f) abs(p.(f)) > 1e-9, nz));
    fprintf('  active gates: %s\n', strjoin(strrep(nz,'pAtp_',''), ', '));

    % ---- slack, all five cycles ------------------------------------------
    ps = p; ps.RunSlack = true; ps.RunSlackSegments = 'All';
    ps.RunStairs = false; ps.RunKtr = false; ps.RunForceVelocity = false;
    ps.RunForceVelocityTime = false; ps.RunSlackPassive = false;
    t0 = tic; [Eslack, out_slack, fmodel, fdata] = runSlackExperiment(ps);
    fprintf('  slack     %.0f s\n', toc(t0));

    % ---- staircase: a PREDICTION, there is no 2 mM staircase data ---------
    out_stairs = []; E_stairs = NaN;
    pst = p; pst.RunSlack = false; pst.RunStairs = true;
    t0 = tic;
    try
        [E_stairs, out_stairs] = runStairsExperiment(pst);
        fprintf('  staircase %.0f s (cost vs 8 mM data %.4f -- NOT a fit)\n', toc(t0), E_stairs);
    catch ME
        fprintf(2,'  stairs FAILED: %s\n', ME.message);
    end

    % ---- export ----------------------------------------------------------
    % Trimmed to what the figures need: the two traces, the model's own
    % features, and enough provenance that a stale export cannot be mistaken
    % for a fresh one.
    S = struct();
    S.snap = SNAP; S.rung = rung; S.core = CORE; S.atp_mM = ATPMM;
    S.Pi_mM = p.Pi; S.MgADP_mM = p.MgADP;
    S.gates = struct(); for i = 1:numel(nz); S.gates.(nz{i}) = p.(nz{i}); end
    S.slack  = struct('t',out_slack.t,'Force',out_slack.Force,'SL',out_slack.SL, ...
                      'datatable',out_slack.datatable);
    S.velocitytableonfile        = p.velocitytableonfile;
    S.stairs_velocitytableonfile = p.stairs_velocitytableonfile;
    if ~isempty(out_stairs)
        S.stairs = struct('t',out_stairs.t,'Force',out_stairs.Force,'SL',out_stairs.SL);
    end
    S.fmodel = fmodel; S.fdata = fdata; S.Eslack = Eslack; S.E_stairs = E_stairs;
    S.built  = datestr(now,'yyyy-mm-dd HH:MM:SS');                    %#ok<TNOW1,DATST>

    OUTS = { fullfile(resDir, sprintf('bothatp_lowatp_runs_%s.mat', rung))
             fullfile(FIGROOT, rung, 'lowatp_runs.mat') };
    for i = 1:numel(OUTS)
        d = fileparts(OUTS{i});
        if ~isfolder(d)
            if contains(d, 'ATP-depletion-and-heart-failure') && isfolder(fileparts(d))
                mkdir(d);                       % per-rung subfolder in the fig tree
            else
                fprintf(2,'  skipped (no folder): %s\n', OUTS{i}); continue;
            end
        end
        save(OUTS{i},'-struct','S'); fprintf('  wrote %s\n', OUTS{i});
    end
end
fprintf('\nDONE\n');
