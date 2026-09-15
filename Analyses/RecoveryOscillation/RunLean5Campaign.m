% RunLean5Campaign.m -- the 15-hour refit on the lean objective, with slack
% cycle 5 scored again.
%
% Engine: ../FeatureWeighting/RunOptimSpreadW.m. This file sets the seed, the
% budget, the reduced draw pool, and ONE of four search PROFILES, so several
% profiles can run side by side as independent MATLAB sessions.
%
% =========================================================================
% WHAT IS BEING FITTED, AND WHAT CHANGED
% =========================================================================
% The 3-preparation pooled 8 mM target (Model/pooledTarget8mM), regenerated so
% that EVERY feature's spread -- rsK80 included -- is measured over all five
% slack cycles. Seed built by BuildLean5Seed.m, which documents the four
% changes from `floored_seed`. In brief:
%
%   1. rsK80[1:4] -> rsK80. Cycle 5 is the protocol's only re-stretch to
%      ~2.0 um SL rather than 2.2, at the smallest amplitude and slowest ramp,
%      so it is the only probe of how the re-stretch recovery depends on the
%      length it recovers at. It also carries a real ATP signal: measured
%      alone its between-prep spread is 19.5 % with z = 3.01, and its 2/8 mM
%      ratio (0.390) sits inside the range of cycles 1-4 (0.369-0.404), which
%      transducer noise cannot do. Cost: cycles 1-4 lose ~13 % of their weight
%      (s 15.8 % -> 16.9 %). Paid deliberately.
%   2. Three measured-free mechanisms deleted outright.
%   3. The four lattice OBSERVATIONS fixed (six lattice parameters carry two
%      identifiable combinations).
%   4. coolDownLS retired, rsK80 x2, F-V pin 0.40 -> 0.25.
%
% THE POOL is MechanismAblation sec. 7.5 applied: twelve parameters out, one
% in. Out as INERT (they have no effect on the trajectory and have been drawn
% and "optimised" for this project's whole campaign history): sigma1,
% sigma_srd2, LXBpivot. Out as DELETED with their mechanism: kA2re, kA2hop,
% sA2hop, x_M_slack. Out as REDUNDANT: d10_ref, SL_ref_lattice, R_thick,
% R_thin (four of the six lattice parameters, which enter only through two
% ratios) and L_hbare (a pure gain that ka already spans).
%
% In as the one place a redundant coordinate HELPS: `xrate`, which scales the
% whole core cycle at once (ka, kd, k1, k_1, k2, kah, kamh -- verified: x1.5 on
% the modifier gives exactly x1.5 on each, and leaves ksr0 alone). "All rates
% together" is the sloppiest direction in any rate model and without a common
% handle the simplex has to discover it by coordinating seven parameters.
% WRITE THIS DOWN WHEN REPORTING: xrate's VALUE is meaningless in isolation --
% it can always be traded back into the individual rates. It is a search
% coordinate, not a result.
%
% WHAT TO WATCH. Cycle 5 is back precisely so it can be fitted, and the
% mechanism that can fit it is the lattice pair (d_optimal, sigma_lattice):
% f_lattice is 1.00 at SL 2.2 and ~0.95 at SL 2.0, so cycle 5 is the only
% manoeuvre in the objective that the lattice attenuation touches at all. If
% the campaign improves rsK80 by moving something else entirely -- kstiff2, the
% strain-dependence shape -- that is worth knowing before the 15 h is spent.
%
% Watch also for the failure the last two campaigns had: a solved, smooth term
% absorbing the whole gain. F-V is pinned at 0.25 AND floored, so it should now
% be structurally unable to; if the consolidation report shows F-V as the
% biggest mover again, the floor is not doing its job and that is the finding.
%
% =========================================================================
% HOW TO RUN
% =========================================================================
% Pre-flight first, always. It costs one evaluation and it is the difference
% between finding a broken seed now and finding it in fifteen hours:
%
%   matlab -batch "PROFILE='a'; DRYRUN=true; run('Analyses/RecoveryOscillation/RunLean5Campaign.m')"
%
% One profile, 15 h, resumable (resume is automatic if the state file exists):
%
%   matlab -batch "PROFILE='a'; TIME_HRS=15; run('Analyses/RecoveryOscillation/RunLean5Campaign.m')"
%
% TWO profiles, not four. Analyses/RecoveryOscillation sec. 7 records three of
% four sessions dying with 24 MATLAB processes on one node; the engine
% saturates at ~5-6 workers per instance, so two instances x 5 workers = 10.
%
% Then consolidate:  RunRigAwareReport.m with CAMPAIGN='lean5'
%
% SIZING. ~38 s per evaluation at 5 workers; a round is SIMPLEX_EVALS of them.
% Profile a is ~40 min/round -> ~22 rounds in 15 h. Profile c is ~25 min/round
% -> ~35 rounds. The budget is checked BETWEEN rounds and never inside one, so
% the real stop time overshoots by up to one round.
%
% THE PROFILES -- each varies HOW the search moves, not what it scores.
%   a  the pooled8 setting. The control: the configuration that was still
%      descending monotonically when its budget ran out.
%   b  wider draw (16) and more polish (90).
%   c  narrow (6) and fast (40) with a hard kick (0.25). Many more rounds and
%      the strongest escape from a flat ridge -- which is how estiff125
%      stalled (Analyses/StairsPrediction: the incumbent was where the search
%      STOPPED, not a minimum).
%   d  profile a with the global surrogate on (120 evals). The only member that
%      can leave the basin -- and the one that wrote a 19.8 MB log of
%      integration-failure traces last time, because surrogateopt samples the
%      raw parameter box. Watch its log size.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('PROFILE','var');  PROFILE  = 'a'; end
if ~exist('TIME_HRS','var'); TIME_HRS = 15;  end
if ~exist('NWORKERS','var'); NWORKERS = 5;   end
if ~exist('CAMPAIGN','var'); CAMPAIGN = 'lean5'; end

SPEC = campaignSpec(CAMPAIGN);
if ~exist('SEEDSNAP','var');    SEEDSNAP    = SPEC.seedSnap;   end
if ~exist('SEED_EXPECT','var'); SEED_EXPECT = SPEC.seedExpect; end
if ~exist('TAGSTEM','var');     TAGSTEM     = SPEC.stem;       end

switch lower(PROFILE)
    case 'a'; N_DRAW=10; SIMPLEX_EVALS=60; KICK_FRAC=0.10; SURR_EVALS=0;
    case 'b'; N_DRAW=16; SIMPLEX_EVALS=90; KICK_FRAC=0.10; SURR_EVALS=0;
    case 'c'; N_DRAW= 6; SIMPLEX_EVALS=40; KICK_FRAC=0.25; SURR_EVALS=0;
    case 'd'; N_DRAW=10; SIMPLEX_EVALS=60; KICK_FRAC=0.10; SURR_EVALS=120;
    otherwise
        error('RunLean5Campaign:profile', ...
              'PROFILE must be a, b, c or d (got "%s").', PROFILE);
end
TAG          = sprintf('%s_%s', TAGSTEM, lower(PROFILE));
MAX_RUN_TIME = 300;

%% ---- the objective must match the target it is scored against ----------
% A chi-2 is only a chi-2 if the residual and the spread come from the same
% manoeuvre. This campaign exists to score cycle 5, so refuse to start against
% a target whose rsK80 spread was measured on four cycles.
PT = pooledTarget8mM();
kR = find(strcmp(PT.name,'rsK80'), 1);
if isempty(kR) || ~isequal(sort(PT.seg{kR}(:)'), 1:5)
    error('RunLean5Campaign:targetSegments', ...
        ['Model/pooledTarget8mM.m carries rsK80 over %s, not 1:5 -- this is the\n' ...
         '`floored`-era target. Regenerate and rebuild the seed:\n' ...
         '    clear OPTS; run(''Analyses/DataPooling/WritePooledTarget.m'')\n' ...
         '    run(''Analyses/RecoveryOscillation/BuildLean5Seed.m'')'], ...
        mat2str(PT.seg{kR}));
end

if ~isfile(SEEDSNAP)
    error('RunLean5Campaign:noSeed', ...
        ['%s not found. Build it first:\n' ...
         '    clear OPTS; run(''Analyses/DataPooling/WritePooledTarget.m'')\n' ...
         '    run(''Analyses/RecoveryOscillation/BuildLean5Seed.m'')'], SEEDSNAP);
end

%% ---- the lean pool (MechanismAblation sec. 7.5) ------------------------
POOL_OVERRIDE = { ...
    ... % core cycle, plus the common-mode handle that spans it (sec. 7.4)
    'xrate', ...
    'ka','kd','k1','k_1','k2','k_2','kah','kamh', ...
    ... % stiffness / mechanics   (x_M_slack deleted with UseMaxwellTensionOnly)
    'kstiff1','kstiff2','kstiff1_n','kstiff2_n','kSE','estiff','ekSE', ...
    'kSE_M','eta_M','mu_neg','mu', ...
    'dr','dr1','dr2', ...
    ... % SRX                     (sigma1, sigma_srd2 inert: their gates are 0)
    'ksr0','kmsrd','ksr2srd','ksrd2sr','sigma2','sigma_srd1', ...
    ... % attachment gating / FV shoulder
    'A0','v_ref_reg','tau_reg','P_bound_max', ...
    ... % the realignment mechanism + the SRX routing
    'k_cr','tau_cr','F_cr','SRXFromR2HighStrain','sSRXrip', ...
    ... % geometry               (L_hbare degenerate with ka; LXBpivot reset by getParams)
    'L_thick','L_thin', ...
    ... % lattice: the two that are unmeasured. The other four are observations
    ... % and are now fixed -- and cycle 5, at SL ~2.0, is the only manoeuvre in
    ... % the objective these two can actually be identified from.
    'd_optimal','sigma_lattice', ...
    'Lsc0','gamma','k_pas', ...
    ... % strain-dependence shape (this is what sets FV)
    'PieceWiseStrainDep2X_logOffset','PieceWiseStrainDepR1DX_logOffset', ...
    'PieceWiseStrainDepR21X_logOffset','PieceWiseStrainDepX_logOffset', ...
    'PieceWiseStrainDepX__2','PieceWiseStrainDepX__3','PieceWiseStrainDepX__4', ...
    'PieceWiseStrainDep2X__2','PieceWiseStrainDep2X__3','PieceWiseStrainDep2X__4','PieceWiseStrainDep2X__5', ...
    'PieceWiseStrainDepR1DX__1','PieceWiseStrainDepR1DX__2','PieceWiseStrainDepR1DX__3','PieceWiseStrainDepR1DX__4', ...
    'PieceWiseStrainDepR21X__2','PieceWiseStrainDepR21X__3','PieceWiseStrainDepR21X__4', ...
    'PieceWiseStrainDepParams__2','PieceWiseStrainDepParams__3','PieceWiseStrainDepParams__4', ...
    'PieceWiseStrainDep2Params__2','PieceWiseStrainDep2Params__3','PieceWiseStrainDep2Params__4', ...
    'PieceWiseStrainDep2Params__5','PieceWiseStrainDep2Params__6','PieceWiseStrainDep2Params__7', ...
    'PieceWiseStrainDepR1DParams__1','PieceWiseStrainDepR1DParams__2', ...
    'PieceWiseStrainDepR1DParams__3','PieceWiseStrainDepR1DParams__4', ...
    'PieceWiseStrainDepR21Params__2','PieceWiseStrainDepR21Params__3','PieceWiseStrainDepR21Params__4'};

% The core trio, the F-V shoulder, the realignment gain, the passive lever, and
% -- now that rsK80 carries the whole post-re-stretch window on its own -- the
% SRX drain that sets recovery speed.
COMPULSORY_OVERRIDE = {'k2','kstiff2','ka','v_ref_reg','k_cr','eta_M','kmsrd'};

fprintf(['RunLean5Campaign  campaign %s, profile %s -> tag %s\n' ...
         '  %s\n' ...
         '  seed %s (expect %.4f)\n' ...
         '  rsK80 scored on cycles %s, s = %.1f %%\n' ...
         '  budget %.1f h | N_DRAW %d | SIMPLEX %d | KICK %.2f | SURR %d | workers %d\n\n'], ...
        SPEC.name, upper(PROFILE), TAG, SPEC.note, SEEDSNAP, SEED_EXPECT, ...
        mat2str(PT.seg{kR}), 100*PT.s(kR), ...
        TIME_HRS, N_DRAW, SIMPLEX_EVALS, KICK_FRAC, SURR_EVALS, NWORKERS);

run(fullfile(root, 'Analyses', 'FeatureWeighting', 'RunOptimSpreadW.m'));
