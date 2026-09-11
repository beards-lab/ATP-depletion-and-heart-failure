% RunOptimRigAware.m -- the 30-hour refit on the rig-aware objective.
%
% Engine: ../FeatureWeighting/RunOptimSpreadW.m (same parameter pool, same
% optimizeFeatures configuration). This file sets the seed, the budget, and
% ONE of several search PROFILES, so a strong machine can run the profiles
% side by side as independent MATLAB sessions.
%
% =========================================================================
% WHAT IS BEING FITTED
% =========================================================================
% The 3-preparation pooled 8 mM target (Model/pooledTarget8mM), with the
% re-stretch rate scored on slack cycles 1-4 only -- cycle 5 re-stretches to a
% different final length by the smallest amplitude, and injecting each run's
% own measured noise puts its rsK80 at 26-50 % against 4-14 % for cycles 1-4
% (Analyses/RecoveryOscillation). Built by BuildRigAwareSeed.m.
%
% SEED: params/rigaware_seed.m = pooled8_opt reweighted, cost-neutral at
%       2.6858. The 0.0141 below pooled8_opt's 2.6999 is exactly cycle 5's
%       rsK80 residual leaving the objective -- it is deterministic, not the
%       search noise floor.
%
% WHERE THE RESIDUAL IS  (chi-2 = residual in units of the between-prep
% spread; 1.0 = one preparation-width, per slack segment in brackets)
%
%       peak1_y     20.8  [4.2]   <- the lengthening overshoot
%       peak1_dSL   18.8  [3.8]   <- and where that peak sits after the release
%       vall_t      13.9  [2.8]
%       vall2_dy     8.0  [1.6]
%       A            6.7  [1.3]   steady 6.5 [1.3]   ktr 5.2 [1.0]
%       peak2        4.5  [0.9]   restretchSlopeStart 2.7 [0.5]
%       vall_y       2.0  [0.4]   rsK80 1.5 [0.3]    t0 1.5 [0.3]
%       FV           0.06 [0.01]
%
% Nine of seventeen terms are already at or below one preparation-width. The
% whole remaining residual is the SHAPE and TIMING of the re-stretch
% transient, and those four terms are coupled (they are four readings of one
% overshoot), which is why profile B widens the draw: a 10-parameter subset
% may not span the direction that moves all four together.
%
% WATCH FOR: rsK80 = 1.5 and ktr = 5.2 are the ATP-carrying rates and are
% nearly solved. If either climbs while peak1_y falls, the search is trading a
% solved feature for the broken one and the result is not progress. featCost
% per iteration is in params/<TAG>_state.mat -> state.optIter; the snapshot
% gallery is params/<TAG>_iter/. Both are what analyzeOptIterUniqueness reads.
%
% =========================================================================
% HOW TO RUN
% =========================================================================
% One instance, 30 h, resumable (resume is automatic if the state file exists):
%
%   matlab -batch "PROFILE='a'; TIME_HRS=30; run('Analyses/RecoveryOscillation/RunOptimRigAware.m')"
%
% Pre-flight only -- validates the seed, the objective coverage and the pool,
% then stops WITHOUT starting the search. Always do this first: it costs one
% evaluation and it is the difference between finding a broken seed now and
% finding it in 30 hours.
%
%   matlab -batch "PROFILE='a'; DRYRUN=true; run('Analyses/RecoveryOscillation/RunOptimRigAware.m')"
%
% All four profiles at once on a big box (launch_rigaware.ps1 / .sh do this):
%
%   for p in a b c d; do matlab -batch "PROFILE='$p'; TIME_HRS=30; \
%       run('Analyses/RecoveryOscillation/RunOptimRigAware.m')" & sleep 20; done
%
% Then consolidate:  RunRigAwareReport.m
%
% SIZING. One objective evaluation is ~38 s at 5 workers, and a round is
% SIMPLEX_EVALS of them, so profile A is ~40 min/round -> ~45 rounds in 30 h.
% pooled8 improved on every one of its 9 rounds and stopped only because its
% 4-hour budget ran out (stall = 0 at the end), so this budget is buying
% descent that was already happening, not a search for a new basin.
%
% CORES. The engine runs one evaluation at a time and parallelises INSIDE it
% over the five slack segments, so a single instance saturates at ~5-6
% workers. Extra cores are spent by running more PROFILES, not by raising
% NWORKERS. Four instances x 5 workers = 20 cores.
%
% =========================================================================
% THE PROFILES -- each varies HOW the search moves, not what it scores
% =========================================================================
%   a  the pooled8 setting, unchanged. The control: it is the configuration
%      that was still descending monotonically when its budget ran out.
%   b  wider draw (16) and more polish (90). The four open terms are four
%      readings of one overshoot, so the direction that fixes them is
%      probably higher-dimensional than a 10-parameter draw can span.
%   c  narrow (6) and fast (40) with a hard kick (0.25). Many more rounds,
%      and the strongest escape from a flat ridge -- which is how estiff125
%      stalled (Analyses/StairsPrediction: the incumbent was where the search
%      STOPPED, not a minimum).
%   d  profile a with the global surrogate turned back on (120 evals). The
%      other three are local from the incumbent; this one is the only member
%      that can leave the basin.

clc;
% TAGSTEM/SEEDSNAP/SEED_EXPECT are guarded below so one campaign file serves
% both the rigaware and the floored runs without a copy.
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('PROFILE','var');  PROFILE  = 'a'; end
if ~exist('TIME_HRS','var'); TIME_HRS = 30;  end
if ~exist('NWORKERS','var'); NWORKERS = 5;   end

% Seed is a knob: the floored campaign reuses this whole file and only changes
% what it starts from and what that scores.
%   rigaware_seed  2.6858  BuildRigAwareSeed.m   -- the rsK80[1:4] objective
%   floored_seed   0.9081  BuildFlooredSeed.m    -- the same, with every data
%                                                   term given a one-width
%                                                   deadband (sec. 7: the
%                                                   unfloored run put 110 % of
%                                                   its gain into F-V)
if ~exist('SEEDSNAP','var');    SEEDSNAP    = fullfile('params','rigaware_seed.m'); end
if ~exist('SEED_EXPECT','var'); SEED_EXPECT = 2.6858; end
%                       BuildRigAwareSeed.m, cost-neutral at pooled8_opt.
                        % The engine rejects a seed more than 0.25 away -- the
                        % objective's own peak-to-peak roughness, not a
                        % tolerance to tighten (Docs objective_noise_floor).

switch lower(PROFILE)
    case 'a'; N_DRAW=10; SIMPLEX_EVALS=60; KICK_FRAC=0.10; SURR_EVALS=0;
    case 'b'; N_DRAW=16; SIMPLEX_EVALS=90; KICK_FRAC=0.10; SURR_EVALS=0;
    case 'c'; N_DRAW= 6; SIMPLEX_EVALS=40; KICK_FRAC=0.25; SURR_EVALS=0;
    case 'd'; N_DRAW=10; SIMPLEX_EVALS=60; KICK_FRAC=0.10; SURR_EVALS=120;
    otherwise
        error('RunOptimRigAware:profile', ...
              'PROFILE must be a, b, c or d (got "%s").', PROFILE);
end
if ~exist('TAGSTEM','var'); TAGSTEM = 'rigaware'; end
TAG          = sprintf('%s_%s', TAGSTEM, lower(PROFILE));
MAX_RUN_TIME = 300;

if ~isfile(SEEDSNAP)
    error('RunOptimRigAware:noSeed', ...
        ['%s not found. Build it first:\n' ...
         '    OPTS = struct(''segments'', struct(''rsK80'', 1:4));\n' ...
         '    run(''Analyses/DataPooling/WritePooledTarget.m'')\n' ...
         '    run(''Analyses/RecoveryOscillation/BuildRigAwareSeed.m'')'], SEEDSNAP);
end

fprintf(['RunOptimRigAware  profile %s -> tag %s\n' ...
         '  budget %.1f h | N_DRAW %d | SIMPLEX %d | KICK %.2f | SURR %d | workers %d\n\n'], ...
        upper(PROFILE), TAG, TIME_HRS, N_DRAW, SIMPLEX_EVALS, KICK_FRAC, SURR_EVALS, NWORKERS);

run(fullfile(root, 'Analyses', 'FeatureWeighting', 'RunOptimSpreadW.m'));
