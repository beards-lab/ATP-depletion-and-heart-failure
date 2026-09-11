% BuildFlooredSeed.m -- take the reward away from the terms that are already
% solved, so the search has nothing left to buy but the overshoot.
%
% WHY. The rigaware campaign (conclusions.md sec. 7) spent 14 hours and put
% 110 % of its gain into force-velocity: F-V held 40 % of the objective while
% sitting at chi-2 0.10, a third of one preparation-width and well inside the
% disagreement between its own two source curves. Strip F-V and the fit got
% 1.8 % WORSE -- peak1_dSL, vall_t, ktr, rsK80 and coolDownLS all went
% backwards. That is not a bad search; it is a search doing exactly what the
% objective asked.
%
% Re-pinning FV_SHARE cannot fix it. 0.70 -> 0.40 was already tried and the
% campaign above is the result: ANY nonzero weight on a smooth, nearly-solved
% term still offers the cheapest descent in the problem. What has to change is
% not how much the term is worth but whether it pays at all once it is inside
% the noise.
%
% SO: FLOOR_CHI2 = 1. Every data term gets a deadband one preparation-width
% wide. Inside it the term contributes nothing; outside it, it pays as before.
% It is a deadband and NOT a deletion -- a floored term still defends itself
% against regression, which is the property that lets F-V keep its pin without
% the pin being exploitable.
%
% WHAT THIS SHOULD DO TO THE OBJECTIVE, at the seed:
%   FV       chi-2 0.07  -> pays 0
%   rsK80    chi-2 0.42  -> pays 0
%   ktr      chi-2 0.97  -> pays ~0
%   peak1_y  chi-2 3.9   -> pays ~74 % of what it did
%   vall_t   chi-2 2.9   -> pays ~65 %
% i.e. the objective stops being 40 % F-V and becomes almost entirely the
% re-stretch transient, which is the only thing left that is actually wrong.
%
% BASE is rigaware_a_opt -- the campaign's best point, re-evaluated HERE at
% 2.4088 (its cluster header says 2.4153; the cluster reads this objective
% +0.015 high throughout). Starting from it keeps the F-V gain that was
% genuinely made, and stops paying for more of it.
%
% Output -> params/floored_seed.m
% Then:     RunOptimRigAware.m with SEEDSNAP/SEED_EXPECT pointed at it.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

BASE          = 'rigaware_a_opt';
SEEDTAG       = 'floored_seed';
POOLED_TARGET = true;
SPREAD_SOURCE = 'pooled';
FLOOR_CHI2    = 1.0;         % one preparation-width of deadband
FN_RENAME     = cell(0,2);   % rsK80[1:4] is already in rigaware_a_opt
ADD_FEATS     = {};
% The FV pin stays where it is. With a floor in place the pin is no longer
% exploitable -- it buys nothing until F-V leaves its own noise band -- so
% there is no reason to argue about the number as well as the mechanism.
FV_SHARE      = 0.40;
SHARE_CAP     = 0.30;
ANCHOR_GUARD  = true;
ANCHOR_ENTRY  = 'steady[1]|52.5-87.7|5';

fprintf(['BuildFlooredSeed: %s -> params/%s.m\n' ...
         '  FLOOR_CHI2 = %.2f preparation-width(s) of deadband on every data term\n\n'], ...
        BASE, SEEDTAG, FLOOR_CHI2);

run(fullfile(root, 'Analyses', 'FeatureWeighting', 'BuildWeightedObjective.m'));
