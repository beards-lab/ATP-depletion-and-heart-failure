% BuildPooledSeed.m -- move the fit target from ONE preparation to the pooled
% three, and rebuild the objective weights from the pooled spreads.
%
% This is a knob-setting wrapper. All the machinery is
% ../FeatureWeighting/BuildWeightedObjective.m; the point of the file is that
% the exact configuration behind params/pooled8_seed.m is in version control
% rather than typed at a prompt before a multi-hour run.
%
% THREE THINGS CHANGE AT ONCE, and it matters which is which:
%
%   1. THE TARGET VALUES.  POOLED_TARGET = true routes the slack block through
%      Model/pooledTargetFeatures, so `features_data` becomes the 3-preparation
%      pooled mean at the pooled level Fref = 69.35 kPa instead of 03/27's own
%      79.76 kPa. Net effect on the force features: a uniform -13.1 % level
%      step PLUS a shape difference. The shape difference is the interesting
%      part -- the pool's first re-stretch peak sits 3 % BELOW its isometric
%      where 03/27's sits 12 % above, so the pooled target pushes directly
%      against the model's known lengthening overshoot.
%
%   2. THE WEIGHTS.  SPREAD_SOURCE = 'pooled' takes s_f from the same pooled
%      analysis: Fss22-anchored, rundown-corrected, variance-moderated. The
%      moderation is the substantive change -- it compresses the 1/s^2 dynamic
%      range 446x -> 13x, so `steady` (4 % raw, a within-run repeatability
%      masquerading as a between-prep spread) stops taking the largest weight,
%      and `vall_t` / `peak1_dSL` (84 % / 76 % raw) stop being unscored.
%
%   3. THE FV PIN.  FV_SHARE 0.70 -> 0.40. Deliberate, and it does confound
%      attribution, so it is called out here rather than buried: at 0.70 the
%      fourteen non-FV data terms held 6.6 % of the objective between them, so
%      the slack shape could not be bought at all. The pooled target's whole
%      contribution is in that slack block; leaving FV at 0.70 would have made
%      the campaign another FV refit. FV's own residual is already inside the
%      disagreement between its two source curves.
%
% AND ONE ESTIMATOR SWAP.  FN_RENAME {'rsK','rsK80'}. Not a new term and not a
% reweighting: the same quantity, measured properly. Analyses/RestretchMetricAudit
% showed rsK reads its amplitude off the tail of the fit window (a 60 ms window
% shift moved it +77 %); rsK80 = ln4/(t80-t20) takes the amplitude from `steady`.
% The pooled target is built on rsK80, so scoring rsK against a pooled value
% would pair an audited target with an unaudited estimator. Its pooled spread
% is also better: 16.9 % against rsK's 20.8 %.
%
% NOT changed: rsT20 is NOT added to the objective. It is the audited
% re-stretch DELAY and the model's worst single feature (rsT0 z = -4.1), but
% adding a term while the target moves would make the result unreadable. It is
% reported as a candidate at the bottom of the weight table; add it in the NEXT
% campaign with ADD_FEATS = {'rsT20'}.
%
% Output -> params/pooled8_seed.m
%           ../FeatureWeighting/results/weight_table.txt, feature_weighting.png
%
% Then: RunOptimPooled.m

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

BASE          = 'estiff125_opt';      % current best on the single-prep objective
SEEDTAG       = 'pooled8_seed';
POOLED_TARGET = true;                 % -> params0.UsePooledTarget
SPREAD_SOURCE = 'pooled';
FN_RENAME     = {'rsK', 'rsK80'};
ADD_FEATS     = {};
FV_SHARE      = 0.40;
SHARE_CAP     = 0.30;

% The force level is now the POOLED one, and the guard's range is the span the
% three preparations actually cover (58.3 - 79.8 kPa corrected), widened 10 %
% at each end so it is a wall and not a target. With FV at 0.40 the level is
% carried mostly by the force features themselves, so this is insurance only.
ANCHOR_GUARD  = true;
ANCHOR_ENTRY  = 'steady[1]|52.5-87.7|5';

fprintf(['BuildPooledSeed: %s -> params/%s.m\n' ...
         '  target = pooled (3 preps) | spreads = pooled/moderated | FV pin %.2f\n' ...
         '  rename %s -> %s\n\n'], BASE, SEEDTAG, FV_SHARE, FN_RENAME{1}, FN_RENAME{2});

run(fullfile(root, 'Analyses', 'FeatureWeighting', 'BuildWeightedObjective.m'));
