% BuildRigAwareSeed.m -- rebuild the objective now that the ripple on every
% force recovery has been shown to be the rig.
%
% A knob-setting wrapper, exactly like ../DataPooling/BuildPooledSeed.m; all
% the machinery is ../FeatureWeighting/BuildWeightedObjective.m. The point of
% the file is that the configuration behind params/rigaware_seed.m is in
% version control rather than typed at a prompt before a 30-hour run.
%
% WHAT THE MEASUREMENT SAID  (conclusions.md; RunNoiseBudget, RunFeatureNoise)
%
%   The recoveries ARE single exponentials. The residual of a one-exponential
%   fit is 90 % white transducer noise (sigma 1.2-1.6 kPa, identical in the
%   PNB+Mava control with the cross-bridges chemically switched off), 6 % a
%   slow ramp-locked wander whose period does not move with nucleotide, and
%   4 % anything else. A noise-free model cannot exceed R2 ~ 0.83 against this
%   data, where the data's own exponential reaches 0.795 and the model 0.999.
%
% ONE THING CHANGES, and it is not a weight.
%
%   rsK80  ->  rsK80[1:4].  The re-stretch rate is read off the 20-80 % band
%   of the rise, so its precision is set by how well those two crossings are
%   located, i.e. by sigma / amplitude. On slack cycle 5 the re-stretch is to a
%   DIFFERENT final length by the SMALLEST amplitude (0.061 ML against
%   0.080-0.140; Analyses/RestretchMetricAudit sec. 3), and injecting each
%   run's own measured noise 200x puts cycle 5 at 26-50 % while cycles 1-4 sit
%   at 4-14 % -- in every run, at both ATP concentrations. Averaged in, cycle 5
%   takes the feature's noise to 10.7 %, against a between-preparation spread
%   of 10.2 %: with cycle 5 included, rsK80's ENTIRE apparent spread is the
%   transducer and the feature carries no resolvable between-preparation
%   information at all. Dropping it: noise 5.7 %, spread 7.7 %. There is signal
%   again, and the weight is legitimate.
%
%   This is a scoring-window change, not a de-weighting. The weight follows
%   from the tighter spread and goes UP (w 35.0 -> 40.1), which is the right
%   direction: cycle 5 was adding scatter, not information.
%
% NOT changed, deliberately:
%   * the three trace-least-squares terms (coolDownLS, PS_rampupRMSE,
%     PS_holdDecayRMSE). They have a genuine noise floor -- a perfect model
%     still scores sigma -- but RunObjectiveFloor measures that floor at 0.7 %
%     of the objective against the 12.8 % those terms carry. They are 2-6x
%     above their floor, so they are fitting signal, not rig.
%   * rsR2. It is a BOUNDARY (0.4-1.0), the model sits at 0.999 inside it, and
%     it costs exactly 0. Leave it a boundary: turning it into a fit-to-data
%     term would ask the model to reproduce 0.795, which is 90 % noise, and
%     the ceiling for a noise-free model is 0.834. Never score rsR2 as data.
%   * ktr_rmse's 0-0.5 bound, which asks the model's ktr rise to be well
%     described by ONE exponential. The measurement just validated that as the
%     right shape, so this guardrail is now earned rather than assumed.
%
% ADD_RST20 -- an optional, argued addition, default OFF.
%   If the recovery is a pure exponential then it IS three numbers: amplitude,
%   rate, delay. The objective scores the amplitude (vall2_dy, steady) and the
%   rate (rsK80[1:4]) but not the delay, so a third of a description now known
%   to be sufficient is unconstrained. rsT20 is the audited delay, its ATP
%   effect is resolved (z = 2.4) and its noise (13 %) is well inside its spread
%   (21 %). ../DataPooling/BuildPooledSeed.m deferred it to "the next
%   campaign", which is this one. It is OFF by default because adding a term
%   while testing a scoring change would confound the attribution: turn it on
%   deliberately, and read the campaign as two changes rather than one.
%
% Output -> params/rigaware_seed.m
%           ../FeatureWeighting/results/weight_table.txt, feature_weighting.png
%
% PREREQUISITE: Model/pooledTarget8mM.m must have been regenerated with
%   rsK80's spread measured on cycles 1-4, i.e.
%     OPTS = struct('segments', struct('rsK80', 1:4));
%     run('Analyses/DataPooling/WritePooledTarget.m')
%   The check below refuses to build a seed against a target that disagrees.
%
% Then: RunOptimRigAware.m

clearvars -except ADD_RST20; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('ADD_RST20','var'); ADD_RST20 = false; end

%% ---- refuse to build against a stale target ----------------------------
PT = pooledTarget8mM();
k  = find(strcmp(PT.name, 'rsK80'), 1);
if ~isfield(PT, 'seg') || ~isequal(PT.seg{k}, 1:4)
    error('BuildRigAwareSeed:staleTarget', ...
        ['Model/pooledTarget8mM.m has rsK80''s spread on segments %s, not 1:4.\n' ...
         'Regenerate it first:\n' ...
         '    OPTS = struct(''segments'', struct(''rsK80'', 1:4));\n' ...
         '    run(''Analyses/DataPooling/WritePooledTarget.m'')'], ...
        mat2str(PT.seg{k}));
end

BASE          = 'pooled8_opt';        % current best on the pooled objective
SEEDTAG       = 'rigaware_seed';
POOLED_TARGET = true;                 % -> params0.UsePooledTarget
SPREAD_SOURCE = 'pooled';
FN_RENAME     = {'rsK80', 'rsK80[1:4]'};
ADD_FEATS     = {};
if ADD_RST20; ADD_FEATS = {'rsT20'}; end

% Held at the pooled8 values so the campaign is a test of the scoring change
% and not of a re-pinned FV share.
FV_SHARE      = 0.40;
SHARE_CAP     = 0.30;
ANCHOR_GUARD  = true;
ANCHOR_ENTRY  = 'steady[1]|52.5-87.7|5';

fprintf(['BuildRigAwareSeed: %s -> params/%s.m\n' ...
         '  target = pooled (3 preps) | spreads = pooled/moderated | FV pin %.2f\n' ...
         '  scoring window %s -> %s   (cycle 5 is a different manoeuvre at 26-50 %% noise)\n' ...
         '  rsT20 added: %d\n\n'], BASE, SEEDTAG, FV_SHARE, FN_RENAME{1}, FN_RENAME{2}, ADD_RST20);

run(fullfile(root, 'Analyses', 'FeatureWeighting', 'BuildWeightedObjective.m'));
