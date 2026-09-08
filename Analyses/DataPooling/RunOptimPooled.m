% RunOptimPooled.m -- refit against the 3-preparation POOLED 8 mM target.
%
% The engine is ../FeatureWeighting/RunOptimSpreadW.m (same parameter pool,
% same optimizeFeatures configuration); this file only points it at the pooled
% seed. Keeping one engine means the pool list is not duplicated and the two
% campaigns are searching the same space -- which is the only way the totals
% are comparable at all.
%
% WHAT IS BEING FITTED, AND WHAT IT COSTS AT THE SEED
%
% The target is Model/pooledTarget8mM: three 2026 protocol days, each force
% feature divided by that preparation's own pre-slack isometric, the run
% recorded second in each session lifted for rundown, features pooled per
% preparation, spreads moderated, and the whole force block put back on the
% ONE pooled level Fref = 69.35 kPa. Built by BuildPooledSeed.m; regenerate
% the target itself with WritePooledTarget.m when a new preparation lands.
%
% estiff125_opt, the best fit to the single preparation, scores
%       0.9547 against 03/27 alone      ->      4.4309 against the pool
% and the jump is NOT mostly the -13 % level step. In chi-square units (cost
% divided by the measured between-prep spread, per slack segment):
%
%       peak1_y      6.6   <- the lengthening overshoot. THE defect.
%       peak1_dSL    5.2   <- and where that peak sits relative to the release
%       vall_t       4.4
%       A            2.3
%       vall_y       1.7
%       peak2        1.3     vall2_dy 1.3     steady 1.3
%       restretchSlopeStart 0.4   t0 0.3   ktr 0.15   rsK80 0.06   FV 0.02
%
% Eight of sixteen terms are already at or below one preparation-width, so
% there is no broad gain to be had: the whole residual is the shape and timing
% of the re-stretch transient. The pooled target is what made this legible --
% the pool's first peak sits 3 % BELOW its own isometric where 03/27's sits
% 12 % above, so 03/27 was quietly rewarding the overshoot.
%
% WATCH FOR: rsK80 = 0.06 and ktr = 0.15 are essentially exact. If either
% climbs while peak1_y falls, the search is trading a solved feature for the
% broken one and the result is not progress. featCost per iteration is in
% params/pooled8_state.mat -> state.optIter; the per-iteration snapshot
% gallery is params/pooled8_iter/. Both are what analyzeOptIterUniqueness
% reads.
%
% Usage:
%   matlab -batch "TIME_HRS=4; run('Analyses/DataPooling/RunOptimPooled.m')"
% Resume is automatic if params/pooled8_state.mat exists.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('TIME_HRS','var'); TIME_HRS = 4; end
if ~exist('NWORKERS','var'); NWORKERS = 4; end

TAG      = 'pooled8';
SEEDSNAP = fullfile('params', 'pooled8_seed.m');
% Seed total from BuildPooledSeed.m. The engine rejects a seed that differs by
% more than 0.25 -- the objective's own peak-to-peak numerical roughness, not a
% tolerance to be tightened (Docs objective_noise_floor).
SEED_EXPECT = 4.4309;

run(fullfile(root, 'Analyses', 'FeatureWeighting', 'RunOptimSpreadW.m'));
