% AnalyzeUniquenessExample.m
% =========================================================================
% Worked example for analyzeOptIterUniqueness -- inspect the param_opt_iter
% cloud that optimizeFeatures accumulates (state.optIter + params/<tag>_iter/)
% and decide whether the solutions clustered around the BEST cost are UNIQUE
% (parameters reconverge each time the optimiser -- incl. a random kick --
% settles) or DEGENERATE (different parameter vectors reach the same cost, i.e.
% a sloppy/unidentifiable manifold).
%
% Run from the repo root after (or during) a RunOptimFull run:
%   cd('C:\home\git\ATP-depletion-and-heart-failure'); addpath(genpath('.'));
%
% See also: analyzeOptIterUniqueness, optimizeFeatures, comparePUParams
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));

tag = 'optfull_FourthSlackTimebase';   % <<< the cfg.tag used in RunOptimFull

%% 1) Full readout: radar + identifiability + feature-influence + verdict.
%    The feature panel is drawn from the PER-FEATURE costs optimizeFeatures
%    stores per record (rec.featCost) -- no re-simulation, seconds to render.
%    It shows which features got better/worse across the near-best cloud.
out = analyzeOptIterUniqueness(tag);

fprintf('\nVerdict: %s\n', out.verdict);
fprintf('near-best solutions: %d | median pairwise param spread: %.1f%%\n', ...
    numel(out.cost), 100*out.medPairSpread);
if ~isempty(out.sloppy);     fprintf('sloppy directions : %s\n', strjoin(out.sloppy, ', ')); end
if ~isempty(out.identified); fprintf('identified (pinned): %d params\n', numel(out.identified)); end

%% 2) Legacy states with no stored featCost: runFeatures=false disables the
%    re-simulation FALLBACK, giving a params-only readout (radar +
%    identifiability). New runs store featCost, so their feature panel renders
%    either way.
% out = analyzeOptIterUniqueness(tag, struct('runFeatures', false));

%% 3) Tighten/loosen the near-optimal band, or use only improved incumbents.
% out = analyzeOptIterUniqueness(tag, struct('relTol', 0.05));            % within 5% of best
% out = analyzeOptIterUniqueness(tag, struct('absTol', 0.5));             % within +0.5 cost
% out = analyzeOptIterUniqueness(tag, struct('improvedOnly', true));      % converged incumbents only

%% 4) Restrict the fingerprint to a hand-picked set of key levers.
% out = analyzeOptIterUniqueness(tag, struct('runFeatures',false, ...
%     'paramList', {{'ka','k1','k2','kstiff2','ksr0','kmsrd','A0','tau_reg'}}));
