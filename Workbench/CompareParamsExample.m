    % CompareParamsExample  Circular ("radar") comparison of active dxdt parameters.
%
%   Demonstrates comparePUParams across several paramsets in the available
%   normalisations. Run from the repo root after:
%       cd('C:\home\git\ATP-depletion-and-heart-failure'); addpath(genpath('.'));
%
%   Normalisation modes:
%     'minmax' - per-axis min-max across sets; spotlights which set is largest
%                on each param (per-axis self-scaled).
%     'log'    - shared radial log10 scale; magnitude-aware fingerprint.
%     'first'  - ratio to the FIRST set: set 1 is the unit reference circle,
%                radius = value/value_set1, rings marked 0.5x/1x/2x..., and the
%                first set's actual value is printed next to each axis label.
%
%   See also: comparePUParams, loadParams, getParams

%% Pick paramsets to compare -- the BASELINE goes first (matters for 'first').
snapshots = {'opt2state_v2_opt', 'optfull4_opt', 'optfull_FourthSlackTimebase_opt', 'passiveCoupling_opt'};
sets  = cellfun(@loadParams, snapshots, 'UniformOutput', false);
labs  = {'v2opt (base)', 'restretch', 'reseeded'};
labs = snapshots;

%% 1) Default: per-axis min-max (who is biggest on each param).
comparePUParams(sets, [], struct('labels',{labs}, 'title','minmax'));

%% 2) Magnitude-aware: shared radial log scale.
comparePUParams(sets, [], struct('norm','log', 'labels',{labs}, ...
    'title','log (magnitude-aware)'));

%% 3) NEW -- ratios vs the first (baseline) set.
%   The baseline is the dashed 1x circle; each other set bulges out (>1x) or
%   pulls in (<1x). The baseline's absolute value is printed next to each label.
comparePUParams(sets, [], struct('norm','first', 'labels',{labs}, ...
    'title','ratio to baseline'));

%% 4) Ratios vs baseline, restricted to a hand-picked set of key levers.
comparePUParams(sets, {'ka','k1','k2','kstiff1','kstiff2','ksr0','A0','tau_reg'}, ...
    struct('norm','first', 'labels',{labs}, 'title','key levers vs baseline'));
