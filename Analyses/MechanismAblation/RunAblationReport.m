% RunAblationReport.m -- consolidate the ablation arms into one verdict table.
%
% The statistic is the arm's re-optimised cost against the CONTROL arms'
% re-optimised cost, never against the un-re-optimised 0.9081. sigma is taken
% across the control replicates: it is what one round of the same budget buys
% for no mechanism change at all, and nothing smaller than it is a result.
%
% Also prints the raw chi-2 movers per arm, because the floored objective pays
% nothing inside one preparation-width: a mechanism can be deleted at zero
% scalar cost while a solved term visibly degrades (plan.md trap T4).
%
%   matlab -batch "run('Analyses/MechanismAblation/RunAblationReport.m')"
%
% See also: RunAblationArm, ablationSpec

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

D = dir(fullfile(here, 'results', 'arms', 'abl_*.mat'));
if isempty(D); error('RunAblationReport:noResults','no arm results in results/arms'); end

R = struct('id',{},'name',{},'note',{},'craw',{},'c0',{},'c1',{},'dN',{},'ct0',{},'ct1',{},'nm',{});
for i = 1:numel(D)
    S = load(fullfile(D(i).folder, D(i).name), 'res'); r = S.res;
    R(end+1) = struct('id',r.id,'name',r.name,'note',r.note,'craw',r.cost_raw, ...
        'c0',r.cost_start,'c1',r.cost_end,'dN',r.dNParams, ...
        'ct0',r.featCost_start,'ct1',r.featCost_end,'nm',{r.featNames}); %#ok<SAGROW>
end

isCtrl = startsWith({R.id}, 'CTRL');
ctrl = [R(isCtrl).c1];
% THE THRESHOLD. Two things set it, and the larger wins.
%
% (1) Spread across the control replicates. This came out at ZERO: every
%     control re-optimised the incumbent to 0.9081 exactly. At this budget a
%     round buys the incumbent nothing, which is a result in itself -- it says
%     the arms are not being charged for a restart.
% (2) The objective's own irreproducibility, which is NOT zero and is the real
%     floor. optimizeFeatures re-evaluates its best point before writing a
%     snapshot and refuses if it has drifted; that guard fired on control 3
%     (found 0.8561, re-evaluated 0.9081) and on C2 (0.8488 -> 0.9042). Those
%     are phantom improvements of 0.052 and 0.055 -- measured, on this machine,
%     on this objective. An arm claiming less than that has claimed nothing.
%
% So sigma is floored at 0.055. Anything inside +-2 sigma of the control is
% indistinguishable from the incumbent.
PHANTOM = 0.055;   % largest non-reproducible "improvement" seen in this campaign
if isempty(ctrl)
    warning('RunAblationReport:noControls', ...
        'no control replicate completed -- thresholds fall back to the phantom floor alone');
    ctrlMean = 0.9081; sigRaw = 0;
else
    ctrlMean = mean(ctrl); sigRaw = std(ctrl);
end
sig = max(sigRaw, PHANTOM);
fprintf(['control: n=%d, mean %.4f, spread %.4f\n' ...
         'phantom floor (non-reproducible gains seen in this campaign): %.4f\n' ...
         'threshold = 2 sigma = %.4f\n\n'], numel(ctrl), ctrlMean, sigRaw, PHANTOM, 2*sig);

A = R(~isCtrl);
d = [A.c1] - ctrlMean;
[~, ord] = sort(d, 'ascend');

fprintf('%-4s %-24s %10s %8s %8s %8s %6s  %s\n', ...
        'id','arm','raw','calib','reopt','vs ctrl','dPar','verdict');
fprintf('%s\n', repmat('-', 1, 108));
for i = ord
    v = iVerdict(d(i), 2*sig, A(i).dN);
    fprintf('%-4s %-24s %10.4g %8.4f %8.4f %+8.4f %+6d  %s\n', ...
        A(i).id, A(i).name, A(i).craw, A(i).c0, A(i).c1, d(i), A(i).dN, v);
end

fprintf('\n=== raw chi-2 movers (unfloored), arm vs its own frozen start ===\n');
for i = ord
    if isempty(A(i).ct0) || isempty(A(i).ct1) || numel(A(i).ct0) ~= numel(A(i).ct1); continue; end
    dc = A(i).ct1 - A(i).ct0;
    [~, o] = sort(dc, 'descend');
    worse = o(1:min(3,numel(o))); better = flip(o(max(1,numel(o)-2):end));
    fprintf('%-4s %-22s  worse: %s   better: %s\n', A(i).id, A(i).name, ...
        iList(A(i).nm, dc, worse), iList(A(i).nm, dc, better));
end

save(fullfile(here,'results','report.mat'), 'R','sig','ctrlMean');
fprintf('\nwritten: %s\n', fullfile(here,'results','report.mat'));

% -------------------------------------------------------------------------
function v = iVerdict(d, thr, dN)
if d < -thr
    v = 'IMPROVES -- promote to a full run';
elseif d > thr
    v = 'LOAD-BEARING -- keep';
elseif dN < 0
    v = sprintf('FREE -- delete (%d params for nothing)', -dN);
elseif dN > 0
    v = sprintf('NOT WORTH IT -- %d params buy nothing', dN);
else
    v = 'NEUTRAL';
end
end

function s = iList(nm, dc, ix)
p = cell(1, numel(ix));
for i = 1:numel(ix); p{i} = sprintf('%s %+.2f', nm{ix(i)}, dc(ix(i))); end
s = strjoin(p, ', ');
end
