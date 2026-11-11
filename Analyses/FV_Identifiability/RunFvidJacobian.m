% RunFvidJacobian.m -- the model evaluations behind the two feature Jacobians.
%
%   TIER = 'beard'    every evaluation is three simulations, one per [ATP], of a
%                     Beard2022 isovelocity record (~2 min each); the model trace
%                     is saved whole, so features are chosen afterwards
%   TIER = 'battery'  every evaluation is the slack battery (slack + ktr +
%                     passive + F-V), and the raw feature structs are saved
%
%   matlab -batch "TIER='beard'; NW=4; run('Analyses/FV_Identifiability/RunFvidJacobian.m')"
%
% One checkpoint file per task in results/jac_<TIER>/, so a run can be stopped
% and resumed. ONLY = [..] restricts to those evaluations (timing, debugging).

here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
cd(root); addpath(genpath(root));
% never let a worktree copy shadow the model (genpath picks up .claude/worktrees)
wt = strsplit(genpath(fullfile(root, '.claude')), pathsep); wt = wt(~cellfun(@isempty, wt));
if ~isempty(wt), rmpath(wt{:}); end
assert(startsWith(which('getParams'), root), 'getParams resolves outside this checkout: %s', which('getParams'));
if ~exist('TIER', 'var'), TIER = 'beard'; end
if ~exist('NW', 'var'), NW = 4; end
if ~exist('ONLY', 'var'), ONLY = []; end

[pBat, pBeard, pool, recs] = fvidBase();
E = fvidEvalList(pool, TIER);
outP = fullfile(here, 'results', ['jac_' TIER]); if ~isfolder(outP), mkdir(outP); end
save(fullfile(outP, 'evals.mat'), 'E', 'pool');

% task list: (evaluation, record) for beard, (evaluation, 1) for the battery
nR = 1; if strcmp(TIER, 'beard'), nR = numel(recs); end
[ie, ir] = ndgrid(1:numel(E), 1:nR); T = [ie(:), ir(:)];
if ~isempty(ONLY), T = T(ismember(T(:, 1), ONLY), :); end
ckf = @(t) fullfile(outP, sprintf('e%03d_r%d.mat', t(1), t(2)));
files = arrayfun(@(i) ckf(T(i, :)), (1:size(T, 1))', 'uni', 0);
todo = find(~cellfun(@isfile, files));
fprintf('%s: %d evaluations, %d tasks, %d to run on %d workers\n', TIER, numel(E), size(T, 1), numel(todo), NW);
if isempty(todo), return; end

pl = gcp('nocreate'); if isempty(pl) || pl.NumWorkers ~= NW, delete(pl); parpool('local', NW); end
Tt = T(todo, :); Ft = files(todo); Et = E(Tt(:, 1)); isBeard = strcmp(TIER, 'beard');
tAll = tic;
parfor ii = 1:numel(todo)
    e = Et(ii); r = struct('e', e, 'rec', Tt(ii, 2), 'ok', false, 'msg', '', 'tSim', NaN, 'F', [], 'fm', [], 'fd', []);
    t0 = tic;
    try
        if isBeard
            r.F = simBeard(pBeard, e, recs(Tt(ii, 2)));
        else
            p = pBat;
            if ~isempty(e.name), p.mods = cellstr(e.name); p.g = e.f; end
            [r.fm, r.fd] = batteryFeatures(p);
        end
        r.ok = true;
    catch err
        r.msg = err.message;
    end
    r.tSim = toc(t0);
    parsaveR(Ft{ii}, r);
    lab = e.name; if iscell(lab), lab = sprintf('synth(%d)', numel(lab)); end
    fprintf('%s %s x%.3f rec %d: ok=%d %.0f s %s\n', datestr(now, 'HH:MM:SS'), lab, e.f(1), Tt(ii, 2), r.ok, r.tSim, r.msg); %#ok<TNOW1,DATST>
end
fprintf('%s done: %d tasks in %.1f min\n', TIER, numel(todo), toc(tAll) / 60);

function parsaveR(f, r) %#ok<INUSD>
save(f, 'r');
end
