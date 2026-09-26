%% RunObjectiveNoise.m -- the full objective's numerical noise under a setting
%
% Same design as RunSlackNoise, but each point is a full evaluateBakersExp
% (slack, ktr protocol, F-V, passive) scored on PARAMSET's own objective --
% the number the refit campaigns see. Set OVR (struct of parameter overrides,
% e.g. struct('odeRelTol',1e-3)) and TAG before calling; OVR = struct() is the
% current behaviour.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit8_opt'; end
if ~exist('NREP', 'var'), NREP = 8; end
if ~exist('OVR', 'var'), OVR = struct(); end
if ~exist('TAG', 'var'), TAG = 'baseline'; end
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
addpath(genpath(root));
res  = fullfile(here, 'results');

P = getParams(loadParams(PARAMSET), [], true, false);
P.RunForceVelocity = true; P.RunSlack = true; P.RunStairs = false; P.RunKtr = true;
P.RunForceVelocityTime = false; P.RunSlackPassive = true; P.EvalFeatures = true;
P.BreakOnODEUnstable = false; P.PlotEachSeparately = 0; P.PlotFeatureFitting = 0;
P.RunSlackSegments = 'All'; P.MaxRunTime = 900; P.OptimizeOn = 'Feats';
P.ghostSave = ''; P.ghostLoad = ''; P.SaveBest = 0; P.mods = {}; P.g = ones(1,0);
f = fieldnames(OVR); for k = 1:numel(f), P.(f{k}) = OVR.(f{k}); end
P = getParams(P, P.g, true, false);

rng(11);
pool = P.leanPool; D = randn(NREP, numel(pool));
if isempty(gcp('nocreate')), parpool(5); end
cost = nan(NREP+1, 1); FC = cell(NREP+1, 1); wall = nan(NREP+1, 1); names = {};
tic;
parfor r = 1:NREP+1
    pv = P;
    if r > 1
        for i = 1:numel(pool)
            nm = pool{i};
            if ~contains(nm, '__') && isfield(pv, nm) && isscalar(pv.(nm)) && pv.(nm) ~= 0
                pv.(nm) = pv.(nm) * (1 + 1e-6*D(r-1, i));
            end
        end
        pv = getParams(pv, pv.g, true, false);
    end
    t0 = tic;
    [c, ~, fc] = evaluateBakersExp(1, pv, true, false);
    cost(r) = c; FC{r} = fc(:)'; wall(r) = toc(t0);
end
[~, ~, ~, names] = evaluateBakersExp(1, P, true, false);
C = vertcat(FC{:});
lab = cellfun(@(s) regexprep(regexprep(s, '\|.*$', ''), '\[.*$', ''), names, 'uni', 0);
fprintf('\n=== full objective, %s, %s (%.0f s, %.0f s/eval) ===\n', PARAMSET, TAG, toc, mean(wall));
fprintf('  base %.4f | mean %.4f | SD %.4f (%.1f %%) | range %.4f..%.4f\n', cost(1), mean(cost), std(cost), ...
        100*std(cost)/mean(cost), min(cost), max(cost));
for j = 1:size(C, 2)
    if mean(C(:,j)) < 1e-6 && std(C(:,j)) < 1e-6, continue; end
    fprintf('  %-22s mean %.4f  SD %.4f (%5.1f %%)\n', lab{j}, mean(C(:,j)), std(C(:,j)), 100*std(C(:,j))/mean(C(:,j)));
end
save(fullfile(res, sprintf('objnoise_%s_%s.mat', PARAMSET, TAG)), 'cost', 'C', 'lab', 'wall', 'OVR');
