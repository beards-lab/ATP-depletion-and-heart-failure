%% BuildSlipSeeds.m -- seeds and objective for the A1-slip campaign (RunSlipCampaign.m).
%
% OBJECTIVE = the yield campaigns' common fn (timing trio floored at nEl*s^2/3,
% BuildYieldSeeds) with the two sample-quantised timing terms replaced by their
% sub-sample versions (sec. 10a): peak1_dSL -> peak1_dSLq, vall_t -> vall_tq,
% same rule fl = nEl*s^2/3, w = LAMBDA/s^2, s = the pooled moderated spread.
%
% ARMS
%   slipA_seed  yieldB_opt + A1 slip to the neighbouring actin site
%               (YieldHopShift = d_actin = 5.5 nm), YieldHop1 = HOP_A (the rest
%               unbinds to PD); YieldHop1 joins the pool, s_yield is already in it
%   slipB_seed  slipA + RegAvailLengtheningGain = GAIN_B (joins the pool)
%   slipC_seed  yieldA_opt (no gate) -- control
% Each seed is scored here (AllPar) and the cost stored in results/slip_seeds.mat;
% RunSlipCampaign refuses to start if a machine re-scores it more than 5 % away.
%   matlab -batch "run('Analyses/OverstretchedHeads/BuildSlipSeeds.m')"
if ~exist('HOP_A', 'var'),   HOP_A   = 0.98;   end   % just below 1, so the split can move both ways
if ~exist('SY_A', 'var'),    SY_A    = 0.012;  end   % slip 5.5 nm + s_yield 12 nm = 0.501 with no re-fit (sec. 10b); [] keeps yieldB_opt's 14.1 nm
if ~exist('GAIN_B', 'var'),  GAIN_B  = 0.005;  end   % half of the no-refit valley-matching dose (sec. 8)
LAMBDA = 0.0106294;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
res = fullfile(here, 'results');

PT = pooledTarget8mM();
assert(all(ismember({'peak1_dSLq', 'vall_tq'}, PT.name)), 'run UpdatePooledTargetVertex.m first');
sq = @(nm) PT.s(strcmp(PT.name, nm));
fn = loadParams('yieldA_seed'); fn = fn.fn;
base = @(e) regexprep(regexprep(e, '\|.*$', ''), '\[.*$', '');
SWAP = {'peak1_dSL', 'peak1_dSLq'; 'vall_t', 'vall_tq'};
for k = 1:size(SWAP, 1)
    i = find(cellfun(@(e) ~contains(e, ',') && strcmp(base(e), SWAP{k, 1}), fn));
    assert(numel(i) == 1, 'expected one fn entry for %s', SWAP{k, 1});
    s = sq(SWAP{k, 2}); nEl = 5;
    new = sprintf('%s|fl:%.6g|%.6g', SWAP{k, 2}, nEl*s^2/3, LAMBDA/s^2);
    fprintf('%-40s -> %s   (s %.4f)\n', fn{i}, new, s);
    fn{i} = new;
end

pA = getParams(loadParams('yieldB_opt'), [], true, false);
pA.fn = fn; pA.YieldHopShift = pA.d_actin; pA.YieldHop1 = HOP_A;
if ~isempty(SY_A), pA.s_yield = SY_A; end
pA.leanPool = unique([pA.leanPool, {'YieldHop1', 'YieldHopShift'}], 'stable');   % slip distance free: 5.5 and 11 nm both scored well (sec. 10b)
pB = pA; pB.RegAvailLengtheningGain = GAIN_B;
pB.leanPool = unique([pB.leanPool, {'RegAvailLengtheningGain'}], 'stable');
pC = getParams(loadParams('yieldA_opt'), [], true, false); pC.fn = fn;

NOTE = 'Analyses/OverstretchedHeads sec. 10; objective: yield trio floored + sub-sample peak1_dSLq / vall_tq';
seeds = {'slipA_seed', 'slipB_seed', 'slipC_seed'};
P = {pA, pB, pC};
DESC = {sprintf('yieldB_opt + A1 slip %.1f nm, YieldHop1 %.2f', 1e3*pA.YieldHopShift, HOP_A), ...
        sprintf('slipA + RegAvailLengtheningGain %.4g', GAIN_B), 'yieldA_opt, no gate (control)'};
for k = 1:3
    writeParamsToMFile(fullfile(root, 'params', [seeds{k} '.m']), P{k}, [], [DESC{k} '. ' NOTE]);
end

if isempty(gcp('nocreate')), refreshPool(5); end
cost = nan(1, 3);
for k = 1:3
    p = getParams(loadParams(seeds{k}), [], true, false);
    p.RunForceVelocity = true; p.RunSlack = true; p.RunStairs = false; p.RunKtr = true;
    p.RunForceVelocityTime = false; p.RunSlackPassive = true; p.EvalFeatures = true;
    p.BreakOnODEUnstable = false; p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0;
    p.RunSlackSegments = 'AllPar'; p.MaxRunTime = 300; p.OptimizeOn = 'Feats';
    p.ghostSave = ''; p.ghostLoad = ''; p.SaveBest = 0; p.mods = {}; p.g = ones(1,0);
    [cost(k), ~, ct, nm] = evaluateBakersExp(1, p, true, false);
    [~, o] = sort(ct(:), 'descend');
    fprintf('\n%s = %.4f  (%s)\n', seeds{k}, cost(k), DESC{k});
    for i = o(1:6)', if ct(i) > 0, fprintf('    %-22s %.4f\n', nm{i}, ct(i)); end, end
end
save(fullfile(res, 'slip_seeds.mat'), 'seeds', 'cost', 'DESC', 'fn', 'HOP_A', 'GAIN_B');
