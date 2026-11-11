%% RunYieldRefit.m -- one campaign of the controlled pair (BuildYieldSeeds).
% Same optimiser, draw sizes and numerics as RunRefit9Optim; only the seed and
% the compulsory set differ. The compulsory set is the re-stretch machinery:
% attachment (ka, A0), weak-state detachment and reversal (kd, k_1), the two
% stiffnesses, the load-realign deficit (k_cr, tau_cr) -- plus s_yield in B.
%   matlab -batch "TAG='yieldA'; run('Analyses/OverstretchedHeads/RunYieldRefit.m')"
%   matlab -batch "TAG='yieldB'; run('Analyses/OverstretchedHeads/RunYieldRefit.m')"
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
if ~exist('TAG', 'var');      TAG      = 'yieldB';           end
if ~exist('SEEDSNAP', 'var'); SEEDSNAP = [TAG '_seed'];      end
if ~exist('TIME_HRS', 'var'); TIME_HRS = 3.0;                end
if ~exist('NWORKERS', 'var'); NWORKERS = 5;                  end

params0 = getParams(); s = loadParams(SEEDSNAP);
f = fieldnames(s); for i = 1:numel(f); params0.(f{i}) = s.(f{i}); end
params0.RunForceVelocity = true; params0.RunSlack = true; params0.RunStairs = false; params0.RunKtr = true;
params0.RunForceVelocityTime = false; params0.RunSlackPassive = true;
params0.EvalFeatures = true; params0.BreakOnODEUnstable = false;
params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0;
params0.RunSlackSegments = 'AllPar'; params0.MaxRunTime = 300;
params0.OptimizeOn = 'Feats'; params0.ghostSave = ''; params0.ghostLoad = '';
params0.SaveBest = 0; params0.mods = {}; params0.g = ones(1,0);

pool = params0.leanPool;
keep = true(1, numel(pool));
for i = 1:numel(pool)
    nm = pool{i}; us = strfind(nm, '__');
    if isempty(us)
        keep(i) = isfield(params0, nm) && isscalar(params0.(nm)) && params0.(nm) ~= 0;
    else
        b = nm(1:us(1)-1); ix = str2double(nm(us(1)+2:end));
        keep(i) = isfield(params0, b) && numel(params0.(b)) >= ix && params0.(b)(ix) ~= 0;
    end
end
pool = pool(keep);
COMP = {'ka', 'kd', 'k_1', 'kstiff1', 'kstiff2', 'k_cr', 'tau_cr', 'A0'};
if params0.UseYieldGate, COMP{end+1} = 's_yield'; end
assert(all(ismember(COMP, pool)), 'compulsory params absent from the pool');

cfg = struct();
cfg.baseSnap = fullfile('params', [SEEDSNAP '.m']);
cfg.params0 = params0; cfg.tag = TAG; cfg.fn = params0.fn;
cfg.pool = pool; cfg.compulsory = COMP;
cfg.N_DRAW = 12; cfg.SIMPLEX_EVALS = 30; cfg.SURR_EVALS = 0; cfg.KICK_FRAC = 0.10;
cfg.CONFIRM_N = 0; cfg.CONFIRM_EPS = 1e-6; cfg.MaxRunTime = 300;
cfg.N_ROUNDS = 1000; cfg.TIME_BUDGET_HRS = TIME_HRS; cfg.DEBUG = false;
cfg.RESUME = isfile(fullfile(root, 'params', [TAG '_state.mat']));

if isempty(gcp('nocreate'))
    cl = parcluster('local');
    jsl = fullfile(tempdir, sprintf('mjob_%s_%d', TAG, feature('getpid')));
    if ~isfolder(jsl); mkdir(jsl); end
    cl.JobStorageLocation = jsl;
    if cl.NumWorkers < NWORKERS; cl.NumWorkers = NWORKERS; end
    parpool(cl, NWORKERS);
end
pchk = params0; pchk.mods = {}; pchk.g = ones(1,0);
[c0, ~, ct0, nm0] = evaluateBakersExp(1, pchk, true, false);
fprintf('\n[%s] seed %s = %.4f | pool %d | compulsory %s | budget %.1f h | RESUME %d\n', ...
        TAG, SEEDSNAP, c0, numel(pool), strjoin(COMP, ','), TIME_HRS, cfg.RESUME);
[~, o] = sort(ct0, 'descend');
for i = o(1:min(10, numel(o)))
    if ct0(i) <= 0; continue; end
    fprintf('        %-24s %7.4f  (%4.1f %%)\n', nm0{i}, ct0(i), 100*ct0(i)/c0);
end
optimizeFeatures(cfg);
fprintf('\nCAMPAIGN DONE\n');
