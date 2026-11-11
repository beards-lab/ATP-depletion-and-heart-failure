% RunLeanOptim.m -- 2 h local refit of the reduced model on the lean objective.
%
% Seed:      params/lean_seed.m  (BuildLeanSeed.m)
% Mechanisms: A2Reattaching, MaxwellTensionOnly, DynamicPassive deleted
% Objective:  coolDownLS dropped, rsK80 doubled, FV share 0.40 -> 0.25
% Pool:       params0.leanPool -- inert, deleted and redundant parameters removed
%
% WHY ROUNDS OF 30 AND NOT ONE LONG ONE. optimizeFeatures checks
% TIME_BUDGET_HRS between rounds and never inside one, so a single long round
% cannot be stopped on time (this is what let one ablation arm run 5.7 h).
% Short rounds make the 2 h budget the binding constraint, as intended, and
% each round redraws the parameter subset, which on a 70-parameter pool is the
% point.
%
%   matlab -batch "run('Analyses/MechanismAblation/RunLeanOptim.m')"
%   matlab -batch "TIME_HRS=2; TAG='lean'; run('Analyses/MechanismAblation/RunLeanOptim.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('TAG','var');       TAG      = 'lean';       end
if ~exist('SEEDSNAP','var');  SEEDSNAP = 'lean_seed';  end
if ~exist('TIME_HRS','var');  TIME_HRS = 2;            end
if ~exist('NWORKERS','var');  NWORKERS = 5;            end

params0 = getParams(); s = loadParams(SEEDSNAP);
f = fieldnames(s); for i = 1:numel(f); params0.(f{i}) = s.(f{i}); end
params0.RunForceVelocity=true; params0.RunKtr=false; params0.RunSlack=true;
params0.RunStairs=false; params0.RunForceVelocityTime=false;
params0.RunSlackPassive=true; params0.EvalFeatures=true;
params0.BreakOnODEUnstable=false; params0.PlotEachSeparately=0;
params0.PlotFeatureFitting=0; params0.RunSlackSegments='AllPar';
params0.MaxRunTime=150; params0.OptimizeOn='Feats';
params0.mods = {}; params0.g = [];

pool = params0.leanPool;
keep = true(1,numel(pool));
for i = 1:numel(pool)
    nm = pool{i}; us = strfind(nm,'__');
    if isempty(us)
        keep(i) = isfield(params0,nm) && isscalar(params0.(nm)) && params0.(nm) ~= 0;
    else
        b = nm(1:us(1)-1); ix = str2double(nm(us(1)+2:end));
        keep(i) = isfield(params0,b) && numel(params0.(b)) >= ix && params0.(b)(ix) ~= 0;
    end
end
if any(~keep); fprintf('dropping %d zero-base params: %s\n', nnz(~keep), strjoin(pool(~keep),', ')); end
pool = pool(keep);

cfg = struct();
cfg.baseSnap   = fullfile('params',[SEEDSNAP '.m']);
cfg.params0    = params0;
cfg.tag        = TAG;
cfg.fn         = params0.fn;
cfg.pool       = pool;
% The levers the lean objective actually needs on the table every round: the
% core trio, the F-V shoulder, the realignment gain, the passive lever, and --
% now that rsK80 carries the re-stretch on its own -- the SRX drain that sets
% recovery speed.
cfg.compulsory = {'k2','kstiff2','ka','v_ref_reg','k_cr','eta_M','kmsrd'};
cfg.N_DRAW          = 12;
cfg.SIMPLEX_EVALS   = 30;
cfg.SURR_EVALS      = 0;
cfg.KICK_FRAC       = 0.10;
cfg.MaxRunTime      = 150;
cfg.N_ROUNDS        = 1000;          % the clock stops it, not this
cfg.TIME_BUDGET_HRS = TIME_HRS;
cfg.DEBUG           = false;
cfg.RESUME          = isfile(fullfile(root,'params',[TAG '_state.mat']));

if isempty(gcp('nocreate'))
    try
        cl = parcluster('local');
        jsl = fullfile(tempdir, sprintf('mjob_%s_%d', TAG, feature('getpid')));
        if ~isfolder(jsl); mkdir(jsl); end
        cl.JobStorageLocation = jsl;
        if cl.NumWorkers < NWORKERS; cl.NumWorkers = NWORKERS; end
        parpool(cl, NWORKERS);
    catch
        if isempty(gcp('nocreate')); parpool('local', NWORKERS); end
    end
end

rep = checkFeatureCoverage(params0);
if ~rep.ok
    fprintf(2,'UNSCORABLE fn entries:\n');
    if ~isempty(rep.missingModel); fprintf(2,'  model: %s\n', strjoin(rep.missingModel,', ')); end
    if ~isempty(rep.missingData);  fprintf(2,'  data : %s\n', strjoin(rep.missingData ,', ')); end
    error('RunLeanOptim:coverage','objective has unscorable terms');
end

pchk = params0; pchk.mods = {}; pchk.g = ones(1,0);
[c0,~,ct0,nm0] = evaluateBakersExp(1, pchk, true, false);
fprintf('\n[%s] seed %s = %.4f | pool %d | budget %.1f h\n', TAG, SEEDSNAP, c0, numel(pool), TIME_HRS);
[~,o] = sort(ct0,'descend');
for i = o(1:min(10,numel(o)))
    fprintf('        %-22s %7.4f  (%4.1f %%)\n', nm0{i}, ct0(i), 100*ct0(i)/c0);
end

optimizeFeatures(cfg);
