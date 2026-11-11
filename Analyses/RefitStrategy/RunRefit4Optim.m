% RunRefit4Optim.m -- fourth refit, 2026-09-21, on the RE-CURATED pooled target.
%   The peak-prominence fix (Analyses/RestretchPeakAudit) cut peak1_dSL's and
%   vall_t's between-prep spread by 3x; floors and weights re-derived from the
%   new spreads. See Analyses/RefitStrategy/BuildRefit4Seed.m.
%   restart from ktrfence_opt (refit2 stalled 9 rounds from ktrfence_seed),
%   objective = refit2's + the PS_* block floored at seed x1.2 so passive can
%   be neither bought nor sold. See Analyses/RefitStrategy/BuildRefit3Seed.m.
%   restart from ktrfence_seed (NOT ktrfence_opt), + A_over_steady term,
%   + slack floors x0.75, FV floor x0.5, ktr fence retargeted on the 3-prep mean.
%   See Analyses/RefitStrategy/BuildRefit2Seed.m for the rationale of each.
%
% Seed:      params/ktrfence_seed.m  (BuildKtrFenceSeed.m), cost 0.6997
% Model:     lean_opt + UseMaxwellTensionOnly restored + tau_reg 0.2 ms
% Objective: lean (floored, POOLED, rsK80, no coolDownLS, FV share 0.25) with
%            every floor rescaled to s^2/n, and an inert ktr fence.
% Pool:      lean_seed's leanPool + x_M_slack (the mechanism is live again).
%
% WHY ROUNDS AND NOT ONE LONG ONE. optimizeFeatures checks TIME_BUDGET_HRS
% between rounds and never inside one, so short rounds make the clock the
% binding constraint, and each round redraws the subset -- the point on a
% 77-parameter pool.
%
% x_M_slack IS COMPULSORY. It is the only parameter in the pool whose current
% value is transplanted from a different parametrisation (fastgate_best) and
% has never been fitted on this objective: BuildLeanSeed zeroed it when it
% deleted the mechanism, and this seed brings the mechanism back.
%
%   matlab -batch "run('Analyses/RefitStrategy/RunKtrFenceOptim.m')"
%   matlab -batch "TIME_HRS=2; run('Analyses/RefitStrategy/RunKtrFenceOptim.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));

if ~exist('TAG','var');       TAG      = 'refit4';      end
if ~exist('SEEDSNAP','var');  SEEDSNAP = 'refit4_seed'; end
if ~exist('TIME_HRS','var');  TIME_HRS = 6.0;             end
if ~exist('NWORKERS','var');  NWORKERS = 5;               end

params0 = getParams(); s = loadParams(SEEDSNAP);
f = fieldnames(s); for i = 1:numel(f); params0.(f{i}) = s.(f{i}); end
params0.RunForceVelocity=true; params0.RunSlack=true; params0.RunStairs=false;
params0.RunKtr=true;                       % the fence needs the protocol simulated
params0.RunForceVelocityTime=false; params0.RunSlackPassive=true;
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false;
params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0;
params0.RunSlackSegments='AllPar'; params0.MaxRunTime=150;
params0.OptimizeOn='Feats'; params0.ghostSave=''; params0.ghostLoad='';
params0.SaveBest=0; params0.mods = {}; params0.g = ones(1,0);

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
if any(~keep)
    fprintf('dropping %d zero-base params: %s\n', nnz(~keep), strjoin(pool(~keep),', '));
end
pool = pool(keep);

cfg = struct();
cfg.baseSnap   = fullfile('params',[SEEDSNAP '.m']);
cfg.params0    = params0;
cfg.tag        = TAG;
cfg.fn         = params0.fn;
cfg.pool       = pool;
cfg.compulsory = {'k2','kstiff2','ka','v_ref_reg','k_cr','eta_M','kmsrd','x_M_slack'};
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
    error('RunKtrFenceOptim:coverage','objective has unscorable terms');
end

pchk = params0; pchk.mods = {}; pchk.g = ones(1,0);
[c0,~,ct0,nm0] = evaluateBakersExp(1, pchk, true, false);
fprintf('\n[%s] seed %s = %.4f | pool %d | compulsory %d | budget %.1f h | RESUME %d\n', ...
        TAG, SEEDSNAP, c0, numel(pool), numel(cfg.compulsory), TIME_HRS, cfg.RESUME);
[~,o] = sort(ct0,'descend');
for i = o(1:min(10,numel(o)))
    if ct0(i) <= 0; continue; end
    fprintf('        %-24s %7.4f  (%4.1f %%)\n', nm0{i}, ct0(i), 100*ct0(i)/c0);
end
iF = find(strcmp(nm0,'ktrP_k63'),1);
fprintf('        %-24s %7.4f  <- the fence, inert at the seed\n', 'ktrP_k63', ct0(iF));

optimizeFeatures(cfg);
fprintf('\nCAMPAIGN DONE\n');
