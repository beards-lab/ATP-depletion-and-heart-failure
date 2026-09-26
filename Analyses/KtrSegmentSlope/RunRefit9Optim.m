% RunRefit9Optim.m -- refit9, 2026-09-24 02:30: refit8's ATP-carrying objective
% under the FIXED numerics (Analyses/StrainWindowNoise): SiLU rectifier,
% RelTol 1e-3, dense 0.5 ms output, extrapolation clamp -- objective noise
% 19.5 % -> 0.3 %, so CONFIRM_N is back to 0. Floors nEl*s^2/3 on each ATP
% term's own raw between-prep spread; force shape, vall_t, peak1_dSL as guards.
% Seed refit9_seed (refit8_opt re-scored: 0.2422). MaxRunTime 300: dense output
% at RelTol 1e-3 costs ~1.5x per evaluation.
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunRefit9Optim.m')"
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunRefit9Optim.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));

if ~exist('TAG','var');       TAG      = 'refit9';      end
if ~exist('SEEDSNAP','var');  SEEDSNAP = 'refit9_seed'; end
if ~exist('TIME_HRS','var');  TIME_HRS = 12.0;          end
if ~exist('NWORKERS','var');  NWORKERS = 5;             end

params0 = getParams(); s = loadParams(SEEDSNAP);
f = fieldnames(s); for i = 1:numel(f); params0.(f{i}) = s.(f{i}); end
params0.RunForceVelocity=true; params0.RunSlack=true; params0.RunStairs=false;
params0.RunKtr=true;                       % the k63 fence needs the protocol
params0.RunForceVelocityTime=false; params0.RunSlackPassive=true;
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false;
params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0;
params0.RunSlackSegments='AllPar'; params0.MaxRunTime=300;
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

COMP = {'k2','ka','ksrd','sigma_srd2','PieceWiseStrainDep2X_logOffset','ekSE','dr2','L_thin','k_cr','tau_cr'};
missing = COMP(~ismember(COMP, pool));
if ~isempty(missing)
    error('RunRefit9Optim:compulsory', ...
          'compulsory params absent from the pool: %s', strjoin(missing,', '));
end

cfg = struct();
cfg.baseSnap   = fullfile('params',[SEEDSNAP '.m']);
cfg.params0    = params0;
cfg.tag        = TAG;
cfg.fn         = params0.fn;
cfg.pool       = pool;
cfg.compulsory = COMP;
cfg.N_DRAW          = 12;
cfg.SIMPLEX_EVALS   = 30;
cfg.SURR_EVALS      = 0;
cfg.KICK_FRAC       = 0.10;
cfg.CONFIRM_N       = 0;      % the objective is reproducible to ~0.3 % now; no confirmation needed
cfg.CONFIRM_EPS     = 1e-6;
cfg.MaxRunTime      = 300;
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

% Coverage check on the SERIAL slack path. With 'AllPar' this is the first
% parfeval call of the session, so the workers pay JIT compilation inside the
% same MaxRunTime clock; when that tripped the watchdog on 2026-09-23, fetching
% the failed future crashed MATLAB outright with an access violation in
% parallel::plr::parfeval::handleWorkerAborted. The check needs one evaluation,
% not a fast one -- run it serially and leave 'AllPar' to the campaign, whose
% pool is warm by then.
pcov = params0; pcov.RunSlackSegments = 'All';
rep = checkFeatureCoverage(pcov);
if ~rep.ok
    fprintf(2,'UNSCORABLE fn entries:\n');
    if ~isempty(rep.missingModel); fprintf(2,'  model: %s\n', strjoin(rep.missingModel,', ')); end
    if ~isempty(rep.missingData);  fprintf(2,'  data : %s\n', strjoin(rep.missingData ,', ')); end
    error('RunRefit9Optim:coverage','objective has unscorable terms');
end

pchk = params0; pchk.mods = {}; pchk.g = ones(1,0);
[c0,~,ct0,nm0] = evaluateBakersExp(1, pchk, true, false);
fprintf('\n[%s] seed %s = %.4f | pool %d | compulsory %d | budget %.1f h | RESUME %d\n', ...
        TAG, SEEDSNAP, c0, numel(pool), numel(cfg.compulsory), TIME_HRS, cfg.RESUME);
[~,o] = sort(ct0,'descend');
for i = o(1:min(12,numel(o)))
    if ct0(i) <= 0; continue; end
    fprintf('        %-24s %7.4f  (%4.1f %%)\n', nm0{i}, ct0(i), 100*ct0(i)/c0);
end
iK = find(strcmp(nm0,'ktr_shape'),1);
if ~isempty(iK)
    fprintf('        %-24s %7.4f  <- the ladder slope, scored for the first time\n', ...
            'ktr_shape', ct0(iK));
end

optimizeFeatures(cfg);
fprintf('\nCAMPAIGN DONE\n');
