% RunRefit5Optim.m -- fifth refit, 2026-09-22, on the SHAPE-floored objective.
%
% WHAT CHANGED UNDER IT (Analyses/KtrSegmentSlope/BuildRefit5Seed.m):
%   * every level-dominated force feature re-floored on `s := 2 * s_shape`.
%     The old floors were built from the 21-26 % difference in absolute force
%     between the three hearts (80.6 / 68.9 / 50.6 kPa) while the shape each
%     feature describes agrees to 1.6-2.6 %; since fl goes as s^2 they were
%     ~100x too wide and refit4_opt sat inside all of them.
%   * new level-free term `ktr_shape[2:5] = ktr/ktr(1)`, so the slack ladder's
%     SLOPE is scored for the first time. All three preparations fall across
%     it (-4.9 / -23.3 / -15.2 %); the model rises.
%   Seed cost 0.7988 (refit4_opt's parameters, unchanged, re-scored).
%
% WHAT THE CAMPAIGN IS POINTED AT. The objective is now dominated by the force
% block's shape -- vall_y 24 %, vall_t 15 %, peak2 12 %, peak1_y 11 %,
% vall2_dy 9 %, steady 8 % -- and by ktr_shape 4 %.
%
% WHY P_bound_max IS NOT PRE-SET. It is the biggest single lever on the ktr
% slope (-3.00 at x0.25) but it is welded to F-V: the scan in conclusions.md
% sec. 9 shows every value that buys slope costs 10-500 units of F-V, against a
% ktr_shape term worth 0.032. It stays in the pool and stays free; the slope has
% to come from the cheap levers instead, which is why three of them are
% compulsory below.
%
% COMPULSORY SET. k2 / kstiff2 / ka carry the force level and so the whole
% re-floored block; PieceWiseStrainDep2X_logOffset (-3.77 at +15 %), ekSE
% (-1.71) and dr2 (-1.43) are the three cheap ktr-slope levers; kmsrd (-1.31)
% is the surviving force-gated SRX exit; k_cr is the re-stretch realignment that
% vall_y and peak2 run through.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunRefit5Optim.m')"
%   matlab -batch "TIME_HRS=2; run('Analyses/KtrSegmentSlope/RunRefit5Optim.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));

if ~exist('TAG','var');       TAG      = 'refit5';      end
if ~exist('SEEDSNAP','var');  SEEDSNAP = 'refit5_seed'; end
if ~exist('TIME_HRS','var');  TIME_HRS = 6.0;           end
if ~exist('NWORKERS','var');  NWORKERS = 5;             end

params0 = getParams(); s = loadParams(SEEDSNAP);
f = fieldnames(s); for i = 1:numel(f); params0.(f{i}) = s.(f{i}); end
params0.RunForceVelocity=true; params0.RunSlack=true; params0.RunStairs=false;
params0.RunKtr=true;                       % the k63 fence needs the protocol
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

COMP = {'k2','kstiff2','ka','PieceWiseStrainDep2X_logOffset','ekSE','dr2','kmsrd','k_cr'};
missing = COMP(~ismember(COMP, pool));
if ~isempty(missing)
    error('RunRefit5Optim:compulsory', ...
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
    error('RunRefit5Optim:coverage','objective has unscorable terms');
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
