%% RunSlipCampaign.m -- one arm of the A1-slip campaign (Analyses/OverstretchedHeads sec. 10).
%
% ARMS (PROFILE), all on ONE objective: the campaigns' common fn with the
% re-stretch timing read BETWEEN samples (peak1_dSLq, vall_tq; sec. 10a):
%   a  yield gate + A1 slip-rebinding to the neighbouring actin site
%      (YieldHopShift = d_actin) with the rest unbinding to PD; the split
%      YieldHop1 and the yield strain s_yield are fitted
%   b  a + lengthening recruitment (RegAvailLengtheningGain fitted)
%   c  control: no gate, from yieldA_opt -- can the existing parameters move
%      the timing now that the objective can SEE sub-sample changes?
% Seeds: params/slip{A,B,C}_seed.m, written by BuildSlipSeeds.m; their costs are
% recorded below (SEED_EXPECT) for the pre-flight.
%
% USAGE (repo root; one MATLAB session per arm, ~5 workers each):
%   matlab -batch "PROFILE='a'; DRYRUN=true; run('Analyses/OverstretchedHeads/RunSlipCampaign.m')"
%   matlab -batch "PROFILE='a'; TIME_HRS=20; run('Analyses/OverstretchedHeads/RunSlipCampaign.m')"
%   matlab -batch "PROFILE='a'; SMOKE=true;  run('Analyses/OverstretchedHeads/RunSlipCampaign.m')"
% or all arms at once with launch_slip.sh / launch_slip.ps1 / slip.sbatch.
% Resume is automatic: re-running an arm continues from params/slip_<p>_state.mat.
%
% DRYRUN  pre-flight only: data files present, every fn term scorable, and the
%         seed re-scores within SEED_TOL of its recorded cost on this machine.
% SMOKE   one tiny round (3 draws, 4 simplex evaluations) under tag
%         slip_<p>_smoke, to prove the loop, the state file and the snapshot
%         work end to end. Delete params/slip_<p>_smoke* afterwards.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cd(root); addpath(genpath(root));
if ~exist('PROFILE', 'var');  PROFILE  = 'a';   end
if ~exist('TIME_HRS', 'var'); TIME_HRS = 20;    end
if ~exist('NWORKERS', 'var'); NWORKERS = 5;     end
if ~exist('DRYRUN', 'var');   DRYRUN   = false; end
if ~exist('SMOKE', 'var');    SMOKE    = false; end
if ~exist('SEED_TOL', 'var'); SEED_TOL = 0.05;  end   % relative; cross-machine offsets are ~1-2 %
if ~exist('MAXRUN', 'var');   MAXRUN   = 600;   end   % per-solve guard (s): the gated arms ran 210-380 s per AllPar evaluation under load; a false timeout feeds the search a fake penalty

COMP0 = {'ka', 'kd', 'k_1', 'kstiff1', 'kstiff2', 'k_cr', 'tau_cr', 'A0'};
switch lower(PROFILE)
    case 'a', SEEDSNAP = 'slipA_seed'; SEED_EXPECT = 0.5848; COMP = [COMP0, {'s_yield', 'YieldHop1'}];
    case 'b', SEEDSNAP = 'slipB_seed'; SEED_EXPECT = 0.7678; COMP = [COMP0, {'s_yield', 'YieldHop1', 'RegAvailLengtheningGain'}];
    case 'c', SEEDSNAP = 'slipC_seed'; SEED_EXPECT = 0.5457; COMP = COMP0;
    otherwise, error('RunSlipCampaign:profile', 'PROFILE must be a, b or c (got "%s")', PROFILE);
end
TAG = ['slip_' lower(PROFILE)]; if SMOKE, TAG = [TAG '_smoke']; end
fprintf('[%s] arm %s | seed %s | %s | budget %.1f h | %d workers | release %s\n', TAG, upper(PROFILE), ...
        SEEDSNAP, ternary(DRYRUN, 'DRY RUN', ternary(SMOKE, 'SMOKE', 'PRODUCTION')), TIME_HRS, NWORKERS, version('-release'));

params0 = getParams(); s = loadParams(SEEDSNAP);
f = fieldnames(s); for i = 1:numel(f); params0.(f{i}) = s.(f{i}); end
params0.RunForceVelocity = true; params0.RunSlack = true; params0.RunStairs = false; params0.RunKtr = true;
params0.RunForceVelocityTime = false; params0.RunSlackPassive = true;
params0.EvalFeatures = true; params0.BreakOnODEUnstable = false;
params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0;
params0.RunSlackSegments = 'AllPar'; params0.MaxRunTime = MAXRUN;
params0.OptimizeOn = 'Feats'; params0.ghostSave = ''; params0.ghostLoad = '';
params0.SaveBest = 0; params0.mods = {}; params0.g = ones(1,0);

%% ---- pre-flight 1: data files (data/ is not in git) ----------------------
need = {'velocitytableonfile', 'ktr_velocitytableonfile', 'passive_velocitytableonfile'};
miss = {};
for q = need
    if isfield(params0, q{1}) && ~isempty(params0.(q{1})) && ~isfile(fullfile(root, 'data', params0.(q{1})))
        miss{end+1} = params0.(q{1}); %#ok<SAGROW>
    end
end
if ~isempty(miss), error('RunSlipCampaign:data', 'missing in data/: %s', strjoin(miss, ', ')); end
PT = pooledTarget8mM();
assert(all(ismember({'peak1_dSLq', 'vall_tq'}, PT.name)), ...
    'Model/pooledTarget8mM.m has no vertex features -- run Analyses/OverstretchedHeads/UpdatePooledTargetVertex.m');

%% ---- the pool ---------------------------------------------------------------
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
assert(all(ismember(COMP, pool)), 'compulsory params absent from the pool: %s', strjoin(COMP(~ismember(COMP, pool)), ', '));

if isempty(gcp('nocreate'))
    cl = parcluster('local');
    jsl = fullfile(tempdir, sprintf('mjob_%s_%d', TAG, feature('getpid')));
    if ~isfolder(jsl); mkdir(jsl); end
    cl.JobStorageLocation = jsl;             % one job store per session: parallel sessions must not share it
    if cl.NumWorkers < NWORKERS; cl.NumWorkers = NWORKERS; end
    parpool(cl, NWORKERS);
end

%% ---- pre-flight 2: coverage (serial, so worker JIT does not hit MaxRunTime) --
pcov = params0; pcov.RunSlackSegments = 'All'; pcov.MaxRunTime = 1800;   % one serial pass; the gated arms exceed 300 s here under load
rep = checkFeatureCoverage(pcov);
if ~rep.ok
    error('RunSlipCampaign:coverage', 'unscorable fn terms: model {%s} data {%s}', ...
          strjoin(rep.missingModel, ', '), strjoin(rep.missingData, ', '));
end

%% ---- pre-flight 3: the seed scores what it scored where it was built -------
pchk = params0; pchk.mods = {}; pchk.g = ones(1,0);
[c0, ~, ct0, nm0] = evaluateBakersExp(1, pchk, true, false);
% Expected seed costs are written into the switch above (BuildSlipSeeds, R2023a,
% 2026-09-25), because results/*.mat is gitignored and would not reach another
% machine; a local results/slip_seeds.mat from a rebuild takes precedence.
sf = fullfile(root, 'Analyses', 'OverstretchedHeads', 'results', 'slip_seeds.mat');
if isfile(sf)
    E = load(sf); j = find(strcmp(E.seeds, SEEDSNAP)); if ~isempty(j), SEED_EXPECT = E.cost(j); end
end
fprintf('[%s] seed %.4f here vs %.4f recorded (%+.1f %%)\n', TAG, c0, SEED_EXPECT, 100*(c0/SEED_EXPECT - 1));
if abs(c0/SEED_EXPECT - 1) > SEED_TOL
    error('RunSlipCampaign:seed', 'seed cost %.4f is more than %.0f %% off the recorded %.4f -- code or data differ', ...
          c0, 100*SEED_TOL, SEED_EXPECT);
end
[~, o] = sort(ct0(:), 'descend');
for i = o(1:min(8, numel(o)))'
    if ct0(i) > 0, fprintf('        %-24s %7.4f  (%4.1f %%)\n', nm0{i}, ct0(i), 100*ct0(i)/c0); end
end
fprintf('[%s] pool %d | compulsory %s\n', TAG, numel(pool), strjoin(COMP, ','));
if DRYRUN, fprintf('[%s] PRE-FLIGHT OK\n', TAG); return; end

%% ---- the search -------------------------------------------------------------
cfg = struct();
cfg.baseSnap = fullfile('params', [SEEDSNAP '.m']);
cfg.params0 = params0; cfg.tag = TAG; cfg.fn = params0.fn;
cfg.pool = pool; cfg.compulsory = COMP;
cfg.N_DRAW = 12; cfg.SIMPLEX_EVALS = 30; cfg.SURR_EVALS = 0; cfg.KICK_FRAC = 0.10;
cfg.CONFIRM_N = 0; cfg.CONFIRM_EPS = 1e-6; cfg.MaxRunTime = MAXRUN;
cfg.N_ROUNDS = 1000; cfg.TIME_BUDGET_HRS = TIME_HRS; cfg.DEBUG = false;
if SMOKE
    cfg.N_DRAW = numel(COMP) + 1; cfg.SIMPLEX_EVALS = 4; cfg.N_ROUNDS = 1; cfg.TIME_BUDGET_HRS = 1;
end
cfg.RESUME = isfile(fullfile(root, 'params', [TAG '_state.mat']));
optimizeFeatures(cfg);
fprintf('\n[%s] CAMPAIGN DONE\n', TAG);

function r = ternary(c, a, b)
if c, r = a; else, r = b; end
end
