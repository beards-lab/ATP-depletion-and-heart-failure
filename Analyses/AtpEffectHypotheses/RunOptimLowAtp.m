% RunOptimLowAtp.m -- refit the low-ATP gates, using all the workers.
%
%   matlab -batch "NWORKERS=25; TIME_HRS=6; run('Analyses/AtpEffectHypotheses/RunOptimLowAtp.m')"
%
% BRANCH: this needs `floored-campaign` (worktree C:\home\git\ATP-wt-floored).
% It is the only branch carrying BOTH the ATP gate machinery and the pooled
% 2 mM target; `lowatp-2state` has the gates but no pooled targets and an older
% core. The pre-flight below refuses to start anywhere else.
%
% =========================================================================
% WHY THIS IS NOT RunOptimSpreadW WITH A DIFFERENT SEED
% =========================================================================
% The 8 mM engine parallelises INSIDE one objective evaluation, over the five
% slack segments, so it saturates at ~5-6 workers and extra cores are spent by
% running more profiles as separate MATLAB sessions. That is the wrong shape
% here, twice over: the ask is to keep it in ONE MATLAB for progress tracking,
% and the low-ATP problem has only a handful of free parameters.
%
% So the parallelism moves outward. This is a PARALLEL PATTERN SEARCH: each
% iteration builds a batch of candidate gate vectors -- both coordinate steps
% and random multi-gate perturbations -- and evaluates the whole batch in one
% parfor. Inside a worker the segment loop runs serially (RunSlackSegments
% 'All', not 'AllPar'), because a nested parfor degrades to a serial loop
% anyway and pretending otherwise just costs pool overhead.
%
% Arithmetic at 25 workers: 25 evaluations per round in roughly the wall time
% of ONE serial evaluation (~2-4 min), against 25 x 40 s = 17 min if they were
% run one at a time on a 5-worker pool. The batch is sized to the pool, so
% every worker gets exactly one candidate and none idles.
%
% =========================================================================
% WHAT IS BEING FITTED
% =========================================================================
% ONLY the ATP gate exponents. The core is frozen and 8 mM is invariant BY
% CONSTRUCTION -- every gate is normalised to 1 at MgATP_ref, so no value of
% any exponent can move the high-ATP fit. The pre-flight asserts that at
% |delta| = 0 rather than trusting it.
%
% Objective: the pooled 2 mM target (Model/pooledTarget2mM), with the weights
% and floors rescaled to the 2 mM between-preparation spread
% (rescaleFnToCondition). Passive PS_* terms are dropped -- they are not
% pooled at 2 mM.
%
% READ THE DIAGNOSIS FIRST (RunLowAtpDiagnosis.m). If most terms already sit
% inside their deadband, extra gates have nothing to buy and a campaign will
% confirm that expensively. The reason C0, C4 and P1 came out identical is a
% question about the OBJECTIVE, not about the optimiser.

clearvars -except NWORKERS TIME_HRS CORE TAG GATES SEED_EXPS; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

if ~exist('NWORKERS','var');  NWORKERS  = 25;            end
if ~exist('TIME_HRS','var');  TIME_HRS  = 6;             end
if ~exist('CORE','var');      CORE      = 'lean5_a_opt'; end
if ~exist('TAG','var');       TAG       = 'lowatp25';    end
% The chemically defensible steps plus the two Pi proxies. k2 is the canonical
% ATP step and stays compulsory; the rest are on trial.
if ~exist('GATES','var')
    GATES = {'k2','k_1','k1','ksrd2sr','kmsrd','ksr0','kah'};
end
if ~exist('SEED_EXPS','var'); SEED_EXPS = struct('k2', 0.44); end

STEP0   = 0.25;    % initial pattern step, in exponent units
STEP_MIN= 0.02;    % stop refining below this
LB      = -1.5;    % exponent bounds. |p| = 1.5 is a x2.8 factor at 2 mM,
UB      =  1.5;    %   already beyond anything physiological
STATE   = fullfile(resDir, [TAG '_state.mat']);

%% ---- branch / machinery pre-flight -------------------------------------
need = {'Model/pooledTarget2mM.m','Model/rescaleFnToCondition.m'};
for i = 1:numel(need)
    assert(isfile(fullfile(root,need{i})), ...
        ['%s is missing -- you are on the wrong branch. This needs ' ...
         'floored-campaign (worktree ATP-wt-floored).'], need{i});
end

if isempty(gcp('nocreate')); parpool('local', NWORKERS); end
pool = gcp; fprintf('[%s] pool of %d workers\n', TAG, pool.NumWorkers);
NBATCH = pool.NumWorkers;

LoadData;

%% ---- the objective -----------------------------------------------------
p0 = getParams(loadParams(CORE), [], true, false);
p0.EvalFeatures=true; p0.OptimizeOn='Feats'; p0.BreakOnODEUnstable=false;
p0.RunKtr=false; p0.RunStairs=false; p0.RunForceVelocityTime=false;
p0.RunSlackPassive=false; p0.RunSlack=true; p0.RunForceVelocity=true;
p0.MaxRunTime=400; p0.PlotEachSeparately=0; p0.PlotFeatureFitting=0;
p0.UsePooledTarget=true; p0.mods={}; p0.g=ones(1,0);
p0.UseAtpPowerGates=true; p0.MgATP_ref=8;
p0.UseAtpDetachR2=0; p0.UseAtpKmsrd=0; p0.UseStrainDepAtpGate=0;
% Explicitly OFF. It cannot re-extract slack features (that branch is disabled)
% but it DOES reach runFVExperiment, where re-extracting the 2 mM F-V data
% takes the FV term from 0.0000 to 10.8651 -- 84 % of the objective. A search
% run with it on optimises that artefact.
p0.recalculateDataFeats = false;
% SERIAL segments: the parallelism is across candidates now. A nested parfor
% would run serially regardless, so asking for 'AllPar' only adds overhead.
p0.RunSlackSegments='All';

prim = @(c) cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), c, 'uni', 0);
fn = p0.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
fn = fn(~startsWith(prim(fn),'PS_'));
fn = rescaleFnToCondition(fn, 2);

p2 = p0; p2.MgATP=2; p2.MgADP=0; p2.Pi=0.8; p2.fn = fn;
p2.velocitytableonfile = 'protocol_03_27_2026_2mM_slack.mat';
for i=1:numel(GATES); p2.(['pAtp_' GATES{i}]) = 0; end

rep = checkFeatureCoverage(p2);
if ~rep.ok
    bad = regexprep(unique([rep.missingModel(:);rep.missingData(:)])','\s*\(in\s.*$','');
    fprintf(2,'[%s] dropping unscorable: %s\n', TAG, strjoin(bad,', '));
    fn = fn(~ismember(prim(fn), bad)); p2.fn = fn;
    assert(checkFeatureCoverage(p2).ok, 'objective still has unscorable terms');
end
fprintf('[%s] %d scored terms at 2 mM\n', TAG, numel(fn));

% Does 'All' score the same as the 'AllPar' the 8 mM work used? If not, this
% optimiser is minimising a different function from the one already reported.
pChk = p2; pChk.RunSlackSegments='AllPar';
cSer = iCost(p2, GATES, zeros(1,numel(GATES)), fn);
cPar = iCost(pChk, GATES, zeros(1,numel(GATES)), fn);
% Threshold is the objective's OWN numerical roughness, not bit-identity. The
% two modes take different execution paths (see the LXBpivot handling in
% runSlackExperiment), so they will not agree to the last bit and demanding
% that just fails on noise. sigma ~ 0.076, peak-to-peak 0.25
% (Docs objective_noise_floor): warn above sigma, refuse above peak-to-peak.
dSeg = abs(cSer - cPar);
fprintf('[%s] segment-mode check: All %.6f vs AllPar %.6f (|d| %.4f)\n', ...
        TAG, cSer, cPar, dSeg);
if dSeg > 0.25
    error('RunOptimLowAtp:segMode', ...
        ['All and AllPar differ by %.4f, beyond the objective''s peak-to-peak ' ...
         'roughness -- the parallel search would minimise a different function'], dSeg);
elseif dSeg > 0.076
    warning('RunOptimLowAtp:segMode', ...
        ['All and AllPar differ by %.4f, above sigma. Inside peak-to-peak, so ' ...
         'the search proceeds, but the final point is re-evaluated under ' ...
         'AllPar before it is quoted.'], dSeg);
end
SEG_OFFSET = cSer - cPar;   % reported again at the end

% 8 mM invariance, asserted not assumed.
p8 = p2; p8.MgATP=8; p8.Pi=0;
p8.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
c8a = iCost(p8, GATES, zeros(1,numel(GATES)), fn);
c8b = iCost(p8, GATES, 0.5*ones(1,numel(GATES)), fn);
fprintf('[%s] 8 mM invariance |d| = %.3e  (MUST be 0)\n', TAG, abs(c8a-c8b));
assert(abs(c8a-c8b) < 1e-9, 'gates leak into the 8 mM fit');

%% ---- start point -------------------------------------------------------
nG = numel(GATES);
x  = zeros(1,nG);
fs = fieldnames(SEED_EXPS);
for i=1:numel(fs)
    j = find(strcmp(GATES, fs{i}),1);
    if ~isempty(j); x(j) = SEED_EXPS.(fs{i}); end
end
if isfile(STATE)
    S = load(STATE); x = S.x; fBest = S.fBest; iter0 = S.iter; step = S.step;
    fprintf('[%s] resumed at iter %d, best %.4f\n', TAG, iter0, fBest);
else
    fBest = iCost(p2, GATES, x, fn); iter0 = 0; step = STEP0;
    fprintf('[%s] seed %.4f   (%s)\n', TAG, fBest, iStr(GATES,x));
end

%% ---- parallel pattern search ------------------------------------------
tStart = tic; iter = iter0; hist = [];
while toc(tStart) < TIME_HRS*3600 && step >= STEP_MIN
    iter = iter + 1;

    % Batch: 2 coordinate steps per gate, remainder random multi-gate moves.
    %
    % ROTATE the coordinate order by iteration. The batch is capped at the pool
    % size, so with 7 gates (14 coordinate directions) and fewer than 14
    % workers a fixed order would test the same leading gates every round and
    % never move the rest -- the search would silently be over a subset of the
    % parameters it claims to fit. Rotating means every direction comes up.
    C = [];
    rot = mod((0:nG-1) + iter - 1, nG) + 1;
    for j = rot
        for sgn = [1 -1]
            c = x; c(j) = c(j) + sgn*step; C(end+1,:) = c;              %#ok<AGROW>
        end
    end
    nCoord = size(C,1);
    while size(C,1) < NBATCH
        C(end+1,:) = x + step*randn(1,nG);                              %#ok<AGROW>
    end
    C = min(max(C(1:NBATCH,:), LB), UB);
    if NBATCH < nCoord && iter == 1
        fprintf(['[%s] NOTE pool (%d) < coordinate directions (%d): the batch ' ...
                 'rotates, so a full sweep takes %d iterations\n'], ...
                TAG, NBATCH, nCoord, ceil(nCoord/NBATCH));
    end

    f = inf(NBATCH,1);
    parfor b = 1:NBATCH
        f(b) = iCost(p2, GATES, C(b,:), fn);
    end

    [fMin, iMin] = min(f);
    if fMin < fBest - 1e-6
        x = C(iMin,:); fBest = fMin;
        fprintf('[%s] iter %3d  step %.3f  BEST %.4f   %s\n', ...
                TAG, iter, step, fBest, iStr(GATES,x));
        iWrite(resDir, TAG, GATES, x, fBest, iter);
    else
        step = step/2;
        fprintf('[%s] iter %3d  no gain (best cand %.4f) -> step %.3f\n', ...
                TAG, iter, fMin, step);
    end
    hist(end+1,:) = [iter, fBest, step];                                %#ok<AGROW>
    save(STATE, 'x','fBest','iter','step','GATES','CORE','hist','fn');
end

%% ---- re-evaluate the winner on the reported objective ------------------
% The search minimised the SERIAL segment mode. Everything else in the project
% quotes 'AllPar'. Re-evaluate there before the number leaves this script --
% same discipline as never quoting a snapshot header.
pFinal = p2; pFinal.RunSlackSegments = 'AllPar';
fFinal = iCost(pFinal, GATES, x, fn);
fprintf('\n[%s] DONE  %.2f h, %d iters\n', TAG, toc(tStart)/3600, iter);
fprintf('[%s] best (All)    %.4f\n', TAG, fBest);
fprintf('[%s] best (AllPar) %.4f   <- QUOTE THIS. pre-flight offset was %+.4f\n', ...
        TAG, fFinal, SEG_OFFSET);
fprintf('[%s] %s\n', TAG, iStr(GATES,x));
save(STATE, 'x','fBest','fFinal','iter','step','GATES','CORE','hist','fn','SEG_OFFSET');

%% ---- helpers -----------------------------------------------------------
function c = iCost(pBase, gates, xv, fn)
% One objective evaluation at a gate vector. A FUNCTION so the scripts it
% calls (LoadData, RunBakersExp) get a workspace of their own -- a script that
% creates variables cannot run directly inside parfor (transparency), and
% RunBakersExp clobbers whatever it touches.
p = pBase;
for i = 1:numel(gates); p.(['pAtp_' gates{i}]) = xv(i); end
params0 = p;                                                    %#ok<NASGU>
try
    LoadData; RunBakersExp;
    c = sum(evalFeatureCost(features_data, features_model, fn));
    if ~isfinite(c); c = 1e6; end
catch
    c = 1e6;      % a failed integration is a bad point, not a crashed run
end
end

function s = iStr(gates, xv)
s = '';
for i = 1:numel(gates)
    if xv(i) ~= 0; s = [s sprintf('%s %+0.3f  ', gates{i}, xv(i))]; end %#ok<AGROW>
end
if isempty(s); s = '(all gates zero)'; end
end

function iWrite(resDir, tag, gates, xv, c, iter)
fid = fopen(fullfile(resDir,[tag '_best.txt']),'w');
fprintf(fid, '%% %s best, iter %d, 2 mM cost %.6f\n', tag, iter, c);
for i = 1:numel(gates)
    fprintf(fid, 'params0.pAtp_%s = %.17g;\n', gates{i}, xv(i));
end
fclose(fid);
end
