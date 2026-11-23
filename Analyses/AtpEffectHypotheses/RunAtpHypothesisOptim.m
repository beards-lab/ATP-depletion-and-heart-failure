% RunAtpHypothesisOptim.m -- PRODUCTION optimiser for the ATP-pathway hypotheses.
%
% Fits the ATP-gate exponents for H1..H4 on a frozen core, one rung at a time,
% and writes a snapshot + state per rung. Designed to be dropped on another
% machine and run unattended.
%
% ---------------------------------------------------------------------------
% QUICK START (another machine, from the repo root)
%
%   matlab -batch "run('Analyses/AtpEffectHypotheses/RunAtpHypothesisOptim.m')"
%
% or with options:
%
%   matlab -batch "CORE='spreadw_iter_018'; RUNGS={'H2'}; BUDGET_H=6; ...
%                  run('Analyses/AtpEffectHypotheses/RunAtpHypothesisOptim.m')"
%
% Dry run (no fitting, ~2 min: validates data, coverage, invariance, cost):
%
%   matlab -batch "DRYRUN=true; run('Analyses/AtpEffectHypotheses/RunAtpHypothesisOptim.m')"
%
% ---------------------------------------------------------------------------
% OPTIONS (set before running; every one has a working default)
%
%   CORE       char   frozen 8 mM parameter snapshot, by name (params/<CORE>.m)
%   RUNGS      cell   which rungs, e.g. {'H1','H2'}       default all four
%   BUDGET_H   double wall-clock hours for the WHOLE run  default 8
%   NWORKERS   double parallel pool size                  default 5
%   ACCEPT     double accept threshold on the cost        default 3x the
%                     measured noise peak-to-peak (0.106) = 0.32
%   TAG        char   output prefix                       default 'atphyp'
%   RESUME     logical continue from a previous state     default true
%   DRYRUN     logical validate only, fit nothing         default false
%
% ---------------------------------------------------------------------------
% WHAT IT DOES, AND WHY IT IS CHEAP
%
% Every ATP gate is exactly 1 at MgATP_ref, so the 8 mM fit is invariant to all
% of them BY CONSTRUCTION. That is asserted numerically before any fitting (it
% must read 0.000e+00) and it is what makes this affordable: only the 2 mM leg
% is ever re-evaluated. One evaluation = 2 mM slack + 2 mM force-velocity.
%
% Coordinate descent on ABSOLUTE exponent values, not multipliers, so a gate can
% leave exactly zero and can change sign -- both of which happened in the pilot
% and neither of which a multiplier-based optimiser can do.
%
% The accept threshold matters. The 2 mM objective carries numerical scatter of
% sigma ~0.036 / peak-to-peak ~0.106 (measured, non-monotone in the perturbation
% => scatter, not sensitivity). ACCEPT defaults to 3x that. Do NOT lower it
% below ~0.11: below the floor the descent chases solver noise, which shows up
% as a parameter oscillating between two values while the cost "improves".
%
% ---------------------------------------------------------------------------
% OUTPUTS (all under Analyses/AtpEffectHypotheses/results/)
%
%   <TAG>_state.mat        resumable state, written after every rung
%   <TAG>_<rung>_opt.m     parameter snapshot per rung (params/)
%   <TAG>_summary.txt      the ladder table, human-readable
%
% ---------------------------------------------------------------------------
% REQUIREMENTS
%   MATLAB with Parallel Computing Toolbox (falls back to serial if absent).
%   data/ must contain, at minimum:
%       protocol_03_27_2026_2mM_slack.mat
%       protocol_03_27_2026_8mM_slack.mat
%   Force-velocity data is compiled into Model/LoadData.m -- no file needed.
%   Checked explicitly at startup; the script stops with a clear message if a
%   file is missing rather than failing 40 minutes in.

% NOTE: deliberately NO `clear` here. This is a script, so `clear` would wipe the
% very options the caller set on the command line -- which silently reverted
% DRYRUN=true to a live 8-hour fit the first time it was tested.
clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root));

%% ---- options --------------------------------------------------------------
if ~exist('CORE','var');     CORE     = 'spreadw_iter_018'; end
if ~exist('RUNGS','var');    RUNGS    = {'H1','H2','H3','H4'}; end
if ~exist('BUDGET_H','var'); BUDGET_H = 8;        end
if ~exist('NWORKERS','var'); NWORKERS = 5;        end
if ~exist('ACCEPT','var');   ACCEPT   = 0.32;     end
if ~exist('TAG','var');      TAG      = 'atphyp'; end
if ~exist('RESUME','var');   RESUME   = true;     end
if ~exist('DRYRUN','var');   DRYRUN   = false;    end

resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end
STATE  = fullfile(resDir, [TAG '_state.mat']);

fprintf('=== ATP hypothesis optimiser ===\n');
fprintf('core     : %s\n', CORE);
fprintf('rungs    : %s\n', strjoin(RUNGS,', '));
fprintf('budget   : %.1f h   accept %.3f   workers %d\n', BUDGET_H, ACCEPT, NWORKERS);
fprintf('dry run  : %d\n\n', DRYRUN);

%% ---- preflight: everything that can fail cheaply, fails now ---------------
NEED = {'protocol_03_27_2026_2mM_slack.mat','protocol_03_27_2026_8mM_slack.mat'};
miss = NEED(~cellfun(@(f) isfile(fullfile('data',f)), NEED));
if ~isempty(miss)
    error('RunAtpHypothesisOptim:missingData', ...
          'missing data file(s) in data/:\n  %s', strjoin(miss, sprintf('\n  ')));
end
if ~isfile(fullfile('params',[CORE '.m']))
    error('RunAtpHypothesisOptim:missingCore','core snapshot params/%s.m not found', CORE);
end
LoadData;
if exist('parpool','file')==2 && isempty(gcp('nocreate'))
    parpool('local', NWORKERS);
end

%% ---- gate classification (decided by chemistry, before fitting) -----------
CHEM  = {'k2','kah','kamh'};                       % nucleotide is a substrate
SRX   = {'kmsrd','ksr0','ksr2srd','ksrd2sr'};      % nucleotide-state transitions
MECH3 = {'eta_M','kSE','mu_neg'};                  % the 3 in the historical fit
REST  = {'kstiff2','kSE_M','dr','dr2','P_bound_max','ka','kd','k1','k_1','tau_reg','A0'};
ALL21 = [CHEM SRX MECH3 REST];

HDEF = struct( ...
  'H1', {{CHEM,                      'ATP chemistry only'}}, ...
  'H2', {{[CHEM SRX],                '+ SRX nucleotide states'}}, ...
  'H3', {{[CHEM SRX MECH3],          '+ mechanical'}}, ...
  'H4', {{ALL21,                     'all 21 (any rate)'}});

%% ---- the 2 mM objective ---------------------------------------------------
p0 = getParams(loadParams(CORE), [], true, false);
p0.UseAtpJointFit = false; p0.UseAtpPowerGates = true; p0.MgATP_ref = 8;
p0.UseAtpDetachR2 = 0; p0.UseAtpKmsrd = 0; p0.UseStrainDepAtpGate = 0;
p0.UseR2Ceiling = 0; p0.UseD0State = 0;
p0.EvalFeatures = true; p0.BreakOnODEUnstable = false; p0.OptimizeOn = 'Feats';
p0.RunSlackSegments = 'AllPar'; p0.MaxRunTime = 400;
p0.PlotEachSeparately = 0; p0.PlotFeatureFitting = 0; p0.PlotFullscreen = 0;
p0.RunKtr = false; p0.RunStairs = false; p0.RunForceVelocityTime = false;
p0.RunSlackPassive = false; p0.RunSlack = true; p0.RunForceVelocity = true;
p0.mods = {}; p0.g = ones(1,0);

% 2 mM has slack AND force-velocity, but no passive protocol -> drop PS_* only.
prim = @(c) cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), c, 'uni', 0);
fn2 = p0.fn; [~,iu] = unique(fn2,'stable'); fn2 = fn2(iu);
fn2 = fn2(~startsWith(prim(fn2),'PS_'));

p2 = p0; p2.MgATP = 2; p2.MgADP = 0; p2.Pi = 0.8;
p2.velocitytableonfile = 'protocol_03_27_2026_2mM_slack.mat';
p2.fn = fn2;

% Any fn entry the 2 mM leg cannot produce costs a flat MISSING_FEATURE_COST=100
% every evaluation and silently dominates. Drop them, loudly.
rep = checkFeatureCoverage(p2);
if ~rep.ok
    bad = unique([rep.missingModel(:); rep.missingData(:)])';
    fprintf(2,'dropping %d unscorable feature(s): %s\n', numel(bad), strjoin(bad,', '));
    fn2 = fn2(~ismember(prim(fn2), bad)); p2.fn = fn2;
    rep = checkFeatureCoverage(p2);
    assert(rep.ok, 'objective still unscorable after dropping');
end
fprintf('2 mM objective: %d features (force-velocity included: %d)\n', ...
        numel(fn2), any(startsWith(prim(fn2),'FV')));

cost = @(q) iCost(q, fn2, ATP_c, Data_ATP);

%% ---- preflight 2: seed cost + the 8 mM invariance that makes this cheap ---
q0 = p2; for g = 1:numel(ALL21); q0.(['pAtp_' ALL21{g}]) = 0; end
tic; cSeed = cost(q0); tEval = toc;
fprintf('seed (all gates 0, i.e. H0): %.4f   [%.0f s/evaluation]\n', cSeed, tEval);

b8 = q0; b8.MgATP = 8; b8.Pi = 0;
b8.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
c8a = cost(b8);
b8g = b8; for g = 1:numel(ALL21); b8g.(['pAtp_' ALL21{g}]) = 0.3; end
c8b = cost(b8g);
fprintf('8 mM invariance |d| = %.3e  (MUST be 0 -- gates are 1 at MgATP_ref)\n', abs(c8b-c8a));
assert(abs(c8b-c8a) < 1e-9, 'ATP gates are leaking into the 8 mM condition');

nEvalTotal = numel(RUNGS) * 4 * 8;   % rough: 4 steps x ~8 gates x 1 sweep
fprintf('\nrough scale: ~%.0f min per full sweep of H2 (7 gates x 4 steps)\n', 7*4*tEval/60);

if DRYRUN
    fprintf('\nDRY RUN OK -- data present, objective scorable, invariance exact.\n');
    fprintf('Remove DRYRUN to fit.\n');
    return
end

%% ---- resume ---------------------------------------------------------------
OUT = struct('rung',{},'dof',{},'c0',{},'c',{},'exps',{},'sweeps',{},'evals',{});
if RESUME && isfile(STATE)
    S = load(STATE); OUT = S.OUT;
    fprintf('resumed: %d rung(s) already done (%s)\n', numel(OUT), strjoin({OUT.rung},', '));
end

%% ---- fit each rung --------------------------------------------------------
todo = RUNGS(~ismember(RUNGS, {OUT.rung}));
if isempty(todo); fprintf('nothing left to do.\n'); else
budgetPer = BUDGET_H*3600 / numel(todo);
tAll = tic; nEval = 0;

for r = 1:numel(todo)
    rung = todo{r};
    gates = HDEF.(rung){1}; lbl = HDEF.(rung){2};
    q = p2; for g = 1:numel(ALL21); q.(['pAtp_' ALL21{g}]) = 0; end
    q.UseAtpPowerGates = true;

    % warm start from the previous (nested) rung if we have it
    prev = find(ismember({OUT.rung}, {'H1','H2','H3'}), 1, 'last');
    if ~isempty(prev)
        for g = 1:numel(ALL21)
            f = ['pAtp_' ALL21{g}];
            if ismember(ALL21{g}, gates) && isfield(OUT(prev).exps, f)
                q.(f) = OUT(prev).exps.(f);
            end
        end
    end

    c = cost(q); nEval = nEval + 1; c0 = c;
    fprintf('\n=== %s: %s (%d gates) seed %.4f ===\n', rung, lbl, numel(gates), c);

    tR = tic; improved = true; sweep = 0; fullSweep = false;
    while improved && toc(tR) < budgetPer
        improved = false; sweep = sweep + 1; fullSweep = true;
        steps = [0.15 -0.15 0.30 -0.30]; if sweep > 1; steps = [0.07 -0.07]; end
        for gi = 1:numel(gates)
            % Running out of budget mid-sweep is NOT convergence: the gates after
            % this one were never tested. Record that, or a truncated rung gets
            % reported as "converged" and its untested gates are wrongly called
            % worthless -- which is exactly what happened to H3/eta_M in the pilot.
            if toc(tR) > budgetPer; fullSweep = false; break; end
            f = ['pAtp_' gates{gi}]; base = q.(f);
            bestC = c; bestV = base;
            for s = steps
                cand = q; cand.(f) = base + s;
                cc = cost(cand); nEval = nEval + 1;
                if cc < bestC - ACCEPT; bestC = cc; bestV = cand.(f); end
                if toc(tR) > budgetPer; break; end
            end
            if bestV ~= base
                fprintf('  [%5.1f min|%4d] %-12s %+.3f -> %+.3f   %.4f (%+.4f)\n', ...
                        toc(tAll)/60, nEval, gates{gi}, base, bestV, bestC, bestC-c);
                q.(f) = bestV; c = bestC; improved = true;
            end
        end
    end
    % CONVERGED means a FULL sweep completed with nothing accepted. Without the
    % fullSweep term, a rung that ran out of budget partway through its gate list
    % is reported as converged and its UNTESTED gates are wrongly concluded to be
    % worthless.
    conv = ~improved && fullSweep;
    fprintf('  %s: %.4f -> %.4f  (%d sweeps, %.1f min, %s)\n', rung, c0, c, sweep, ...
            toc(tR)/60, iIf(conv,'CONVERGED','TRUNCATED - not all gates tested'));

    exps = struct();
    for g = 1:numel(ALL21); exps.(['pAtp_' ALL21{g}]) = q.(['pAtp_' ALL21{g}]); end
    OUT(end+1) = struct('rung',rung,'dof',numel(gates),'c0',c0,'c',c, ...
                        'exps',exps,'sweeps',sweep,'evals',nEval); %#ok<SAGROW>
    save(STATE,'OUT','CORE','fn2','ACCEPT');

    snap = q; snap.MgATP = 8; snap.Pi = 0;
    snap.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
    snap.mods = {}; snap.g = [];
    writeParamsToMFile(fullfile(root,'params',sprintf('%s_%s_opt.m',TAG,rung)), snap, [], ...
        sprintf('%s %s on core %s: 2 mM cost %.4f (%d gates, %s)', ...
                rung, lbl, CORE, c, numel(gates), iIf(conv,'converged','budget')), true);
end
end

%% ---- summary --------------------------------------------------------------
txt = fullfile(resDir, [TAG '_summary.txt']);
fid = fopen(txt,'w');
for f = [1 fid]
    fprintf(f,'\n================ ATP HYPOTHESIS LADDER ================\n');
    fprintf(f,'core: %s   accept: %.3f   features: %d\n\n', CORE, ACCEPT, numel(fn2));
    fprintf(f,'%-6s %-26s %5s %10s %10s %9s\n','rung','gates','DOF','seed','fitted','per DOF');
    for k = 1:numel(OUT)
        fprintf(f,'%-6s %-26s %5d %10.3f %10.3f %9.3f\n', OUT(k).rung, ...
                HDEF.(OUT(k).rung){2}, OUT(k).dof, OUT(k).c0, OUT(k).c, ...
                OUT(k).c/max(OUT(k).dof,1));
    end
    fprintf(f,'\nfitted exponents:\n');
    for g = 1:numel(ALL21)
        v = nan(1,numel(OUT));
        for k = 1:numel(OUT); v(k) = OUT(k).exps.(['pAtp_' ALL21{g}]); end
        if any(abs(v) > 1e-9)
            fprintf(f,'  %-14s', ALL21{g}); fprintf(f,' %8.3f', v); fprintf(f,'\n');
        end
    end
end
fclose(fid);
fprintf('\nwrote %s\n', txt);
fprintf('ATP HYPOTHESIS OPTIM DONE\n');

%% ==========================================================================
function c = iCost(q, fn2, ATP_c, Data_ATP)
% One evaluation = 2 mM slack + 2 mM force-velocity, scored on fn2.
%
% A single stiff/slow parameter combination must NOT kill a multi-hour
% unattended run: evaluateModel throws when it exceeds MaxRunTime, and an
% uncaught throw here loses every rung still to come. Return a large finite
% penalty instead, so the descent simply rejects that step and carries on.
try
    [~,~,fm,fd] = runSlackExperiment(q);
    qf = q; qf.RunSlack = false;
    [~,~,fmv,fdv] = runFVExperiment(qf, ATP_c, Data_ATP);
    f = fieldnames(fmv); for i=1:numel(f); fm.(f{i}) = fmv.(f{i}); end
    f = fieldnames(fdv); for i=1:numel(f); fd.(f{i}) = fdv.(f{i}); end
    c = sum(evalFeatureCost(fd, fm, fn2, 2));
    if ~isfinite(c); c = 1e6; end
catch ME
    fprintf(2,'    [eval failed: %s -> penalty]\n', ME.message);
    c = 1e6;
end
end

function s = iIf(c,a,b); if c; s=a; else; s=b; end; end
