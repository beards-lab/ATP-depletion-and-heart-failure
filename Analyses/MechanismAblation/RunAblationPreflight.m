% RunAblationPreflight.m -- stage 0 of the mechanism ablation scan.
%
% Flip every arm's switches at FROZEN parameters and record what happens to the
% objective. One evaluation per arm, ~38 s each, ~25 min for the whole list.
%
% WHY THIS RUNS FIRST. It is the cheapest possible filter and it answers three
% questions that would otherwise be answered after hours of optimisation:
%
%   delta == 0      the switch is dead code, or the mechanism is on but its
%                   parameter is 0 and therefore inert (plan.md trap T1). No
%                   re-optimisation can rescue it -- drop the arm, or find the
%                   parameter that was supposed to make it bite.
%   delta finite    the arm is real. Its size is also the "uncompensated" cost,
%                   which is the number stage 1 gets to buy back.
%   penalty / NaN   the flip broke the model. B4 is the known case: turning off
%                   the SRD pool with kmsr = 0 leaves SRX with no exit.
%
% It also checks the pool surgery BEFORE the search sees it: every `add` must be
% a field that is non-zero after the flip (or RunOptimSpreadW's own filter drops
% it and the mechanism is never drawn -- traps T1/T2), and every `drop` must
% actually be in the base pool.
%
% Usage:
%   matlab -batch "run('Analyses/MechanismAblation/RunAblationPreflight.m')"
%   matlab -batch "ARMS={'B4','C1'}; run('Analyses/MechanismAblation/RunAblationPreflight.m')"
%
% Writes: Analyses/MechanismAblation/results/preflight.mat
%
% See also: ablationSpec, plan.md

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('SEEDSNAP','var');  SEEDSNAP  = 'floored_seed'; end
if ~exist('SEED_EXPECT','var'); SEED_EXPECT = 0.9081;     end
if ~exist('NWORKERS','var');  NWORKERS  = 5;              end
if ~exist('ARMS','var');      ARMS      = {};             end   % {} = all

OUTDIR = fullfile(here, 'results');
if ~isfolder(OUTDIR); mkdir(OUTDIR); end

%% ---- base params0: defaults, then the snapshot pasted over --------------
% Exactly RunOptimSpreadW's `params0 = getParams(); run(SEEDSNAP)`, but through
% loadParams so a stale params0 cannot leak fields in from the base workspace.
base = getParams();
snap = loadParams(SEEDSNAP);
f = fieldnames(snap);
for i = 1:numel(f); base.(f{i}) = snap.(f{i}); end

% Run configuration -- identical to the campaign engine, or the costs are not
% comparable with 0.9081.
base.RunForceVelocity = true;  base.RunKtr = false;   base.RunSlack = true;
base.RunStairs = false;        base.RunForceVelocityTime = false;
base.RunSlackPassive = true;   base.EvalFeatures = true;
base.BreakOnODEUnstable = false;
base.PlotEachSeparately = 0;   base.PlotFeatureFitting = 0;
base.RunSlackSegments = 'AllPar';
base.MaxRunTime = 300;         base.OptimizeOn = 'Feats';
base.mods = {}; base.g = ones(1,0);

%% ---- the base pool, for checking each arm's `drop` entries --------------
% Mirrors the pool in Analyses/FeatureWeighting/RunOptimSpreadW.m.
basePool = iBasePool();

%% ---- pool ---------------------------------------------------------------
if isempty(gcp('nocreate'))
    try
        refreshPool(NWORKERS);
    catch ME
        warning('RunAblationPreflight:pool','refreshPool failed (%s) -- falling back', ME.message);
        if isempty(gcp('nocreate')); parpool('local', NWORKERS); end
    end
end

%% ---- baseline -----------------------------------------------------------
fprintf('\n=== ablation preflight: baseline ===\n');
t0 = tic;
[c0, ~, ct0, nm0] = evaluateBakersExp(1, base, true, false);
fprintf('base %s = %.4f   (expected %.4f, %.0f s)\n', SEEDSNAP, c0, SEED_EXPECT, toc(t0));
if ~isempty(SEED_EXPECT) && abs(c0 - SEED_EXPECT) > 0.25
    % 0.25 = the objective's own peak-to-peak roughness, not a tolerance
    error('RunAblationPreflight:seedDrift', ...
        'baseline %.4f is more than the noise floor from the expected %.4f', c0, SEED_EXPECT);
end

%% ---- arms ---------------------------------------------------------------
A = ablationSpec();
if ~isempty(ARMS)
    A = A(ismember({A.id}, ARMS) | ismember({A.name}, ARMS));
end
n = numel(A);
fprintf('\n=== %d arms, one frozen-parameter evaluation each ===\n\n', n);

R = struct('id',{},'name',{},'cost',{},'delta',{},'status',{},'poolNote',{}, ...
           'featCost',{},'featNames',{},'msg',{},'secs',{});

for k = 1:n
    a = A(k);
    p = base;
    for j = 1:2:numel(a.set); p.(a.set{j}) = a.set{j+1}; end

    % --- pool surgery check, before anything is spent on the arm ---------
    notes = {};
    for j = 1:numel(a.add)
        nmj = a.add{j};
        v = iResolve(p, nmj);
        if isnan(v)
            notes{end+1} = sprintf('add %s: ABSENT', nmj);                    %#ok<SAGROW>
        elseif v == 0
            notes{end+1} = sprintf('add %s: ZERO -> would be dropped', nmj);  %#ok<SAGROW>
        end
    end
    for j = 1:numel(a.drop)
        if ~ismember(a.drop{j}, basePool)
            notes{end+1} = sprintf('drop %s: not in base pool', a.drop{j});   %#ok<SAGROW>
        end
    end
    poolNote = strjoin(notes, '; ');

    % --- the evaluation --------------------------------------------------
    tA = tic; msg = '';
    try
        [c, ~, ct, nm] = evaluateBakersExp(1, p, true, false);
        if ~isfinite(c)
            status = 'NONFINITE';
        elseif c > 50
            % A big FINITE cost is not a broken arm: it is an uncompensated
            % one, which is exactly what stage 1 gets to buy back. Only
            % NO-OP / NONFINITE / ERROR disqualify an arm from the search.
            status = 'UNCOMP';
        elseif abs(c - c0) < 1e-9
            status = 'NO-OP';
        else
            status = 'live';
        end
    catch ME
        c = NaN; ct = []; nm = {}; status = 'ERROR'; msg = ME.message;
    end
    secs = toc(tA);

    R(k).id = a.id; R(k).name = a.name; R(k).cost = c; R(k).delta = c - c0;
    R(k).status = status; R(k).poolNote = poolNote;
    R(k).featCost = ct; R(k).featNames = nm; R(k).msg = msg; R(k).secs = secs;

    fprintf('%-4s %-24s %9s  d=%+8.4f  %-9s %3.0fs %s\n', a.id, a.name, ...
        iFmt(c), c - c0, status, secs, poolNote);
    if ~isempty(msg); fprintf('       ! %s\n', msg); end

    save(fullfile(OUTDIR,'preflight.mat'), 'R','c0','ct0','nm0','SEEDSNAP','A');
end

%% ---- report -------------------------------------------------------------
fprintf('\n=== ranked by |delta| at frozen parameters ===\n');
d = [R.delta]; d(~isfinite(d)) = Inf;
[~, ord] = sort(abs(d), 'descend');
fprintf('%-4s %-24s %9s %9s  %s\n','id','arm','cost','delta','status');
for i = ord
    fprintf('%-4s %-24s %9s %+9.4f  %s\n', R(i).id, R(i).name, iFmt(R(i).cost), R(i).delta, R(i).status);
end

nop = find(strcmp({R.status},'NO-OP'));
if ~isempty(nop)
    fprintf('\nNO-OP (dead switch or inert parameter -- drop, or find the missing value):\n');
    for i = nop; fprintf('  %-4s %s  |  %s\n', R(i).id, R(i).name, A(i).note); end
end
bad = find(ismember({R.status},{'ERROR','NONFINITE'}));
if ~isempty(bad)
    fprintf('\nBROKEN (needs a companion revival or a different seed value):\n');
    for i = bad; fprintf('  %-4s %-22s %s %s\n', R(i).id, R(i).name, R(i).status, R(i).msg); end
end

fprintf('\nwritten: %s\n', fullfile(OUTDIR,'preflight.mat'));

%% ---- helpers ------------------------------------------------------------
function s = iFmt(c)
if isfinite(c); s = sprintf('%.4f', c); else; s = 'NaN'; end
end

function v = iResolve(p, nm)
% Scalar value of a pool entry, including the fieldname__index form.
us = strfind(nm,'__');
if isempty(us)
    if isfield(p,nm) && isscalar(p.(nm)) && isnumeric(p.(nm)); v = p.(nm); else; v = NaN; end
else
    b = nm(1:us(1)-1); ix = str2double(nm(us(1)+2:end));
    if isfield(p,b) && numel(p.(b)) >= ix; v = p.(b)(ix); else; v = NaN; end
end
end

function pool = iBasePool()
pool = { ...
    'ka','kd','k1','k_1','k2','k_2','kah','kamh', ...
    'kstiff1','kstiff2','kstiff1_n','kstiff2_n','kSE','estiff','ekSE', ...
    'kSE_M','eta_M','mu_neg','mu','x_M_slack', ...
    'dr','dr1','dr2', ...
    'ksr0','kmsrd','ksr2srd','ksrd2sr','sigma1','sigma2','sigma_srd1','sigma_srd2', ...
    'A0','v_ref_reg','tau_reg','P_bound_max', ...
    'k_cr','tau_cr','F_cr','SRXFromR2HighStrain','sSRXrip', ...
    'kA2re','kA2hop','sA2hop', ...
    'L_thick','L_hbare','L_thin','LXBpivot', ...
    'd10_ref','SL_ref_lattice','R_thick','R_thin','d_optimal','sigma_lattice', ...
    'Lsc0','gamma','k_pas', ...
    'PieceWiseStrainDep2X_logOffset','PieceWiseStrainDepR1DX_logOffset', ...
    'PieceWiseStrainDepR21X_logOffset','PieceWiseStrainDepX_logOffset', ...
    'PieceWiseStrainDepX__2','PieceWiseStrainDepX__3','PieceWiseStrainDepX__4', ...
    'PieceWiseStrainDep2X__2','PieceWiseStrainDep2X__3','PieceWiseStrainDep2X__4','PieceWiseStrainDep2X__5', ...
    'PieceWiseStrainDepR1DX__1','PieceWiseStrainDepR1DX__2','PieceWiseStrainDepR1DX__3','PieceWiseStrainDepR1DX__4', ...
    'PieceWiseStrainDepR21X__2','PieceWiseStrainDepR21X__3','PieceWiseStrainDepR21X__4', ...
    'PieceWiseStrainDepParams__2','PieceWiseStrainDepParams__3','PieceWiseStrainDepParams__4', ...
    'PieceWiseStrainDep2Params__2','PieceWiseStrainDep2Params__3','PieceWiseStrainDep2Params__4', ...
    'PieceWiseStrainDep2Params__5','PieceWiseStrainDep2Params__6','PieceWiseStrainDep2Params__7', ...
    'PieceWiseStrainDepR1DParams__1','PieceWiseStrainDepR1DParams__2', ...
    'PieceWiseStrainDepR1DParams__3','PieceWiseStrainDepR1DParams__4', ...
    'PieceWiseStrainDepR21Params__2','PieceWiseStrainDepR21Params__3','PieceWiseStrainDepR21Params__4'};
end
