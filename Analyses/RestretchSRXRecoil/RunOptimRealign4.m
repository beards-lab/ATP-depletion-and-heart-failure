% RunOptimRealign4.m — large-scale refit on the POWER objective.
%
% Objective change from realign3, both cost-neutral at the seed:
%
%   1. FV_normpowerAvg added, weighted per point by the inverse variance between
%      the two FV source datasets ('w:FV_normpowerVar'). Power is what the ATP
%      effect actually shows up in, and unlike FV_fnorm it is not anchored by an
%      isometric point that is 1 by construction.
%
%   2. FV_fnorm REMOVED and its cost folded into the power term (weight 20 ->
%      62.2779). Power = fnorm x v, so scoring both double-counts the same four
%      residuals. At the seed the power term now carries 2.3107, exactly what
%      FV_fnorm (1.5686) + power@20 (0.7421) carried before.
%
%   3. Passive weights x0.4466 (done in realign3's successor objective), moving
%      passive from 26.4% to 11.8% of the objective. Deliberately NOT zero:
%      eta_M / kSE_M / mu_neg are shared with the re-stretch transient, so with
%      no passive term the dashpot would be unconstrained.
%
% Net: seed total is UNCHANGED at 5.0819, so any further drop is real progress
% rather than a change of yardstick. (Neutrality is exact only at the seed; the
% groups scale differently as the optimiser moves. That is intended.)
%
% Seed: params/realign3_opt.m — already carries UseMaxwellTensionOnly = 1,
% UseLoadRealign = 1 (k_cr = 1.590), x_M_slack = 0.0103.
%
% LARGE-SCALE pool: the comprehensive candidate set (as RunOptimRsK) plus the
% realignment mechanism, the SRX routing, the Maxwell slack length and the FV
% shoulder. Roughly 70 parameters; block-coordinate draws of 10.
%
% Usage:
%   matlab -batch "TIME_HRS=12; run('Analyses/RestretchSRXRecoil/RunOptimRealign4.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('TIME_HRS','var'); TIME_HRS = 12; end
if ~exist('NWORKERS','var'); NWORKERS = 4;  end   % 4, not 5 — realign3 lost a worker at 5
if ~exist('TAG',     'var'); TAG      = 'realign4'; end
if ~exist('W_POW',   'var'); W_POW    = 62.2779; end
if ~exist('PS_SCALE','var'); PS_SCALE = 0.446621; end
if ~exist('DROP',    'var'); DROP     = {'ovrsht_dy', 'FV_fnorm'}; end

params0 = getParams(); run('params/realign3_opt.m');

%% ---- build the objective -------------------------------------------------
fn   = params0.fn;
prim = cellfun(@(s) regexprep(strtok(s,'|'), '\[.*$',''), fn, 'uni', 0);
drop = ismember(prim, DROP);
fprintf('dropping %d term(s): %s\n', nnz(drop), strjoin(fn(drop), ' , '));
fn = fn(~drop); prim = prim(~drop);

nPS = 0;
for i = find(startsWith(prim, 'PS_'))
    tk = strsplit(fn{i}, '|');
    if numel(tk) == 1; wOld = 1; tk{end+1} = ''; else; wOld = str2double(tk{end}); end %#ok<AGROW>
    tk{end} = num2str(wOld*PS_SCALE, '%.6g');
    fn{i} = strjoin(tk, '|'); nPS = nPS + 1;
end
fprintf('scaled %d passive weights by %.6g\n', nPS, PS_SCALE);

if any(strcmp(prim,'FV_normpowerAvg'))
    error('RunOptimRealign4:dupPower','FV_normpowerAvg already present in fn.');
end
fn{end+1} = sprintf('FV_normpowerAvg|w:FV_normpowerVar|%.6g', W_POW);
fprintf('added: %s\n', fn{end});

params0.fn = fn; params0.mods = {}; params0.g = [];

SEEDSNAP = fullfile('params', ['realign_seed_' TAG '.m']);
writeParamsToMFile(SEEDSNAP, params0, [], ...
    sprintf('seed for RunOptimRealign4: realign3_opt + power objective (W=%g, PS x%g)', W_POW, PS_SCALE));

%% ---- large-scale pool ----------------------------------------------------
cfg = struct();
cfg.baseSnap = SEEDSNAP; cfg.tag = TAG; cfg.fn = fn;
cfg.pool = { ...
    ... % core cycle
    'ka','kd','k1','k_1','k2','k_2','kah','kamh', ...
    ... % stiffness / mechanics
    'kstiff1','kstiff2','kstiff1_n','kstiff2_n','kSE','estiff','ekSE', ...
    'kSE_M','eta_M','mu_neg','mu','x_M_slack', ...
    'dr','dr1','dr2', ...
    ... % SRX
    'ksr0','kmsrd','ksr2srd','ksrd2sr','sigma1','sigma2','sigma_srd1','sigma_srd2', ...
    ... % attachment gating / FV shoulder
    'A0','v_ref_reg','tau_reg','P_bound_max', ...
    ... % the realignment mechanism + the SRX routing
    'k_cr','tau_cr','F_cr','SRXFromR2HighStrain','sSRXrip', ...
    ... % restretch transient
    'kA2re','kA2hop','sA2hop', ...
    ... % geometry / lattice
    'L_thick','L_hbare','L_thin','LXBpivot', ...
    'd10_ref','SL_ref_lattice','R_thick','R_thin','d_optimal','sigma_lattice', ...
    'Lsc0','gamma','k_pas', ...
    ... % strain-dependence shape (this is what sets FV)
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
% k2/kstiff2/ka are the standard trio; v_ref_reg keeps FV tradeable every round;
% k_cr keeps the mechanism itself always on the table (it can be switched off).
cfg.compulsory = {'k2','kstiff2','ka','v_ref_reg','k_cr'};

keep = true(1, numel(cfg.pool));
for i = 1:numel(cfg.pool)
    f = cfg.pool{i}; us = strfind(f,'__');
    if isempty(us)
        keep(i) = isfield(params0,f) && isscalar(params0.(f)) && params0.(f) ~= 0;
    else
        b = f(1:us(1)-1); ix = str2double(f(us(1)+2:end));
        keep(i) = isfield(params0,b) && numel(params0.(b)) >= ix && params0.(b)(ix) ~= 0;
    end
end
if any(~keep); fprintf('dropping %d zero-base params: %s\n', nnz(~keep), strjoin(cfg.pool(~keep),', ')); end
cfg.pool = cfg.pool(keep);

cfg.N_DRAW = 10; cfg.SIMPLEX_EVALS = 60; cfg.MaxRunTime = 300;
cfg.SURR_EVALS = 0; cfg.KICK_FRAC = 0.1; cfg.DEBUG = false;
cfg.TIME_BUDGET_HRS = TIME_HRS;
cfg.RESUME = isfile(fullfile(root,'params',[cfg.tag '_state.mat']));

params0.RunForceVelocity = true; params0.RunKtr = false; params0.RunSlack = true;
params0.RunStairs = false; params0.RunForceVelocityTime = false;
params0.RunSlackPassive = true;
params0.EvalFeatures = true; params0.BreakOnODEUnstable = false;
params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0;
params0.RunSlackSegments = 'AllPar';
params0.MaxRunTime = cfg.MaxRunTime; params0.OptimizeOn = 'Feats';
cfg.params0 = params0;

if isempty(gcp('nocreate')); parpool('local', NWORKERS); end

pchk = cfg.params0; pchk.mods = {}; pchk.g = ones(1,0);
rep = checkFeatureCoverage(pchk);
if ~rep.ok
    fprintf(2, 'UNSCORABLE fn entries — each costs a flat MISSING_FEATURE_COST=100:\n');
    if ~isempty(rep.missingModel); fprintf(2,'  model: %s\n', strjoin(rep.missingModel,', ')); end
    if ~isempty(rep.missingData);  fprintf(2,'  data : %s\n', strjoin(rep.missingData ,', ')); end
    error('RunOptimRealign4:coverage','objective has unscorable terms');
end
fprintf('[%s] coverage OK (%d fn entries)\n', cfg.tag, numel(fn));

[c0,~,ct0,nm0] = evaluateBakersExp(1, pchk, true, false);
iP = find(contains(nm0,'FV_normpowerAvg'),1);
fprintf('[%s] SEED = %.4f   (realign3 scored 5.0819 on the OLD objective)\n', cfg.tag, c0);
fprintf('[%s] FV power term = %.4f\n', cfg.tag, ct0(iP));
if ~isfinite(c0) || c0 > 1e4; error('RunOptimRealign4:badSeed','seed %.4g looks like a penalty', c0); end
fprintf('[%s] pool=%d | N_DRAW=%d | budget %.1f h | resume=%d\n', ...
        cfg.tag, numel(cfg.pool), cfg.N_DRAW, cfg.TIME_BUDGET_HRS, cfg.RESUME);

optimizeFeatures(cfg);
