% RunOptimSpreadW.m — large-scale refit on the SPREAD-WEIGHTED objective.
%
% Objective change from realign4 — cost-neutral at the seed (3.8513), so any
% drop is real progress rather than a change of yardstick:
%
%   Every data-fit weight is now  w = lambda / s^2, where
%     s   = between-preparation relative SD of that feature at 8 mM, measured
%           on the THREE 2026 protocol days (Baker excluded: different
%           protocol, ATP-dependently truncated amplitudes).
%     lambda = one global constant fixing cost-neutrality.
%   Because evalFeatureCost sums SQUARED RELATIVE residuals, 1/s^2 makes each
%   term a chi-square: 1.0 = "one preparation-width away from the data".
%   ATP relevance is measured and PLOTTED but does not enter the weights — see
%   ../FeatureWeighting/conclusions.md Sec 4.
%
%   Two shares are pinned BY HAND, not by the data:
%     FV_SHARE  = 0.70  *** ARBITRARY ***  FV is a separate experiment and the
%                 whole muscle work profile is in it. Its measured spread would
%                 put it at 1.4 % (chi2 1.67 over 7 pts — already inside the
%                 disagreement between its own two source curves).
%     SHARE_CAP = 0.30  on every other term (inert at FV_SHARE >= ~0.55).
%
%   Built by ../FeatureWeighting/BuildWeightedObjective.m from
%   ../FeatureWeighting/RunFeatureSpread.m. Seed snapshot: params/spreadw_seed.m
%   (it already CARRIES the reweighted fn — no fn surgery happens here).
%
% What moved, and why it matters for the search:
%   FV       34.8 % -> 70.0 %  (pinned)
%   rsK       9.9 % -> 12.5 %  (chi2 33.4 — the one genuinely broken feature)
%   peak1_dSL 8.5 % ->  0.2 %  (76 % between-prep spread — barely measurable)
%   t0_cross  6.8 % ->  0.3 %  (38 % spread)
%   ktr       3.7 % ->  0.4 %  (chi2 1.11 — already inside one prep-width)
%   ...the fourteen non-FV data terms hold 6.6 % between them, so the slack
%   SHAPE is largely along for the ride at this FV_SHARE. Watch it: the seed
%   composition printed below is the yardstick, and the run's own featCost
%   history (params/spreadw_state.mat -> state.optIter) records what was paid.
%
% NOTE on the absolute force level: extractForceVelocityAttributes normalises
%   the model's power by the MODEL's OWN isometric force, so 70 % of this
%   objective is scale-blind. The weight backing absolute force falls 13x
%   (125 -> 9.4). A uniform +20 % drift still costs ~1.9 (half the objective),
%   but if that is not enough, set ANCHOR_GUARD = true in
%   BuildWeightedObjective.m and rebuild the seed.
%
% Usage:
%   matlab -batch "TIME_HRS=12; run('Analyses/FeatureWeighting/RunOptimSpreadW.m')"
% Resume is automatic if params/<TAG>_state.mat exists.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('TIME_HRS','var'); TIME_HRS = 12;             end
if ~exist('NWORKERS','var'); NWORKERS = 4;              end   % 5 lost a worker on realign3
if ~exist('TAG',     'var'); TAG      = 'spreadw';      end
if ~exist('SEEDSNAP','var'); SEEDSNAP = fullfile('params','spreadw_seed.m'); end
if ~exist('SEED_EXPECT','var'); SEED_EXPECT = 3.8513;   end   % [] to skip the check

if ~isfile(SEEDSNAP)
    error('RunOptimSpreadW:noSeed', ['%s not found — run ' ...
        'Analyses/FeatureWeighting/RunFeatureSpread.m then BuildWeightedObjective.m first.'], SEEDSNAP);
end

params0 = getParams(); run(SEEDSNAP);
fn = params0.fn;
params0.mods = {}; params0.g = [];

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

% k2/kstiff2/ka are the standard trio; v_ref_reg keeps FV tradeable every round
% (it is now only 5 % of the objective, so it must stay explicitly drawable
% rather than rely on weight to defend it); k_cr keeps the realignment
% mechanism on the table; eta_M/kSE_M are the passive<->restretch shared
% levers, and passive is now 24 % of the objective.
cfg.compulsory = {'k2','kstiff2','ka','v_ref_reg','k_cr','eta_M'};

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
if any(~keep)
    fprintf('dropping %d zero-base params: %s\n', nnz(~keep), strjoin(cfg.pool(~keep),', '));
end
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

%% ---- pre-flight ----------------------------------------------------------
pchk = cfg.params0; pchk.mods = {}; pchk.g = ones(1,0);
rep = checkFeatureCoverage(pchk);
if ~rep.ok
    fprintf(2, 'UNSCORABLE fn entries — each costs a flat MISSING_FEATURE_COST=100:\n');
    if ~isempty(rep.missingModel); fprintf(2,'  model: %s\n', strjoin(rep.missingModel,', ')); end
    if ~isempty(rep.missingData);  fprintf(2,'  data : %s\n', strjoin(rep.missingData ,', ')); end
    error('RunOptimSpreadW:coverage','objective has unscorable terms');
end
fprintf('[%s] coverage OK (%d fn entries)\n', cfg.tag, numel(fn));

[c0,~,ct0,nm0] = evaluateBakersExp(1, pchk, true, false);
fprintf('[%s] SEED = %.4f   (realign4_opt scored 3.8513 on the OLD objective)\n', cfg.tag, c0);
if ~isfinite(c0) || c0 > 1e4
    error('RunOptimSpreadW:badSeed','seed %.4g looks like a penalty', c0);
end
if ~isempty(SEED_EXPECT) && abs(c0 - SEED_EXPECT) > 0.25
    % 0.25 = the objective's own peak-to-peak numerical roughness
    % (Docs objective_noise_floor); a bigger gap means the seed is not the
    % cost-neutral one BuildWeightedObjective wrote.
    error('RunOptimSpreadW:seedDrift', ...
          'seed %.4f differs from the expected %.4f by more than the noise floor', c0, SEED_EXPECT);
end
[~, ord] = sort(ct0,'descend');
fprintf('[%s] seed composition:\n', cfg.tag);
for i = ord(1:min(10,numel(ord)))
    fprintf('        %-22s %7.4f  (%4.1f %%)\n', nm0{i}, ct0(i), 100*ct0(i)/c0);
end
fprintf('[%s] pool=%d | N_DRAW=%d | budget %.1f h | resume=%d\n', ...
        cfg.tag, numel(cfg.pool), cfg.N_DRAW, cfg.TIME_BUDGET_HRS, cfg.RESUME);

% DRYRUN=true validates seed / coverage / pool and stops before the search.
if exist('DRYRUN','var') && DRYRUN
    fprintf('[%s] DRYRUN — pre-flight passed, not starting the search.\n', cfg.tag);
    return;
end

optimizeFeatures(cfg);
