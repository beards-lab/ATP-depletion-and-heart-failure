% RunJointRefit.m -- a REAL refit: core model parameters AND the ATP gates,
% scored jointly on 8 mM (FV + slack + passive) and 2 mM (slack).
%
% Why joint, and why core parameters. The previous round refitted only the 21
% ATP-gate exponents. Those are exactly 1 at MgATP_ref by construction, so they
% cannot touch 8 mM -- and the figures showed most of the residual is a defect
% the two conditions SHARE (rsK is 53 % off at 8 mM, FV is -20 % mid-range at
% 8 mM). A gate-only fit is structurally incapable of addressing that; only the
% core parameters can, and they must be scored at both [ATP] so the fix is not
% bought by wrecking the other condition.
%
% Objective:  Et = w8*E(8 mM, full battery) + w2*E(2 mM, slack)
% via params0.UseAtpJointFit, which recurses through evaluateBakersExp, so every
% existing path (featCost decomposition, iter gallery, resume) is reused.
%
% Seed: params/jointrefit_seed.m, written and VALIDATED by SetupJointRefit.m
% (both legs fully scorable -- no MISSING_FEATURE_COST terms -- and 8 mM proven
% invariant to the gates to 0.0e+00).
%
% Usage:
%   matlab -batch "TIME_HRS=6; run('Analyses/LowATP_Realign4Transfer/RunJointRefit.m')"

clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;

if ~exist('TIME_HRS','var'); TIME_HRS = 6;         end
if ~exist('NWORKERS','var'); NWORKERS = 5;         end
if ~exist('TAG','var');      TAG      = 'joint82'; end

params0 = getParams(); run(fullfile(root,'params','jointrefit_seed.m'));
params0 = getParams(params0, [], true, false);
params0.mods = {}; params0.g = [];

%% ---- pool ---------------------------------------------------------------
% realign4's core pool (it produced the current 8 mM best) PLUS the ATP gates
% that are actually carrying an effect. The gates are multiplier-tuned like
% everything else; only the 10 non-zero ones can be, since g*0 = 0.
CORE = { ...
    'ka','kd','k1','k_1','k2','k_2','kah','kamh', ...
    'kstiff1','kstiff2','kstiff1_n','kstiff2_n','kSE','estiff','ekSE', ...
    'kSE_M','eta_M','mu_neg','mu','x_M_slack','dr','dr1','dr2', ...
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

GATEP = {'pAtp_k2','pAtp_kmsrd','pAtp_kah','pAtp_kamh','pAtp_ksr0', ...
         'pAtp_ksr2srd','pAtp_ksrd2sr','pAtp_eta_M','pAtp_kSE','pAtp_mu_neg'};

cfg = struct();
cfg.tag = TAG; cfg.fn = params0.fn;
cfg.baseSnap = fullfile(root,'params','jointrefit_seed.m');
cfg.pool = [CORE, GATEP];

% Always on the table: the standard trio, the FV lever, the realignment gain
% that owns rsK, and the strongest ATP lever -- so every round can trade the
% shared defect against the ATP-specific one rather than only one of them.
cfg.compulsory = {'k2','kstiff2','ka','v_ref_reg','k_cr','pAtp_kmsrd'};

% drop zero-base params: a multiplier cannot move them off zero
keep = true(1,numel(cfg.pool));
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

cfg.N_DRAW = 12; cfg.SIMPLEX_EVALS = 50; cfg.SURR_EVALS = 0;
cfg.MaxRunTime = 300; cfg.KICK_FRAC = 0.1; cfg.DEBUG = false;
cfg.TIME_BUDGET_HRS = TIME_HRS;
cfg.RESUME = isfile(fullfile(root,'params',[cfg.tag '_state.mat']));

params0.MaxRunTime = cfg.MaxRunTime;
cfg.params0 = params0;

if isempty(gcp('nocreate')); parpool('local', NWORKERS); end

%% ---- validate before spending hours -------------------------------------
pchk = params0; pchk.mods = {}; pchk.g = ones(1,0);
p8 = pchk; p8.UseAtpJointFit = false;
r8 = checkFeatureCoverage(p8);
p2 = p8; p2.MgATP = params0.AtpJoint_MgATP2; p2.Pi = params0.AtpJoint_Pi2;
p2.velocitytableonfile = params0.AtpJoint_vt2; p2.fn = params0.AtpJoint_fn2;
p2.RunForceVelocity = false; p2.RunForceVelocityTime = false; p2.RunSlackPassive = false;
r2 = checkFeatureCoverage(p2);
if ~r8.ok || ~r2.ok
    fprintf(2,'UNSCORABLE fn entries (flat cost 100 each, every eval):\n');
    fprintf(2,'  8 mM: %s\n', strjoin(unique([r8.missingModel(:);r8.missingData(:)])',', '));
    fprintf(2,'  2 mM: %s\n', strjoin(unique([r2.missingModel(:);r2.missingData(:)])',', '));
    error('RunJointRefit:coverage','objective has unscorable terms');
end
fprintf('[%s] coverage OK  (8 mM %d feats | 2 mM %d feats)\n', TAG, numel(p8.fn), numel(p2.fn));

[c0,~,ct0,nm0] = evaluateBakersExp(1, pchk, true, false);
s8 = sum(ct0(startsWith(nm0,'hi:'))); s2 = sum(ct0(startsWith(nm0,'lo:')));
fprintf('[%s] SEED = %.4f   [ 8 mM %.4f | 2 mM %.4f ]\n', TAG, c0, s8, s2);
if ~isfinite(c0) || c0 > 1e4
    error('RunJointRefit:badSeed','seed %.4g looks like a penalty', c0);
end
fprintf('[%s] pool=%d (%d gates) | N_DRAW=%d | budget %.1f h | resume=%d\n', ...
        TAG, numel(cfg.pool), sum(startsWith(cfg.pool,'pAtp_')), cfg.N_DRAW, ...
        cfg.TIME_BUDGET_HRS, cfg.RESUME);

optimizeFeatures(cfg);
fprintf('\nJOINT REFIT DONE\n');
