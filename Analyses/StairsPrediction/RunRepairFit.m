% RunRepairFit.m — "is estiff > 1 repairable by the other parameters?"
%
% estiff is FROZEN at 1.25 (the mildest genuine flip past 1, iso-force
% compensated: kstiff2 x3.662) and every other parameter in the spreadw pool is
% free, including kSE. The objective is the one the incumbent was actually
% bought with (params/spreadw_opt.m fn); the staircase is held OUT, so whatever
% the search reaches is scored on the staircase afterwards as a prediction.
%
%   incumbent  e=0.965 : Et = 0.979,  staircase feature cost 5.42
%   seed       e=1.250 : Et = 5.480,  staircase feature cost 18.30
%
% The question this run answers: how far back toward ~1 can Et be brought with
% estiff pinned at 1.25, and does the staircase follow or stay broken.
%
% Usage:  matlab -batch "TIME_HRS=1.5; run('Analyses/StairsPrediction/RunRepairFit.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));

if ~exist('TIME_HRS','var'); TIME_HRS = 1.5;                 end
if ~exist('NWORKERS','var'); NWORKERS = 3;                   end  % share the box
if ~exist('TAG',     'var'); TAG      = 'estiff125';         end
if ~exist('SEEDSNAP','var'); SEEDSNAP = fullfile('params','estiff125_seed.m'); end


if ~isfile(SEEDSNAP); error('RunRepairFit:noSeed','%s not found', SEEDSNAP); end

params0 = getParams(); run(SEEDSNAP);
fn = params0.fn;
params0.mods = {}; params0.g = [];

cfg = struct();
cfg.baseSnap = SEEDSNAP; cfg.tag = TAG; cfg.fn = fn;

% Everything the spreadw run was allowed to move, MINUS estiff (the treatment).
% kSE and ekSE stay in: the 2-D grid showed series compliance is the parameter
% that trades against the exponent, so excluding it would rig the answer.
cfg.pool = { ...
    'ka','kd','k1','k_1','k2','k_2','kah','kamh', ...
    'kstiff1','kstiff2','kstiff1_n','kstiff2_n','kSE','ekSE', ...
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

% FV took the worst damage from the flip (0.18 -> 3.40), so its levers are drawn
% every round; kSE is compulsory because it is the direct counterweight to estiff.
cfg.compulsory = {'k2','kstiff2','ka','v_ref_reg','kSE'};

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
cfg.pool = cfg.pool(keep);
assert(~any(strcmp(cfg.pool,'estiff')), 'estiff must stay frozen');

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
    fprintf(2,'UNSCORABLE fn entries:\n');
    if ~isempty(rep.missingModel); fprintf(2,'  model: %s\n', strjoin(rep.missingModel,', ')); end
    if ~isempty(rep.missingData);  fprintf(2,'  data : %s\n', strjoin(rep.missingData ,', ')); end
    error('RunRepairFit:coverage','objective has unscorable terms');
end
[c0,~,ct0,nm0] = evaluateBakersExp(1, pchk, true, false);
fprintf('[%s] estiff = %.4g (frozen) | SEED Et = %.4f  (incumbent at e=0.965 was 0.979)\n', ...
        cfg.tag, pchk.estiff, c0);
[~,ord] = sort(ct0,'descend');
for i = ord(1:min(8,numel(ord)))
    fprintf('        %-22s %7.4f (%4.1f %%)\n', nm0{i}, ct0(i), 100*ct0(i)/c0);
end
fprintf('[%s] pool=%d  budget %.2f h  resume=%d\n', cfg.tag, numel(cfg.pool), TIME_HRS, cfg.RESUME);

optimizeFeatures(cfg);

%% ---- score the result on the staircase, out of sample -------------------
optSnap = fullfile('params',[TAG '_opt.m']);
if isfile(optSnap)
    % MATLAB caches parsed script files; the optimiser just rewrote this one.
    clear functions; rehash; %#ok<CLFUNC>
    pOpt = getParams(loadParams([TAG '_opt']), [], true, false);
    d  = load('data/protocol_03_27_2026_velocitytable_stairs.mat');
    fD = extractStairsAttributes(d.datatable(:,1), d.datatable(:,3), d.velocitytable, 'st_', struct(), true);
    ps = pOpt; ps.RunStairs = 1; ps.PlotEachSeparately = false; ps.ghostSave=''; ps.ghostLoad='';
    [Elsq, os] = runStairsExperiment(ps);
    tMov = d.velocitytable(find(d.velocitytable(:,2)~=0,1),1);
    Fm = interp1(os.t, os.Force, tMov,'linear','extrap');
    [Cf, parts] = stairsCost(os.t, os.Force/Fm, d.velocitytable, fD);
    [EtOpt,~,~] = evaluateBakersExp(pOpt.g, pOpt, true, false);
    fprintf('\n=== RESULT (%s, estiff frozen at %.4g) ===\n', TAG, pOpt.estiff);
    fprintf('  objective Et  : seed %.4f -> %.4f   (incumbent 0.975)\n', c0, EtOpt);
    fprintf('  staircase feat: %.3f              (incumbent 5.415)\n', Cf);
    fprintf('  staircase LSQ : %.5f  [peak %.2f nadir %.2f plat %.3f rec %.2f]\n', ...
            Elsq, parts.st_peak, parts.st_nadir, parts.st_plat, parts.st_rec);
    fprintf('  kSE %.5g (seed %.5g)   kstiff2 %.5g (seed %.5g)\n', ...
            pOpt.kSE, pchk.kSE, pOpt.kstiff2, pchk.kstiff2);
    save(fullfile(here,'results',[TAG '_fit.mat']), 'EtOpt','Cf','Elsq','parts','c0','pOpt');
end
