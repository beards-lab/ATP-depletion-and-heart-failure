% RunOptim3State.m
% =========================================================================
% Optimize the 3-state model (Design B: strain-dependent post-stroke p3) from
% the non-collapsing seed params/params_3state_seedB.m (hand-calibrated to
% cost ~4.60: NoS=3, kstiff3=12000, k3=2000, alpha3=120, drp3=0.02, s3=0).
%
% Question this answers: can a 3rd state beat the 2-state floor (~3.22), or does
% the optimizer just drive p3 toward transience (k3 large => 2-state-like)?
%
% Pool = the v3 2-state pool + the p3 levers {k3, kstiff3, alpha3, drp3, s3}.
% (alpha3/s3 were dead until this campaign wired them into R3's strain-dependent
% detachment: R3 = k3*p3*max(1,exp(-alpha3*(s-s3))).) compulsory forces the p3
% force/turnover balance (k3, kstiff3, alpha3) into every draw so p3 is always
% exercised. R2 knots PieceWiseStrainDep2Params__5/6/7 shape the p2->p3 feed.
%
% NOTE 3-state is ~3-4x slower per eval than 2-state and some param regions time
% out (evaluateModel THROWS on MaxRunTime / instability). optimizeFeatures'
% safeCost / safeCostClamp / final verify are now try/catch-guarded, so a timeout
% is a high finite cost the optimizer routes around, not a crash. SIMPLEX_EVALS
% halved (30) for throughput; MaxRunTime 75 s to cut spurious timeouts.
%
% Isolation: cfg.tag='opt3state' -> params/opt3state_opt.m. Never touches
% opt2state*.  Run:  cd(root); addpath(genpath('.')); RunOptim3State
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));

cfg = struct();
cfg.baseSnap = 'params/params_3state_seedB.m';
cfg.tag      = 'opt3state';

cfg.fn = {'FV_fnorm|FV_v|10', 'ktr|2', 'A|50', 'ktr_rmse|0-.5|.1', ...
    'XTOR[1]|3-15|15,XTOR_vmax[1]|5-25,SRX_ss[1]|0.01-0.40|5,attached_ss[1]|0.10-0.45|10,PT_ss[1]|0.01-0.70',...
    't0_crossing|SLdiff|2', 'restretchSlopeStart|1', 'peak1_y|10', 'peak1_dSL|1', 'vall_y|10', ...
    'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|1.0', 'ovrsht_dy|1', ...
    'k2|100-1500|0.01,ksrd|0.-20|0.001,kmsrd|0.0-20|0.001'};

cfg.pool = {'k2', 'A0', 'ka', 'kstiff2', 'v_ref_reg', 'kstiff1', 'kSE', ...
    'k_2', 'k1', 'PieceWiseStrainDepR1DParams__4', 'PieceWiseStrainDepR1DParams__3', ...
    'ksr0', 'dr', 'tau_reg', 'kmsrd', 'ksr2srd', 'ksrd2sr', ...
    'PieceWiseStrainDep2Params__5', 'PieceWiseStrainDep2Params__6', 'PieceWiseStrainDep2Params__7', ...
    'kstiff2_n', ...
    'k3', 'kstiff3', 'alpha3', 'drp3', 's3'};                 % <<< p3 levers (Design B)
cfg.compulsory = {'k2', 'ka', 'kstiff2', 'k3', 'kstiff3', 'alpha3'};  % p3 force/turnover always drawn

cfg.SURR_EVALS    = 0;
cfg.KICK_FRAC     = 0;
cfg.SIMPLEX_EVALS = 30;     % halved for 3-state throughput
cfg.MaxRunTime    = 75;     % fewer spurious timeouts than 60
cfg.DEBUG           = false;
cfg.TIME_BUDGET_HRS = 1.5;

% Auto-resume: 3-state runs are long and background launches can be killed at
% session teardown. If a prior (interrupted) run left a state file, CONTINUE
% from its incumbent instead of restarting from the seed. Delete
% params/opt3state_state.mat to force a fresh start.
cfg.RESUME = isfile(fullfile(root, 'params', 'opt3state_state.mat'));

optimizeFeatures(cfg);
