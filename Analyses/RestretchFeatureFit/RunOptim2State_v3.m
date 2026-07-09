% RunOptim2State_v3.m
% =========================================================================
% Establish the 2-state FLOOR by optimizing ALL identified levers together.
% Seed = params/params_2state_v3_seed.m (= opt2state_v2_opt.m, cost 3.22).
% Pool = v2 pool (regavail + freed SRX-return kmsrd/ksr2srd/ksrd2sr + R2
% shortening knots) PLUS kstiff2_n (negative-strain stiffness).
%
% WHY kstiff2_n: probePeakLevers showed it is a genuine FV-mid-tail lever
% (compressed p2 heads over-contribute resistive force; lowering kstiff2_n
% lifts the tail) but it is COUPLED to ktr (lowering it re-inflates ktr) --
% the same ktr<->FV-tail frontier as kmsrd. Freeing both lets the optimizer
% pick the best joint point on that frontier. (kstiff1_n was inert; catch
% bond couples valley<->peak the wrong way -- both excluded.)
%
% If this run cannot beat ~3.2, the 2-state ceiling is proven and the 3rd
% state is justified (two residual clusters -- restretch valley/peak and
% FV-tail/ktr -- both need an independent population/timescale).
%
% Isolation: cfg.tag='opt2state_v3' -> params/opt2state_v3_opt.m.
% Run:  cd(root); addpath(genpath('.')); RunOptim2State_v3
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));

cfg = struct();
cfg.baseSnap = 'params/params_2state_v3_seed.m';
cfg.tag      = 'opt2state_v3';

cfg.fn = {'FV_fnorm|FV_v|10', 'ktr|2', 'A|50', 'ktr_rmse|0-.5|.1', ...
    'XTOR[1]|3-15|15,XTOR_vmax[1]|5-25,SRX_ss[1]|0.01-0.40|5,attached_ss[1]|0.10-0.45|10,PT_ss[1]|0.01-0.70',...
    't0_crossing|SLdiff|2', 'restretchSlopeStart|1', 'peak1_y|10', 'peak1_dSL|1', 'vall_y|10', ...
    'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|1.0', 'ovrsht_dy|1', ...
    'k2|100-1500|0.01,ksrd|0.-20|0.001,kmsrd|0.0-20|0.001'};

cfg.pool = {'k2', 'A0', 'ka', 'kstiff2', 'v_ref_reg', 'kstiff1', 'kSE', ...
    'k_2', 'k1', 'PieceWiseStrainDepR1DParams__4', 'PieceWiseStrainDepR1DParams__3', ...
    'ksr0', 'dr', 'tau_reg', 'kmsrd', 'ksr2srd', 'ksrd2sr', ...
    'PieceWiseStrainDep2Params__5', 'PieceWiseStrainDep2Params__6', 'PieceWiseStrainDep2Params__7', ...
    'kstiff2_n'};                                              % <<< negative-strain stiffness (FV-tail frontier)
cfg.compulsory = {'k2', 'A0', 'ka', 'v_ref_reg', 'kmsrd', 'kstiff2_n'};  % both ktr<->FV-tail knobs always drawn

cfg.SURR_EVALS = 0;
cfg.KICK_FRAC  = 0;
cfg.DEBUG           = false;
cfg.TIME_BUDGET_HRS = 0.7;

optimizeFeatures(cfg);
