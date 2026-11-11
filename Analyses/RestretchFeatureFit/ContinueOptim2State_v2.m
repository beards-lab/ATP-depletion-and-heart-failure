% ContinueOptim2State_v2.m
% One-shot RESUME of the opt2state_v2 optimizer from its saved incumbent
% (params/opt2state_v2_state.mat, cost ~3.49) with the CORRECTED pool that now
% includes the R2 shortening-branch knots (FV-tail levers). Tests whether the
% optimizer can recover the FV tail (regressed to .53/.21 in the first run when
% ktr was driven down via kmsrd) while HOLDING ktr near data (~49).
% If cost keeps dropping and FV recovers -> ktr and FV coexist in 2-state.
% If FV stays collapsed -> a genuine ktr<->FV-tail residual weld -> motivates 3-state.
%
% Run:  cd(root); addpath(genpath('.')); ContinueOptim2State_v2
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));

cfg = struct();
cfg.baseSnap = 'params/params_2state_v2_seed.m';   % ignored on RESUME except for identity
cfg.tag      = 'opt2state_v2';
cfg.fn = {'FV_fnorm|FV_v|10', 'ktr|2', 'A|50', 'ktr_rmse|0-.5|.1', ...
    'XTOR[1]|3-15|15,XTOR_vmax[1]|5-25,SRX_ss[1]|0.01-0.40|5,attached_ss[1]|0.10-0.45|10,PT_ss[1]|0.01-0.70',...
    't0_crossing|SLdiff|2', 'restretchSlopeStart|1', 'peak1_y|10', 'peak1_dSL|1', 'vall_y|10', ...
    'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|1.0', 'ovrsht_dy|1', ...
    'k2|100-1500|0.01,ksrd|0.-20|0.001,kmsrd|0.0-20|0.001'};
cfg.pool = {'k2', 'A0', 'ka', 'kstiff2', 'v_ref_reg', 'kstiff1', 'kSE', ...
    'k_2', 'k1', 'PieceWiseStrainDepR1DParams__4', 'PieceWiseStrainDepR1DParams__3', ...
    'ksr0', 'dr', 'tau_reg', 'kmsrd', 'ksr2srd', 'ksrd2sr', ...
    'PieceWiseStrainDep2Params__5', 'PieceWiseStrainDep2Params__6', 'PieceWiseStrainDep2Params__7'};
cfg.compulsory = {'k2', 'A0', 'ka', 'v_ref_reg', 'kmsrd', 'PieceWiseStrainDep2Params__6'};  % force FV-tail knot in every draw
cfg.SURR_EVALS = 0;
cfg.KICK_FRAC  = 0;
cfg.DEBUG           = false;
cfg.RESUME          = true;    % continue from the 3.49 incumbent
cfg.TIME_BUDGET_HRS = 1.0;

optimizeFeatures(cfg);
