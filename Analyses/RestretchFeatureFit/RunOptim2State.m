% RunOptim2State.m
% =========================================================================
% Thin driver for the 2-state regavail feature fit. Seeds from
% params/params_2state_seed.m (a bare snapshot with NO OptimizeParams /
% OptimizeCompulsory fields -- the tunable pool is therefore defined here,
% not read off the snapshot) and hands everything to optimizeFeatures.m.
%
% Isolation: cfg.tag='opt2state' routes all output to
% params/opt2state_opt.m and params/opt2state_state.mat. A second instance
% run concurrently with a DIFFERENT tag (e.g. RunOptim3State.m's
% 'opt3state') never touches these files -- see optimizeFeatures.m header
% for the accept/save invariant that keeps the incumbent snapshot correct.
%
% Run:  cd(root); addpath(genpath('.')); RunOptim2State
% =========================================================================

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));

cfg = struct();
cfg.baseSnap = 'params/params_2state_seed.m';
cfg.tag      = 'opt2state';

% --- feature list: hardcoded per spec (verbatim from params_restretch_best) ---
cfg.fn = {'FV_fnorm|FV_v|10', 'ktr|2', 'A|50', 'ktr_rmse|0-.5|.1', ...
    'XTOR[1]|3-15|15,XTOR_vmax[1]|5-25,SRX_ss[1]|0.01-0.40|5,attached_ss[1]|0.10-0.45|10,PT_ss[1]|0.01-0.70',...
    't0_crossing|SLdiff|2', 'restretchSlopeStart|1', 'peak1_y|10', 'peak1_dSL|1', 'vall_y|10', ...
    'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|1.0', 'ovrsht_dy|1', ...
    'k2|100-1500|0.01,ksrd|0.-20|0.001,kmsrd|0.0-20|0.001'};

% --- tunable pool: 2-state regavail configuration ---
cfg.pool = {'k2', 'A0', 'ka', 'kstiff2', 'v_ref_reg', 'kstiff1', 'kSE', ...
    'k_2', 'k1', 'PieceWiseStrainDepR1DParams__4', 'PieceWiseStrainDepR1DParams__3', ...
    'ksr0', 'dr', 'tau_reg'};
cfg.compulsory = {'k2', 'A0', 'ka', 'kstiff2', 'v_ref_reg'};

% --- optimizer mode: clean fminsearch on random-drawn subsets ---
cfg.SURR_EVALS = 0;   % surrogateopt DISABLED (start with pure fminsearch)
cfg.KICK_FRAC  = 0;   % no stall-kick; a non-improving round just draws a fresh random subset

% <<< true = ~10 min smoke test; false = production run (set DEBUG=false
%     for production) >>>
cfg.DEBUG = true;

optimizeFeatures(cfg);
