% RunOptim2State_v2.m
% =========================================================================
% Expanded-pool 2-state feature fit. Same machinery as RunOptim2State.m, but
% seeds from params/params_2state_v2_seed.m (= current best opt2state_opt.m
% with the frozen SRX-return rate kmsrd reset from 60 -> 18, i.e. inside its
% own parameterBounds range [0.1,20]) and ADDS the previously-frozen
% SRX-return timescale levers to the tunable pool:
%
%   kmsrd, ksr2srd, ksrd2sr   (the SRD-manifold return kinetics)
%
% WHY: coupling-probe (Analyses/RestretchFeatureFit/probeCoupling_report.txt)
% showed the ktr residual (1.70, the single largest cost term at 4.55) is NOT
% a structural iron wall. Core-cycle levers (k2, global cycle) lower ktr only
% by breaking FV. But kmsrd -- frozen at 60 (3x over its [0.1,20] bound and
% absent from every prior pool) -- lowers ktr 69->58 while HOLDING steady
% force, peaks, FV shoulder and all physiology bounds; single-lever cost
% 4.55 -> 3.67. kmsrd sets the slow force-redevelopment timescale via the
% SR->SRD->PD return path (mechanism-evaluation.md sec 1.4 / 3.2): a fast
% return adds a fast component to redevelopment that inflates the fitted ktr.
%
% This run frees that timescale (+ its two partner rates) so the optimizer
% can trade ktr down to ~49 and recover the small FV-tail dip via the
% existing pool (ka, kSE, k2, ...).
%
% Isolation: cfg.tag='opt2state_v2' -> params/opt2state_v2_opt.m and
% params/opt2state_v2_state.mat. Never touches opt2state_* or opt3state_*.
%
% Run:  cd(root); addpath(genpath('.')); RunOptim2State_v2
% =========================================================================

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));

cfg = struct();
cfg.baseSnap = 'params/params_2state_v2_seed.m';
cfg.tag      = 'opt2state_v2';

% --- feature list: identical to RunOptim2State (verbatim USERFN) ---
cfg.fn = {'FV_fnorm|FV_v|10', 'ktr|2', 'A|50', 'ktr_rmse|0-.5|.1', ...
    'XTOR[1]|3-15|15,XTOR_vmax[1]|5-25,SRX_ss[1]|0.01-0.40|5,attached_ss[1]|0.10-0.45|10,PT_ss[1]|0.01-0.70',...
    't0_crossing|SLdiff|2', 'restretchSlopeStart|1', 'peak1_y|10', 'peak1_dSL|1', 'vall_y|10', ...
    'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|1.0', 'ovrsht_dy|1', ...
    'k2|100-1500|0.01,ksrd|0.-20|0.001,kmsrd|0.0-20|0.001'};

% --- tunable pool: 2-state regavail + freed SRX-return timescale + FV-tail knots ---
% NOTE: the R2 negative-strain (shortening-branch) knots
% PieceWiseStrainDep2Params__5/6/7 shape the FV tail. They MUST be in the pool:
% a first bounded run without them lowered ktr via kmsrd (69->~50, ktr solved)
% but let the FV tail collapse (.66->.53, .32->.21) with no lever to restore it.
cfg.pool = {'k2', 'A0', 'ka', 'kstiff2', 'v_ref_reg', 'kstiff1', 'kSE', ...
    'k_2', 'k1', 'PieceWiseStrainDepR1DParams__4', 'PieceWiseStrainDepR1DParams__3', ...
    'ksr0', 'dr', 'tau_reg', ...
    'kmsrd', 'ksr2srd', 'ksrd2sr', ...                                  % freed SRX-return timescale (ktr)
    'PieceWiseStrainDep2Params__5', 'PieceWiseStrainDep2Params__6', 'PieceWiseStrainDep2Params__7'};  % R2 shortening knots (FV tail)
cfg.compulsory = {'k2', 'A0', 'ka', 'kstiff2', 'v_ref_reg', 'kmsrd'};  % kmsrd always drawn

% --- optimizer mode: clean fminsearch on random-drawn subsets ---
cfg.SURR_EVALS = 0;   % surrogateopt DISABLED
cfg.KICK_FRAC  = 0;   % (kick only perturbs a throwaway working copy; incumbent is protected)

% --- BOUNDED confirmation run: ~36 min, then stop and write best snapshot.
%     Raise TIME_BUDGET_HRS (or delete it -> default 28h) for a full production run. ---
cfg.DEBUG           = false;
cfg.TIME_BUDGET_HRS = 0.6;

optimizeFeatures(cfg);
