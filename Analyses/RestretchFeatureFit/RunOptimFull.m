% RunOptimFull.m
% =========================================================================
% COMPREHENSIVE feature-fit optimizer. One driver, model-agnostic: it reads the
% seed, includes the p3 levers only when NumberOfStates==3, and auto-drops any
% candidate whose seed value is 0 (a multiplicative modifier g can never move a
% zero base). Hands the assembled pool to optimizeFeatures.m (random-subset
% block-coordinate descent: each round draws N_DRAW params and runs fminsearch on
% just that block; over many rounds every param -- and every co-drawn coupling --
% is sampled).
%
% ---- SEED (pick one) -------------------------------------------------------
%   2-state, current best + A2 hop (cost ~2.95):  params/params_2state_a2hop.m
%   3-state (Design B) + A2 hop:                   params/params_full3_seed.m
% Auto-resume: re-running continues from params/optfull_state.mat (delete it to
% start fresh). Safe to leave running for hours; incumbent is written on every
% improvement.
%
% ---- EFFECTIVE NUMBER (N_DRAW) --------------------------------------------
% fminsearch is Nelder-Mead: reliable to ~10 dimensions, degrades sharply beyond
% (N+1 simplex vertices, ill-conditioned contraction in high-D). So each round
% optimizes an EFFECTIVE block of N_DRAW=10 params, NOT the whole ~40-param pool.
% The block captures the pairwise couplings that dominate this fit (ktr<->FV via
% k2+kmsrd, force<->peak via kstiff2+ka, peak1_dSL<->ktr via kSE+kmsrd) whenever
% those params co-draw; the multi-round structure samples the rest. SIMPLEX_EVALS
% is set so the simplex gets ~10-12 evals/dim (enough to descend a 10-D block).
%
% ---- ALWAYS INCLUDED (compulsory) -----------------------------------------
% Every drawn subset perturbs some shape/rate/timescale param with side effects on
% the two globally-dominant observables -- isometric FORCE LEVEL and cycling RATE.
% If a subset can't compensate those, the round can't improve and is wasted. So we
% ALWAYS include the "rebalancers":
%   k2       -- master cycling rate (ktr / FV / ATPase)
%   kstiff2  -- master force scale (stiffness; also sets the peaks)
%   ka       -- master attachment: compensates force via bridge NUMBER (peak-neutral),
%               distinct from kstiff2 (peak-inflating). Having both lets each round
%               pick peak-neutral vs peak-affecting force compensation.
% 3-state adds k3 + kstiff3 so the new p3 turnover/force balance -- which everything
% now trades against -- is always co-tuned (else a cycle change leaves p3 stranded).
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));

cfg = struct();
cfg.baseSnap = 'params/params_2state_a2hop.m';   % <<< 2-state ~2.95. Switch to params/params_full3_seed.m for 3-state.
cfg.baseSnap = 'params/optfull3_opt.m';
cfg.tag      = 'optfull4';

cfg.fn = {'FV_fnorm|FV_v|10', 'ktr|4', 'A|50', 'ktr_rmse|0-.5|.1', ...
    'XTOR[1]|3-15|15,XTOR_vmax[1]|5-25,SRX_ss[1]|0.01-0.40|5,attached_ss[1]|0.10-0.45|10,PT_ss[1]|0.01-0.70',...
    't0_crossing|SLdiff|2', 'restretchSlopeStart|1', 'peak1_y|10', 'peak1_dSL|1', 'vall_y|10', ...
    'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|1.0', 'ovrsht_dy|1', ...
    'k2|100-1500|0.01,ksrd|0.-20|0.001,kmsrd|0.0-20|0.001', 'doublePeak|10', 'coolDownLS|0.05'};

% ---- read the seed so the pool matches the model ----
params0 = getParams(); run(cfg.baseSnap); p0 = params0;
NoS = p0.NumberOfStates;

% ---- COMPREHENSIVE candidate pool (everything with demonstrated leverage) ----
cand = { ...
    ... % core cross-bridge cycle (force, ktr, FV, ATPase)
    'ka','kd','k1','k_1','k2','k_2','kah', ...
    ... % stiffness / force scale
    'kstiff1','kstiff2','kstiff1_n','kstiff2_n','kSE', ...
    ... % series viscoelastic + restretch transient (peak1_dSL region)
    'kSE_M','eta_M','mu_neg', ...
    ... % power-stroke strain offsets
    'dr','dr2', ...
    ... % SRX / SRD manifold (ktr timescale + slow recovery)
    'ksr0','kmsrd','ksr2srd','ksrd2sr','sigma1','sigma2','sigma_srd1','sigma_srd2', ...
    ... % registration availability (FV shoulder + restretch undershoot)
    'A0','v_ref_reg','tau_reg', ...
    ... % occupancy ceiling + passive titin (restretch force)
    'P_bound_max','k_pas', ...
    ... % A2 restretch hop (valley / peak2)
    'kA2hop','sA2hop', ...
    ... % strain-dependent rate knots (PCHIP interior knots; boundary clamps excluded)
    'PieceWiseStrainDep2X_logOffset', 'PieceWiseStrainDepR1DX_logOffset', 'PieceWiseStrainDepR21X_logOffset', ...
    'PieceWiseStrainDepX__2','PieceWiseStrainDepX__3','PieceWiseStrainDepX__4', ...
    'PieceWiseStrainDep2X__2','PieceWiseStrainDep2X__3','PieceWiseStrainDep2X__4', ...
    'PieceWiseStrainDep2X__2','PieceWiseStrainDep2X__3','PieceWiseStrainDep2X__4', 'PieceWiseStrainDep2X__5',...
    'PieceWiseStrainDepR1DX__1','PieceWiseStrainDepR1DX__2','PieceWiseStrainDepR1DX__3','PieceWiseStrainDepR1DX__4', ...
    'PieceWiseStrainDepR21X__2','PieceWiseStrainDepR21X__3','PieceWiseStrainDepR21X__4',...
    'PieceWiseStrainDepParams__2','PieceWiseStrainDepParams__3','PieceWiseStrainDepParams__4', ...
    'PieceWiseStrainDep2Params__2','PieceWiseStrainDep2Params__3','PieceWiseStrainDep2Params__4', ...
    'PieceWiseStrainDep2Params__2','PieceWiseStrainDep2Params__3','PieceWiseStrainDep2Params__4', 'PieceWiseStrainDep2Params__5',...
    'PieceWiseStrainDepR1DParams__1','PieceWiseStrainDepR1DParams__2','PieceWiseStrainDepR1DParams__3','PieceWiseStrainDepR1DParams__4', ...
    'PieceWiseStrainDepR21Params__2','PieceWiseStrainDepR21Params__3','PieceWiseStrainDepR21Params__4'};

if NoS == 3
    cand = [cand, {'k3','k3m','kstiff3','drp3','alpha3','s3'}];   % 3rd-state levers
end

% ---- drop degenerate candidates (0-valued in seed -> immovable by g) ----
keep = true(1, numel(cand));
for i = 1:numel(cand); keep(i) = (seedval(p0, cand{i}) ~= 0); end
dropped = cand(~keep);
cfg.pool = cand(keep);
if ~isempty(dropped)
    fprintf('[RunOptimFull] dropped %d zero-base params (untunable): %s\n', numel(dropped), strjoin(dropped, ', '));
end

% ---- always-included rebalancers ----
cfg.compulsory = {'k2','kstiff2','ka'};
if NoS == 3; cfg.compulsory = [cfg.compulsory, {'k3','kstiff3'}]; end

% ---- effective number + budgets ----
cfg.N_DRAW        = 10;                          % Nelder-Mead-tractable block size
cfg.SIMPLEX_EVALS = (NoS==3)*70 + (NoS==2)*80;  % ~7-12 evals/dim
cfg.MaxRunTime    = (NoS==3)*90 + (NoS==2)*60;   % generous ceiling: the A2-hop config is stiffer and flaked at 45s under load
cfg.SURR_EVALS    = 0;                           % pure fminsearch (no surrogate)
cfg.KICK_FRAC     = 0;
cfg.DEBUG           = false;
cfg.TIME_BUDGET_HRS = 48;                        % long haul; raise/lower as desired
cfg.RESUME          = false;%isfile(fullfile(root, 'params', 'optfull_state.mat'));

fprintf('[RunOptimFull] model=%d-state | pool=%d | N_DRAW=%d | compulsory={%s}\n', ...
    NoS, numel(cfg.pool), cfg.N_DRAW, strjoin(cfg.compulsory, ','));

optimizeFeatures(cfg);

% ---- local: read a scalar/array-element param value by name (handles Field__N) ----
function v = seedval(p, name)
    us = strfind(name, '__');
    if isempty(us)
        v = p.(name);
    else
        f = name(1:us(1)-1); idx = str2double(name(us(1)+2:end));
        arr = p.(f); v = arr(idx);
    end
end
