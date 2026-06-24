% SRX_FastGentle.m  —  Fast, gently force-sensitive SRX/DRX reparametrization
% =========================================================================
% Overlay: run AFTER ModelParamsSlackKtrOpt.m (or any base), then getParams.
%
%   params0 = getParams(); run('ModelOptParams\ModelParamsSlackKtrOpt.m');
%   run('ModelOptParams\SRX_FastGentle.m');  params0 = getParams(params0);
%
% DESIGN GOALS (and how each parameter delivers them)
% ---------------------------------------------------
%   1. FAST: DRX<->SRX equilibrates in ~10 ms at baseline (was ~940 ms).
%   2. GENTLE force sensitivity: ~2-3x SRX swing over 0->85 kPa (was ~50x,
%      i.e. SRX no longer collapses to ~0 under load).
%   3. VISIBLE: active-isometric SRX ~0.11, resting ~0.20-0.25 (was 0.005).
%   4. NO baseline skew: isometric & force-velocity forces unchanged, because
%      the parked heads are drawn from the non-force-generating free pool.
%
% THE KINETICS (model variable names in parentheses)
% --------------------------------------------------
%   SRX has two sub-pools that behave identically here (symmetric design):
%     P_SR  <-> PT   mobilize=kmsr  (SR->PT),  park=ksr0  (PT->SR)
%     P_SRD <-> PD   mobilize=kmsrd (SRD->PD), park=ksrd  (PD->SRD)
%   plus fast P_SR<->P_SRD mixing (ksr2srd, ksrd2sr).
%
%   * SPEED is set by the MOBILIZATION rate. The slow relaxation eigenvalue of
%     the symmetric subsystem equals the mobilization rate, so
%         tau_rest ~ 1000/kmsr ms.  kmsr=100 -> tau=10 ms.  (verified in sim)
%     Force accelerates mobilization, so the active working point is even
%     faster (~3.4 ms).
%   * POOL SIZE is set by the park/mobilize RATIO. ksr0/kmsr = 40/100 = 0.4.
%     Scale ksr0 & ksrd together to grow/shrink the whole SRX pool without
%     changing speed or sensitivity.
%   * FORCE SENSITIVITY is carried by sigma1/sigma_srd1 only (mobilization);
%     parking is made force-independent (sigma2/sigma_srd2 huge), restoring
%     the original getParams design intent (default sigma2 was 1e6).
%     Smaller sigma1 = MORE sensitive:  sigma1=60 -> 3.0x swing
%                                       sigma1=80 -> 2.3x swing  (recommended)
%                                       sigma1=120-> 1.7x swing
%
% VALIDATION (activated isometric, 8 mM ATP, vs current ModelParamsSlackKtrOpt)
% ---------------------------------------------------------------------------
%                       current     this file (sigma1=80)
%   tau equilibration   ~938 ms     10.0 ms (rest) / 3.4 ms (active)
%   SRX @ active force   0.005       0.111
%   SRX @ ~zero force    0.27        0.25
%   force swing          ~50x        2.3x
%   isometric force      85.1        85.2   (unchanged)
%   force-velocity       --          unchanged at v=0,-1,-2,-3
% =========================================================================

params0.UseSuperRelaxed    = 1;
params0.UseSuperRelaxedADP = 1;

% --- Mobilization (out of SRX): sets the ~10 ms equilibration speed ---
params0.kmsr  = 100;   % SR  -> PT
params0.kmsrd = 100;   % SRD -> PD

% --- Parking (into SRX): sets the pool size via ratio to mobilization ---
params0.ksr0  = 40;    % PT  -> SR
params0.ksrd  = 40;    % PD  -> SRD

% --- Force sensitivity: gentle, mobilization-only ---
params0.sigma1     = 80;    % exp(F/80): ~2.9x mobilization accel at F=85
params0.sigma_srd1 = 80;
params0.sigma2     = 1e4;   % parking ~ force-independent
params0.sigma_srd2 = 1e4;

% --- Fast, symmetric SR <-> SRD mixing (keeps the two sub-pools balanced) ---
params0.ksr2srd = 50;
params0.ksrd2sr = 50;

% --- Initial SRX fractions ~ active steady state (optional: fast tau makes
%     any startup transient decay in ~10 ms, so exact values are not critical) ---
params0.SRXT_0 = 0.06;   % initial P_SR
params0.SRXD_0 = 0.05;   % initial P_SRD
