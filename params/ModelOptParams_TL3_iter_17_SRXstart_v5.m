% SRX+SRXD overlay v5 — run on top of ModelOptParams_TL3_iter_17.m
% Target: SRXT+SRXD = 0.200 at STEADY STATE during peak active contraction.
%
% Design:
%   ksr0=62.5      — HIGH entry rate: pumps heads into SRX despite force-dependent suppression
%   sigma_srd1=80  — SRXD is force-sensitive (20-80 kPa range) but less than SRXT
%   sigma_srd2=35.4 — SRXD entry suppressed by force (like SRXT)
%   kmsrd=1.893    — empirically tuned in full simulation for SRX_ss=0.200 at F~74 kPa
%
% Steady-state SRX:
%   F=0   (rest):    SRX ≈ 0.676  (SRXT≈0.079 + SRXD≈0.597)
%   F~74  (active):  SRX ≈ 0.200  (SRXT≈0.052 + SRXD≈0.148)
%
% Physiological rationale:
%   Higher PT→SR entry rate (ksr0 ÷12.5) overcomes force-dependent mobilization,
%   allowing SRX to plateau at 0.2 instead of depleting to zero at high force.
%   This tests whether imbalance in the SRX pool affects cycling kinetics.

params0.UseSuperRelaxed    = 1;
params0.UseSuperRelaxedADP = 1;

% PT ⇌ SR: increased entry rate
params0.ksr0 = 62.5;   % ← INCREASED from 5 to overcome force suppression
params0.kmsr = 4.0;

% SR ⇌ SRD (inter-SRX, force-independent)
params0.ksr2srd = 20.1111;
params0.ksrd2sr = 2.82418;

% PD ⇌ SRD: force-sensitive SRXD; kmsrd tuned for SRX=0.200 at active F
params0.ksrd    = 2.82418;
params0.kmsrd   = 1.893;      % empirically tuned: SRX_ss=0.200 at F~74 kPa

% Force-sensitivity exponents
params0.sigma1     = 24.7961;   % SRXT exit accelerated by force
params0.sigma2     = 42.0324;   % SRXT entry suppressed by force
params0.sigma_srd1 = 80;        % SRXD exit accelerated by force (less sensitive than SRXT)
params0.sigma_srd2 = 35.4197;   % SRXD entry suppressed by force (like SRXT)

% Initial fractions = zero-force 4-pool analytical steady state
params0.SRXT_0 = 0.0791;   % P_SR  at F=0
params0.SRXD_0 = 0.5971;   % P_SRD at F=0
