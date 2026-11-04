% SRX+SRXD overlay — run on top of ModelOptParams_TL3_iter_17.m
% Tuned for total SRX ~ 0.2 at steady-state (resting / low-force phases).
% Full ATP-cycle analogy: ksr2srd=kah, ksrd2sr=kamh, kmsrd=kah, ksrd=kamh.
% This keeps SRD drainable so it doesn't accumulate as a sink.
% Steady-state testing (zero-velocity, 3s) confirms SRXT+SRXD ≈ 0.2.
% At F=80 kPa SRX drops well below 0.2 — sigma1/sigma2 are the optimizer's
% knobs for adjusting the force-dependent dynamic range.

params0.UseSuperRelaxed = 1;
params0.UseSuperRelaxedADP = 1;

% PT ⇌ SR: ksr0 (forward), kmsr (backward)
% Optimized kmsr=4.0 to balance SR→PT exit rate and maintain ~0.2 total SRX
params0.ksr0 = 5;
params0.kmsr = 4.0;

% ATP-cycle analogy for SRD transitions
% SR ⇌ SRD: ksr2srd=kah (ATP-driven), ksrd2sr=kamh (reverse)
% PD ⇌ SRD: ksrd=kamh, kmsrd=kah (prevents SRD accumulation)
params0.ksr2srd = 20.1111;  % = kah
params0.ksrd2sr = 2.82418;  % = kamh
params0.ksrd = 2.82418;     % = kamh
params0.kmsrd = 20.1111;    % = kah

% Force-sensitive exponents
params0.sigma1 = 24.7961;
params0.sigma2 = 42.0324;
params0.sigma_srd1 = 24.573;
params0.sigma_srd2 = 35.4197;

% Steady-state fractions at zero force (resting)
% These correspond to steady-state after ~3 seconds
params0.SRXT_0 = 0.094;
params0.SRXD_0 = 0.106;
