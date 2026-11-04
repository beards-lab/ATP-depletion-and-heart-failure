% SRX+SRXD overlay — run on top of ModelOptParams_TL3_iter_17.m
% Version 2: Adjusted to reach steady-state SRXT+SRXD ≈ 0.2
% After testing v1, found that decay was too rapid. Increased kmsr to slow SR→PT exit.

params0.UseSuperRelaxed = 1;
params0.UseSuperRelaxedADP = 1;

% Kinetic rates for SRX states
% Increased kmsr from 2 → 2.7 to reduce SR→PT exit flux and maintain higher SRX pool
params0.ksr0 = 5;
params0.kmsr = 2.7;  % increased from 2.0 to slow SRX decay

% ATP-cycle analogy: prevent SRD from accumulating as a sink
% ksr2srd = kah (SR loses ATP to become SRD)
% ksrd2sr = kamh (SRD regains ADP, becomes SR)
% ksrd = kamh (PD loses ADP to become SRD)
% kmsrd = kah (SRD regains ATP, becomes PD)
params0.ksr2srd = 20.1111;  % = kah
params0.ksrd2sr = 2.82418;  % = kamh
params0.ksrd = 2.82418;     % = kamh
params0.kmsrd = 20.1111;    % = kah

% Initial SRX fractions (zero-force SS, resting phase)
% These will decay slightly during high-force simulation; target ~0.2 at steady state
params0.SRXT_0 = 0.094;
params0.SRXD_0 = 0.106;
params0.sigma1 = 24.7961;
params0.sigma2 = 42.0324;
params0.sigma_srd1 = 24.573;
params0.sigma_srd2 = 35.4197;

% Optimizer knobs: all 8 SRX kinetic parameters
params0.mods = {'ksr0', 'kmsr', 'ksr2srd', 'ksrd2sr', 'ksrd', 'kmsrd', 'sigma1', 'sigma2'};
params0.g = [1  1  1  1  1  1  1  1];
