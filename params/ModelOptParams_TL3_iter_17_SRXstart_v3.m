% SRX+SRXD overlay — run on top of ModelOptParams_TL3_iter_17.m
% Version 3: Further optimized to reach steady-state SRXT+SRXD ≈ 0.2
% v2 achieved 0.1742; v3 increases kmsr further to 3.5

params0.UseSuperRelaxed = 1;
params0.UseSuperRelaxedADP = 1;

% Kinetic rates for SRX states
% v2 had kmsr=2.7 → achieved 0.1742, error=0.0258
% v3 increases to kmsr=3.5 (estimate based on linear extrapolation)
params0.ksr0 = 5;
params0.kmsr = 3.5;  % increased from 2.7 to further slow SR→PT exit flux

% ATP-cycle analogy: prevent SRD from accumulating as a sink
params0.ksr2srd = 20.1111;  % = kah
params0.ksrd2sr = 2.82418;  % = kamh
params0.ksrd = 2.82418;     % = kamh
params0.kmsrd = 20.1111;    % = kah

% Initial SRX fractions
params0.SRXT_0 = 0.094;
params0.SRXD_0 = 0.106;
params0.sigma1 = 24.7961;
params0.sigma2 = 42.0324;
params0.sigma_srd1 = 24.573;
params0.sigma_srd2 = 35.4197;

% Optimizer knobs: all 8 SRX kinetic parameters
params0.mods = {'ksr0', 'kmsr', 'ksr2srd', 'ksrd2sr', 'ksrd', 'kmsrd', 'sigma1', 'sigma2'};
params0.g = [1  1  1  1  1  1  1  1];
