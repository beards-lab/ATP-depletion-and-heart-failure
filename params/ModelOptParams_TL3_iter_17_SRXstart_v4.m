% SRX+SRXD overlay v4 — run on top of ModelOptParams_TL3_iter_17.m
% Analytically derived: kmsrd=14.4 achieves SRXT+SRXD=0.200 at zero force.
%
% Root cause of v1-v3 failure: with kmsrd=20.11 the PD pool (≈70% of free
% heads) drains P_SRD back to PD too fast. Maximum achievable SRX ≈ 17%
% regardless of kmsr. Reducing kmsrd to 14.4 solves the drain.
%
% 4-pool balance at zero force (PT,PD,P_SR,P_SRD), TL3 kah/kamh=20.11/2.82:
%   PD/PT ≈ 7.1   →  PD ≈ 0.70 of free pool
%   kmsrd = 14.4  →  P_SR = 0.039, P_SRD = 0.161, total SRX = 0.200
%
% Physiological context:
%   SRXD exit lifetime  = 1/(ksrd2sr + kmsrd) = 1/17.2  ≈  58 ms
%   SRXT exit lifetime  = 1/(kmsr   + ksr2srd) = 1/24.1  ≈  41 ms
%   Both ~40-60 ms — consistent with heads transiently parking in IHM.

params0.UseSuperRelaxed    = 1;
params0.UseSuperRelaxedADP = 1;

% PT ⇌ SR
params0.ksr0 = 5;
params0.kmsr = 4.0;

% SR ⇌ SRD  (ATP-cycle analogy; ksr2srd = kah, ksrd2sr = kamh)
params0.ksr2srd = 20.1111;
params0.ksrd2sr = 2.82418;

% PD ⇌ SRD
params0.ksrd  = 2.82418;
params0.kmsrd = 14.4;      % KEY: reduced from 20.11 to allow SRX pool to fill

% Force-sensitivity exponents (unchanged from v1-v3)
params0.sigma1     = 24.7961;
params0.sigma2     = 42.0324;
params0.sigma_srd1 = 24.573;
params0.sigma_srd2 = 35.4197;

% Initial fractions = analytical zero-force steady state
params0.SRXT_0 = 0.039;   % P_SR  at F=0
params0.SRXD_0 = 0.161;   % P_SRD at F=0
