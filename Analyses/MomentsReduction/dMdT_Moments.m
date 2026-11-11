function [f, outputs] = dMdT_Moments(t, PU, params)
% DMDT_MOMENTS  Moment-closure twin of dPUdT_CombinedTransitions (prototype).
%
%   Same tail layout and same scalar physics as the grid model; only the two
%   strain distributions p1(s), p2(s) are replaced by K scaled raw moments each,
%       mu_k = sum over heads of (s/s0)^k,   k = 0..K-1,
%   so with ss := K the tail indices Ns*ss+1.. are byte-identical to the dxdt.
%
%   PU = [mu1(0..K-1), mu2(0..K-1), P_SR, NP, SL, LSE, PD, P_SRD, X_visc, A_reg, A_def]
%
%   Closure (params.MomClosure):
%     'gauss'  K = 3. Gaussian from (mu0, mu1, mu2), kernels by Gauss-Hermite.
%     'qmom'   K = 2N. N-node Gauss rule by Wheeler's algorithm, falling back
%              to fewer nodes if the moment set stops being realizable.
%
%   Strain kernels come from xbKernels (shared with the grid model's parity
%   test); the scalar periphery below is copied from dPUdT_CombinedTransitions
%   for the branches active in fastgate_best and must be kept in sync by hand
%   until the dxdt is refactored (see conclusions.md).
%
%   There is no strain grid, no window, no attachment pointer: the attachment
%   source enters as the exact moments of the continuous triangle of
%   half-width A1AttachmentWidth, and sliding enters as the exact drift term
%   d mu_k/dt += k * u * mu_{k-1}, u = (ds/dt)/s0 = velHS/(2 s0).

vel = params.Vums;
K = params.ss; s0 = params.MomScale; Ns = 2;
MAX_RATE = 1e4;

mu1 = PU(1:K); mu2 = PU(K+1:2*K);
P_SR  = max(0, PU(Ns*K+1));
NP    = 0;
SL    = PU(Ns*K+3);
LSE   = PU(Ns*K+4);
PD    = PU(Ns*K+5);
P_SRD = max(0, PU(Ns*K+6));
X_visc = PU(Ns*K+7);
A_reg = min(1, max(params.A0, PU(Ns*K+8)));
A_def = min(params.Adef_max, max(0, PU(Ns*K+9)));

%% closure: nodes (in s) and masses for p1 and p2
[x1, w1] = closeMoments(mu1, params); [x2, w2] = closeMoments(mu2, params);
y1 = x1/s0; y2 = x2/s0;                  % scaled node strains
p1_0 = max(0, mu1(1)); p2_0 = max(0, mu2(1));

%% geometry (UseOverlap, UseLatticeSpacing) -- verbatim
L_T_HS1 = min(params.L_thick*0.5, SL*0.5);
L_T_HS2 = max((SL-LSE)*0.5 - ((SL-LSE)-params.L_thin), params.L_hbare*0.5);
N_overlap = (L_T_HS1 - L_T_HS2)*2/(params.L_thick - params.L_hbare);
d10 = params.d10_ref * sqrt(params.SL_ref_lattice / SL);
excess = max(0, d10 - params.R_thick - params.R_thin - params.d_optimal);
f_lattice = exp(-(excess^2) / (2 * params.sigma_lattice^2));

%% force
K1 = xbKernels(x1, params, vel); K2 = xbKernels(x2, params, vel);
F_active = sum(w1.*K1.F1) + sum(w2.*K2.F2);
PT = max(0, 1 - (p1_0 + p2_0 + PD + P_SR + P_SRD));
N_overlap = (N_overlap - (p1_0 + p2_0))/(P_SR + P_SRD + PT + PD);        % UseOverlapFactor

F_passive = params.k_pas*max(SL-LSE - params.Lsc0, 0)^params.gamma;
if params.UseDynamicPassive && vel > 0
    F_passive = F_passive + params.c_titin_visc * vel;
end
F_passive = F_passive + X_visc;                                           % UseMaxwellDashpot
F_total = F_active + F_passive;

Force = sign(LSE)*abs(max(params.MaxSlackNegativeForce, params.kSE*LSE)).^params.ekSE;
if (Force - F_total) > 0
    velHS = (Force - F_total)/params.mu;
else
    velHS = (Force - F_total)/params.mu_neg;
end
dLSEdt = vel - velHS;
Force = Force + params.mu2*vel;

dSLc = vel - dLSEdt;                                                      % UseMaxwellTensionOnly, eps = 0
dx_dash_dt = params.kSE_M * max(0, dSLc) - X_visc / params.eta_M - X_visc * max(0, -dSLc) / params.x_M_slack;

%% scalar transitions -- verbatim for the active branches
RTD = params.kah*PT;
RDT = params.kamh*PD;
ka_eff = params.ka * A_reg * (1 - A_def);
RD1 = ka_eff*PD*N_overlap*f_lattice;
P_bound = p1_0 + p2_0;
if strcmpi(params.OccupancyForm, 'langmuir')
    RD1 = RD1 / (1 + P_bound / params.P_bound_max);
else
    RD1 = RD1 * max(0, 1 - P_bound / params.P_bound_max);
end
F_SR = F_total;

RSRD2PD = params.kmsrd*exp(F_SR/params.sigma_srd1)*P_SRD;
RPD2SRD = params.ksrd*exp(-F_SR/params.sigma_srd2)*PD;
RSR2PT  = params.kmsr*exp(F_SR/params.sigma1)*P_SR;
RPT2SR  = params.ksr0*exp(-F_SR/params.sigma2)*PT;
RSRD2SR = params.ksrd2sr*P_SRD;
RSR2SRD = params.ksr2srd*P_SR;

%% node fluxes (per node, mass/s)
J1D = w1.*K1.k1D; J12 = w1.*K1.k12;
J21 = w2.*K2.k21; J2 = w2.*K2.k2; Jhop = w2.*K2.kHop;
RdetSRX = sum(K2.fSRX .* J2);

dU_SRD = RSR2SRD - RSRD2SR - RSRD2PD + RPD2SRD;
dU_SR  = -RSR2PT + RPT2SR + RSRD2SR - RSR2SRD + RdetSRX;   % SRXRecoilToD = 0 -> lands in P_SR

RT2 = 0;
if params.UseA2Reattaching && vel > 0
    RT2 = params.kA2re * PT;
end

dPD = RSRD2PD - RPD2SRD + RTD - RDT - RD1 + sum(J1D);

%% moment equations
k = (0:K-1)';
P1 = y1(:)'.^k; P2 = y2(:)'.^k;           % K x nodes, row k = y^k
u = velHS / (2*s0);
drift1 = [0; k(2:end).*mu1(1:end-1)] * u;
drift2 = [0; k(2:end).*mu2(1:end-1)] * u;
aAtt = triMoments(params.A1AttachmentWidth/s0, K);
yLand = (params.dA2re/s0).^k; if params.dA2re == 0, yLand = double(k == 0); end
dmu1 = RD1*aAtt + drift1 - P1*(J1D + J12) + P2*J21;
dmu2 =            drift2 + P1*J12 - P2*(J21 + J2 + Jhop) + (RT2 + sum(Jhop))*yLand;

%% availability states -- verbatim
vreg = max(0, -vel) + params.RegAvailLengtheningGain * max(0, vel);
A_inf  = params.A0 + (1 - params.A0) * vreg / (vreg + params.v_ref_reg);
dA_reg = (A_inf - A_reg) / params.tau_reg;

dBound = dmu1(1) + dmu2(1);
g_load = max(0, F_total) / (max(0, F_total) + params.F_cr);
uL = -dBound; loss = uL / (1 + exp(-uL/0.5));
dA_def = params.k_cr * loss * g_load - A_def/params.tau_cr;
if A_def >= params.Adef_max; dA_def = min(0, dA_def); end

f = [dmu1; dmu2; dU_SR; 0; vel; dLSEdt; dPD; dU_SRD; dx_dash_dt; dA_reg; dA_def];
outputs = [Force, F_active, F_passive, N_overlap, p1_0, p2_0, PT, numel(x1), numel(x2)];
end

function a = triMoments(h, K)
% moments of the symmetric triangle on [-h, h], unit mass
k = (0:K-1)';
a = 2*h.^k ./ ((k+1).*(k+2)); a(mod(k,2) == 1) = 0;
end

function [x, w] = closeMoments(mu, params)
s0 = params.MomScale;
m0 = mu(1);
if m0 <= 1e-12, x = 0; w = 0; return; end
switch params.MomClosure
    case 'gauss'
        m = mu(2)/m0; v = max(0, mu(3)/m0 - m^2);
        x = s0*(m + sqrt(v)*params.cachedGH.x); w = m0*params.cachedGH.w;
    case 'qmom'
        [y, w] = wheeler(mu, numel(mu)/2);
        x = s0*y;
end
end

function [x, w] = wheeler(m, N)
% Wheeler's algorithm (Marchisio & Fox 2013), with realizability fallback:
% stops at the last node count whose recurrence coefficient b_k stayed positive.
sig = zeros(N+1, 2*N);                  % row r <-> k = r-2 (row 1 is k = -1)
sig(2,:) = m(1:2*N)';
a = zeros(N,1); b = zeros(N,1);
a(1) = m(2)/m(1); b(1) = m(1);
Nok = 1;
for kk = 1:N-1
    r = kk + 2;
    for l = kk:2*N-kk-1
        sig(r, l+1) = sig(r-1, l+2) - a(kk)*sig(r-1, l+1) - b(kk)*sig(r-2, l+1)*(kk > 1);
    end
    bk = sig(r, kk+1)/sig(r-1, kk);
    if ~(bk > 1e-10*max(1, abs(a(1)))^2) || ~isfinite(bk), break; end
    a(kk+1) = sig(r, kk+2)/sig(r, kk+1) - sig(r-1, kk+1)/sig(r-1, kk);
    b(kk+1) = bk; Nok = kk + 1;
end
a = a(1:Nok); b = b(1:Nok);
J = diag(a) + diag(sqrt(b(2:end)), 1) + diag(sqrt(b(2:end)), -1);
[V, D] = eig(J);
x = diag(D); w = m(1) * V(1,:)'.^2;
end
