function [f, outputs, rates] = dPUdT_CombinedTransitions_mex_flat(t, PU, params)
% DPUDT_COMBINEDTRANSITIONS_MEX_FLAT  Pure-MATLAB fixed-code-path version.
%
%   Implements exactly the same computation as dPUdT_CombinedTransitions but
%   with all boolean switches hard-coded to their values in ModelParamsSlackKtrOpt:
%     UseOverlap=1, UseNegativeKstiff=1, UsePassive=1 (not Titin),
%     UsePieceWiseStrainDep=1, UseA2AttachmentShift=1,
%     UseSuperRelaxed=1, UseSuperRelaxedADP=1,
%     UseDirectSRXTransition=0, UsePassiveForSR=0, useHalfActiveForSR=0,
%     UseA1AttachmentKernel=0, UseTargetZoneSaturation=0,
%     FudgeVmax=0, UseViscoelasticSE=0, UseMaxwellDashpot=0,
%     LegacyStrainFlipping=0, NumberOfStates=2.
%
%   Reads all scalars from params.mex_vals (packed by packMexVals).
%   Reads params.Vums and params.LXBpivot directly (change between ode15s calls).
%
%   Purpose: debugging step before C MEX translation. Outputs must match
%   dPUdT_CombinedTransitions to within 1e-12 on the same inputs.
%
%   See also: packMexVals, dPUdT_CombinedTransitions_mex, dPUdT_core_mex

% --- Index constants (must match packMexVals.m and dPUdT_core_mex.c) ---
N_SCALARS = 43;
N_SS = 60;   % params.ss — fixed for this config (MaxStrainArraySize=60)

% Pchip layout offsets (1-based)
N_R1D   = 3;   % R1D has 3 segments → 4 breaks, 12 coefs
N_PWSD  = 4;   % pwsd has 4 segments → 5 breaks, 16 coefs
N_R21   = 3;   % R21 has 3 segments → 4 breaks, 12 coefs
N_PWSD2 = 6;   % pwsd2 has 6 segments → 7 breaks, 24 coefs

OFF_R1D_BRK  = N_SCALARS + N_SS + 1;          % 104
OFF_R1D_COF  = OFF_R1D_BRK + (N_R1D+1);       % 108
OFF_PWSD_BRK = OFF_R1D_COF + N_R1D*4;         % 120
OFF_PWSD_COF = OFF_PWSD_BRK + (N_PWSD+1);     % 125
OFF_R21_BRK  = OFF_PWSD_COF + N_PWSD*4;       % 141
OFF_R21_COF  = OFF_R21_BRK + (N_R21+1);       % 145
OFF_PWSD2_BRK = OFF_R21_COF + N_R21*4;        % 157
OFF_PWSD2_COF = OFF_PWSD2_BRK + (N_PWSD2+1);  % 164

% --- unpack flat array ---
mv = params.mex_vals;
vel      = params.Vums;
LXBpivot = params.LXBpivot;

kd     = mv( 1);  k1     = mv( 2);  k_1    = mv( 3);  k2     = mv( 4);
kstiff1   = mv( 5);  kstiff2   = mv( 6);
kstiff1_n = mv( 7);  kstiff2_n = mv( 8);
ka     = mv( 9);  kah    = mv(10);  kamh   = mv(11);
dS     = mv(12);
% ss  = mv(13);   % = N_SS = 60
dr     = mv(14);  dr2    = mv(15);
L_thick = mv(16); L_hbare = mv(17); L_thin  = mv(18);
k_pas   = mv(19); Lsc0    = mv(20); gam     = mv(21);
kSE     = mv(22); ekSE    = mv(23);
MaxSlackNegativeForce = mv(24);
mu      = mv(25); mu_neg  = mv(26); mu2     = mv(27);
ksr0    = mv(28); sigma1  = mv(29); sigma2  = mv(30);
kmsr    = mv(31); kmsrd   = mv(32); sigma_srd1 = mv(33);
ksrd    = mv(34); sigma_srd2 = mv(35);
ksrd2sr = mv(36); ksr2srd = mv(37);
slope   = mv(38); d_actin = mv(39);
s_threshold_R = mv(40); s_threshold_L = mv(41);
Ax      = mv(42); estiff  = mv(43);

% Strain grid (column vector)
s_grid = mv(N_SCALARS+1 : N_SCALARS+N_SS);

% Pchip R1D
b_R1D  = mv(OFF_R1D_BRK  : OFF_R1D_BRK+N_R1D);
c_R1D  = reshape(mv(OFF_R1D_COF  : OFF_R1D_COF +N_R1D *4-1), N_R1D,  4);
% Pchip pwsd (R12)
b_pwsd = mv(OFF_PWSD_BRK : OFF_PWSD_BRK+N_PWSD);
c_pwsd = reshape(mv(OFF_PWSD_COF : OFF_PWSD_COF+N_PWSD*4-1), N_PWSD, 4);
% Pchip R21
b_R21  = mv(OFF_R21_BRK  : OFF_R21_BRK+N_R21);
c_R21  = reshape(mv(OFF_R21_COF  : OFF_R21_COF +N_R21 *4-1), N_R21,  4);
% Pchip pwsd2 (R2)
b_pwsd2 = mv(OFF_PWSD2_BRK : OFF_PWSD2_BRK+N_PWSD2);
c_pwsd2 = reshape(mv(OFF_PWSD2_COF : OFF_PWSD2_COF+N_PWSD2*4-1), N_PWSD2, 4);

%% Unpack state vector (Ns=2, ss=60)
p1 = PU(1:N_SS); p1(p1<0) = 0;
p2 = PU(N_SS+1:2*N_SS); p2(p2<0) = 0;

P_SR  = PU(2*N_SS+1);
% NP = 0 — not a state we track
SL    = PU(2*N_SS+3);
LSE   = PU(2*N_SS+4);
PD    = PU(2*N_SS+5);
P_SRD = PU(2*N_SS+6);
% x_dash = 0, F_Maxwell = 0, dx_dash_dt = 0

%% Overlap (UseOverlap=1, UseOverlapFactor=0)
L_T_HS1 = min(L_thick*0.5, SL*0.5);
L_T_HS2 = max((SL-LSE)*0.5 - ((SL-LSE)-L_thin), L_hbare*0.5);
L_ov    = L_T_HS1 - L_T_HS2;
N_overlap = (L_ov*2/(L_thick - L_hbare));

% f_lattice = 1 (UseLatticeSpacing=0)

%% Strain shift (LegacyStrainFlipping=0)
s = s_grid - (-(SL - LSE) + LXBpivot)/2;

%% Force: UseNegativeKstiff=1
p1_0 = dS*sum(p1);
p2_0 = dS*sum(p2);

sdr = sign(s+dr).*abs(s+dr).^estiff;   % p2 strain weighting

p1_1_pos = dS*sum(s.*p1.*(s>=0));
p1_1_neg = dS*sum(s.*p1.*(s<=0));
p2_1_pos = dS*sum(sdr.*p2.*(s >= -dr));
p2_1_neg = dS*sum(sdr.*p2.*(s < -dr));

F_active = kstiff2*p2_1_pos + kstiff2_n*p2_1_neg + kstiff1*p1_1_pos + kstiff1_n*p1_1_neg;
p1_1 = p1_1_neg + p1_1_pos;
p2_1 = p2_1_pos + p2_1_neg;

%% PT
PT = max(0, 1 - (p1_0 + p2_0 + PD + P_SR + P_SRD));

%% Passive (UsePassive=1, UseTitinIdentifiedPassive=0)
F_passive = k_pas * max(SL-LSE - Lsc0, 0)^gam;

F_total = F_active + F_passive;

%% Serial stiffness / force (FudgeVmax=0, UseViscoelasticSE=0)
Force = sign(LSE)*abs(max(MaxSlackNegativeForce, kSE*LSE))^ekSE;
if (Force - F_total) > 0
    velHS = (Force - F_total)/mu;
else
    velHS = (Force - F_total)/mu_neg;
end
dLSEdt = vel - velHS;
Force = Force + mu2*vel;   % mu2 = 1e-9 (negligible)

%% Transition rates: UsePieceWiseStrainDep=1 (inline polynomial eval)
R1D = kd  * p1 .* peval(b_R1D,  c_R1D,  s);
R12 = k1  * p1 .* peval(b_pwsd, c_pwsd, s);
R21 = k_1 * p2 .* peval(b_R21,  c_R21,  s);
s2  = s + dr2 - dr;
R2  = k2  * p2 .* peval(b_pwsd2, c_pwsd2, s2);

%% UseA2AttachmentShift=1
slope_over_dr = slope / dr;
RA_K = max(0, slope_over_dr*(s - s_threshold_R)) ...
     + max(0, -slope_over_dr*(s + s_threshold_L));
dp2_RAm = p2 .* RA_K;

s_target = (s > s_threshold_R).*(s - d_actin) ...
         + (s < -s_threshold_L).*(s + d_actin);

L_idx = max(1, min(N_SS-1, 1 + floor((s_target - s(1)) / dS)));
R_idx = L_idx + 1;
w_R = (s_target - s(L_idx)) ./ (s(R_idx) - s(L_idx));
w_L = 1 - w_R;

dp2_RAL = accumarray(L_idx, dp2_RAm .* w_L, [N_SS, 1]);
dp2_RAR = accumarray(R_idx, dp2_RAm .* w_R, [N_SS, 1]);

%% SRX dynamics (UseSuperRelaxed=1, UseSuperRelaxedADP=1, UsePassiveForSR=0)
F_SR = F_total;

RTD = kah  * PT;
RDT = kamh * PD;
RD1 = ka   * PD * N_overlap;   % UseStretchActivation=0, f_lattice=1

RSRD2PD = kmsrd  * exp(F_SR/sigma_srd1) * max(0, P_SRD);
RPD2SRD = ksrd   * exp(-F_SR/sigma_srd2) * PD;
RSR2PT  = kmsr   * exp(F_SR/sigma1) * max(0, P_SR);
RPT2SR  = ksr0   * exp(-F_SR/sigma2) * max(0, PT);
RSRD2SR = ksrd2sr * max(0, P_SRD);
RSR2SRD = ksr2srd * max(0, P_SR);

dU_SRD = RSR2SRD - RSRD2SR - RSRD2PD + RPD2SRD;
dU_SR  = -RSR2PT + RPT2SR + RSRD2SR - RSR2SRD;   % UseDirectSRXTransition=0

%% Governing flows
dPD = RSRD2PD - RPD2SRD + RTD - RDT - RD1 + sum(R1D)*dS;
dp1 = -R1D - R12 + R21;
dp2 = R12 - R21 - R2 + dp2_RAR + dp2_RAL - dp2_RAm;

%% Attachment: UseA1AttachmentKernel=0, Ax=0.006 → triangular window
halfspan = ceil(Ax/dS);   % = 3
s_0 = 1 + round(-s(1)/dS, 6);
att = zeros(N_SS, 1);
for jj = floor(max(1, s_0-halfspan)):ceil(min(N_SS, s_0+halfspan))
    dist = abs(s(jj));
    w = max(0, 1 - dist/Ax);
    att(jj) = w;
end
att_norm = att / max(1, sum(att)) / dS;
dp1 = dp1 + att_norm * RD1;   % f_sat = 1

dNP  = 0;
dSL  = vel;
% dx_dash_dt = 0

%% Assemble (Ns=2, dx_dash_dt=0)
f = [dp1; dp2; dU_SR; dNP; dSL; dLSEdt; dPD; dU_SRD; 0];

%% Outputs and rates (only when requested)
if nargout > 1
    F_Maxwell    = 0;
    f_lattice    = 1;
    f_saturation = 1;
    outputs = [Force, F_active, F_passive, N_overlap, ...
               p1_0, p2_0, p1_1, p2_1, PT, F_Maxwell, f_lattice, f_saturation];
    % rates layout matches dPUdT_CombinedTransitions (16 elements):
    % [RTD, RDT, RD1, sum(R1D)*dS, sum(R12)*dS, sum(R21)*dS, sum(R2)*dS,
    %  sum(XB_Ripped)*dS, RSR2PT, RPT2SR, RSRD2PD, RPD2SRD, RSR2SRD,
    %  RSRD2SR, RT2, sum(R2D)*dS]
    RT2 = 0;
    rates = [RTD, RDT, RD1, sum(R1D)*dS, sum(R12)*dS, sum(R21)*dS, sum(R2)*dS, ...
             0, RSR2PT, RPT2SR, RSRD2PD, RPD2SRD, RSR2SRD, RSRD2SR, RT2, 0];
end

end

%% Local functions
function v = peval(b, c, xi)
% PEVAL  Evaluate piecewise cubic polynomial.
%   B is the break vector (column or row, length n+1).
%   C is the coefficient matrix (n × 4), highest degree first.
%   XI is the query vector (column). Returns column of same length.
%   Equivalent to ppval but without binary search overhead.
    brow = b(1:end-1)';      % (1 × n) row for broadcast
    si   = max(1, min(size(c,1), sum(xi >= brow, 2)));  % (ss × 1)
    dx   = xi - b(si);
    v    = ((c(si,1).*dx + c(si,2)).*dx + c(si,3)).*dx + c(si,4);
end
