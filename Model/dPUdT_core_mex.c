/*
 * dPUdT_core_mex.c  —  C MEX for the cardiac cross-bridge ODE RHS.
 *
 * Implements the fixed code path corresponding to ModelParamsSlackKtrOpt:
 *   UseOverlap=1, UseNegativeKstiff=1, UsePassive=1 (k_pas branch),
 *   UsePieceWiseStrainDep=1 (inline polynomial), UseA2AttachmentShift=1,
 *   UseSuperRelaxed=1, UseSuperRelaxedADP=1,
 *   UseDirectSRXTransition=0, UsePassiveForSR=0, useHalfActiveForSR=0,
 *   UseA1AttachmentKernel=0, FudgeVmax=0, UseViscoelasticSE=0,
 *   UseMaxwellDashpot=0, LegacyStrainFlipping=0, NumberOfStates=2.
 *
 * MATLAB signature (via dPUdT_CombinedTransitions_mex.m wrapper):
 *   [f, outputs, rates] = dPUdT_core_mex(t, PU, mex_vals, Vums, LXBpivot)
 *
 * Inputs:
 *   prhs[0] = t          (1×1 double, ignored — no explicit t-dependence)
 *   prhs[1] = PU         ((2*ss+7)×1 double, state vector; ss from mv[12])
 *   prhs[2] = mex_vals   (187×1 double, packed by packMexVals.m; s_grid padded to 60)
 *   prhs[3] = Vums       (1×1 double, current sarcomere shortening velocity, um/s)
 *   prhs[4] = LXBpivot   (1×1 double, current pivot position, um)
 *
 * Outputs (nlhs controls how many are returned):
 *   plhs[0] = f          ((2*ss+7)×1 double, dPU/dt; ss from mv[12])
 *   plhs[1] = outputs    (1×12 double, Force etc.)
 *   plhs[2] = rates      (1×13 double, scalar rates)
 *
 * Index layout of mex_vals (1-based in MATLAB, 0-based here):
 *   0-42  : scalar params (see IDX_* below)
 *   43-102: strain grid s[60]
 *   103-106: pchip R1D breaks (4 values, n=3 segs)
 *   107-118: pchip R1D coefs (3×4, col-major = 12 values)
 *   119-123: pchip pwsd breaks (5 values, n=4 segs)
 *   124-139: pchip pwsd coefs (4×4, col-major = 16 values)
 *   140-143: pchip R21 breaks (4 values, n=3 segs)
 *   144-155: pchip R21 coefs (3×4, col-major = 12 values)
 *   156-162: pchip pwsd2 breaks (7 values, n=6 segs)
 *   163-186: pchip pwsd2 coefs (6×4, col-major = 24 values)
 *
 * IMPORTANT: this index layout must match packMexVals.m exactly.
 *
 * Compile with:  compileMex.m  (runs mex -O dPUdT_core_mex.c)
 *
 * See also: packMexVals, dPUdT_CombinedTransitions_mex, compileMex
 */

#include "mex.h"
#include <math.h>
#include <string.h>

/* =========================================================
 * Scalar parameter indices into mex_vals (0-based)
 * ========================================================= */
#define IDX_KD          0
#define IDX_K1          1
#define IDX_K_1         2
#define IDX_K2          3
#define IDX_KSTIFF1     4
#define IDX_KSTIFF2     5
#define IDX_KSTIFF1_N   6
#define IDX_KSTIFF2_N   7
#define IDX_KA          8
#define IDX_KAH         9
#define IDX_KAMH        10
#define IDX_DS          11
#define IDX_SS          12
#define IDX_DR          13
#define IDX_DR2         14
#define IDX_L_THICK     15
#define IDX_L_HBARE     16
#define IDX_L_THIN      17
#define IDX_K_PAS       18
#define IDX_LSC0        19
#define IDX_GAMMA       20
#define IDX_KSE         21
#define IDX_EKSE        22
#define IDX_MAX_SLACK   23
#define IDX_MU          24
#define IDX_MU_NEG      25
#define IDX_MU2         26
#define IDX_KSR0        27
#define IDX_SIGMA1      28
#define IDX_SIGMA2      29
#define IDX_KMSR        30
#define IDX_KMSRD       31
#define IDX_SIGMA_SRD1  32
#define IDX_KSRD        33
#define IDX_SIGMA_SRD2  34
#define IDX_KSRD2SR     35
#define IDX_KSR2SRD     36
#define IDX_SLOPE       37
#define IDX_D_ACTIN     38
#define IDX_S_THR_R     39
#define IDX_S_THR_L     40
#define IDX_A1_WIDTH    41
#define IDX_ESTIFF      42

#define N_SCALARS  43
#define N_SS       60   /* fixed strain bins */

/* Pchip segment counts */
#define N_R1D    3
#define N_PWSD   4
#define N_R21    3
#define N_PWSD2  6

/* Pchip layout offsets (0-based) */
#define OFF_S          N_SCALARS                          /* 43 */
#define OFF_R1D_BRK   (OFF_S + N_SS)                     /* 103 */
#define OFF_R1D_COF   (OFF_R1D_BRK + N_R1D + 1)          /* 107 */
#define OFF_PWSD_BRK  (OFF_R1D_COF + N_R1D * 4)          /* 119 */
#define OFF_PWSD_COF  (OFF_PWSD_BRK + N_PWSD + 1)        /* 124 */
#define OFF_R21_BRK   (OFF_PWSD_COF + N_PWSD * 4)        /* 140 */
#define OFF_R21_COF   (OFF_R21_BRK + N_R21 + 1)          /* 144 */
#define OFF_PWSD2_BRK (OFF_R21_COF + N_R21 * 4)          /* 156 */
#define OFF_PWSD2_COF (OFF_PWSD2_BRK + N_PWSD2 + 1)      /* 163 */

#define TOTAL_MV_LEN  (OFF_PWSD2_COF + N_PWSD2 * 4)      /* 187 */

/* Max state vector length (ss=N_SS): 2*N_SS + 7 scalars = 127.
 * Actual runtime length is 2*ss+7 where ss = (int)mv[IDX_SS]. */
#define N_STATE  (2*N_SS + 7)

/* =========================================================
 * Utility: clamp integer to [lo, hi]
 * ========================================================= */
static int clamp_int(int x, int lo, int hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

/* =========================================================
 * peval: evaluate pchip polynomial at one point xi.
 *   brk  — break vector (n_seg+1 values)
 *   cof  — coefficient matrix stored col-major (n_seg × 4)
 *           cof[col*n_seg + row] = C[row][col]
 *   n_seg — number of segments
 * Returns C[si][0]*dx^3 + C[si][1]*dx^2 + C[si][2]*dx + C[si][3]
 * ========================================================= */
static double peval_scalar(const double *brk, const double *cof, int n_seg, double xi)
{
    /* Find segment index (0-based) by linear scan — n_seg <= 6, so trivial */
    int si = 0;
    while (si < n_seg - 1 && xi >= brk[si + 1]) ++si;
    double dx = xi - brk[si];
    /* cof stored col-major: column col, row si → cof[col*n_seg + si] */
    double c0 = cof[0*n_seg + si];
    double c1 = cof[1*n_seg + si];
    double c2 = cof[2*n_seg + si];
    double c3 = cof[3*n_seg + si];
    return ((c0*dx + c1)*dx + c2)*dx + c3;
}

/* =========================================================
 * mexFunction — entry point
 * ========================================================= */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* --- input validation --- */
    if (nrhs != 5)
        mexErrMsgIdAndTxt("dPUdT_core_mex:nrhs", "Expected 5 inputs: t, PU, mex_vals, Vums, LXBpivot");
    if (mxGetNumberOfElements(prhs[2]) != TOTAL_MV_LEN)
        mexErrMsgIdAndTxt("dPUdT_core_mex:mex_vals",
            "mex_vals must have %d elements (got %zu)", TOTAL_MV_LEN,
            mxGetNumberOfElements(prhs[2]));

    /* Runtime ss — may be < N_SS when FV conditions use a narrow Slim range.
     * s_grid in mex_vals is always padded to N_SS entries; only [0..ss-1] are used.
     * Stack arrays (p1, p2, s, …) are allocated for N_SS; loops run only ss times. */
    const double *mv_early = mxGetPr(prhs[2]);
    int ss = (int)mv_early[IDX_SS];
    if (ss <= 0 || ss > N_SS)
        mexErrMsgIdAndTxt("dPUdT_core_mex:ss",
            "ss=%d out of range [1,%d]. Check packMexVals.", ss, N_SS);

    /* --- unpack inputs --- */
    const double *PU      = mxGetPr(prhs[1]);
    const double *mv      = mxGetPr(prhs[2]);
    double vel            = mxGetScalar(prhs[3]);   /* Vums */
    double LXBpivot       = mxGetScalar(prhs[4]);

    /* Scalar params */
    double kd      = mv[IDX_KD];
    double k1      = mv[IDX_K1];
    double k_1     = mv[IDX_K_1];
    double k2      = mv[IDX_K2];
    double kstiff1 = mv[IDX_KSTIFF1];
    double kstiff2 = mv[IDX_KSTIFF2];
    double kstiff1_n = mv[IDX_KSTIFF1_N];
    double kstiff2_n = mv[IDX_KSTIFF2_N];
    double ka      = mv[IDX_KA];
    double kah     = mv[IDX_KAH];
    double kamh    = mv[IDX_KAMH];
    double dS      = mv[IDX_DS];
    double L_thick = mv[IDX_L_THICK];
    double L_hbare = mv[IDX_L_HBARE];
    double L_thin  = mv[IDX_L_THIN];
    double k_pas   = mv[IDX_K_PAS];
    double Lsc0    = mv[IDX_LSC0];
    double gam     = mv[IDX_GAMMA];
    double kSE     = mv[IDX_KSE];
    double ekSE    = mv[IDX_EKSE];
    double MaxSlackNeg = mv[IDX_MAX_SLACK];
    double mu      = mv[IDX_MU];
    double mu_neg  = mv[IDX_MU_NEG];
    double mu2     = mv[IDX_MU2];
    double ksr0    = mv[IDX_KSR0];
    double sigma1  = mv[IDX_SIGMA1];
    double sigma2  = mv[IDX_SIGMA2];
    double kmsr    = mv[IDX_KMSR];
    double kmsrd   = mv[IDX_KMSRD];
    double sigma_srd1 = mv[IDX_SIGMA_SRD1];
    double ksrd    = mv[IDX_KSRD];
    double sigma_srd2 = mv[IDX_SIGMA_SRD2];
    double ksrd2sr = mv[IDX_KSRD2SR];
    double ksr2srd = mv[IDX_KSR2SRD];
    double slope   = mv[IDX_SLOPE];
    double dr      = mv[IDX_DR];
    double dr2     = mv[IDX_DR2];
    double d_actin = mv[IDX_D_ACTIN];
    double s_thr_R = mv[IDX_S_THR_R];
    double s_thr_L = mv[IDX_S_THR_L];
    double Ax      = mv[IDX_A1_WIDTH];
    double estiff  = mv[IDX_ESTIFF];

    /* Strain grid */
    const double *s_grid = mv + OFF_S;

    /* Pchip pointers */
    const double *b_R1D   = mv + OFF_R1D_BRK;
    const double *c_R1D   = mv + OFF_R1D_COF;
    const double *b_pwsd  = mv + OFF_PWSD_BRK;
    const double *c_pwsd  = mv + OFF_PWSD_COF;
    const double *b_R21   = mv + OFF_R21_BRK;
    const double *c_R21   = mv + OFF_R21_COF;
    const double *b_pwsd2 = mv + OFF_PWSD2_BRK;
    const double *c_pwsd2 = mv + OFF_PWSD2_COF;

    /* --- unpack state vector (0-based: p1[0..ss-1], p2[ss..2ss-1], scalars after) --- */
    double p1[N_SS], p2[N_SS];
    int k;
    for (k = 0; k < ss; ++k) {
        p1[k] = PU[k]      < 0.0 ? 0.0 : PU[k];
        p2[k] = PU[ss + k] < 0.0 ? 0.0 : PU[ss + k];
    }
    double P_SR  = PU[2*ss];      /* MATLAB: 2*ss+1 */
    /* NP = PU[2*ss+1] — always 0, skip */
    double SL    = PU[2*ss + 2];
    double LSE   = PU[2*ss + 3];
    double PD    = PU[2*ss + 4];
    double P_SRD = PU[2*ss + 5];
    /* x_dash = PU[2*ss+6] — UseMaxwellDashpot=0, dx_dash_dt=0 */

    /* --- Overlap --- */
    double L_T_HS1 = (L_thick*0.5 < SL*0.5) ? L_thick*0.5 : SL*0.5;
    double tmp1    = (SL-LSE)*0.5 - ((SL-LSE) - L_thin);
    double tmp2    = L_hbare*0.5;
    double L_T_HS2 = (tmp1 > tmp2) ? tmp1 : tmp2;
    double L_ov    = L_T_HS1 - L_T_HS2;
    double N_overlap = L_ov * 2.0 / (L_thick - L_hbare);

    /* --- Strain shift (LegacyStrainFlipping=0) --- */
    double shift = (-(SL - LSE) + LXBpivot) * 0.5;
    double s[N_SS];
    for (k = 0; k < ss; ++k)
        s[k] = s_grid[k] - shift;

    /* --- Force (UseNegativeKstiff=1) --- */
    double p1_0 = 0.0, p2_0 = 0.0;
    double p1_1_pos = 0.0, p1_1_neg = 0.0;
    double p2_1_pos = 0.0, p2_1_neg = 0.0;
    double sdr[N_SS];
    for (k = 0; k < ss; ++k) {
        double sk = s[k];
        double sdk = sk + dr;
        sdr[k] = (sdk >= 0.0 ? 1.0 : -1.0) * pow(fabs(sdk), estiff);
        p1_0 += p1[k];
        p2_0 += p2[k];
        p1_1_pos += sk * p1[k] * (sk >= 0.0 ? 1.0 : 0.0);
        p1_1_neg += sk * p1[k] * (sk <= 0.0 ? 1.0 : 0.0);
        p2_1_pos += sdr[k] * p2[k] * (sk >= -dr ? 1.0 : 0.0);
        p2_1_neg += sdr[k] * p2[k] * (sk <  -dr ? 1.0 : 0.0);
    }
    p1_0 *= dS; p2_0 *= dS;
    p1_1_pos *= dS; p1_1_neg *= dS;
    p2_1_pos *= dS; p2_1_neg *= dS;
    double p1_1 = p1_1_neg + p1_1_pos;
    double p2_1 = p2_1_pos + p2_1_neg;

    double F_active = kstiff2*p2_1_pos + kstiff2_n*p2_1_neg
                    + kstiff1*p1_1_pos + kstiff1_n*p1_1_neg;

    /* PT */
    double PT = 1.0 - (p1_0 + p2_0 + PD + P_SR + P_SRD);
    if (PT < 0.0) PT = 0.0;

    /* Passive */
    double arg_pas = (SL - LSE) - Lsc0;
    double F_passive = (arg_pas > 0.0) ? k_pas * pow(arg_pas, gam) : 0.0;

    double F_total = F_active + F_passive;

    /* Serial stiffness */
    double kSE_LSE = kSE * LSE;
    double kSE_LSE_clamped = (kSE_LSE > MaxSlackNeg) ? kSE_LSE : MaxSlackNeg;
    double Force = (LSE >= 0.0 ? 1.0 : -1.0) * pow(fabs(kSE_LSE_clamped), ekSE);
    double velHS = (Force - F_total) / ((Force - F_total > 0.0) ? mu : mu_neg);
    double dLSEdt = vel - velHS;
    Force += mu2 * vel;

    /* --- Transition rates: inline pchip eval --- */
    double R1D[N_SS], R12[N_SS], R21[N_SS], R2[N_SS];
    double sum_R1D = 0.0, sum_R12 = 0.0, sum_R21 = 0.0, sum_R2 = 0.0;
    for (k = 0; k < ss; ++k) {
        double sk = s[k];
        double s2k = sk + dr2 - dr;
        R1D[k] = kd  * p1[k] * peval_scalar(b_R1D,  c_R1D,  N_R1D,  sk);
        R12[k] = k1  * p1[k] * peval_scalar(b_pwsd,  c_pwsd, N_PWSD, sk);
        R21[k] = k_1 * p2[k] * peval_scalar(b_R21,   c_R21,  N_R21,  sk);
        R2[k]  = k2  * p2[k] * peval_scalar(b_pwsd2, c_pwsd2, N_PWSD2, s2k);
        sum_R1D += R1D[k];
        sum_R12 += R12[k];
        sum_R21 += R21[k];
        sum_R2  += R2[k];
    }

    /* --- UseA2AttachmentShift --- */
    double slope_over_dr = slope / dr;
    double dp2_RAm[N_SS];
    double dp2_RAL[N_SS], dp2_RAR[N_SS];
    for (k = 0; k < ss; ++k) { dp2_RAL[k] = 0.0; dp2_RAR[k] = 0.0; }

    for (k = 0; k < ss; ++k) {
        double sk = s[k];
        double ra_k_pos = slope_over_dr * (sk - s_thr_R);
        double ra_k_neg = -slope_over_dr * (sk + s_thr_L);
        double RA_K = (ra_k_pos > 0.0 ? ra_k_pos : 0.0)
                    + (ra_k_neg > 0.0 ? ra_k_neg : 0.0);
        dp2_RAm[k] = p2[k] * RA_K;

        if (dp2_RAm[k] == 0.0) continue;  /* skip zero contribution */

        double s_tgt = 0.0;
        if (sk > s_thr_R)       s_tgt = sk - d_actin;
        else if (sk < -s_thr_L) s_tgt = sk + d_actin;
        /* else s_tgt = 0 and dp2_RAm[k] should be 0 (already handled above) */

        /* Bin index (0-based): use shifted s[0], matching MATLAB's s(1) */
        double frac = (s_tgt - s[0]) / dS;
        int L_idx = (int)floor(frac);
        if (L_idx < 0) L_idx = 0;
        if (L_idx > ss - 2) L_idx = ss - 2;
        int R_idx = L_idx + 1;

        double w_R = (s_tgt - s[L_idx]) / (s[R_idx] - s[L_idx]);
        double w_L = 1.0 - w_R;

        dp2_RAL[L_idx] += dp2_RAm[k] * w_L;
        dp2_RAR[R_idx] += dp2_RAm[k] * w_R;
    }

    /* --- SRX dynamics --- */
    double F_SR = F_total;  /* UsePassiveForSR=0, useHalfActiveForSR=0 */
    double P_SR_pos  = (P_SR  > 0.0) ? P_SR  : 0.0;
    double P_SRD_pos = (P_SRD > 0.0) ? P_SRD : 0.0;

    double RTD    = kah   * PT;
    double RDT    = kamh  * PD;
    double RD1    = ka    * PD * N_overlap;  /* f_lattice=1 */

    double RSRD2PD = kmsrd  * exp(F_SR / sigma_srd1) * P_SRD_pos;
    double RPD2SRD = ksrd   * exp(-F_SR / sigma_srd2) * PD;
    double RSR2PT  = kmsr   * exp(F_SR / sigma1) * P_SR_pos;
    double RPT2SR  = ksr0   * exp(-F_SR / sigma2) * PT;
    double RSRD2SR = ksrd2sr * P_SRD_pos;
    double RSR2SRD = ksr2srd * P_SR_pos;

    double dU_SRD = RSR2SRD - RSRD2SR - RSRD2PD + RPD2SRD;
    double dU_SR  = -RSR2PT + RPT2SR + RSRD2SR - RSR2SRD;

    /* --- Governing flows --- */
    double dPD = RSRD2PD - RPD2SRD + RTD - RDT - RD1 + sum_R1D * dS;
    double dp1[N_SS], dp2[N_SS];
    for (k = 0; k < ss; ++k) {
        dp1[k] = -R1D[k] - R12[k] + R21[k];
        dp2[k] = R12[k] - R21[k] - R2[k] + dp2_RAR[k] + dp2_RAL[k] - dp2_RAm[k];
    }

    /* --- Attachment: triangular window, halfspan=3 --- */
    int halfspan = (int)ceil(Ax / dS);
    double s0_frac = 1.0 + (-s[0] / dS);  /* fractional index of zero strain (1-based) */
    double att[N_SS];
    double att_sum = 0.0;
    for (k = 0; k < ss; ++k) att[k] = 0.0;

    int jlo = (int)floor(s0_frac - halfspan) - 1;   /* convert to 0-based */
    int jhi = (int)ceil(s0_frac + halfspan) - 1;
    if (jlo < 0) jlo = 0;
    if (jhi >= ss) jhi = ss - 1;

    for (k = jlo; k <= jhi; ++k) {
        double dist = fabs(s[k]);
        double w = 1.0 - dist / Ax;
        if (w > 0.0) { att[k] = w; att_sum += w; }
    }
    if (att_sum < 1.0) att_sum = 1.0;
    double att_scale = RD1 / (att_sum * dS);
    for (k = jlo; k <= jhi; ++k)
        dp1[k] += att[k] * att_scale;

    /* --- Assemble output f ((2*ss+7)×1) --- */
    plhs[0] = mxCreateDoubleMatrix(2*ss + 7, 1, mxREAL);
    double *f = mxGetPr(plhs[0]);
    for (k = 0; k < ss; ++k) f[k]      = dp1[k];
    for (k = 0; k < ss; ++k) f[ss + k] = dp2[k];
    f[2*ss]   = dU_SR;
    f[2*ss+1] = 0.0;    /* dNP = 0 */
    f[2*ss+2] = vel;     /* dSL */
    f[2*ss+3] = dLSEdt;
    f[2*ss+4] = dPD;
    f[2*ss+5] = dU_SRD;
    f[2*ss+6] = 0.0;    /* dx_dash_dt = 0 */

    if (nlhs < 2) return;

    /* --- outputs (1×12) --- */
    plhs[1] = mxCreateDoubleMatrix(1, 12, mxREAL);
    double *outputs = mxGetPr(plhs[1]);
    outputs[0]  = Force;
    outputs[1]  = F_active;
    outputs[2]  = F_passive;
    outputs[3]  = N_overlap;
    outputs[4]  = p1_0;
    outputs[5]  = p2_0;
    outputs[6]  = p1_1;
    outputs[7]  = p2_1;
    outputs[8]  = PT;
    outputs[9]  = 0.0;   /* F_Maxwell = 0 */
    outputs[10] = 1.0;   /* f_lattice = 1 */
    outputs[11] = 1.0;   /* f_saturation = 1 */

    if (nlhs < 3) return;

    /* --- rates (1×16) --- matches MATLAB layout:
         [RTD, RDT, RD1, sum(R1D)*dS, sum(R12)*dS, sum(R21)*dS, sum(R2)*dS,
          sum(XB_Ripped)*dS, RSR2PT, RPT2SR, RSRD2PD, RPD2SRD, RSR2SRD,
          RSRD2SR, RT2, sum(R2D)*dS]                                        */
    plhs[2] = mxCreateDoubleMatrix(1, 16, mxREAL);
    double *rates = mxGetPr(plhs[2]);
    rates[0]  = RTD;
    rates[1]  = RDT;
    rates[2]  = RD1;
    rates[3]  = sum_R1D  * dS;   /* sum(R1D)*dS  */
    rates[4]  = sum_R12  * dS;   /* sum(R12)*dS  */
    rates[5]  = sum_R21  * dS;   /* sum(R21)*dS  */
    rates[6]  = sum_R2   * dS;   /* sum(R2)*dS   */
    rates[7]  = 0.0;              /* sum(XB_Ripped)*dS = 0 (k2rip=0) */
    rates[8]  = RSR2PT;
    rates[9]  = RPT2SR;
    rates[10] = RSRD2PD;
    rates[11] = RPD2SRD;
    rates[12] = RSR2SRD;
    rates[13] = RSRD2SR;
    rates[14] = 0.0;              /* RT2 = 0 (UseA2Reattaching=0) */
    rates[15] = 0.0;              /* sum(R2D)*dS = 0 (UseA2MechanicalRecocking=0) */
}
