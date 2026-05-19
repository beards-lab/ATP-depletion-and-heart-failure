% DriverSlackSimplePW.m
% Compare full vs simplified piecewise strain-dependent rates (R1D 7→5, R2 10→5).
% Uses justPlotStateTransitionsFlag to visualise the curves side-by-side,
% then runs the slack experiment with the simplified parameterisation.

clear; close all; clc;

%% Load latest optimised params
params0 = getParams();
ModelOptParams_SRXTD_iter_9

params0.RunForceVelocity    = false;
params0.RunKtr              = false;
params0.RunStairs           = false;
params0.RunSlack            = true;
params0.RunForceLengthEstim = false;
params0.RunSlackSegments    = 'All';
params0.EvalFeatures        = true;
params0.PlotEachSeparately  = true;
params0.BreakOnODEUnstable  = true;
params0.MaxRunTime          = 600;

%% Simplified piecewise definitions
%
% R1D (kd multiplier, p1→PD): 7→5 pts
%   dropped: X=-0.015 (P=4000, redundant with far-left clamp)
%            X=+0.00417 (P=8.6, smooth PCHIP covers 0→0.01 adequately)
R1D_X_full  = [-0.05  -0.015  -0.00366  0       0.00417  0.01048  0.1];
R1D_P_full  = [4000    4000    43.891   4.187   8.578   135.308  2000];
R1D_X_simp  = [-0.05          -0.00366  0                0.01048  0.1];
R1D_P_simp  = [4000             43.891  4.187           135.308  2000];

% R2 (k2 multiplier, p2→PD): 10→5 pts
%   collapsed left plateau (X=-0.01,-0.0105,-0.0085, all P≈50) into single far-left clamp
%   kept left corner, center, right corner, far-right clamp
R2_X_full   = [-0.01  -0.0105  -0.0085  -0.00555  -0.00427  0       0.00618  0.01  0.0233  0.1];
R2_P_full   = [50      50       50.608   1.133     0.968     0.467   1        20    50      50];
R2_X_simp   = [-0.01                    -0.00555            0       0.00618  0.01          0.1];
R2_P_simp   = [50                        1.133              0.467   1        20            50];

%% Figure: compare full vs simplified curves
s_eval = linspace(-0.06, 0.12, 500);

figure(200); clf;
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(s_eval, ppval(pchip(R1D_X_full, R1D_P_full), s_eval), 'b-',  'LineWidth', 2, 'DisplayName', 'full (7 pts)');
plot(s_eval, ppval(pchip(R1D_X_simp, R1D_P_simp), s_eval), 'r--', 'LineWidth', 2, 'DisplayName', 'simplified (5 pts)');
plot(R1D_X_full, R1D_P_full, 'bo', 'MarkerFaceColor', 'b');
plot(R1D_X_simp, R1D_P_simp, 'rs', 'MarkerFaceColor', 'r');
set(gca, 'YScale', 'log'); xlim([-0.06 0.12]); grid on; box on;
xlabel('Strain'); ylabel('R1D multiplier'); title('R1D  (kd, p1→PD)');
legend('Location', 'best');

nexttile; hold on;
plot(s_eval, ppval(pchip(R2_X_full, R2_P_full), s_eval), 'b-',  'LineWidth', 2, 'DisplayName', 'full (10 pts)');
plot(s_eval, ppval(pchip(R2_X_simp, R2_P_simp), s_eval), 'r--', 'LineWidth', 2, 'DisplayName', 'simplified (5 pts)');
plot(R2_X_full, R2_P_full, 'bo', 'MarkerFaceColor', 'b');
plot(R2_X_simp, R2_P_simp, 'rs', 'MarkerFaceColor', 'r');
set(gca, 'YScale', 'log'); xlim([-0.06 0.12]); grid on; box on;
xlabel('Strain'); ylabel('R2 multiplier'); title('R2  (k2, p2→PD)');
legend('Location', 'best');

sgtitle('Full vs simplified piecewise — shape comparison');

%% Skip simplified run (unstable); go straight to v2
% The simplified shapes (5-pt versions) cause ODE instability.
% We will build v2 (rescaled) and v3 (with ATP reduction + FV shoulder) directly from full shapes instead.

%% ============================================================
%  Physiologically grounded shapes v2 + v3
%  (rescaled to min~1, + experimental R21, + ATP reduction)
% ============================================================
%
% Build v2 and v3 directly from full original shapes (more stable than simplified)
%
% STRAIN CONVENTION (for all four transitions):
%   s = 0        : attachment point — P1 binds actin here
%   s > 0        : spring stretched, force-generating
%   s < 0        : spring compressed, resistive / shortening regime
%   P2 eff strain: s + dr  (dr = 10.5 nm = power stroke distance)
%   P2 nat. len. : s = -dr ≈ -10.5 nm (zero P2 tension)
%
% LITERATURE BASIS (mouse cardiac trabeculae, α-MHC, ~20 °C):
%   Greenberg 2016 PNAS  https://doi.org/10.1073/pnas.1516598113
%     δ = 1.1 nm for force-dependent ADP release in cardiac myosin
%   Piazzesi 2007 Cell   https://doi.org/10.1016/j.cell.2007.09.029
%     Power stroke rate ~6000 s⁻¹ at low load (x-ray + mechanics, frog fibre)
%   Capitanio 2012 NatStructMolBiol https://doi.org/10.1038/nsmb.1820
%     Myosin-V reverse power stroke: d ≈ 1.2 nm sensitivity distance
%   Huxley & Simmons 1971 Nature
%     Power stroke rate increases with decreasing load (negative strain)
%
% R1D_v2 — pure reparametrisation, IDENTICAL effective rates, better interpretability.
%   min=1 at s=0, kd_v2 = kd_orig × R1D(s=0).
%   kd_v2 ≈ 122 s⁻¹ matches literature for cardiac α-MHC (Greenberg 2016: 88–121 s⁻¹).
%
% R2_v2 — reparametrisation, effectively same rates.
%   min=1 at s=0, k2_v2 = k2_orig × R2(s≈0).
%   k2_v2 ≈ 60 s⁻¹ = ADP release at optimal loaded P2 (eff strain ≈ 11 nm).
%
% R12_v2 — UNCHANGED from current optimized shape.
%   Peak at s = −8.3 nm: thermodynamic bias (compressed P1 closer to P2 natural length).
%   R12_max = 543 × 10.8 = 5864 s⁻¹ ≈ Piazzesi 2007 ~6000 s⁻¹ at low load.
%
% R21_v2 — NEW strain dependence (EXPERIMENTAL).
%   Mechanical basis: reverse power stroke works against load F = kstiff2 × (s + dr).
%   At s = −dr = −10.5 nm : P2 at natural length, no load     → 3× easier reversal
%   At s = 0               : P2 loaded by dr, reference        → 1×
%   At s = +10 nm          : P2 highly loaded                  → 0.2× (suppressed)

% --- compute normalization factors from CURRENT (original) shapes in params0 ---
R1D_at_0   = ppval(pchip(params0.PieceWiseStrainDepR1DX,    params0.PieceWiseStrainDepR1DParams), 0);
R2_at_eval = ppval(pchip(params0.PieceWiseStrainDep2X,      params0.PieceWiseStrainDep2Params),   params0.dr2 - params0.dr);

kd_v2 = params0.kd * R1D_at_0;    % absorb R1D(s=0) multiplier into base rate
k2_v2 = params0.k2 * R2_at_eval;  % absorb R2(s≈0) multiplier into base rate

fprintf('v2 base rates: kd_v2 = %.1f s-1 (was %.1f), k2_v2 = %.1f s-1 (was %.1f)\n', ...
    kd_v2, params0.kd, k2_v2, params0.k2);

% R1D_v2: rescaled breakpoints (same X as original full, P divided by R1D_at_0)
R1D_X_v2 = params0.PieceWiseStrainDepR1DX;
R1D_P_v2 = params0.PieceWiseStrainDepR1DParams / R1D_at_0;

% R2_v2: rescaled breakpoints (same X as original full, P divided by R2_at_eval)
R2_X_v2 = params0.PieceWiseStrainDep2X;
R2_P_v2 = params0.PieceWiseStrainDep2Params / R2_at_eval;

% R12_v2: unchanged (literature validates peak at negative strain)
R12_X_v2 = params0.PieceWiseStrainDepX;
R12_P_v2 = params0.PieceWiseStrainDepParams;

% R21_v2: new strain-dependent shape (EXPERIMENTAL)
R21_X_v2 = [-0.050,  -0.010,  0,    0.010,  0.050];
R21_P_v2 = [  4.0,     3.0,   1.0,  0.2,    0.05];

% --- 4-panel comparison plot (v2 vs original) ---
s_ev  = linspace(-0.015, 0.020, 1000)';
s_nm  = s_ev * 1000;
dr_s  = params0.dr2 - params0.dr;
dr_nm = params0.dr * 1000;

figure(202); clf;
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(s_nm, params0.kd * ppval(pchip(params0.PieceWiseStrainDepR1DX, params0.PieceWiseStrainDepR1DParams), s_ev), ...
    'b-',  'LineWidth', 1.5, 'DisplayName', 'original (full)');
plot(s_nm, kd_v2      * ppval(pchip(R1D_X_v2, R1D_P_v2), s_ev), ...
    'g-',  'LineWidth', 2,   'DisplayName', sprintf('v2 (min=1, kd=%.0f)',kd_v2));
xline(0, 'k--'); set(gca, 'YScale', 'log'); xlim([-15 20]); grid on; box on;
xlabel('Strain s (nm)'); ylabel('Effective kd rate (s⁻¹)');
title('R1D  (P1→PD detachment)');
legend('Location', 'southeast', 'FontSize', 7);

nexttile; hold on;
plot(s_nm, ppval(pchip(params0.PieceWiseStrainDepX, params0.PieceWiseStrainDepParams), s_ev), ...
    'b-', 'LineWidth', 2, 'DisplayName', 'unchanged (same as v2)');
xline(0,'k--'); xline(-dr_nm,'b:','P2 nat. len.');
set(gca, 'YScale', 'log'); xlim([-15 20]); grid on; box on;
xlabel('Strain s (nm)'); ylabel('R12 multiplier');
title('R12  (P1→P2 power stroke, UNCHANGED)');
legend('Location', 'northeast', 'FontSize', 7);

nexttile; hold on;
plot(s_nm, params0.k2 * ppval(pchip(params0.PieceWiseStrainDep2X, params0.PieceWiseStrainDep2Params), s_ev + dr_s), ...
    'b-',  'LineWidth', 1.5, 'DisplayName', 'original (full)');
plot(s_nm, k2_v2      * ppval(pchip(R2_X_v2, R2_P_v2), s_ev + dr_s), ...
    'g-',  'LineWidth', 2,   'DisplayName', sprintf('v2 (min=1, k2=%.0f)',k2_v2));
xline(0, 'k--'); set(gca, 'YScale', 'log'); xlim([-15 20]); grid on; box on;
xlabel('Strain s (nm)'); ylabel('Effective k2 rate (s⁻¹)');
title('R2  (P2→PD, ADP release)');
legend('Location', 'northeast', 'FontSize', 7);

nexttile; hold on;
plot(s_nm, ppval(pchip(params0.PieceWiseStrainDepR21X, params0.PieceWiseStrainDepR21Params), s_ev), ...
    'b-', 'LineWidth', 1.5, 'DisplayName', 'original (flat ≈ 1×)');
plot(s_nm, ppval(pchip(R21_X_v2, R21_P_v2), s_ev), ...
    'g-', 'LineWidth', 2, 'DisplayName', 'v2 (strain-dep, NEW)');
xline(0,'k--'); xline(-dr_nm,'b:','P2 nat. len.'); set(gca, 'YScale', 'log');
xlim([-15 20]); ylim([0.01 10]); grid on; box on;
xlabel('Strain s (nm)'); ylabel('R21 multiplier');
title('R21  (P2→P1 reverse stroke, NEW)');
legend('Location', 'northeast', 'FontSize', 7);

sgtitle('v2 shapes: rescaled + new R21 (no ATP/FV changes yet)', 'FontSize', 11, 'FontWeight', 'bold');

%% Apply v2 shapes and run slack
params0.kd                          = kd_v2;
params0.PieceWiseStrainDepR1DX      = R1D_X_v2;
params0.PieceWiseStrainDepR1DParams = R1D_P_v2;
params0.k2                          = k2_v2;
params0.PieceWiseStrainDep2X        = R2_X_v2;
params0.PieceWiseStrainDep2Params   = R2_P_v2;
params0.PieceWiseStrainDepR21X      = R21_X_v2;
params0.PieceWiseStrainDepR21Params = R21_P_v2;

figure(203); clf;
params0.RunForceVelocity = true;
params0.FV_velocities = [0 -0.5000 -1.0000 -2.0000 -4.0000];
disp('=== v2 shapes: slack + FV experiment ===');
tic; RunBakersExp; toc

if params0.EvalFeatures
    plotFeatures(features_data, features_model, [], params0.fn);
end

%% ============================================================
%  v4: freeze SRX as fixed pool + deepen R2 well at force-generating zone
% ============================================================
%
% Post-critique strategy (see plan file):
%   1. Disable SRX dynamics → no force-dependent SRX↔DRX confounding while
%      we probe the R2-well effect. With both flags off, all SRX fluxes
%      (RSR2PT, RPT2SR, RSRD2SR, RSR2SRD, RSRD2PD, RPD2SRD) are zero, so
%      P_SR and P_SRD stay at their initial values for the whole run.
%   2. Deepen R2 well at X=0 (= s2=0 in P-list frame ≈ s=+0.5 nm physical,
%      i.e. where bulk-P2 mass sits at isometric and force per bridge =
%      kstiff·(s+dr) is at its peak ≈ +11 nm spring extension).
%   3. Add flanking shoulders at X=±0.002 (P=0.5) to widen the well and
%      avoid PCHIP pinching.
%   4. ka, kah, kamh, k1, k_1, kd, k2 base — all UNTOUCHED.
%      (ka÷5 would re-introduce biphasic ktr at ~63 vs ~12 s⁻¹;
%       kah÷5 moves hydrolysis into the ktr band. Both broke v3.)
%
% This is a probe to isolate the R2-well effect, not a final calibration.

% --- 1. Freeze SRX as fixed parked pool
params0.UseSuperRelaxed    = 0;
params0.UseSuperRelaxedADP = 0;
% Split matches v2 internal SRX ratio P_SRD/P_SR = ksr2srd/ksrd2sr ≈ 14.2
params0.SRXT_0 = 0.033;   % ATP-bound parked head (minor)
params0.SRXD_0 = 0.467;   % ADP·Pi-bound parked head (bulk of OFF state)
% Total OFF = 0.5 → cycling pool = 0.5 throughout the run.

% --- 2. R2 well: deepen at X=0, add ±2 nm shoulders
R2_X_v4 = sort([R2_X_v2, -0.002, +0.002]);            % insert shoulder breakpoints
R2_P_v4 = ppval(pchip(R2_X_v2, R2_P_v2), R2_X_v4);    % interp existing values
R2_P_v4(R2_X_v4 == -0.002) = 0.5;   % left shoulder
R2_P_v4(R2_X_v4 ==  0)     = 0.2;   % well bottom (was ≈ 1.0 → 5× deeper)
R2_P_v4(R2_X_v4 ==  0.002) = 0.5;   % right shoulder

params0.PieceWiseStrainDep2X      = R2_X_v4;
params0.PieceWiseStrainDep2Params = R2_P_v4;

% --- 3. Comparison plot
s_ev4  = linspace(-0.015, 0.020, 1000)';
s_nm4  = s_ev4 * 1000;
dr_s4  = params0.dr2 - params0.dr;

figure(206); clf; hold on;
plot(s_nm4, k2_v2 * ppval(pchip(R2_X_v2, R2_P_v2), s_ev4 + dr_s4), ...
    'g-', 'LineWidth', 2, 'DisplayName', sprintf('v2 (k2=%.0f, well min≈1.0)', k2_v2));
plot(s_nm4, k2_v2 * ppval(pchip(R2_X_v4, R2_P_v4), s_ev4 + dr_s4), ...
    'm-', 'LineWidth', 2, 'DisplayName', 'v4 (well min=0.2, shoulders=0.5)');
xline(0, 'k--', 'attach'); xline(-2, 'k:'); xline(+2, 'k:');
set(gca, 'YScale', 'log'); xlim([-15 20]); grid on; box on;
xlabel('Strain s (nm)'); ylabel('Effective k2 rate (s⁻¹)');
title('R2: v2 vs v4 — deepened well at force-generating zone');
legend('Location', 'southeast', 'FontSize', 9);

% --- 4. Run slack + FV with v4
params0.RunForceVelocity = true;
params0.FV_velocities = [0 -0.5 -1 -2 -4];
figure(207); clf;
disp('=== v4: frozen SRX (50% OFF) + R2 well deepened ===');
tic; RunBakersExp; toc

if params0.EvalFeatures
    total_cost = plotFeatures(features_data, features_model, [], params0.fn);
end

sum(total_cost)

%% ============================================================
%  v5: tackle multimodal force redevelopment after slack
% ============================================================
%
% Diagnosis (v3/v4 results were dominated by failure to fit a single-
% exponential force redevelopment): force redevelopment after slack rises
% with at least two distinct time scales. Three candidate causes ranked:
%
%   (A) R12 peak is OFFSET from where new bridges land.
%       Current R12 peak at s = -8.3 nm × 10.8.
%       At s = 0 (where fresh P1 lives after slack-induced reattachment),
%       R12 multiplier ≈ 1-2.
%       → fresh bridges power-stroke slowly (k1 × R12(0) ≈ 800 s⁻¹);
%         bridges that drift to negative s power-stroke fast (≈5860 s⁻¹).
%       Two cohorts → bimodal rise.
%
%   (B) Strain-dep functions are TOO NARROW.
%       R12 transition spans only ~5 nm; R2 minimum is ~1-2 nm wide.
%       PCHIP through clustered breakpoints amplifies peakiness.
%       Cohort behaviour from narrow rates → distinct time scales.
%
%   (C) A1AttachmentWidth = 4 nm — newly attached bridges in tight strain
%       window → cohort dynamics.
%
% Strategy: layered probes. v5a tests SRX initial fractions alone.
% v5b adds R12 broadening. v5c adds A1 kernel broadening.
% We start from the v2 shapes (re-applied) so v4's frozen-SRX / deep-R2
% changes do not contaminate.

% --- restore v2 baseline (undo v4 freezes / well deepening) ---
params0.UseSuperRelaxed    = 1;
params0.UseSuperRelaxedADP = 1;
params0.PieceWiseStrainDep2X      = R2_X_v2;
params0.PieceWiseStrainDep2Params = R2_P_v2;
params0.k2                        = k2_v2;

% --- v5a: SRX initial conditions ---
% Park 20% in SRX-T (ATP-bound parked) and 20% in SRX-D (ADP·Pi parked).
% Total OFF ≈ 40% at t=0; force-dependent kinetics still active.
params0.SRXT_0 = 0.2;
params0.SRXD_0 = 0.2;

% --- v5b: broaden R12 (P1→P2 power stroke) ---
% Original: 5 breakpoints, peak 10.8× at -8.3 nm, falls to ~1× by s=0.
% v5b   : plateau ≈ 6× across s ∈ [-10, 0] nm; smooth decay positive side.
% Rationale: any P1 bridge in the cycling-relevant strain range should be
% able to power-stroke at a similar rate, instead of only those near the
% Huxley-Simmons optimum. This collapses two cohorts into one.
R12_X_v5 = [-1,    -0.010, -0.005, -0.002,  0,      0.005,  0.021];
R12_P_v5 = [10,     6.0,    6.0,    6.0,    6.0,    2.0,    0.015];

% --- v5c: broaden A1 attachment kernel ---
A1_width_v5 = 0.012;   % 4 nm → 12 nm (3× wider)

% --- comparison plot: R12 v2 vs v5b ---
figure(208); clf; hold on;
plot(s_nm, params0.k1 * ppval(pchip(R12_X_v2,   R12_P_v2),   s_ev), ...
    'g-', 'LineWidth', 2, 'DisplayName', sprintf('v2 (peak at -8.3 nm, k1=%.0f)', params0.k1));
plot(s_nm, params0.k1 * ppval(pchip(R12_X_v5,   R12_P_v5),   s_ev), ...
    'm-', 'LineWidth', 2, 'DisplayName', sprintf('v5b (plateau across [-10,0] nm)'));
xline(0, 'k--', 'attach'); xline(-dr_nm, 'b:', 'P2 nat. len.');
set(gca, 'YScale', 'log'); xlim([-15 20]); grid on; box on;
xlabel('Strain s (nm)'); ylabel('Effective k1 rate (s⁻¹)');
title('R12: v2 vs v5b — broadened power-stroke plateau');
legend('Location', 'southeast', 'FontSize', 9);

%% v5a — SRX initial conditions only
disp('=== v5a: SRXT_0 = SRXD_0 = 0.2 (v2 shapes otherwise) ===');
params0.SkipParametersInBounds = true;
params0.PieceWiseStrainDepX      = R12_X_v2;
params0.PieceWiseStrainDepParams = R12_P_v2;
params0.A1AttachmentWidth        = 0.004;   % unchanged
params0.SRXD_0 = 0.1;
params0.SRXT_0 = 0.1;
params0.RunForceVelocity = true;
params0.FV_velocities = -[0, 0.5, 1, 2, 4];
params0.FV_dataset = 'Baker2022';
params0.UseSuperRelaxed = false;
params0.UseSuperRelaxedADP = false;


params0.fn = {
    'FV_fnorm|FV_v|10', 'ktr|2', 'A|50', ...
    'ktr_rmse|0-0.2|.1', ...
    'XTOR[1]|2-8,XTOR_vmax[1]|4-30,SRX_ss[1]|0.1-0.8,attached_ss[1]|0.2-0.6', ...
    't0_crossing|SLdiff|2', ...
    'restretchSlopeStart|0.1',  ...
    'peak1_y|10', 'peak1_dSL|0.2', ...
    'vall_y|10', 'vall_t|0.2', 'peak2|5', ...
    'steady|50', 'vall2_dy|0.1', 'ovrsht_dy|0.1', ...
    'AssertParams|1'...
};

params0.RunSlack = false;
params0.RunForceVelocity = false;
params0.SkipParametersInBounds = true;
params0.ka = 317*100;
params0.k1 = 543*100;
figure(209); clf;
tic; RunBakersExp; toc
%

if params0.EvalFeatures
    cost_v5a = plotFeatures(features_data, features_model, [], params0.fn);
    fprintf('v5a total cost = %.3f\n', sum(cost_v5a));
end

%% v5b — v5a + broadened R12
disp('=== v5b: v5a + R12 broadened plateau ===');
params0.PieceWiseStrainDepX      = R12_X_v5;
params0.PieceWiseStrainDepParams = R12_P_v5;
figure(210); clf;
tic; RunBakersExp; toc
if params0.EvalFeatures
    cost_v5b = plotFeatures(features_data, features_model, [], params0.fn);
    fprintf('v5b total cost = %.3f\n', sum(cost_v5b));
end

%% v5c — v5b + broadened A1 attachment kernel
disp('=== v5c: v5b + A1AttachmentWidth = 12 nm ===');
params0.A1AttachmentWidth = A1_width_v5;
figure(211); clf;
tic; RunBakersExp; toc
if params0.EvalFeatures
    cost_v5c = plotFeatures(features_data, features_model, [], params0.fn);
    fprintf('v5c total cost = %.3f\n', sum(cost_v5c));
end

fprintf('\n=== v5 cost comparison ===\n');
fprintf('  v5a (SRX init only)             : %.3f\n', sum(cost_v5a));
fprintf('  v5b (+R12 broadened)            : %.3f\n', sum(cost_v5b));
fprintf('  v5c (+A1 kernel broadened)      : %.3f\n', sum(cost_v5c));