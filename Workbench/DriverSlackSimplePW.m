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

%% Apply simplified curves
params0.PieceWiseStrainDepR1DX      = R1D_X_simp;
params0.PieceWiseStrainDepR1DParams = R1D_P_simp;
params0.PieceWiseStrainDep2X        = R2_X_simp;
params0.PieceWiseStrainDep2Params   = R2_P_simp;

%% Run slack experiment
figure(201); clf;
tic
RunBakersExp;
toc

if params0.EvalFeatures
    plotFeatures(features_data, features_model, [], params0.fn);
end

%% ============================================================
%  Physiologically grounded shapes v2
%  (rescaled to min~1, simplified left plateau, + experimental R21)
% ============================================================
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
%   min=1 at s=0, kd_v2 = kd_orig × R1D(s=0).  Drops spurious bump at +4.2 nm.
%   kd_v2 ≈ 122 s⁻¹ matches literature for cardiac α-MHC (Greenberg 2016: 88–121 s⁻¹).
%
% R2_v2 — reparametrisation + left-plateau simplification, effectively same rates.
%   min=1 at s=0, k2_v2 = k2_orig × R2(s≈0).
%   k2_v2 ≈ 60 s⁻¹ = ADP release at optimal loaded P2 (eff strain ≈ 11 nm).
%   Left flank: δ_left = 4.45/ln(44) ≈ 1.2 nm — consistent with Bell δ = 1.1 nm.
%   Right flank: δ_right ≈ 2.7 nm (implies kstiff_head ≈ 1 pN/nm; within range).
%   Single left clamp replaces 3-pt cluster of full version (PCHIP artifact risk).
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
%   Characteristic length ~6 nm (conservative; Capitanio myosin-V d ≈ 1.2 nm).
%   R2/R21 >> 1 at s < −5 nm → R21 change negligible during fast shortening.
%   Effect near isometric (s ≈ 0 to +5 nm): modest increase in power-stroke commitment.

% --- compute normalization factors from current (simplified) shapes in params0 ---
R1D_at_0   = ppval(pchip(params0.PieceWiseStrainDepR1DX,    params0.PieceWiseStrainDepR1DParams), 0);
R2_at_eval = ppval(pchip(params0.PieceWiseStrainDep2X,      params0.PieceWiseStrainDep2Params),   params0.dr2 - params0.dr);

kd_v2 = params0.kd * R1D_at_0;    % absorb R1D(s=0) multiplier into base rate
k2_v2 = params0.k2 * R2_at_eval;  % absorb R2(s≈0) multiplier into base rate

fprintf('v2 base rates: kd_v2 = %.1f s-1 (was %.1f), k2_v2 = %.1f s-1 (was %.1f)\n', ...
    kd_v2, params0.kd, k2_v2, params0.k2);

% R1D_v2: rescaled breakpoints (same X as simplified, P divided by R1D_at_0)
R1D_X_v2 = params0.PieceWiseStrainDepR1DX;
R1D_P_v2 = params0.PieceWiseStrainDepR1DParams / R1D_at_0;
% Shape: min=1 at s=0; left clamp ~956×; right ~32× at s=+dr (P1 at full power-stroke length)

% R2_v2: rescaled breakpoints (same X as simplified, P divided by R2_at_eval)
R2_X_v2 = params0.PieceWiseStrainDep2X;
R2_P_v2 = params0.PieceWiseStrainDep2Params / R2_at_eval;
% Shape: min=1 at s≈0; left clamp ~107×; right ~43× at s=+10 nm (Bell δ_eff ≈ 2.7 nm)

% R12_v2: unchanged (literature validates peak at negative strain)
R12_X_v2 = params0.PieceWiseStrainDepX;
R12_P_v2 = params0.PieceWiseStrainDepParams;

% R21_v2: new strain-dependent shape (EXPERIMENTAL)
%   At s = −10.5 nm (P2 natural length): easy reversal   → 3–4×
%   At s = 0 (reference):                                 → 1×
%   At s = +10 nm (highly loaded):                        → 0.2×
R21_X_v2 = [-0.050,  -0.010,  0,    0.010,  0.050];
R21_P_v2 = [  4.0,     3.0,   1.0,  0.2,    0.05];

% --- 4-panel comparison plot ---
s_ev  = linspace(-0.015, 0.020, 1000)';
s_nm  = s_ev * 1000;
dr_s  = params0.dr2 - params0.dr;   % R2 evaluation offset (≈ −0.5 nm)
dr_nm = params0.dr * 1000;

figure(202); clf;
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% R1D — show effective rates (kd × multiplier) so v2 overlays simplified exactly
nexttile; hold on;
plot(s_nm, params0.kd   * ppval(pchip(R1D_X_full, R1D_P_full), s_ev), 'b-',  'LineWidth', 1.5, 'DisplayName', 'full (7 pt, kd=29)');
plot(s_nm, params0.kd   * ppval(pchip(R1D_X_simp, R1D_P_simp), s_ev), 'r--', 'LineWidth', 1.5, 'DisplayName', 'simplified (5 pt, kd=29)');
plot(s_nm, kd_v2        * ppval(pchip(R1D_X_v2,   R1D_P_v2),   s_ev), 'g-',  'LineWidth', 2,   'DisplayName', sprintf('v2 (min=1, kd=%.0f)',kd_v2));
xline(0, 'k--', 'attach', 'FontSize', 8);
set(gca, 'YScale', 'log'); xlim([-15 20]); ylim([10 2e5]); grid on; box on;
xlabel('Strain s (nm)'); ylabel('Effective kd rate (s⁻¹)');
title('R1D  (P1→PD detachment)');
legend('Location', 'southeast', 'FontSize', 7);
text(-14, 500, 'v2 overlays simplified — same effective rates', 'FontSize', 7, 'Color', [0 0.5 0]);

% R12 — unchanged; just annotate
nexttile; hold on;
plot(s_nm, ppval(pchip(params0.PieceWiseStrainDepX, params0.PieceWiseStrainDepParams), s_ev), ...
    'b-', 'LineWidth', 2, 'DisplayName', 'current = v2 (unchanged)');
xline(0,    'k--', 'attach',    'FontSize', 8);
xline(-dr_nm,'b:',  'P2 nat. len.','FontSize', 8);
set(gca, 'YScale', 'log'); xlim([-15 20]); grid on; box on;
xlabel('Strain s (nm)'); ylabel('R12 multiplier');
title('R12  (P1→P2 power stroke, UNCHANGED)');
text(-14, 0.1, sprintf('Peak at −8.3 nm = Piazzesi 2007 (6000 s⁻¹ at low load)', []), 'FontSize', 7);
legend('Location', 'northeast', 'FontSize', 7);

% R2 — show effective rates (k2 × multiplier)
nexttile; hold on;
plot(s_nm, params0.k2   * ppval(pchip(R2_X_full, R2_P_full), s_ev + dr_s), 'b-',  'LineWidth', 1.5, 'DisplayName', 'full (10 pt, k2=128)');
plot(s_nm, params0.k2   * ppval(pchip(R2_X_simp, R2_P_simp), s_ev + dr_s), 'r--', 'LineWidth', 1.5, 'DisplayName', 'simplified (5 pt, k2=128)');
plot(s_nm, k2_v2        * ppval(pchip(R2_X_v2,   R2_P_v2),   s_ev + dr_s), 'g-',  'LineWidth', 2,   'DisplayName', sprintf('v2 (min=1, k2=%.0f)',k2_v2));
xline(0, 'k--', 'attach', 'FontSize', 8);
set(gca, 'YScale', 'log'); xlim([-15 20]); grid on; box on;
xlabel('Strain s (nm)  [evaluated at s − 0.5 nm]'); ylabel('Effective k2 rate (s⁻¹)');
title('R2  (P2→PD, ADP release)');
legend('Location', 'northeast', 'FontSize', 7);
text(-14, 200, sprintf('δ_{left}≈1.2 nm (Bell), k2@s=0: %.0f s⁻¹', k2_v2), 'FontSize', 7, 'Color', [0 0.5 0]);

% R21 — multiplier comparison (current flat vs new strain-dependent)
nexttile; hold on;
plot(s_nm, ppval(pchip(params0.PieceWiseStrainDepR21X, params0.PieceWiseStrainDepR21Params), s_ev), ...
    'b-', 'LineWidth', 1.5, 'DisplayName', 'current (flat ≈ 1×)');
plot(s_nm, ppval(pchip(R21_X_v2, R21_P_v2), s_ev), ...
    'g-', 'LineWidth', 2, 'DisplayName', 'v2 (strain-dep, EXPERIMENTAL)');
xline(0,    'k--', 'attach',      'FontSize', 8);
xline(-dr_nm,'b:', 'P2 nat. len.','FontSize', 8);
set(gca, 'YScale', 'log'); xlim([-15 20]); ylim([0.01 10]); grid on; box on;
xlabel('Strain s (nm)'); ylabel('R21 multiplier');
title('R21  (P2→P1 reverse stroke, NEW)');
legend('Location', 'northeast', 'FontSize', 7);
text(-14, 4, 'Capitanio 2012: d≈1.2 nm (myosin-V)', 'FontSize', 7, 'Color', [0 0.5 0]);

sgtitle('Physiologically grounded v2 shapes vs current/simplified', 'FontSize', 11, 'FontWeight', 'bold');

%% Apply v2 shapes and run slack
params0.kd                          = kd_v2;
params0.PieceWiseStrainDepR1DX      = R1D_X_v2;
params0.PieceWiseStrainDepR1DParams = R1D_P_v2;

params0.k2                          = k2_v2;
params0.PieceWiseStrainDep2X        = R2_X_v2;
params0.PieceWiseStrainDep2Params   = R2_P_v2;

% R12 unchanged — no assignment needed
params0.PieceWiseStrainDepR21X      = R21_X_v2;
params0.PieceWiseStrainDepR21Params = R21_P_v2;

figure(203); clf;
params0.RunForceVelocity = true;
params0.FV_velocities = [0 -0.5000   -1.0000   -2.0000   -4.0000];
disp('=== v2 shapes: slack experiment ===');
tic; RunBakersExp; toc

if params0.EvalFeatures
    plotFeatures(features_data, features_model, [], params0.fn);
end
