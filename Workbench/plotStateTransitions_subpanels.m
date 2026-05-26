% Plot all cross-bridge state transitions in a clean tiled layout.
% One transition per panel: forward (solid) and backward (dashed) rates.
% Strain-dependent transitions use s (strain) as x-axis.
% Force-dependent transitions use F_SR (force) as x-axis.

% Create figure and tiledlayout (2 rows × 4 columns = 8 panels)
fig = gcf();
if isempty(fig.Children) || ~isa(fig.Children(1), 'matlab.graphics.layout.TiledChartLayout')
    clf; fig.Position = [100 100 1400 800];
    tlo = tiledlayout(fig, 2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
else
    tlo = fig.Children(1);
    tlo.Children(:).delete();  % clear all tiles
end

% Define force range for force-dependent plots
F_range = linspace(0, 100, 200);
s_range = s;  % strain grid already available

%% Panel 1: Hydrolysis (P_T ↔ P_D)
nexttile; hold on; grid on; box on;
RTD_scalar = params.kah;  % per-unit PT
RDT_scalar = params.kamh; % per-unit PD
plot([s_range(1) s_range(end)], [RTD_scalar RTD_scalar], 'b-', 'LineWidth', 2.5, 'DisplayName', 'P_T → P_D (hydrolysis)');
plot([s_range(1) s_range(end)], [RDT_scalar RDT_scalar], 'r--', 'LineWidth', 2.5, 'DisplayName', 'P_D → P_T (rev. hydrolysis)');
set(gca, 'XLim', [s_range(1) s_range(end)], 'YScale', 'log', 'YLim', [1 1000]);
xlabel('Strain (μm)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Rate (s⁻¹)', 'FontSize', 10, 'FontWeight', 'bold');
title('Hydrolysis (ATP cycle)', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);

%% Panel 2: Attachment (P_D ↔ P1)
nexttile; hold on; grid on; box on;
RD1_mean = mean(RD1);  % nominal attachment rate
R1D_mean = mean(R1D);  % detachment from P1
line_r1d = plot(s_range, interp1(s, R1D, s_range, 'linear'), 'r--', 'LineWidth', 2.5, 'DisplayName', 'P1 → P_D (detach, strain-dep)');
plot([s_range(1) s_range(end)], [RD1_mean RD1_mean], 'b-', 'LineWidth', 2.5, 'DisplayName', 'P_D → P1 (attach)');
% Mark piecewise control points for R1D
if isfield(params, 'PieceWiseStrainDepR1DX') && ~isempty(params.PieceWiseStrainDepR1DX)
    r1d_vals = params.kd * params.PieceWiseStrainDepR1DParams;
    plot(params.PieceWiseStrainDepR1DX, r1d_vals, 'x', 'MarkerSize', 10, 'LineWidth', 2.5, 'Color', line_r1d.Color, 'HandleVisibility', 'off');
end
set(gca, 'XLim', [s_range(1) s_range(end)], 'YScale', 'log', 'YLim', [1 1000]);
xlabel('Strain (μm)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Rate (s⁻¹)', 'FontSize', 10, 'FontWeight', 'bold');
title('Attachment / P1 Detachment', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);

%% Panel 3: Power Stroke (P1 ↔ P2)
nexttile; hold on; grid on; box on;
line_r12 = plot(s_range, interp1(s, R12, s_range, 'linear'), 'b-', 'LineWidth', 2.5, 'DisplayName', 'P1 → P2 (power stroke)');
line_r21 = plot(s_range, interp1(s, R21, s_range, 'linear'), 'r--', 'LineWidth', 2.5, 'DisplayName', 'P2 → P1 (reverse stroke)');
% Mark piecewise control points for R12
if isfield(params, 'PieceWiseStrainDepX') && ~isempty(params.PieceWiseStrainDepX)
    r12_vals = params.k1 * params.PieceWiseStrainDepParams;
    plot(params.PieceWiseStrainDepX, r12_vals, 'x', 'MarkerSize', 10, 'LineWidth', 2.5, 'Color', line_r12.Color, 'HandleVisibility', 'off');
end
% Mark piecewise control points for R21
if isfield(params, 'PieceWiseStrainDepR21X') && ~isempty(params.PieceWiseStrainDepR21X)
    r21_vals = params.k_1 * params.PieceWiseStrainDepR21Params;
    plot(params.PieceWiseStrainDepR21X, r21_vals, 'x', 'MarkerSize', 10, 'LineWidth', 2.5, 'Color', line_r21.Color, 'HandleVisibility', 'off');
end
set(gca, 'XLim', [s_range(1) s_range(end)], 'YScale', 'log', 'YLim', [1 1000]);
xlabel('Strain (μm)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Rate (s⁻¹)', 'FontSize', 10, 'FontWeight', 'bold');
title('Power Stroke (P1 ↔ P2)', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);

%% Panel 4: P2 Detachment
nexttile; hold on; grid on; box on;
line_r2 = plot(s_range, interp1(s, R2, s_range, 'linear'), 'b-', 'LineWidth', 2.5, 'DisplayName', 'P2 → P_D (detach)');
% Backward from P_D to P2 is typically negligible; show baseline
plot([s_range(1) s_range(end)], [0.1 0.1], 'r--', 'LineWidth', 2.5, 'DisplayName', 'P_D → P2 (negligible)');
% Mark piecewise control points for R2 (with strain offset correction)
if isfield(params, 'PieceWiseStrainDep2X') && ~isempty(params.PieceWiseStrainDep2X)
    strain_offset = params.dr2 - params.dr;
    r2_x_corrected = params.PieceWiseStrainDep2X - strain_offset;
    r2_vals = params.k2 * params.PieceWiseStrainDep2Params;
    plot(r2_x_corrected, r2_vals, 'x', 'MarkerSize', 10, 'LineWidth', 2.5, 'Color', line_r2.Color, 'HandleVisibility', 'off');
end
set(gca, 'XLim', [s_range(1) s_range(end)], 'YScale', 'log', 'YLim', [1 1000]);
xlabel('Strain (μm)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Rate (s⁻¹)', 'FontSize', 10, 'FontWeight', 'bold');
title('P2 Detachment', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);

%% Panel 5: SRX-ATP (P_T ↔ P_SR)
nexttile; hold on; grid on; box on;
RPT2SR_vec = params.ksr0 * exp(-F_range / params.sigma2);
RSR2PT_vec = params.kmsr * exp(F_range / params.sigma1);
plot(F_range, RPT2SR_vec, 'b-', 'LineWidth', 2.5, 'DisplayName', 'P_T → P_{SR} (entry)');
plot(F_range, RSR2PT_vec, 'r--', 'LineWidth', 2.5, 'DisplayName', 'P_{SR} → P_T (exit)');
set(gca, 'XLim', [F_range(1) F_range(end)], 'YScale', 'log', 'YLim', [1 1000]);
xlabel('Force (kPa)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Rate (s⁻¹)', 'FontSize', 10, 'FontWeight', 'bold');
title('SRX-ATP (P_T ↔ P_{SR})', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);

%% Panel 6: SRX-ADP (P_D ↔ P_SRD)
nexttile; hold on; grid on; box on;
RPD2SRD_vec = params.ksrd * exp(-F_range / params.sigma_srd2);
RSRD2PD_vec = params.kmsrd * exp(F_range / params.sigma_srd1);
plot(F_range, RPD2SRD_vec, 'b-', 'LineWidth', 2.5, 'DisplayName', 'P_D → P_{SRD} (entry)');
plot(F_range, RSRD2PD_vec, 'r--', 'LineWidth', 2.5, 'DisplayName', 'P_{SRD} → P_D (exit)');
set(gca, 'XLim', [F_range(1) F_range(end)], 'YScale', 'log', 'YLim', [1 1000]);
xlabel('Force (kPa)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Rate (s⁻¹)', 'FontSize', 10, 'FontWeight', 'bold');
title('SRX-ADP (P_D ↔ P_{SRD})', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);

%% Panel 7: Inter-SRX (P_SR ↔ P_SRD)
nexttile; hold on; grid on; box on;
RSR2SRD_scalar = params.ksr2srd;  % constant
RSRD2SR_scalar = params.ksrd2sr;  % constant
plot([F_range(1) F_range(end)], [RSR2SRD_scalar RSR2SRD_scalar], 'b-', 'LineWidth', 2.5, 'DisplayName', 'P_{SR} → P_{SRD}');
plot([F_range(1) F_range(end)], [RSRD2SR_scalar RSRD2SR_scalar], 'r--', 'LineWidth', 2.5, 'DisplayName', 'P_{SRD} → P_{SR}');
set(gca, 'XLim', [F_range(1) F_range(end)], 'YScale', 'log', 'YLim', [1 1000]);
xlabel('Force (kPa)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Rate (s⁻¹)', 'FontSize', 10, 'FontWeight', 'bold');
title('Inter-SRX (P_{SR} ↔ P_{SRD})', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);

%% Panel 8: Summary legend / parameter table
nexttile; axis off;
param_str = sprintf('kah=%.1f, kamh=%.2f, ka=%.1f s^-1', params.kah, params.kamh, params.ka);
sigma_str = sprintf('sigma1=%.1f, sigma2=%.1f kPa', params.sigma1, params.sigma2);
srd_str = sprintf('sigma_srd1=%.0f, sigma_srd2=%.1f kPa', params.sigma_srd1, params.sigma_srd2);

summaryText = ['Cross-Bridge State Transition Rates' char(10) char(10) ...
    'Strain-Dependent (Panels 1-4):' char(10) ...
    '  * Forward (blue) vs Reverse (red dashed)' char(10) ...
    '  * X-axis: Strain s (um)' char(10) ...
    '  * Log scale shows exponential sensitivity' char(10) char(10) ...
    'Force-Dependent (Panels 5-7):' char(10) ...
    '  * X-axis: Force F_SR (kPa)' char(10) ...
    '  * Panel 5: PT/SR (entry suppressed, exit accel)' char(10) ...
    '  * Panel 6: PD/SRD (ADP-bound pool)' char(10) ...
    '  * Panel 7: Inter-SRX (ATP to ADP transition)' char(10) ...
    'Parameters: ' char(10) ...
    '  ' param_str char(10) ...
    '  ' sigma_str char(10) ...
    '  ' srd_str char(10) ...
    'v5: Force-balanced SRX pool'];
text(0.05, 0.95, summaryText, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'FontSize', 9, 'Parent', gca);

% Overall title
% sgtitle('Cross-Bridge State Transitions: Complete Kinetic Scheme', ...
%     'FontSize', 13, 'FontWeight', 'bold');

return;