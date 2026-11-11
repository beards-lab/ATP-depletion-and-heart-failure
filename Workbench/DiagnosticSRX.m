%% DiagnosticSRX — inspect SRX initialization at steady state
% Simple diagnostic to check SRXT and SRXD steady-state values
% Focuses on zero-velocity condition only

clear; clc;
disp('=== SRX Steady-State Diagnostic ===');

% Load base parametrization
disp('Loading base parameters (TL3_iter_17)...');
params0 = getParams();
run('params/ModelOptParams_TL3_iter_17.m');
params0 = getParams(params0);

% Overlay SRX design
disp('Overlaying SRX design (TL3_iter_17_SRXstart)...');
run('params/ModelOptParams_TL3_iter_17_SRXstart.m');
params0 = getParams(params0);

fprintf('\n--- Input SRX Parameters ---\n');
fprintf('UseSuperRelaxed: %d\n', params0.UseSuperRelaxed);
fprintf('UseSuperRelaxedADP: %d\n', params0.UseSuperRelaxedADP);
fprintf('SRXT_0: %.6f (initial)\n', params0.SRXT_0);
fprintf('SRXD_0: %.6f (initial)\n', params0.SRXD_0);
fprintf('SRXT_0 + SRXD_0: %.6f\n', params0.SRXT_0 + params0.SRXD_0);

fprintf('\n--- Kinetic Parameters ---\n');
fprintf('ksr0 (PT→SR): %.6f\n', params0.ksr0);
fprintf('kmsr (SR→PT): %.6f\n', params0.kmsr);
fprintf('ksr2srd (SR→SRD): %.6f\n', params0.ksr2srd);
fprintf('ksrd2sr (SRD→SR): %.6f\n', params0.ksrd2sr);
fprintf('ksrd (PD→SRD): %.6f\n', params0.ksrd);
fprintf('kmsrd (SRD→PD): %.6f\n', params0.kmsrd);

% Run single zero-velocity condition with diagnostics
fprintf('\n--- Running Single Zero-Velocity Condition ---\n');
params0.MgATP = 8;  % 8 mM ATP
params0.Velocity = 0;  % isometric
params0.Ca = 1000;  % activated (high Ca)
params0.ValuesInTime = 1;
params0.PlotEachSeparately = 0;

% Simulate for fixed time to reach steady state
[force, out] = evaluateModel(@dPUdT_CombinedTransitions, [0 3], params0);  % 3 seconds

% Extract SRX fractions from output
if isfield(out, 'PU') && ~isempty(out.PU)
    PU_final = out.PU(:, end);  % last time point
    ss = length(params0.s);

    % SRX state indexing: [p1(1..ss), p2(1..ss), p3(1..ss), P_SR, NP, SL, LSE, PD, P_SRD, ...]
    idx_p_sr = 3*ss + 1;
    idx_p_srd = 3*ss + 5;

    srxt_final = PU_final(idx_p_sr);
    srxd_final = PU_final(idx_p_srd);

    fprintf('\n--- Steady-State Results (t=3s) ---\n');
    fprintf('SRXT: %.6f\n', srxt_final);
    fprintf('SRXD: %.6f\n', srxd_final);
    fprintf('SRXT + SRXD: %.6f\n', srxt_final + srxd_final);
    fprintf('Target: ~0.2\n');
    fprintf('Difference from target: %.6f\n\n', abs(srxt_final + srxd_final - 0.2));
else
    disp('ERROR: Could not extract PU from simulation output');
end

fprintf('--- Assessment ---\n');
if abs(srxt_final + srxd_final - 0.2) < 0.05
    fprintf('✓ PASS: SRX fractions are within acceptable range (±0.05 of 0.2)\n');
else
    fprintf('✗ FAIL: SRX fractions deviate significantly from target\n');
    fprintf('  Recommendation: Adjust kmsrd (currently %.6f) to balance SRD sink\n', params0.kmsrd);
end
