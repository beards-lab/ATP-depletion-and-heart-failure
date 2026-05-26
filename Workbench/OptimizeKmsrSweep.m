%% OptimizeKmsrSweep — find optimal kmsr value for SRX steady-state ≈ 0.2
% Sweeps kmsr values to find best balance

clear; clc;
fprintf('=== Optimizing kmsr for target SRX = 0.2 ===\n\n');

% Sweep values (based on v1=2.0, v2=2.7, v3=3.5, extrapolate further)
kmsr_values = [2.0, 2.7, 3.5, 4.0, 4.5, 5.0];
results = table();

for i = 1:length(kmsr_values)
    kmsr = kmsr_values(i);
    fprintf('Testing kmsr = %.2f ...', kmsr);

    % Initialize base
    params0 = getParams();
    run('params/ModelOptParams_TL3_iter_17.m');
    params0 = getParams(params0);

    % Apply overlay with custom kmsr
    run('params/ModelOptParams_TL3_iter_17_SRXstart.m');  % load base overlay first
    params0.kmsr = kmsr;  % override with test value
    params0 = getParams(params0);

    % Run simulation
    params0.MgATP = 8;
    params0.Velocity = 0;
    params0.Ca = 1000;
    params0.ValuesInTime = 1;
    params0.PlotEachSeparately = 0;

    [force, out] = evaluateModel(@dPUdT_CombinedTransitions, [0 3], params0);

    % Extract SRX fractions
    ss = length(params0.s);
    idx_p_sr = 3*ss + 1;
    idx_p_srd = 3*ss + 5;

    PU_final = out.PU(:, end);
    srxt_final = PU_final(idx_p_sr);
    srxd_final = PU_final(idx_p_srd);
    srx_total = srxt_final + srxd_final;
    error = abs(srx_total - 0.2);

    results = [results; table(kmsr, srxt_final, srxd_final, srx_total, error, 'VariableNames', ...
        {'kmsr', 'SRXT', 'SRXD', 'Total', 'Error'})];

    fprintf(' → Total SRX = %.4f (error = %.6f)\n', srx_total, error);
end

fprintf('\n=== Results Table ===\n');
disp(results);

% Find best
[best_error, best_idx] = min(results.Error);
best_kmsr = results.kmsr(best_idx);
best_total = results.Total(best_idx);

fprintf('\n=== Recommendation ===\n');
fprintf('Best kmsr: %.2f\n', best_kmsr);
fprintf('Achieved SRX: %.4f (target: 0.2, error: %.6f)\n', best_total, best_error);

if best_error <= 0.01
    fprintf('✓ EXCELLENT: Within ±0.01 of target\n');
    fprintf('Recommendation: Use kmsr = %.2f\n', best_kmsr);
elseif best_error <= 0.02
    fprintf('✓ GOOD: Within ±0.02 of target\n');
    fprintf('Recommendation: Use kmsr = %.2f, then fine-tune via optimization\n', best_kmsr);
else
    fprintf('⚠ Fair: Further adjustment needed\n');
    fprintf('Suggest sweeping range [%.2f, %.2f] with finer granularity\n', ...
        max(kmsr_values(1), best_kmsr-1), best_kmsr+1);
end
