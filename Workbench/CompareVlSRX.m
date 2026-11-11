%% CompareVsSRX — test v1 vs v2 SRX parameterization
% Quick test to see if v2 with increased kmsr improves steady-state SRX

clear; clc;

fprintf('=== Comparing SRX Parametrization Versions ===\n\n');

versions = {'v1', 'v2'};
results = struct();

for v = 1:length(versions)
    version = versions{v};

    fprintf('--- Testing %s ---\n', version);

    % Initialize base
    params0 = getParams();
    run('params/ModelOptParams_TL3_iter_17.m');
    params0 = getParams(params0);

    % Load appropriate overlay
    if v == 1
        run('params/ModelOptParams_TL3_iter_17_SRXstart.m');
    else
        run('params/ModelOptParams_TL3_iter_17_SRXstart_v2.m');
    end
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

    results.(version).ksr0 = params0.ksr0;
    results.(version).kmsr = params0.kmsr;
    results.(version).srxt = srxt_final;
    results.(version).srxd = srxd_final;
    results.(version).total = srx_total;
    results.(version).error = abs(srx_total - 0.2);

    fprintf('  kmsr: %.4f\n', params0.kmsr);
    fprintf('  SRXT: %.6f\n', srxt_final);
    fprintf('  SRXD: %.6f\n', srxd_final);
    fprintf('  Total: %.6f (target: 0.2, error: %.6f)\n\n', srx_total, abs(srx_total - 0.2));
end

% Summary table
fprintf('\n=== Summary ===\n');
fprintf('%-6s | %-8s | %-8s | %-8s | %-8s | %-8s\n', 'Ver', 'kmsr', 'SRXT', 'SRXD', 'Total', 'Error');
fprintf('%-6s | %-8s | %-8s | %-8s | %-8s | %-8s\n', '---', '------', '------', '------', '------', '------');
for v = 1:length(versions)
    ver = versions{v};
    fprintf('%-6s | %.4f   | %.4f   | %.4f   | %.4f   | %.4f\n', ...
        ver, results.(ver).kmsr, results.(ver).srxt, results.(ver).srxd, results.(ver).total, results.(ver).error);
end

% Recommendation
[~, best_idx] = min([results.v1.error, results.v2.error]);
best_ver = versions{best_idx};
fprintf('\nBest: %s (error: %.6f)\n', best_ver, results.(best_ver).error);

if best_idx == 2
    fprintf('Recommendation: Use v2 with kmsr = %.4f\n', results.v2.kmsr);
else
    fprintf('Both versions similar. May need further optimization.\n');
end
