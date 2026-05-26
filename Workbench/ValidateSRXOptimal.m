%% ValidateSRXOptimal — validate resting SRX ≈ 0.2 from kmsrd=14.4 parametrization
% Tests that SRXstart_v4 initial conditions give SRXT+SRXD = 0.200 at t=0,
% i.e. at resting/zero-force conditions (physiologically correct: IHM forms at rest).
%
% NOTE: SRX depletes at active force due to force-dependent kinetics — this is
% physiologically correct. sigma1/sigma2 govern force-dependent mobilization.
% Target is resting SRX, not active-force SRX.

clear; clc;
fprintf('=== SRX Parametrization Validation (resting SRX = 0.2) ===\n\n');

% Ensure we run from the project root (params/ is at project root, not in Workbench)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
cd(projectRoot);

%% Initialize and load
params0 = getParams();
run('params/ModelOptParams_TL3_iter_17.m');
params0 = getParams(params0);

% Load SRXstart overlay (kmsrd=14.4 analytically derived for resting SRX=0.2)
run('params/ModelOptParams_TL3_iter_17_SRXstart_v5.m');
% params0 = getParams(params0);
ModelOptParams_TL3_iter_17_SRXstart_v5

fprintf('--- Parameters ---\n');
fprintf('ksr0 (PT→SR):     %.4f  s⁻¹\n', params0.ksr0);
fprintf('kmsr (SR→PT):     %.4f  s⁻¹\n', params0.kmsr);
fprintf('ksr2srd (SR→SRD): %.4f  s⁻¹\n', params0.ksr2srd);
fprintf('ksrd2sr (SRD→SR): %.4f  s⁻¹\n', params0.ksrd2sr);
fprintf('ksrd (PD→SRD):    %.4f  s⁻¹\n', params0.ksrd);
fprintf('kmsrd (SRD→PD):   %.4f  s⁻¹\n\n', params0.kmsrd);

%% Test 1: Initial conditions (resting SRX at zero force)
fprintf('--- Test 1: Resting SRX from initial conditions ---\n');
ss = length(params0.s);
Ns = params0.NumberOfStates;
SRXT_init = params0.PU0(Ns*ss + 1);
SRXD_init = params0.PU0(Ns*ss + 6);
srx_init  = SRXT_init + SRXD_init;

fprintf('SRXT_0 (from PU0): %.6f\n', SRXT_init);
fprintf('SRXD_0 (from PU0): %.6f\n', SRXD_init);
fprintf('TOTAL:             %.6f\n', srx_init);
fprintf('TARGET:            0.200000\n');
fprintf('ERROR:             %.9f\n\n', abs(srx_init - 0.2));

if abs(srx_init - 0.2) < 1e-4
    fprintf('✓ EXCELLENT: Resting SRX = 0.200 (target achieved)\n\n');
else
    fprintf('⚠ Deviation from 0.2 target: %.4f\n\n', abs(srx_init - 0.2));
end

%% Test 2: Dynamic run — show SRX time course and force-dependent depletion
fprintf('--- Test 2: Dynamic SRX time course at isometric Ca=1000 ---\n');
params0.MgATP = 8;
params0.Velocity = 0;
params0.Ca = 1000;
params0.ValuesInTime = 1;
params0.PlotEachSeparately = 0;
params0.SL0 = 2.0;

fprintf('Running 1-second isometric simulation...\n');
[~, out] = evaluateModel(@dPUdT_CombinedTransitions, [0 10], params0);

% Use correct indices from the output params
p  = out.params;
ss_out = length(p.s);
Ns_out = p.NumberOfStates;
iSR  = Ns_out*ss_out + 1;
iSRD = Ns_out*ss_out + 6;

% PU is [n_time × n_states] when ValuesInTime=1
fprintf('\n  t(s)   Force(kPa)   SRX(=SRXT+SRXD)\n');
show_t = [0, 0.01, 0.05, 0.1, 0.5, 1.0];
for ts = show_t
    [~, ti] = min(abs(out.t - ts));
    srx_t = out.PU(ti, iSR) + out.PU(ti, iSRD);
    fprintf('  %.2f   %.1f          %.4f\n', out.t(ti), out.Force(ti), srx_t);
end

% Force at which SRX halves
SRX_t = out.PU(:, iSR) + out.PU(:, iSRD);
half_ti = find(SRX_t < srx_init/2, 1);
if ~isempty(half_ti)
    fprintf('\n  SRX halves (→0.100) at F = %.1f kPa, t = %.3f s\n', ...
        out.Force(half_ti), out.t(half_ti));
end

fprintf('\nNote: SRX depletion at active force is physiologically correct.\n');
fprintf('sigma1=%.1f kPa and sigma2=%.1f kPa govern force-dependent mobilization.\n', ...
    p.sigma1, p.sigma2);
fprintf('Resting parametrization (kmsrd=%.1f) is ready for optimization.\n', p.kmsrd);
