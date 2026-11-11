%% TestSRXParameters — validate SRX+SRXD initialization design
% Test script to verify that the new SRX parametrization keeps total SRX fraction ~0.2
% Runs steady-state simulation (zero velocity, no slack protocol) for fast evaluation

clear;
disp('Starting SRX parameter test...');

% Load experimental data (required for RunBakersExp)
disp('Loading experimental data...');
LoadData;

% Initialize parameters with defaults
disp('Initializing parameters...');
params0 = getParams();

% Load base parametrization from most recent iteration
disp('Loading base parametrization (TL3_iter_17)...');
run('params/ModelOptParams_TL3_iter_17.m');
params0 = getParams(params0);

% Overlay the SRX+SRXD design
disp('Overlaying SRX design (TL3_iter_17_SRXstart)...');
run('params/ModelOptParams_TL3_iter_17_SRXstart.m');
params0 = getParams(params0);

% Configure for steady-state simulation (fast test)
params0.RunSlack = 0;
params0.RunForceVelocity = 1;
params0.FV_velocities = [0];  % only zero-velocity condition
params0.RunKtr = 0;
params0.RunStairs = 0;
params0.RescaleOutputDt = 0;  % keep original time base

% Optional: reduce verbosity and set output options
params0.PlotEachSeparately = 0;
params0.recalculateDataFeats = false;

fprintf('\n=== Testing SRX Parameter Initialization ===\n');
fprintf('Base: ModelOptParams_TL3_iter_17\n');
fprintf('Overlay: ModelOptParams_TL3_iter_17_SRXstart\n');
fprintf('Conditions: Steady-state (zero velocity, no slack)\n');
fprintf('Expected: SRXT + SRXD ≈ 0.2 total\n');
fprintf('UseSuperRelaxed: %d, UseSuperRelaxedADP: %d\n', params0.UseSuperRelaxed, params0.UseSuperRelaxedADP);
fprintf('SRXT_0: %.4f, SRXD_0: %.4f, Sum: %.4f\n\n', params0.SRXT_0, params0.SRXD_0, params0.SRXT_0 + params0.SRXD_0);

% Run the experiment evaluation
disp('Running steady-state simulation...');
tic
E = RunBakersExp;
elapsed = toc;

fprintf('\n=== Results ===\n');
fprintf('Total error: %.4f\n', E);
fprintf('Elapsed time: %.2f seconds\n', elapsed);
fprintf('Check plots for SRXT and SRXD time evolution.\n');
fprintf('If SRXT + SRXD drifts significantly above 0.2, kinetic parameters need adjustment.\n');
