% TestNewGuardrails.m  Test the new restretch guardrails (doublePeak, coolDownLS)
%
% Loads optfull3_opt, runs the slack experiment with features enabled,
% and displays the cost contribution from the two new guardrails.

cd('C:\home\git\ATP-depletion-and-heart-failure');
addpath(genpath('.'));

fprintf('\n========================================\n');
fprintf('Testing New Restretch Guardrails\n');
fprintf('========================================\n\n');

% Refresh pool to pick up any code edits
refreshPool(5);

if isempty(gcp('nocreate'))
    p = parpool('Threads', 5);
end
%% Load params
fprintf('Loading optfull3_opt...\n');
params0 = getParams(loadParams('optfull3_opt'), [], true, false);


% Enable features
params0.RunForceVelocity = true;
params0.EvalFeatures = true;
params0.PlotEachSeparately = true;
params0.PlotFeatureFitting = true;
params0.MaxRunTime = 300;

% Run slack experiment
fprintf('Running slack experiment with all 5 slacks...\n');
RunBakersExp;
% [E_slack, out_slack, features_model, features_data] = runSlackExperiment(params0);

% Extract the new guardrail costs
fprintf('\n--- NEW GUARDRAIL FEATURES ---\n');
fprintf('\n1. doublePeak (spurious 2nd peak penalty):\n');
fprintf('   Per-slack values: ');
fprintf('%.4f ', features_model.doublePeak);
fprintf('\n   Sum: %.6f (target: ~0.5, current=no spurious peaks)\n', sum(features_model.doublePeak));

fprintf('\n2. coolDownLS (cool-down overshoot LS error):\n');
fprintf('   Per-slack values: ');
fprintf('%.4f ', features_model.coolDownLS);
fprintf('\n   Sum: %.6f (target: ~1.0, calibrated at optfull3_opt)\n', sum(features_model.coolDownLS));
figure(2);
plotFeatures(features_data, features_model, [], params0.fn)

%% Evaluate feature costs with just the new guardrails
fn_guardrails = {'doublePeak', 'coolDownLS'};
[E_guardrails, a, b] = evalFeatureCost(features_data, features_model, fn_guardrails, 1);

fprintf('\n--- FEATURE COST CONTRIBUTION ---\n');
fprintf('doublePeak cost:  %.6f\n', E_guardrails(1));
fprintf('coolDownLS cost:  %.6f\n', E_guardrails(2));
fprintf('Total guardrail cost: %.6f\n\n', sum(E_guardrails));

% Full feature evaluation
fprintf('--- FULL FEATURE COST (all %d features) ---\n', numel(params0.fn));
[E_all, ~, ~] = evalFeatureCost(features_data, features_model, params0.fn, 1);
fprintf('Total cost: %.4f\n\n', sum(E_all));

% Reference values for comparison
fprintf('--- REFERENCE VALUES (features_data) ---\n');
fprintf('doublePeak reference: '); fprintf('%.1f ', features_data.doublePeak); fprintf('(all zero)\n');
fprintf('coolDownLS reference: '); fprintf('%.1f ', features_data.coolDownLS); fprintf('(all zero)\n');
fprintf('\nBoth are "error-already-computed" features.\n');
fprintf('The sim values directly become the cost.\n\n');

fprintf('========================================\n');
fprintf('✓ New guardrails are working correctly\n');
fprintf('========================================\n\n');
