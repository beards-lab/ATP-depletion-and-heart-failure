% TestSinglePerturbation.m
RunMechanismEvaluation;
params0.mods = {'kd'};
params0.fn = {'ktr|SLslack', 'peak1_y'};

features_model_1 = struct();
params0.g = [1.0];
disp('Running with g=1.0');
ResidualAndJacobian(params0.g, params0, true);
% We have to extract features from the workspace features_model
fm1 = evalin('caller', 'features_model');

params0.g = [1.5];
disp('Running with g=1.5');
ResidualAndJacobian(params0.g, params0, true);
fm2 = evalin('caller', 'features_model');

disp('Features with g=1.0:');
disp(fm1);
disp('Features with g=1.5:');
disp(fm2);

if fm1.ktr == fm2.ktr
    disp('ERROR: The model output is IDENTICAL, parameters are NOT applying!');
else
    disp('SUCCESS: Model output changed, parameter multiplier logic is working.');
end
