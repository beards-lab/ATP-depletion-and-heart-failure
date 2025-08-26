%% demonstrates basic model usage

% load optimized params
% g = load('gopt.csv');
% Set defaults
params = getParams([], g);

% adjust options
params.ValuesInTime = true;

% adjust params
params.Velocity = 0;
params.MgATP = 8;

% set dxdt func
fcn = @dPUdTCa;

% evaluate
[force out] = evaluateModel(fcn, [0 1], params);

animateStateProbabilities(out, params);