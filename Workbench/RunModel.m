%% demonstrates basic model usage
clear;
LoadData;
ModelParamsSlackKtrOpt;
params0.recalculateDataFeats = false;
params0.ghostSave = '';

% Use 'local' (process-based) pool — MEX is not supported in thread workers.
% Switch back to parpool('threads', 5) if UseCompiledMex is not set.
if isempty(gcp('nocreate'))
    p = parpool('local', 5);
end

params0.UseFastPPval   = true;
params0.UseCompiledMex = true;
params0.modelFcn = '@dPUdT_CombinedTransitions_mex';
params0.RunSlackSegments = 'AllPar';
tic
RunBakersExp;
toc
return;
%%

% load optimized params
% g = load('gopt.csv');
% Set defaults
params = getParams();

% adjust options
params.ValuesInTime = true;

% adjust params
params.Velocity = 0;
params.MgATP = 8;

% set dxdt func
fcn = @dPUdT_CombinedTransitions ;

% evaluate
[force out] = evaluateModel(fcn, [0 1], params);

animate1StateProbabilities(out, params);