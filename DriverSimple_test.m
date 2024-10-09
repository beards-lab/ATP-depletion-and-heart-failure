clear;
figure(2);
clf; 

% initialize parameters
params0 = getParams();
ModelParamsOptim_fudgeSlackVelocities

% new modifiers - to be overwritten
params0.justPlotStateTransitionsFlag = false;

LoadData; 
tic
RunBakersExp;
toc