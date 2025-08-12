%% Explores which parameters are necessary to run the simulation

params0 = getParams();

params0.DryRun = true;
LoadData;

tic
RunBakersExp;
toc

%% Find params and vals using wildcards

matchStructFields(params0, 'Run*', false, true);

%%
params0 = getParams();
params0.DryRun = true;

[passes, fails, messi] = testRequiredFields(params0);

