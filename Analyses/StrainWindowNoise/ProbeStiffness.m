%% ProbeStiffness.m -- where does ode15s spend its work when the tolerance is tightened?
%
% One serial slack run per case, the RHS wrapped by dPUdT_Counted so every
% call's time is recorded. Cases are rows of CASES: {label, overrides}. A run
% that hits BUDGET seconds is stopped by the MaxRunTime event and reported as
% such, with the calls it made up to then.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit8_opt'; end
if ~exist('BUDGET', 'var'), BUDGET = 240; end
if ~exist('CASES', 'var')
    CASES = { 'RelTol 1e-2 (current)',      struct()
              'RelTol 1e-3',                struct('odeRelTol', 1e-3, 'odeAbsTol', 1e-6)
              'RelTol 1e-3, ekSE = 1',      struct('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'ekSE', 1) };
end
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..')));
res = fullfile(here, 'results');

P = getParams(loadParams(PARAMSET), [], true, false);
P.RunForceVelocity = false; P.RunSlack = true; P.RunStairs = false; P.RunKtr = false;
P.RunForceVelocityTime = false; P.RunSlackPassive = true; P.EvalFeatures = true;
P.BreakOnODEUnstable = false; P.PlotEachSeparately = 0; P.PlotFeatureFitting = 0;
P.RunSlackSegments = 'All'; P.MaxRunTime = BUDGET; P.OptimizeOn = 'Feats';
P.ghostSave = ''; P.ghostLoad = ''; P.SaveBest = 0; P.mods = {}; P.g = ones(1,0);
P.modelFcn = '@dPUdT_Counted';

global RHS_T RHS_N %#ok<GVMIS>
edges = [-10 0 70 74.40 74.45 74.6 75.0 75.1 75.35 75.4 75.6 75.7 75.95 76.0 76.25 76.3 76.55 76.6 76.85 76.9 78 80];
ST = struct();
for c = 1:size(CASES, 1)
    Pc = P; ov = CASES{c, 2}; f = fieldnames(ov);
    for k = 1:numel(f), Pc.(f{k}) = ov.(f{k}); end
    Pc = getParams(Pc, Pc.g, true, false);
    RHS_T = zeros(1e6, 1); RHS_N = 0;
    tic; status = 'finished';
    try
        [~, ~, fm] = runSlackExperiment(Pc);
    catch err
        status = ['stopped: ' err.message]; fm = struct();
    end
    w = toc;
    tt = RHS_T(1:RHS_N);
    [cnt, ~] = histcounts(tt, edges);
    ST(c).label = CASES{c,1}; ST(c).wall = w; ST(c).calls = RHS_N; ST(c).cnt = cnt; ST(c).status = status; ST(c).fm = fm;
    [~, worst] = max(cnt);
    fprintf('%-28s %6.0f s  %8d RHS calls  (%s)  busiest window %.2f-%.2f s: %d calls; last t %.4f\n', ...
        CASES{c,1}, w, RHS_N, status, edges(worst), edges(worst+1), cnt(worst), max(tt));
end
save(fullfile(res, ['stiffness_' PARAMSET '.mat']), 'ST', 'edges', 'CASES');
