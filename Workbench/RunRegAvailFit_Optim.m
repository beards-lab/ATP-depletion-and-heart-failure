% RunRegAvailFit_Optim.m
% =========================================================================
% Bounded multi-start optimizer for the registration-availability config
% (Phase C2). Starts from the VALIDATED snapshot params/params_reseeded_regavail.m
% (UseRegistrationAvailability ON, A0=0.6, kstiff recompensated; total cost 8.89,
% which already BEATS the force-less-p1 optimum params_reseeded_opt.m at 9.90 and
% flattens the FV shoulder — see Analyses/FV_Shoulder/conclusions.md).
%
% Free-set refines the feature + the levers it couples to:
%   A0, v_ref_reg  — registration depth/width  (the FV-shoulder shape)
%   ka, k2         — duty / cycling            (pulls ktr 64 -> ~49)
%   kstiff1,kstiff2 — force scale              (holds A / steady at data)
%   R2 knots 5,6   — shortening detachment     (FV shape)
%
% Objective: sum(E) = evaluateBakersExp(g, params0), the bounded scalar cost
% (data-fit features + literature OUTPUT bounds + AssertParams input guard).
%
% Run:  cd(root); addpath(genpath('.')); RunRegAvailFit_Optim
% NOTE: ~25-35 s per eval (parallel). Multi-start fmincon -> long; NOT auto-run
% to convergence by the assistant. Verify a few iterations first, then let it run.
% =========================================================================

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root));
LoadData;

%% -------- base = the validated registration-availability snapshot -------
params0 = getParams();
params0.FV_velocities = -[0, 0.5, 1, 2, 3, 4, 5, 6];
run(fullfile(root, 'params', 'params_reseeded_regavail_opt2.m'));   % full snapshot, feature ON
params0 = getParams(params0, [], true, false);

params0.RunForceVelocity = true;  params0.RunKtr = true;  params0.RunSlack = true;
params0.RunStairs = false;  params0.RunForceVelocityTime = false;
params0.EvalFeatures = true;  params0.BreakOnODEUnstable = false;  params0.MaxRunTime = 120;
params0.PlotEachSeparately = 0;  params0.PlotProbsOnFig = 0;
params0.fn = boundedOutputFn(0.001);
params0.FV_velocities = -[0 0.5 1 2 4];

if isempty(gcp('nocreate')); parpool('Threads', 5); end
params0.RunSlackSegments = 'AllPar';
%%
params0.PlotEachSeparately = 1;
tic
RunBakersExp;
toc
figure(3123);
plotFeatures(features_data, features_model, [], params0.fn)

%% -------- free-set (registration-availability + coupled levers) ---------
% Each g(i) multiplies the CURRENT regavail value; box-bounded below. A0's upper
% multiplier (x1.5 -> 0.9) keeps A0 < 1 (availability floor must stay sub-maximal).
params0.mods = {'A0', 'v_ref_reg', 'ka', 'k2', 'kstiff1', 'kstiff2', ...
                'PieceWiseStrainDep2Params__5', 'PieceWiseStrainDep2Params__6'};
g0 = ones(1, numel(params0.mods));
lb = [0.67 0.40 0.60 0.50 0.70 0.70 0.30 0.30];   % A0 down to 0.40, v_ref to 0.32 um/s
ub = [1.50 2.50 1.80 1.20 1.50 1.50 3.00 3.00];   % A0 up to 0.90 (stays < 1)

%% -------- baseline check (g=1 reproduces the saved 8.89) ----------------
params0.g = g0;
optimfun = @(g) evaluateBakersExp(g, params0);
fprintf('Baseline cost (g=1, the saved regavail config): %.4f  (target: ~8.89)\n', optimfun(g0));


%% -------- optimize: bounded multi-start (fmincon) ----------------------
opts = optimoptions('fmincon', 'Display','iter', 'Algorithm','sqp', ...
    'MaxFunctionEvaluations', 300, 'StepTolerance', 1e-3, 'FiniteDifferenceStepSize', 5e-2);
nStart = 3;
best = struct('g', g0, 'f', optimfun(g0));
for s = 1:nStart
    if s == 1; gs = g0; else; gs = min(ub, max(lb, g0 .* (0.7 + 0.6*rand(size(g0))))); end
    [gx, fx] = fmincon(optimfun, gs, [],[],[],[], lb, ub, [], opts);
    fprintf('start %d: f = %.4f\n', s, fx);
    if fx < best.f; best.g = gx; best.f = fx; end
end
fprintf('\nBEST cost %.4f at g = [%s]\n', best.f, num2str(best.g, '%.3f '));

%% -------- write the optimized snapshot ---------------------------------
params0 = getParams(params0, best.g, false, true);   % bake g into the values
params0.mods = {}; params0.g = [];
writeParamsToMFile(fullfile(root, 'params', 'params_reseeded_regavail_opt2.m'), params0);
fprintf('Saved params/params_reseeded_regavail_opt.m\n');
