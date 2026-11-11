% RunBoundedFit_Optim.m
% =========================================================================
% Bounded multi-start optimizer for the 2-state fit (Phase B/C of the FV<->ktr
% plan). Replicates the DriverBoundedFit basis (iter_17 + inline SRX block +
% params_reseeded), then minimizes the bounded scalar cost
%   sum(E) = evaluateBakersExp(g, params0)
% (data-fit features + literature OUTPUT bounds + the 'AssertParams' input-bound
% guard, all via boundedOutputFn) over g-multipliers on a deliberately chosen
% free-set. Box-bounds each multiplier; the physiological param bounds are also
% enforced softly through the AssertParams term in the cost.
%
% PURPOSE: (a) refine the p1-force (Route B) config — especially the FV R2 shape
% while p1-force holds ktr; (b) the "lock proof" — sweep / confirm no g closes
% FV + ktr + the rest together in 2-state.
%
% Run:  cd(root); addpath(genpath('.')); RunBoundedFit_Optim
% NOTE: ~20-40 s per eval (parallel). fminsearch ~ a few hundred evals -> long;
% start with the short multi-start below, or reduce MaxFunEvals. NOT auto-run to
% convergence by the assistant — verify on a few iterations first.
% =========================================================================

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root));
LoadData;                                                   % data for RunBakersExp

%% -------- replicate the DriverBoundedFit basis --------------------------
params0 = getParams();
run(fullfile(root, 'params', 'ModelOptParams_TL3_iter_17.m'));
params0.mods = {}; params0.g = [];
params0 = getParams(params0, [], true, false);

params0.RunForceVelocity = true;  params0.RunSlack = true;
params0.RunKtr = false;  params0.RunStairs = false;  params0.RunForceVelocityTime = false;
params0.EvalFeatures = true;  params0.BreakOnODEUnstable = false;  params0.MaxRunTime = 120;
params0.RunSlackSegments = 'All';

if isempty(gcp('nocreate')); parpool('Threads', 5); end
params0.RunSlackSegments = 'AllPar';

params0.fn = boundedOutputFn(0.001);

% inline SRX/parking block (mirrors DriverBoundedFit lines 91-126)
params0.ksr2srd = "=kah";  params0.ksrd2sr = "=kah/10";  params0.kah = 100;  params0.kamh = "=kah/10";
params0.ksr0 = 190;  params0.kmsr = 0;  params0.kmsrd = 60;  params0.ksrd = 0;
params0.kstiff2 = 6.9794e+04*0.5;
params0.UseSuperRelaxed = true;  params0.UseSuperRelaxedADP = true;
params0.UseVelGaussAttachment = false;
params0.FV_velocities = -[0 0.5 1 2 4];

run(fullfile(root, 'params', 'params_reseeded.m'));         % current best (p1-force config)

%% -------- free-set (sensitivity-selected; FV/ktr/slack-transient drivers) ----
% Chosen by campaign sensitivity (use wit): ka & R2-shortening knots shape FV;
% k2/P_bound_max set duty (XTOR); dr1 is the Route-B p1-force knob; kstiff1/2 scale
% force. Each g(i) multiplies the CURRENT params_reseeded value; bounded below.
params0.mods = {'ka', 'k2', 'dr1', 'kstiff1', 'kstiff2', 'P_bound_max', ...
                'PieceWiseStrainDep2Params__5', 'PieceWiseStrainDep2Params__6'};
g0 = ones(1, numel(params0.mods));
lb = [0.5 0.7 0.3 0.4 0.4 0.4 0.2 0.2];    % multiplier lower bounds
ub = [4.0 1.5 4.0 2.5 2.5 2.0 5.0 2.0];    % multiplier upper bounds

%%
params0.g = g0;
params0.PlotEachSeparately = 1;
features_ghost = features_model;
params0.k2 = params0.k2/2;
tic
RunBakersExp;
toc
plotFeatures(features_data, features_model, features_ghost, params0.fn)
%%
optimfun = @(g) evaluateBakersExp(g, params0);
fprintf('Baseline cost (g=1): %.4f\n', optimfun(g0));

%% -------- optimize: bounded multi-start (fmincon) -----------------------
opts = optimoptions('fmincon', 'Display','iter', 'Algorithm','sqp', ...
    'MaxFunctionEvaluations', 300, 'StepTolerance', 1e-3, 'FiniteDifferenceStepSize', 5e-2);
nStart = 3;                                   % perturbed restarts
best = struct('g', g0, 'f', optimfun(g0));
for s = 1:nStart
    if s == 1; gs = g0; else; gs = min(ub, max(lb, g0 .* (0.7 + 0.6*rand(size(g0))))); end
    [gx, fx] = fmincon(optimfun, gs, [],[],[],[], lb, ub, [], opts);
    fprintf('start %d: f = %.4f\n', s, fx);
    if fx < best.f; best.g = gx; best.f = fx; end
end

fprintf('\nBEST cost %.4f at g = [%s]\n', best.f, num2str(best.g, '%.3f '));
g0 = gx;
%% -------- write the optimized params snapshot ---------------------------
params0.g = best.g;
params0 = getParams(params0, best.g, false, true);          % bake g into the values
params0.mods = {}; params0.g = [];
writeParamsToMFile(fullfile(root, 'params', 'params_reseeded_opt.m'), params0);
fprintf('Saved params/params_reseeded_opt.m — load it in place of params_reseeded to verify.\n');
