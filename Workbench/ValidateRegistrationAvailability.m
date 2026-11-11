% ValidateRegistrationAvailability.m
% =========================================================================
% Phase C2 validation: does the new UseRegistrationAvailability state fix the
% FV "shoulder" (relative isometric depression) WITHOUT breaking ktr / slack?
%
% Mechanism (Model/dPUdT_CombinedTransitions.m + getParams.m): one extra scalar
% ODE state A_reg in [A0,1] gating attachment (ka_eff *= A_reg), driven by the
% IMPOSED velocity (params.Vums):
%     A_inf  = A0 + (1-A0)*|vel|/(|vel| + v_ref_reg)
%     dA_reg = (A_inf - A_reg)/tau_reg
% Isometric (vel=0) -> A_reg=A0<1 (depressed isometric reference); sliding ->
% A_reg->1 (force maintained -> flat shoulder). tau_reg~1/ktr keeps ktr single-order.
%
% Runs the optimizer-result basis (params/params_reseeded_opt.m) on FV+Ktr+Slack
% with the feature OFF (baseline) then ON, and prints the FV_fnorm shoulder, ktr,
% A and XTOR side-by-side vs data. Does NOT optimize — this is a behavioural probe.
%
% Run:  cd(root); addpath(genpath('.')); ValidateRegistrationAvailability
% =========================================================================

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root));
LoadData;

%% -------- base = the optimizer result (full param snapshot) -------------
params0 = getParams();
run(fullfile(root, 'params', 'params_reseeded_opt.m'));   % overlays the optim result
params0 = getParams(params0, [], true, false);

params0.RunForceVelocity = true;  params0.RunKtr = true;  params0.RunSlack = true;
params0.RunStairs = false;  params0.RunForceVelocityTime = false;
params0.EvalFeatures = true;  params0.BreakOnODEUnstable = false;  params0.MaxRunTime = 120;
params0.PlotEachSeparately = 0;  params0.PlotProbsOnFig = 0;
params0.fn = boundedOutputFn(0.001);
params0.FV_velocities = -[0 0.5 1 2 4];

if isempty(gcp('nocreate')); parpool('Threads', 5); end
params0.RunSlackSegments = 'AllPar';

%% ===== A: feature OFF (baseline = optimizer result) =====================
params0.UseRegistrationAvailability = false;
params0 = getParams(params0, [], true, false);
fprintf('\n################  A) UseRegistrationAvailability = OFF  ################\n');
tic; RunBakersExp; toc
fd  = features_data;
fmA = features_model;
reportFeatureCost(fd, fmA, params0.fn);

%% ===== B: feature ON, no force recompensation ==========================
params0.UseRegistrationAvailability = true;
params0.A0 = 0.6;  params0.v_ref_reg = 0.8;  params0.tau_reg = 0.02;   % first-probe values
params0 = getParams(params0, [], true, false);
fprintf('\n################  B) ON, no recompensation  (A0=%.2f)  ################\n', params0.A0);
tic; RunBakersExp; toc
fmB = features_model;

%% ===== C: feature ON + kstiff force-scale recompensation ================
% A_reg=A0 scales attachment down -> isometric force drops. kstiff scales force
% per head, restoring absolute A/steady WITHOUT changing FV_fnorm (a ratio) ->
% the flattened shoulder survives. ON(0.6) dropped A 67.3->52.7 (x0.783) -> x1.28.
kstiff_recomp = 1.28;
params0.kstiff1 = params0.kstiff1 * kstiff_recomp;
params0.kstiff2 = params0.kstiff2 * kstiff_recomp;
params0 = getParams(params0, [], true, false);
fprintf('\n################  C) ON + kstiff x%.2f recompensation  (A0=%.2f)  ################\n', ...
        kstiff_recomp, params0.A0);
tic; RunBakersExp; toc
fmC = features_model;
reportFeatureCost(fd, fmC, params0.fn);

%% -------- focused FV-shoulder + ktr comparison --------------------------
tA = sum(evalFeatureCost(fd, fmA, params0.fn, 2), 'omitnan');
tB = sum(evalFeatureCost(fd, fmB, params0.fn, 2), 'omitnan');
tC = sum(evalFeatureCost(fd, fmC, params0.fn, 2), 'omitnan');
fprintf('\n==================  FV SHOULDER + ktr: data vs OFF vs ON vs ON+recomp  ==================\n');
fprintf('  FV_v          :  %s\n', sprintf('%7.3f ', [fd.FV_v]));
fprintf('  FV_fnorm data :  %s\n', sprintf('%7.3f ', [fd.FV_fnorm]));
fprintf('  FV_fnorm OFF  :  %s\n', sprintf('%7.3f ', [fmA.FV_fnorm]));
fprintf('  FV_fnorm ON   :  %s\n', sprintf('%7.3f ', [fmB.FV_fnorm]));
fprintf('  FV_fnorm C    :  %s\n', sprintf('%7.3f ', [fmC.FV_fnorm]));
fprintf('  --------------------------------------------------------------\n');
fprintf('  %-12s %8s %8s %8s %8s\n', 'metric', 'data', 'OFF', 'ON', 'ON+rcmp');
fprintf('  %-12s %8.3g %8.3g %8.3g %8.3g\n', 'ktr',      scl(fd,'ktr'),  scl(fmA,'ktr'),  scl(fmB,'ktr'),  scl(fmC,'ktr'));
fprintf('  %-12s %8.3g %8.3g %8.3g %8.3g\n', 'ktr_rmse', scl(fd,'ktr_rmse'), scl(fmA,'ktr_rmse'), scl(fmB,'ktr_rmse'), scl(fmC,'ktr_rmse'));
fprintf('  %-12s %8.3g %8.3g %8.3g %8.3g\n', 'A',        scl(fd,'A'),    scl(fmA,'A'),    scl(fmB,'A'),    scl(fmC,'A'));
fprintf('  %-12s %8.3g %8.3g %8.3g %8.3g\n', 'steady',   scl(fd,'steady'), scl(fmA,'steady'), scl(fmB,'steady'), scl(fmC,'steady'));
fprintf('  %-12s %8s %8.3g %8.3g %8.3g\n',   'XTOR',     '-',            scl(fmA,'XTOR'), scl(fmB,'XTOR'), scl(fmC,'XTOR'));
fprintf('  %-12s %8s %8.3g %8.3g %8.3g\n',   'TOTAL cost','-',           tA, tB, tC);
fprintf('========================================================================================\n');

%% -------- local helpers (operate only on args — no workspace capture) ---
function v = scl(s, f)
    if isfield(s, f); v = double([s.(f)]); v = v(1); else; v = NaN; end
end
