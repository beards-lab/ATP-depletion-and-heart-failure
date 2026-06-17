% RunBoundedFitOccupancyAB
% =====================================================================
% A/B of the global binding-site occupancy term against the BoundedFit basis.
% Mirrors Workbench/DriverBoundedFit.m setup exactly, then re-runs RunBakersExp
% under several occupancy configurations and compares features_model + cost.
%
% Configs:
%   off            : baseline (UseGlobalOccupancySaturation = false)
%   linear 0.15    : global, OccupancyForm='linear',   P_bound_max=0.15
%   linear 0.12    : global, OccupancyForm='linear',   P_bound_max=0.12
%   langmuir 0.15  : global, OccupancyForm='langmuir', P_bound_max=0.15
%
% Output: prints a comparison table and saves results to
%   Analysis/BindingSiteOccupancy/occupancy_AB_results.mat
%
% Run:  cd repo root; addpath(genpath('.')); RunBoundedFitOccupancyAB
% =====================================================================

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));

%% ---- build the BoundedFit base parametrization (copied from DriverBoundedFit) ----
params0 = getParams();
run(fullfile(root, 'params', 'ModelOptParams_TL3_iter_17.m'));
params0.mods = {};
params0.g    = [];
params0 = getParams(params0, [], true, false);

% experiment configuration
params0.RunForceVelocity     = true;
params0.RunSlack             = true;
params0.RunKtr               = false;
params0.RunStairs            = false;
params0.RunForceVelocityTime = false;
params0.RunForceLengthEstim  = false;
params0.EvalFeatures         = true;
params0.BreakOnODEUnstable   = false;
params0.MaxRunTime           = 120;
params0.RunSlackSegments     = 'All';
try
    if isempty(gcp('nocreate')); parpool('Threads', 5); end
    params0.RunSlackSegments = 'AllPar';
catch ME
    warning(ME.identifier, 'Parallel pool unavailable (%s); serial.', ME.message);
end

params0.fn = boundedOutputFn(0.001);

% inline SRX / rate params (from DriverBoundedFit)
params0.ksr2srd = "=kah";
params0.ksrd2sr = "=kah/10";
params0.kah = 100;
params0.kamh = "=kah/10";
params0.ksr0  = 500;
params0.kmsr  = 0;
params0.kmsrd = 40;
params0.ksrd  = 0;
params0.ksr2srd = "=kah";
params0.ksr0    = 190;
params0.kmsrd   = 60;
params0.kstiff2 = 6.9794e+04*0.5;
params0.UseSuperRelaxed = true;
params0.UseSuperRelaxedADP = true;
params0.UseVelGaussAttachment = false;
params0.v_att_amplitude = 0.5;
params0.v_att_center = 0;
params0.v_att_sigma = 0.15;
params0.FV_velocities = - [0, 0.5, 1, 2, 4];

% physiological rate reseed overlay (kah=80, kd=100, k2=350, ka=375 ...)
run(fullfile(root, 'params', 'params_reseeded.m'));

params0.justPlotStateTransitionsFlag = false;
params0.PlotFeatureFitting = false;

params_base = params0;   % frozen baseline configuration

%% ---- A/B configurations ----
cfg(1) = struct('name','off',           'use',false, 'form','linear',   'cap',0.30);
cfg(2) = struct('name','linear 0.30',   'use',true,  'form','linear',   'cap',0.30);
cfg(3) = struct('name','linear 0.20',   'use',true,  'form','linear',   'cap',0.20);
cfg(4) = struct('name','langmuir 0.30', 'use',true,  'form','langmuir', 'cap',0.30);

R = struct([]);
for k = 1:numel(cfg)
    params0 = params_base;
    params0.UseGlobalOccupancySaturation = cfg(k).use;
    params0.OccupancyForm = cfg(k).form;
    params0.P_bound_max   = cfg(k).cap;

    % fresh feature containers each run
    features_data  = struct();
    features_model = struct();

    fprintf('\n===== config: %s =====\n', cfg(k).name); t0 = tic;
    figure(90000+k); clf;
    RunBakersExp;        % populates E, features_model, features_data
    cost_vec = plotFeatures(features_data, features_model, features_model, params0.fn, params0);
    fc = sum(cost_vec, 'omitnan');
    bounds = parameterBounds();
    E_phys = evalPhysiologyCost(params0, bounds);
    fprintf('%s: feature cost=%.4f  E_phys=%.4f  GRAND=%.4f  (%.1fs)\n', ...
        cfg(k).name, fc, E_phys, fc+E_phys, toc(t0));

    g1 = @(s) firstNum(features_model, s);
    R(k).name        = cfg(k).name;
    R(k).feat_cost   = fc;
    R(k).E_phys      = E_phys;
    R(k).grand       = fc + E_phys;
    R(k).E           = E;
    R(k).cost_vec    = cost_vec;
    R(k).fn          = params0.fn;
    R(k).XTOR        = g1('XTOR');
    R(k).XTOR_vmax   = g1('XTOR_vmax');
    R(k).SRX_ss      = g1('SRX_ss');
    R(k).attached_ss = g1('attached_ss');
    R(k).PT_ss       = g1('PT_ss');
    R(k).FV_fnorm    = vecOf(features_model,'FV_fnorm');
    R(k).ktr         = meanOf(features_model,'ktr');
    R(k).peak1_y     = meanOf(features_model,'peak1_y');
    R(k).vall_y      = meanOf(features_model,'vall_y');
    R(k).steady      = meanOf(features_model,'steady');
    R(k).features_model = features_model;
end

%% ---- comparison table ----
fprintf('\n\n================ OCCUPANCY A/B SUMMARY ================\n');
fprintf('%-14s %8s %8s %8s | %7s %9s %7s %8s %7s\n', ...
    'config','featcost','E_phys','GRAND','XTOR','attach_ss','SRX_ss','PT_ss','ktr');
for k = 1:numel(R)
    fprintf('%-14s %8.3f %8.3f %8.3f | %7.3g %9.3f %7.3f %8.3f %7.3g\n', ...
        R(k).name, R(k).feat_cost, R(k).E_phys, R(k).grand, ...
        R(k).XTOR, R(k).attached_ss, R(k).SRX_ss, R(k).PT_ss, R(k).ktr);
end
fprintf('\nFV_fnorm (normalised force at v=[0 .5 1 2 4]):\n');
for k = 1:numel(R)
    fprintf('  %-14s %s\n', R(k).name, num2str(R(k).FV_fnorm, '%6.3f '));
end
% data FV_fnorm for reference (same across configs)
if isfield(R(1).features_model,'FV_fnorm')
    fprintf('  (data target)  see features_data.FV_fnorm in saved .mat\n');
end

save(fullfile(fileparts(mfilename('fullpath')), 'occupancy_AB_results.mat'), 'R', 'cfg');
fprintf('\nSaved occupancy_AB_results.mat\n');

%% ---- local helpers ----
function v = firstNum(fm, name)
    if isfield(fm,name) && ~isempty(fm.(name))
        x = fm.(name); x = x(:); v = x(find(~isnan(x),1));
        if isempty(v), v = NaN; end
    else, v = NaN; end
end
function v = meanOf(fm, name)
    if isfield(fm,name) && ~isempty(fm.(name)), v = mean(fm.(name),'omitnan'); else, v = NaN; end
end
function v = vecOf(fm, name)
    if isfield(fm,name) && ~isempty(fm.(name)), v = fm.(name)(:)'; else, v = NaN; end
end
