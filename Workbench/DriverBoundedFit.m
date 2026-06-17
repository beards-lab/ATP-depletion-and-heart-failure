% DriverBoundedFit.m
% =====================================================================
% Single-shot driver that runs the most recent FULL parametrization and
% scores it against BOTH:
%
%   (1) INPUT boundaries  — physiological bounds on the model PARAMETERS,
%        from Model/parameterBounds.m (re-anchored 2026-06 to mouse/rat
%        alpha-MHC). Evaluated by evalPhysiologyCost and visualised by
%        plotPhysiologyDashboard (figure 80086).
%
%   (2) OUTPUT boundaries — literature ranges on model OUTPUTS, carried in
%        params0.fn as boundary entries 'field|lb-ub|weight'. Evaluated by
%        evalFeatureCost / plotFeatures (figure 80085). The output-bound
%        ranges below were updated 2026-06 from the output-feature audit
%        (see notes next to each entry).
%
% Base parametrization: params/ModelOptParams_TL3_iter_17.m (May 2026, the
% most recent complete params0 written by writeParamsToMFile).
%
% Run from MATLAB with the repo on the path, e.g.:
%   cd('C:\home\git\ATP-depletion-and-heart-failure'); addpath(genpath('.'));
%   DriverBoundedFit
% or non-interactively:
%   matlab -batch "cd('C:\home\git\ATP-depletion-and-heart-failure'); addpath(genpath('.')); DriverBoundedFit"
%
% See also: parameterBounds, evalPhysiologyCost, plotPhysiologyDashboard,
%           evalFeatureCost, plotFeatures, RunBakersExp, extractSlackAttributes
% =====================================================================

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root));

%% -------- options --------------------------------------------------
USE_SRX_OVERLAY = false;   % true => layer the 4-pool SRX overlay (v5) on top
                           % of iter_17. Note: SRX is also enabled inline below
                           % (lines ~138-139), so this flag is for the overlay
                           % rates only (ksr0, sigma etc from v5 file).
USE_RESEEDED    = true;    % true => apply physiological rate reseed overlay
                           % (kah=80, kd=100, k2=350, ka=375) to bring sub-floor
                           % iter_17 rates up to alpha-MHC physiological starting
                           % points before any optimization.
USE_PARPOOL     = true;   % true => parallel slack segments (needs Parallel
                           % Computing Toolbox); false => serial ('All').

%% -------- load most recent full parametrization --------------------
params0 = getParams();                                   % seed defaults
run(fullfile(root, 'params', 'ModelOptParams_TL3_iter_17.m'));   % overlay full set

if USE_SRX_OVERLAY
    run(fullfile(root, 'params', 'ModelOptParams_TL3_iter_17_SRXstart_v5.m'));
end

% Rebuild derived fields (strain grid params.s, PU0, pchip caches) consistently
% with the loaded SL0/switches. mods/g are empty so nothing is re-scaled.
params0.mods = {};
params0.g    = [];
params0 = getParams(params0, [], true, false);

%% -------- experiment configuration ---------------------------------
params0.RunForceVelocity     = true;     % FV curve  -> XTOR/XTOR_vmax not from FV
params0.RunSlack             = true;     % slack/restretch -> XTOR, SRX_ss, attached_ss, PT_ss
params0.RunKtr               = false;
params0.RunStairs            = false;
params0.RunForceVelocityTime = false;
params0.RunForceLengthEstim  = false;
params0.EvalFeatures         = true;     % populate features_model / features_data
params0.BreakOnODEUnstable   = false;
params0.MaxRunTime           = 120;      % s, watchdog per condition
params0.RunSlackSegments     = 'All';    % serial; switch to 'AllPar' with a pool

if USE_PARPOOL
    try
        if isempty(gcp('nocreate')); parpool('Threads', 5); end
        params0.RunSlackSegments = 'AllPar';
    catch ME
        warning(ME.identifier, 'Parallel pool unavailable (%s); falling back to serial.', ME.message);
    end
end

%% -------- OUTPUT + input-guard feature set (params0.fn) -------------
% The fn cell now lives in Model/boundedOutputFn.m — a reusable function that
% returns the data-fit features + literature-grounded OUTPUT bounds (XTOR,
% XTOR_vmax, SRX_ss, attached_ss, PT_ss; mouse/rat alpha-MHC) + the
% 'AssertParams' input guard. All bound choices are cited in that file.
% Edit boundedOutputFn.m to retune the bounds; the argument is the AssertParams
% weight scale (0 to omit the input guard).
params0.fn = boundedOutputFn(0.001);

%% -------- run all configured protocols ------------------------------
params0.ksr2srd = "=kah";
params0.ksrd2sr = "=kah/10";
params0.kah = 100;
params0.kamh = "=kah/10";


% --- Parking (into SRX): sets the pool size via ratio to mobilization ---
params0.ksr0  = 500;    % PT  -> SR
params0.kmsr  = 0;   % SR  -> PT

params0.kmsrd = 40;   % SRD -> PD
params0.ksrd  = 0;    % PD  -> SRD

params0.ksr2srd = "=kah";   % SR→SRD mixing = kah (100)
params0.ksr0    = 190;      % PT→SR parking rate  — sets pool to ~0.20 at max force
params0.kmsrd   = 60;       % SRD→PD outflow      — shapes SR/SRD split

params0.kstiff2 = 6.9794e+04*0.5;

params0.RunForceVelocity = true;
params0.RunSlack = false;

params0.UseSuperRelaxed = true;
params0.UseSuperRelaxedADP = true;


params0.RunSlack             = true;     % slack/restretch -> XTOR, SRX_ss, attached_ss, PT_ss
params0.RunSlackSegments = 'AllPar';
params0.RunForceVelocity = true;

params0.UseVelGaussAttachment = false;
params0.v_att_amplitude = 0.5; % 0.3
params0.v_att_center = 0;
params0.v_att_sigma = 0.15; % 0.5;

params0.FV_velocities = - [0, 0.5, 1, 2, 4];
features_data  = struct();
features_model = struct();

if USE_RESEEDED
    run(fullfile(root, 'params', 'params_reseeded.m'));
end

fprintf('Running protocols (FV + slack)...\n'); tic

params0.justPlotStateTransitionsFlag = false;
params0.PlotFeatureFitting = false;
figure(1);clf;
RunBakersExp;          % populates E, features_model, features_data, out, outs
fprintf('done in %.1f s. Protocol cost terms E = [%s]\n', toc, num2str(E, '%.3g '));

%% -------- OUTPUT-boundary report + figure (fig 80085 + dashboard) ---
% plotFeatures(...,params0) also renders the INPUT-bound physiology dashboard
% in figure 80086, so this single call shows both boundary layers.
if ~exist('features_ghost', 'var')
    features_ghost = features_model;
end
cost_vec = plotFeatures(features_data, features_model, features_ghost, params0.fn, params0);
fprintf('\nTotal OUTPUT feature cost: %.4f\n', sum(cost_vec, 'omitnan'));

% Detailed per-feature cost breakdown + raw data-vs-model physiology dump.
% (cost_vec breakdown sorted by contribution; OUTPUT-bound sub-costs read
%  straight from params0.fn so they can't drift from boundedOutputFn.)
reportFeatureCost(features_data, features_model, params0.fn);

% Conservation check: SRX + attached + PT (+ PuR + NP) should be ~1
if all(isfield(features_model, {'SRX_ss','attached_ss','PT_ss'}))
    s = sum([first(features_model.SRX_ss), first(features_model.attached_ss), ...
             first(features_model.PT_ss)]);
    fprintf('  conservation SRX+attached+PT = %.3f (remainder = PuR+NP)\n', s);
end

%% -------- INPUT-boundary report (parameter physiology cost) ---------
bounds = parameterBounds();
[E_phys, violations] = evalPhysiologyCost(params0, bounds);
fprintf('\nTotal INPUT (parameter) physiology cost: %.4f\n', E_phys);
if isempty(violations)
    fprintf('  all weighted parameter bounds satisfied.\n');
else
    fprintf('  %d parameter bound(s) violated (name: value [lb,ub] w -> penalty):\n', ...
        numel(violations));
    for k = 1:numel(violations)
        v = violations(k);
        fprintf('    %-14s %.4g  [%.4g, %.4g]  w=%g  -> %.3f\n', ...
            v.name, v.value, v.lb, v.ub, v.weight, v.penalty);
    end
end

% Standalone dashboard (in case plotFeatures' embedded one is not wanted)
figure(80087); clf;
plotPhysiologyDashboard(params0, bounds, gca);

fprintf('\nGRAND TOTAL (output features + input physiology) = %.4f\n', ...
    sum(cost_vec, 'omitnan') + E_phys);

%% -------- local helpers --------------------------------------------
function printOutBound(name, fm, lb, ub)
    if ~isfield(fm, name) || isempty(fm.(name))
        fprintf('  %-12s : (missing)\n', name); return;
    end
    val  = first(fm.(name));
    flag = '';
    if val < lb || val > ub; flag = '  <-- OUT OF BOUND'; end
    fprintf('  %-12s = %7.3g   [%.3g, %.3g]%s\n', name, val, lb, ub, flag);
end

function y = first(x)
    x = x(:); y = x(find(~isnan(x), 1));
    if isempty(y); y = NaN; end
end
                         