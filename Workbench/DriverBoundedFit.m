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
                           % of iter_17 so the SRX_ss output bound is exercised.
                           % iter_17 itself has UseSuperRelaxed=0, so with the
                           % overlay OFF, SRX_ss ~ 0 and will read as a low-side
                           % violation of the [0.10,0.50] output bound BY DESIGN.
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

%% -------- OUTPUT boundaries (literature-grounded), in params0.fn ----
% Boundary entries 'field|lb-ub|weight' score the SIMULATED output against a
% literature range (cost 0 inside, quadratic outside). The grouped entry
% (commas, one joint cost / one bar tile) collects the state-occupancy and
% turnover outputs. Ranges updated 2026-06 (mouse/rat alpha-MHC, ~21-25 C,
% maximal Ca, isometric pre-slack steady state):
%
%   XTOR        per-head ATPase, isometric     1.5-10 s^-1  (Rossmanith 1995;
%                                                            Tanner/de Tombe 2007)
%   XTOR_vmax   per-head ATPase at Vmax          4-20 s^-1  (Fenn x2-5; Stienen 2004)
%   SRX_ss      SRX/OFF frac, ACTIVE state    0.10-0.50     (Reconditi 2017;
%               (0.75 was the RESTING value -> lowered)      Brunello 2020; Ma 2022)
%   attached_ss duty ratio (strong+weak)      0.05-0.30     (Pinzauti 2018;
%               (0.2-0.6 was ~2-4x too high -> lowered)      Reconditi 2017)
%   PT_ss       DRX available fraction         0.20-0.70    (residual/conservation)
% [1] selects ATP/segment index 1 (8 mM), matching prior usage.
%
% Plain entries (e.g. 'ktr|2', 'peak1_y|10') and the X-Y entry 'FV_fnorm|FV_v|10'
% are DATA-FIT features (sim vs experiment), not literature bounds — kept as in
% the base parametrization. 'ktr_rmse|0-0.2|.1' is a fit-quality bound.
% 'kstiff1|2000-30000|1' is an INPUT (parameter) bound expressed as a feature;
% the full input-bound set is enforced separately via parameterBounds() below.
params0.fn = {
    'FV_fnorm|FV_v|10', 'ktr|2', 'A|50', ...
    'ktr_rmse|0-0.2|.1', ...
    ['XTOR[1]|1.5-10|15,' ...
     'XTOR_vmax[1]|4-20,' ...
     'SRX_ss[1]|0.05-0.50|5,' ...
     'attached_ss[1]|0.3-0.60|10,' ...
     'PT_ss[1]|0.20-0.70'], ...
    't0_crossing|SLdiff|2', ...
    'restretchSlopeStart|0.1', ...
    'peak1_y|10', 'peak1_dSL|0.2', ...
    'vall_y|10', 'vall_t|0.2', 'peak2|5', ...
    'steady|50', 'vall2_dy|0.1', 'ovrsht_dy|0.1', ...
    'AssertParams|0.001'...
};

%% -------- run all configured protocols ------------------------------
params0.RunSlack             = true;     % slack/restretch -> XTOR, SRX_ss, attached_ss, PT_ss

params0.UseVelGaussAttachment = true;
params0.v_att_amplitude = 0.3;
params0.v_att_center = 0;
params0.v_att_sigma = 0.5;

params0.FV_velocities = - [0, 0.5, 1, 2, 4];
features_data  = struct();
features_model = struct();
fprintf('Running protocols (FV + slack)...\n'); tic

params0.PlotFeatureFitting = false;
figure(1);clf;
RunBakersExp;          % populates E, features_model, features_data, out, outs
fprintf('done in %.1f s. Protocol cost terms E = [%s]\n', toc, num2str(E, '%.3g '));

% params0.fn = {
%     'FV_fnorm|FV_v|10'};
cost_vec = plotFeatures(features_data, features_model, features_ghost, params0.fn, params0)
%% -------- OUTPUT-boundary report + figure (fig 80085 + dashboard) ---
% plotFeatures(...,params0) also renders the INPUT-bound physiology dashboard
% in figure 80086, so this single call shows both boundary layers.
if ~exist("features_ghost", "var")
    features_ghost = features_model;
end
cost_vec = plotFeatures(features_data, features_model, features_ghost, params0.fn);
fprintf('\nTotal OUTPUT feature cost: %.4f\n', sum(cost_vec, 'omitnan'));

printOutBound('XTOR',        features_model, 1.5, 10);
printOutBound('XTOR_vmax',   features_model, 4,   20);
printOutBound('SRX_ss',      features_model, 0.10, 0.50);
printOutBound('attached_ss', features_model, 0.25, 0.30);
printOutBound('PT_ss',       features_model, 0.20, 0.70);

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
