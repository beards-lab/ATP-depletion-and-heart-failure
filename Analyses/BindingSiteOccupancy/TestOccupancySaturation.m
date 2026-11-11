% TestOccupancySaturation  Verify the binding-site occupancy / attachment
% saturation mechanisms added to dPUdT_CombinedTransitions.
%
% Background (see OccupancySaturation_Report.md in this folder):
%   The model is a mean-field strain-distribution model. Its strain axis is the
%   cross-bridge strain (distance from the attachment point), NOT a position
%   along actin and NOT a binding-site index. Therefore a per-strain-bin
%   occupancy factor does not represent site competition. The only faithful
%   occupancy proxy is the GLOBAL bound fraction P_bound = p1_0+p2_0+p3_0.
%
%   New optional params (getParams.m, all default off => behaviour-preserving):
%     UseGlobalOccupancySaturation : RD1 *= f(P_bound)   [recommended]
%     P_bound_max                  : ceiling on bound fraction (~0.12-0.15 rat)
%     OccupancyForm                : 'linear' | 'langmuir'
%     UseTargetZoneSaturation      : [DEPRECATED] per-strain-bin modifier
%     rho_attach_max               : density cap for the deprecated per-bin form
%
% This script runs four checks and prints a table. Sections 1-4 are isometric
% single-condition runs (fast, robust). Section 5 is an optional full A/B
% against RunBakersExp using the user's current fitted basis.
%
% Run:  cd to repo root, addpath(genpath('.')), then run this file.

clearvars; close all force;
here = fileparts(mfilename('fullpath'));
cd(fullfile(here, '..', '..'));            % repo root
addpath(genpath('.'));
warning('off','all');

%% Section 1+2: A/B at the physical ceiling (P_bound_max = 0.15)
fprintf('\n=== A/B: global occupancy clamps the bound fraction ===\n');
p = baseCfg();                                         [F0,Pb0] = runIso(p);
p = baseCfg(); p.UseGlobalOccupancySaturation = true; p.P_bound_max = 0.15; p.OccupancyForm = 'linear';   [F1,Pb1] = runIso(p);
p = baseCfg(); p.UseGlobalOccupancySaturation = true; p.P_bound_max = 0.15; p.OccupancyForm = 'langmuir'; [F2,Pb2] = runIso(p);
fprintf('%-18s Force=%8.4g  P_bound=%.4f\n', 'baseline (off)', F0, Pb0);
fprintf('%-18s Force=%8.4g  P_bound=%.4f  (ceiling 0.15)\n', 'global linear',  F1, Pb1);
fprintf('%-18s Force=%8.4g  P_bound=%.4f  (half-sat 0.15)\n', 'global langmuir', F2, Pb2);

%% Section 3: P_bound_max sweep (global linear)
fprintf('\n=== Sweep P_bound_max (global linear) ===\n');
caps = [0.05 0.10 0.15 0.20 0.30];
for c = caps
    p = baseCfg(); p.UseGlobalOccupancySaturation = true; p.P_bound_max = c; p.OccupancyForm = 'linear';
    [F,Pb] = runIso(p);
    fprintf('P_bound_max=%.2f : Force=%8.4g  P_bound=%.4f\n', c, F, Pb);
end

%% Section 4: dS-invariance — GLOBAL form is grid-robust, per-bin form is not
fprintf('\n=== dS-invariance (fixed strain extent) ===\n');
for dd = [0.004 0.002]
    p = baseCfg(); p.dS = dd; p.UseGlobalOccupancySaturation = true; p.P_bound_max = 0.15; p.OccupancyForm = 'linear';
    [F,Pb,~,ss] = runIso(p);
    fprintf('GLOBAL  dS=%.4f ss=%3d : Force=%8.4g  P_bound=%.5f\n', dd, ss, F, Pb);
end
for dd = [0.004 0.002]
    p = baseCfg(); p.dS = dd; p.UseTargetZoneSaturation = true; p.rho_attach_max = 16;
    [F,Pb,~,ss] = runIso(p);
    fprintf('PER-BIN dS=%.4f ss=%3d : Force=%8.4g  P_bound=%.5f  [deprecated; grid-dependent]\n', dd, ss, F, Pb);
end

%% Section 5 (OPTIONAL): full A/B against the user's fitted basis via RunBakersExp
% Uncomment and point at your working param snapshot to measure the effect on the
% real multi-experiment cost (FV / Ktr / Stairs / Slack). Requires LoadData and a
% configured velocity table, exactly as in Analysis/ExperimentWAttachments.m.
%
% clear params0; params0 = getParams();
% ModelParamsFeats_FVSlackUpdateValley2          % or your current fitted basis
% LoadData;
% params0.RunForceVelocity = true; params0.RunKtr = true; params0.RunSlack = true; params0.RunStairs = true;
% % --- baseline ---
% params0.UseGlobalOccupancySaturation = false;
% RunBakersExp; E_base = E;  fprintf('baseline   E = %.4g\n', sum(E_base));
% % --- global linear, ceiling 0.15 ---
% params0.UseGlobalOccupancySaturation = true; params0.P_bound_max = 0.15; params0.OccupancyForm = 'linear';
% RunBakersExp; E_lin = E;   fprintf('global lin E = %.4g\n', sum(E_lin));
% % --- global langmuir ---
% params0.OccupancyForm = 'langmuir';
% RunBakersExp; E_lang = E;  fprintf('global lng E = %.4g\n', sum(E_lang));

fprintf('\nDone.\n');

%% ---- local functions (must be at end of script file) ----
function p = baseCfg()
    % Unfitted default isometric config used for mechanism verification.
    p = getParams();
    p.SL0 = 2.2; p.Slim_l = 1.6; p.Slim_r = 2.25;
    p.dS = 0.004; p.MaxStrainArraySize = 400;     % large => fixed strain extent
    p.MgATP = 8; p.UseOverlap = true; p.UsePassive = false; p.UseSuperRelaxed = false;
    p.ValuesInTime = true; p.PlotProbsOnFig = false;
    p.BreakOnODEUnstable = false; p.MaxRunTime = 60; p.g = [];
end

function [F, Pb, nt, ss] = runIso(p)
    p = getParams(p, p.g, true, true);
    [~, out] = evaluateModel(@dPUdT_CombinedTransitions, [0 1], p);
    ss = p.ss;
    if isstruct(out) && isfield(out,'p1_0') && ~isempty(out.p1_0)
        % p3_0 is only stored for the 3-state model; treat as 0 otherwise.
        Pb = lastOr0(out.p1_0) + lastOr0(out.p2_0) + lastOr0(out.p3_0);
        F = out.Force(end); nt = numel(out.t);
    else
        F = NaN; Pb = NaN; nt = 0;
    end
end

function v = lastOr0(x)
    if isempty(x), v = 0; else, v = x(end); end
end
