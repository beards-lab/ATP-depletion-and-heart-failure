function [fn_grouped, feat_vals] = assertActiveParams(params0, scale)
%ASSERTACTIVEPARAMS  Build fn assertions for parameters active in the current parametrisation.
%
%   [FN_GROUPED, FEAT_VALS] = ASSERTACTIVEPARAMS(PARAMS0) returns:
%
%     FN_GROUPED  - A single grouped fn string consumable by evalFeatureCost
%                   as a joint-cost tile (one bar in plotFeatures):
%                     'ka|50-500|10,kd|5-200|10,...'
%                   Suitable as a direct replacement for an 'AssertParams|XX'
%                   placeholder in params0.fn.
%     FEAT_VALS   - Struct keyed by param name with current values from params0.
%                   Merge into features_model before calling evalFeatureCost.
%
%   [FN_GROUPED, FEAT_VALS] = ASSERTACTIVEPARAMS(PARAMS0, SCALE) multiplies
%   every weight from parameterBounds by scalar SCALE (default 1).
%   Use SCALE to tune overall physiology-penalty pressure (e.g. 0.1 for a
%   soft hint, 10 for a hard constraint).
%
%   PARAM SELECTION MIRRORS dPUdT_CombinedTransitions switch logic:
%     - When a switch is OFF, the parameters gated by it are excluded.
%     - Core transition rates (ka, kd, k1, k2, ...) are always included.
%     - Only Tier A–C params (weight > 0 in parameterBounds) are included.
%     - Params absent from params0 or non-scalar are silently skipped.
%
%   IMPORTANT: this file is a manual mirror of the switch structure in
%   dPUdT_CombinedTransitions.  Any new switch / parameter added there
%   must be replicated here.
%
%   USAGE — place 'AssertParams|XX' in params0.fn, then call from RunBakersExp:
%
%       params0.fn = {'ktr', 'A', 'AssertParams|1'};
%
%   RunBakersExp expands these placeholders automatically before evalFeatureCost.
%   The expansion block in RunBakersExp (%%ASSERT ACTIVE PARAMS EXPANSION)
%   handles all 'AssertParams|XX' entries in params0.fn.
%
%   If FN_GROUPED is empty (all active params have weight=0 or are absent),
%   the placeholder entry is removed from fn rather than inserted as ''.
%
%   See also: evalFeatureCost, parameterBounds, resolveModelParam, RunBakersExp

if nargin < 2 || isempty(scale)
    scale = 1;
end

bounds    = parameterBounds();
active    = {};  % ordered list of param names to assert

% =========================================================
% ALWAYS-ACTIVE CORE PARAMETERS
% Accessed unconditionally by dPUdT_CombinedTransitions regardless of any switch.
% =========================================================
active = [active, { ...
    'ka', 'kd', 'k1', 'k_1', 'k2', ...  % XB cycle rates
    'kah', 'kamh', ...                    % ATP hydrolysis / reverse
    'kstiff1', 'kstiff2', ...             % XB stiffness
    'kSE', 'ekSE', ...                    % serial elastic element
    'mu', 'mu_neg', ...                   % internal viscous drag
    'dr', 'SL0', ...                      % structural: power stroke, setpoint
    'xrate', 'vmax' ...                   % rate scaler, max velocity
}];

% =========================================================
% THREE-STATE MODEL  (NumberOfStates == 3)
% =========================================================
if isfield(params0, 'NumberOfStates') && params0.NumberOfStates == 3
    active = [active, {'k3', 'k3m'}];
end

% kstiff3 is used as the p3 stiffness even in 2-state via fallback to kstiff2 in ODE.
% Only assert it when explicitly toggled ON (i.e. the user intends kstiff3 ≠ kstiff2).
if isfield(params0, 'UseKstiff3') && params0.UseKstiff3
    active{end+1} = 'kstiff3';
end

% =========================================================
% FILAMENT GEOMETRY  (UseOverlap)
% =========================================================
if isfield(params0, 'UseOverlap') && params0.UseOverlap
    active = [active, {'L_thick', 'L_hbare', 'L_thin'}];
end

% =========================================================
% PASSIVE FORCE  (UsePassive)
% gamma bound is only meaningful with the power-law formulation.
% =========================================================
if isfield(params0, 'UsePassive') && params0.UsePassive
    active{end+1} = 'k_pas';
    titinID = isfield(params0, 'UseTitinIdentifiedPassive') && params0.UseTitinIdentifiedPassive;
    if ~titinID
        active{end+1} = 'gamma';
    end
end

% =========================================================
% NEGATIVE-STRAIN STIFFNESS  (UseNegativeKstiff)
% =========================================================
if isfield(params0, 'UseNegativeKstiff') && params0.UseNegativeKstiff
    active = [active, {'kstiff1_n', 'kstiff2_n'}];
end

% =========================================================
% SUPER-RELAXED (IHM DRX/SRX) STATE  (UseSuperRelaxed)
% =========================================================
if isfield(params0, 'UseSuperRelaxed') && params0.UseSuperRelaxed
    active = [active, {'ksr0', 'kmsr', 'sigma1', 'sigma2'}];
end

% =========================================================
% SRX-ADP STATE  (UseSuperRelaxedADP)
% =========================================================
if isfield(params0, 'UseSuperRelaxedADP') && params0.UseSuperRelaxedADP
    active = [active, {'ksrd', 'kmsrd', 'sigma_srd1', 'sigma_srd2'}];
end

% =========================================================
% SRX ↔ SRD CROSS-TRANSITIONS  (UseSuperRelaxed && UseSuperRelaxedADP)
% =========================================================
if isfield(params0, 'UseSuperRelaxed')   && params0.UseSuperRelaxed && ...
   isfield(params0, 'UseSuperRelaxedADP') && params0.UseSuperRelaxedADP
    active = [active, {'ksrd2sr', 'ksr2srd'}];
end

% =========================================================
% A2 ACTIN HOPPING  (UseA2AttachmentShift)
% d_actin is a Tier-A structural constant but only enters the ODE here.
% =========================================================
if isfield(params0, 'UseA2AttachmentShift') && params0.UseA2AttachmentShift
    active = [active, {'d_actin', 'slope', 's_threshold_R'}];
end

% =========================================================
% MAXWELL DASHPOT (titin viscoelasticity)  (UseMaxwellDashpot)
% =========================================================
if isfield(params0, 'UseMaxwellDashpot') && params0.UseMaxwellDashpot
    active = [active, {'kSE_M', 'eta_M'}];
end

% =========================================================
% LATTICE SPACING CORRECTION  (UseLatticeSpacing)
% =========================================================
if isfield(params0, 'UseLatticeSpacing') && params0.UseLatticeSpacing
    active = [active, {'d10_ref', 'SL_ref_lattice', 'R_thick', 'R_thin', 'd_optimal'}];
end

% =========================================================
% CALCIUM / TROPONIN KINETICS  (UseCa)
% =========================================================
if isfield(params0, 'UseCa') && params0.UseCa
    active = [active, {'K_coop', 'k_on', 'k_off'}];
end

% =========================================================
% BUILD GROUPED STRING
% Resolve each active param, skip absent / weight-0 / NaN.
% =========================================================
assertions = {};
feat_vals  = struct();

skip_ok = isfield(params0, 'SkipParametersInBounds') && params0.SkipParametersInBounds;

for i = 1:numel(active)
    nm = active{i};
    if ~isfield(bounds, nm);  continue; end
    b = bounds.(nm);
    if b.weight == 0;         continue; end          % Tier D — no prior
    val = resolveModelParam(params0, nm);
    if isnan(val);            continue; end          % absent or non-scalar
    feat_vals.(nm) = val;                            % always export value
    in_bounds = (val >= b.lb) && (val <= b.ub);
    if in_bounds && skip_ok;  continue; end          % omit in-bounds from tile
    assertions{end+1} = sprintf('%s|%g-%g|%g', nm, b.lb, b.ub, b.weight * scale); %#ok<AGROW>
end

if isempty(assertions)
    fn_grouped = '';
else
    fn_grouped = strjoin(assertions, ',');
end

end
