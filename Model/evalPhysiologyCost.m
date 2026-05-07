function [E_phys, violations] = evalPhysiologyCost(params, bounds)
%EVALPHYSIOLOGYCOST  Sum of soft-bound penalties across all physiology bounds.
%
%   [E_PHYS, VIOLATIONS] = EVALPHYSIOLOGYCOST(PARAMS, BOUNDS) returns:
%       E_PHYS      scalar penalty cost (sum over all bounded params)
%       VIOLATIONS  struct array, one entry per bound that fired
%                   .name, .value, .lb, .ub, .weight, .penalty
%
%   BOUNDS defaults to PARAMETERBOUNDS() if omitted.
%
%   Penalty form is a quadratic exterior soft bound:
%       p(x; lb, ub) = max(0, (lb-x)/(ub-lb))^2 + max(0, (x-ub)/(ub-lb))^2
%   zero inside [lb, ub], rises as squared fractional violation outside.
%   This keeps fminsearch's simplex geometry well-defined (continuous
%   gradient at the bound boundary).
%
%   Bounds with weight = 0 are skipped (no compute, no penalty).
%   Params not present in PARAMS are silently skipped.
%
%   See also: parameterBounds, plotPhysiologyDashboard

if nargin < 2 || isempty(bounds)
    bounds = parameterBounds();
end

names = fieldnames(bounds);
E_phys = 0;
violations = struct('name', {}, 'value', {}, 'lb', {}, 'ub', {}, ...
                    'weight', {}, 'penalty', {});

for i = 1:numel(names)
    nm = names{i};
    b  = bounds.(nm);

    if b.weight == 0
        continue;
    end

    % Resolve param value: either direct field, or indexed via __N suffix
    val = resolveParam(params, nm);
    if isnan(val)  % param not present
        continue;
    end

    % Soft exterior penalty
    span = b.ub - b.lb;
    if span <= 0  % degenerate bound (Tier D fudge factors with lb=ub=0)
        if abs(val - b.lb) > 0
            p = abs(val - b.lb)^2;
        else
            p = 0;
        end
    else
        p = max(0, (b.lb - val) / span)^2 + max(0, (val - b.ub) / span)^2;
    end

    contribution = b.weight * p;
    E_phys = E_phys + contribution;

    if p > 0
        violations(end+1) = struct( ...
            'name',    nm, ...
            'value',   val, ...
            'lb',      b.lb, ...
            'ub',      b.ub, ...
            'weight',  b.weight, ...
            'penalty', contribution); %#ok<AGROW>
    end
end
end


function val = resolveParam(params, name)
% Return the numeric value for a parameter name. Supports both direct
% fields and the `field__index` notation used by the `mods` machinery.
val = NaN;
if isfield(params, name) && isnumeric(params.(name)) && isscalar(params.(name))
    val = params.(name);
    return;
end

% Try field__index split (e.g. PieceWiseStrainDepParams__3)
sep = strfind(name, '__');
if ~isempty(sep)
    base = name(1:sep(1)-1);
    idx  = str2double(name(sep(1)+2:end));
    if isfield(params, base) && ~isnan(idx) && isnumeric(params.(base)) ...
            && idx >= 1 && idx <= numel(params.(base))
        val = params.(base)(idx);
    end
end
end
