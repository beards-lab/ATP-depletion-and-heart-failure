function v = resolveModelParam(params, name)
%RESOLVEMODELPARAMETER  Extract a scalar from params by name.
%
%   V = RESOLVEMODELPARAMETER(PARAMS, NAME) returns the scalar numeric value
%   for the named parameter.  Supports two naming conventions:
%
%     Direct field:   NAME = 'ka'          → params.ka
%     Indexed field:  NAME = 'arr__2'      → params.arr(2)
%                     (the evalPhysiologyCost / mods machinery convention)
%
%   Returns NaN when the field does not exist, is non-numeric, or is not a
%   scalar (e.g. a string expression not yet resolved, or a vector).
%
%   See also: evalPhysiologyCost, parameterBounds, RunBakersExp

if isfield(params, name)
    v = params.(name);
    if ~isnumeric(v) || numel(v) ~= 1
        v = nan;
    end
    return;
end

% Try field__index split (e.g. PieceWiseStrainDepParams__3)
sep = strfind(name, '__');
if ~isempty(sep)
    base = name(1:sep(1)-1);
    idx  = str2double(name(sep(1)+2:end));
    if isfield(params, base) && ~isnan(idx) && isnumeric(params.(base)) ...
            && idx >= 1 && idx <= numel(params.(base))
        v = params.(base)(idx);
        return;
    end
end

v = nan;
end
