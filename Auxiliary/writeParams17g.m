function writeParams17g(fname, params, comment)
%WRITEPARAMS17G  Serialise a params struct to a .m snapshot at full precision.
%
%   writeParams17g(fname, params, comment)
%
% Same job as writeParamsToMFile, with one difference that matters: scalars are
% written with %.17g rather than %g. %g keeps SIX significant digits, so a
% snapshot round-trip perturbs every parameter by up to ~1e-6 relative and the
% file is therefore NOT the point that was fitted. On this objective that is
% not academic -- %.6g serialisation alone is worth 0.60 of cost
% (Docs objective_noise_floor; Analyses/AtpGateHillForm). %.17g round-trips a
% double exactly.
%
% Factored out of Analyses/MechanismAblation/BuildLeanSeed.m, which carried the
% first copy.
%
% See also writeParamsToMFile, loadParams.

fid = fopen(fname, 'w');
if fid < 0
    error('writeParams17g:open', 'cannot open %s for writing', fname);
end
c = onCleanup(@() fclose(fid));

fprintf(fid, '%% Generated file to manipulate the simulation parameters directly\r\n');
if nargin >= 3 && ~isempty(comment)
    fprintf(fid, '%% %s\r\n', comment);
end

nm = fieldnames(params);
for i = 1:numel(nm)
    v = params.(nm{i});
    if isa(v, 'function_handle') || isstruct(v); continue; end
    if isnumeric(v) && isscalar(v)
        fprintf(fid, 'params0.%s = %.17g;\r\n', nm{i}, v);
    elseif islogical(v) && isscalar(v)
        fprintf(fid, 'params0.%s = %d;\r\n', nm{i}, v);
    elseif ischar(v) && ~isempty(v)
        fprintf(fid, 'params0.%s = ''%s'';\r\n', nm{i}, v);
    elseif iscell(v) && ~isempty(v)
        fprintf(fid, 'params0.%s = {''%s''};\r\n', nm{i}, strjoin(v, ''', '''));
    elseif isnumeric(v) && ~isempty(v)
        fprintf(fid, 'params0.%s = [%s];\r\n', nm{i}, ...
                strjoin(arrayfun(@(x) sprintf('%.17g', x), v(:)', 'uni', 0), ' '));
    end
end
fprintf(fid, 'params0.mods = {};\r\nparams0.g = [];\r\n');
end
