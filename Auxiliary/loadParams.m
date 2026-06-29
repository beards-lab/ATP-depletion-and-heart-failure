function params0 = loadParams(name, base)
%LOADPARAMS Safely load a params/*.m snapshot without clobbering the caller's struct.
%
%   Every snapshot under params/ (written by writeParamsToMFile) assigns to a
%   hardcoded variable named `params0`, e.g.:
%       params0.k2 = 150;
%   Loading such a file with run(file) from a driver script only works if the
%   driver's own struct happens to be named `params0`. If the driver names its
%   struct anything else, run() silently creates a stray `params0` instead, the
%   driver's real struct is left untouched at its defaults, and the overrides
%   are lost with no error or warning.
%
%   LOADPARAMS sidesteps this by run()-ing the snapshot inside its own local
%   workspace (where the local variable is, deliberately, also called
%   params0) and returning the resulting struct under whatever name the
%   caller chooses. The renamed-variable failure mode above cannot occur,
%   because the caller never needs a variable literally named params0:
%
%   Usage:
%       foo = loadParams('params_reseeded_regavail_opt2');
%       p   = getParams(foo, [], true, false);
%
%   or, combined in one line (this is the intended call pattern -- LOADPARAMS
%   itself does NOT call getParams):
%       p = getParams(loadParams('params_reseeded_regavail_opt2'), [], true, false);
%
%   To overlay a snapshot on top of an existing struct (rather than starting
%   from an empty one), pass it as BASE:
%       p = getParams(loadParams('params_reseeded_regavail_opt2', p), [], true, false);
%
%   Inputs:
%     name - bare snapshot name, resolved to <repoRoot>/params/<name>.m, OR a
%            full/relative path to a .m file. '.m' is appended if missing.
%     base - (optional) struct to use as the starting point, so the file's
%            assignments overlay onto it. Default: empty struct (struct()).
%
%   Output:
%     params0 - the struct built by the snapshot file's params0 assignments,
%               returned under whatever variable name the caller chooses.
%
%   Note: LOADPARAMS does NOT call getParams -- callers must do that
%   themselves afterward to resolve expressions, rebuild the strain grid, etc.
%
%   See also: getParams, writeParamsToMFile

    if nargin < 2 || isempty(base)
        base = struct();
    end

    if ~isstruct(base) || numel(base) ~= 1
        error('loadParams:invalidBase', 'BASE must be a scalar struct.');
    end

    if ~contains(name, '.')
        name = [name '.m'];
    end

    % Resolve the snapshot file: try as given (full/relative path) first,
    % then fall back to <repoRoot>/params/<name>.
    if isfile(name)
        theFile = name;
    else
        root = fileparts(fileparts(mfilename('fullpath'))); % Auxiliary is directly under root
        candidate = fullfile(root, 'params', name);
        if isfile(candidate)
            theFile = candidate;
        else
            error('loadParams:fileNotFound', ...
                'Could not find param snapshot "%s" (tried "%s" and "%s").', ...
                name, name, candidate);
        end
    end

    % Run the snapshot in THIS function's local scope. The file's
    % `params0.field = ...` assignments mutate the local variable below;
    % they cannot touch the caller's workspace because run() executes the
    % script as if it were a subfunction of loadParams.
    params0 = base; %#ok<NASGU> (assigned again by run(), kept for clarity)
    run(theFile);

    if ~isstruct(params0) || numel(params0) ~= 1 || isempty(fieldnames(params0))
        error('loadParams:invalidResult', ...
            '"%s" did not produce a non-empty scalar struct named params0.', theFile);
    end
end
