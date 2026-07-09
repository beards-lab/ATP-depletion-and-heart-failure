function params0 = loadParams(name)
%LOADPARAMS  Load a parameter snapshot script into an isolated struct.
%
%   PARAMS0 = LOADPARAMS(NAME) executes the parameter snapshot .m script NAME
%   (a script that populates fields of `params0`, e.g. `params0.ka = 50;`) and
%   returns the resulting struct. The script is run inside this function's own
%   workspace rather than via RUN/EVAL in the base workspace, which prevents two
%   recurring bugs in this project:
%     * a stale `params0` from a previous call leaking fields into the new one;
%     * default-valued fields silently overriding the snapshot.
%
%   NAME may be given with or without the '.m' extension, and with or without a
%   path. Resolution order: an explicit path if given, then the MATLAB path
%   (WHICH), then the repo's `params/` folder.
%
%   Pass the returned struct straight to getParams, e.g.
%     params0 = getParams(loadParams('params_2state_seed'), [], true, false);
%
%   See also: getParams, refreshPool, writeParamsToMFile

fpath = iLocateSnapshot(name);

% Start from a clean struct; the snapshot script assigns into params0, which
% is created in *this* function's workspace because run() executes the file in
% the caller's workspace.
params0 = struct();
run(fpath);

if ~exist('params0', 'var') || ~isstruct(params0) || isempty(fieldnames(params0))
    error('loadParams:noParams0', ...
        'Snapshot "%s" did not define a non-empty params0 struct.', fpath);
end
end


function fpath = iLocateSnapshot(name)
[givenPath, base, ext] = fileparts(name);
if isempty(ext)
    ext = '.m';
end
fname = [base ext];

candidates = {};
if ~isempty(givenPath)
    candidates{end+1} = fullfile(givenPath, fname);
end
w = which(fname);
if ~isempty(w)
    candidates{end+1} = w;
end
% Repo root is the parent of the Auxiliary/ folder that holds this file.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
candidates{end+1} = fullfile(repoRoot, 'params', fname);

for k = 1:numel(candidates)
    if exist(candidates{k}, 'file') == 2
        fpath = candidates{k};
        return;
    end
end

error('loadParams:notFound', ...
    'Could not locate snapshot "%s" (searched the MATLAB path and %s).', ...
    name, fullfile(repoRoot, 'params'));
end
