function twPaths(flavor)
%TWPATHS  Put exactly ONE T-World copy on the path (they share function names).
%   flavor 'twatp'  -> TW_ATP/ of the main checkout (local, git-ignored, READ-ONLY: never edit)
%          'lowatp' -> C:\home\git\TWorld, branch lowAtp (separate repo, READ-ONLY)
%          'plain'  -> ThirdParty/TWorld of this repo
%   Also puts the repo (Model/, Auxiliary/, Analyses/) on the path. Errors if modelRunner,
%   getCurrentsStructure or getStartingState resolve anywhere else.
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
addpath(genpath(root));
rmpath(genpath(fullfile(root, 'ThirdParty', 'TWorld')));
switch flavor
    case 'twatp',  d = 'C:\home\git\ATP-depletion-and-heart-failure\TW_ATP';
    case 'lowatp', d = 'C:\home\git\TWorld';
    case 'plain',  d = fullfile(root, 'ThirdParty', 'TWorld');
    otherwise, error('twPaths:flavor', 'unknown flavor %s', flavor);
end
addpath(d, '-begin');
for f = ["modelRunner" "getCurrentsStructure" "getStartingState" "getDefaultParameters" "model_TWorld"]
    w = which(f);
    if ~startsWith(lower(w), lower(char(java.io.File(d).getCanonicalPath())))
        error('twPaths:shadow', '%s resolves to %s, expected %s', f, w, d);
    end
end
fprintf('twPaths: %s -> %s\n', flavor, d);
end
