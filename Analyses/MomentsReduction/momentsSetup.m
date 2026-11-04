function params = momentsSetup(paramName)
% MOMENTSSETUP  Clean path + base params for the moments analysis.
%
%   genpath('.') from the repo root also adds .claude/worktrees/*, whose
%   copies of Model/ then SHADOW the real ones (which() resolved
%   dPUdT_CombinedTransitions to a worktree copy when this analysis started).
%   So add only what is needed and strip anything under .claude/_archive.
if nargin < 1, paramName = 'fastgate_best'; end
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
root = char(java.io.File(root).getCanonicalPath());
pp = strsplit(path, pathsep);
bad = pp(contains(pp, {[filesep '.claude'], [filesep '_archive']}));
if ~isempty(bad), rmpath(bad{:}); end
addpath(genpath(fullfile(root, 'Model')), genpath(fullfile(root, 'Auxiliary')), ...
        fullfile(root, 'params'), fileparts(mfilename('fullpath')));
assert(startsWith(which('dPUdT_CombinedTransitions'), fullfile(root, 'Model')), 'shadowed dxdt');

params = getParams(loadParams(paramName), [], true, false);
params.odeRelTol = 1e-3;          % the post-StrainWindowNoise standard
params.odeAbsTol = 1e-6;
params.odeOutputDt = 5e-4;
params.MaxRunTime = 900;
params.ShowArrayShiftWarnings = false;
end
