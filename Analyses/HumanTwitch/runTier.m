function runTier(name, nIter)
%RUNTIER  Run one fitting tier by name (see tierConfigs). Usable from `matlab -batch`.
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');  addpath(genpath(root));
cfg = tierConfigs(name);
if nargin < 2, nIter = Inf; end
fitTwitchTier(cfg, nIter);
end
