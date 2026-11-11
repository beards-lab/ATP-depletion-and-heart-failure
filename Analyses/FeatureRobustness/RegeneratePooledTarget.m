%% RegeneratePooledTarget.m -- add rsT20e / rsTfitE / tOn20 to the pooled target, proven non-disruptive.
% Campaigns read Model/pooledTarget8mM.m while this runs, so nothing is written
% there until it is verified:
%   1. poolLoadRuns(true) re-extracts every run; each OLD feature must equal the
%      previous cache bit-for-bit (extractSlackAttributes' data path is unchanged)
%   2. WritePooledTarget writes to STAGE; every existing row of P.val, and Fref,
%      must equal the live file's bit-for-bit (the moderated P.s / P.z of the
%      existing features DO move: the shrinkage prior is fitted over all
%      features -- reported, not used at run time)
% Installing the staged file is a separate, deliberate step (see the log).
if ~exist('STAGE', 'var'), STAGE = 'C:\Users\fjezek\AppData\Local\Temp\featrob_staging'; end
here = fileparts(mfilename('fullpath')); root = fullfile(here, '..', '..');
addpath(genpath(root));
dp = fullfile(root, 'Analyses', 'DataPooling', 'results');
OLD = load(fullfile(dp, 'pool_raw_20260924_pre_robust.mat.bak'), '-mat', 'R'); OLD = OLD.R;

%% 1. re-extract
R = poolLoadRuns(true);
bad = {};
for i = 1:numel(R)
    f = fieldnames(OLD(i).feat);
    for k = 1:numel(f)
        if ~isequaln(OLD(i).feat.(f{k}), R(i).feat.(f{k})), bad{end+1} = sprintf('%s %s:%s', R(i).day, R(i).cond, f{k}); end %#ok<SAGROW>
    end
    if ~isequaln(OLD(i).A, R(i).A), bad{end+1} = sprintf('%s %s:anchors', R(i).day, R(i).cond); end %#ok<SAGROW>
end
fprintf('re-extraction: %d runs, old features differing: %d %s\n', numel(R), numel(bad), strjoin(bad, ', '));
assert(isempty(bad), 'the data path changed -- do not install');

%% 2. write to STAGE and compare with the live target
OUTF = fullfile(STAGE, 'pooledTarget8mM.m'); OPTS = struct(); %#ok<NASGU>
run(fullfile(root, 'Analyses', 'DataPooling', 'WritePooledTarget.m'));   % clears the workspace
STAGE = 'C:\Users\fjezek\AppData\Local\Temp\featrob_staging';
clear pooledTarget8mM; live = pooledTarget8mM();
addpath(STAGE, '-begin'); clear pooledTarget8mM; assert(startsWith(which('pooledTarget8mM'), STAGE));
new = pooledTarget8mM(); rmpath(STAGE); clear pooledTarget8mM;
rowBad = {};
for k = 1:numel(live.name)
    j = find(strcmp(new.name, live.name{k}), 1);
    if isempty(j) || ~isequal(new.val(j,:), live.val(k,:)) || ~strcmp(new.cls{j}, live.cls{k}), rowBad{end+1} = live.name{k}; end %#ok<SAGROW>
end
fprintf('staged target: %d features (live %d); existing rows differing: %d %s; Fref equal: %d\n', ...
        numel(new.name), numel(live.name), numel(rowBad), strjoin(rowBad, ','), isequal(new.Fref, live.Fref));
for q = {'rsT20e','rsTfitE','tOn20'}
    j = strcmp(new.name, q{1});
    fprintf('  %-8s val %s  sraw %.4f  s %.4f  z %.2f\n', q{1}, mat2str(round(new.val(j,:), 5)), new.sraw(j), new.s(j), new.z(j));
end
ds = nan(1, numel(live.name));
for k = 1:numel(live.name), ds(k) = new.s(strcmp(new.name, live.name{k})) / live.s(k) - 1; end
fprintf('moderated s of the existing features moved by %+.1f%% .. %+.1f%% (dof %g -> %g)\n', 100*min(ds), 100*max(ds), live.dof, new.dof);
ok = isempty(rowBad) && isequal(new.Fref, live.Fref);
fprintf('SAFE TO INSTALL: %d\n', ok);
