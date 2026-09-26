%% UpdatePooledTargetVertex.m -- add peak1_tq / vall_tq / peak1_dSLq to the frozen 8 mM target.
%
% The target is re-extracted from the raw preparations by poolLoadRuns and
% frozen by WritePooledTarget. Adding features must not move anything that is
% already in it, so: back up the run cache, re-extract, write the new target to
% a STAGED file, compare every existing field of every existing feature with the
% live Model/pooledTarget8mM.m, and replace the live file only if they agree.
%   matlab -batch "run('Analyses/OverstretchedHeads/UpdatePooledTargetVertex.m')"
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
pdir = fullfile(root, 'Analyses', 'DataPooling', 'results');
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmm'));
cache = fullfile(pdir, 'pool_raw.mat');
if isfile(cache), copyfile(cache, fullfile(pdir, ['pool_raw_' stamp '.mat.bak'])); end
live = fullfile(root, 'Model', 'pooledTarget8mM.m');
copyfile(live, fullfile(pdir, ['pooledTarget8mM_' stamp '.m.bak']));

OLD = pooledTarget8mM();                              % the live target, before anything changes
poolLoadRuns(true);                                   % re-extract with the vertex features
stagedDir = tempname; mkdir(stagedDir);
OUTF = fullfile(stagedDir, 'pooledTarget8mM.m');
save(fullfile(stagedDir, 'state.mat'), 'OLD', 'stamp', 'live', 'pdir', 'root');   % the writer runs clearvars -except OPTS OUTF
run(fullfile(root, 'Analyses', 'DataPooling', 'WritePooledTarget.m'));
stagedDir = fileparts(OUTF); load(fullfile(stagedDir, 'state.mat')); cd(root);
addpath(stagedDir, '-begin'); clear pooledTarget8mM; NEW = pooledTarget8mM(); rmpath(stagedDir); clear pooledTarget8mM;

% The objective reads only VAL (pooledTargetFeatures swaps values; floors and
% weights are numbers in the fn strings), so VAL must not move at all. The
% moderated spreads (s, sraw, ...) share one variance prior across ALL
% features, so adding features shifts them slightly (<~2 %): reported, not fatal.
worst = 0; ds = 0; bad = {};
for k = 1:numel(OLD.name)
    j = find(strcmp(NEW.name, OLD.name{k}));
    if isempty(j), bad{end+1} = [OLD.name{k} ' (missing)']; continue; end %#ok<SAGROW>
    d = max(abs(OLD.val(k, :) - NEW.val(j, :)) ./ max(abs(OLD.val(k, :)), eps), [], 'omitnan');
    worst = max(worst, d);
    if d > 1e-9, bad{end+1} = sprintf('%s.val (rel %.2g)', OLD.name{k}, d); end %#ok<SAGROW>
    ds = max(ds, abs(OLD.s(k) - NEW.s(j)) / OLD.s(k));
end
fprintf('\nmoderated s of existing features: max relative drift %.3g (informational)\n', ds);
new = setdiff(NEW.name, OLD.name, 'stable');
fprintf('\nexisting features: worst relative change %.3g over %d features\n', worst, numel(OLD.name));
fprintf('added: %s\n', strjoin(new, ', '));
for q = new
    j = find(strcmp(NEW.name, q{1}));
    fprintf('   %-12s val %s  s %.3f\n', q{1}, mat2str(NEW.val(j, :), 5), NEW.s(j));
end
if isempty(bad)
    copyfile(OUTF, live);
    fprintf('\nlive target REPLACED (%s); backups *_%s.bak in Analyses/DataPooling/results\n', live, stamp);
else
    fprintf(2, '\nNOT replaced -- existing values moved:\n  %s\nstaged file left at %s\n', strjoin(bad, '\n  '), OUTF);
end
