% RunReproDiagnostic.m -- why does a deterministic model not reproduce a cost?
%
% Three questions, three measurements:
%
%   A. REPEATABILITY. Evaluate the SAME params struct N times in one session on
%      an idle machine. A deterministic objective gives spread 0. Anything else
%      is non-determinism that has nothing to do with storage.
%
%   B. STORAGE ROUND-TRIP. writeParamsToMFile serialises scalars with '%g',
%      which is SIX significant digits. Write the incumbent, read it back,
%      evaluate. This is the user's hypothesis: does the snapshot round-trip
%      alone move the cost?
%
%   C. THE FIX. Same round-trip at '%.17g' (exact double round-trip). If B
%      moves and C does not, the format is the bug and the fix is one line.
%
% Also reports WHICH parameters lose the most precision at %g, since that says
% whether any single parameter is badly conditioned for storage.
%
%   matlab -batch "run('Analyses/MechanismAblation/RunReproDiagnostic.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

NREP     = 4;
SEEDSNAP = 'floored_seed';
OUTDIR   = fullfile(here, 'results'); if ~isfolder(OUTDIR); mkdir(OUTDIR); end
scratch  = tempdir;

base = getParams();
snap = loadParams(SEEDSNAP);
f = fieldnames(snap); for i = 1:numel(f); base.(f{i}) = snap.(f{i}); end
base.RunForceVelocity = true;  base.RunKtr = false;   base.RunSlack = true;
base.RunStairs = false;        base.RunForceVelocityTime = false;
base.RunSlackPassive = true;   base.EvalFeatures = true;
base.BreakOnODEUnstable = false;
base.PlotEachSeparately = 0;   base.PlotFeatureFitting = 0;
base.RunSlackSegments = 'AllPar';
base.MaxRunTime = 150;         base.OptimizeOn = 'Feats';
base.mods = {}; base.g = ones(1,0);

if isempty(gcp('nocreate')); parpool('local', 5); end

%% ---- A. repeatability ---------------------------------------------------
fprintf('\n=== A. same struct, %d evaluations, idle machine ===\n', NREP);
cA = nan(1,NREP);
for i = 1:NREP
    t = tic; cA(i) = evaluateBakersExp(1, base, true, false);
    fprintf('  rep %d: %.10f   (%.0f s)\n', i, cA(i), toc(t));
end
fprintf('  spread: max-min = %.3g,  std = %.3g\n', max(cA)-min(cA), std(cA));

%% ---- B/C. storage round-trip -------------------------------------------
% writeParamsToMFile is not parameterised by format, so write both variants
% here directly from the same struct and compare.
fprintf('\n=== B/C. snapshot round-trip ===\n');
names = fieldnames(base);
cRT = nan(1,2); fmts = {'%g','%.17g'};
for k = 1:2
    fn_out = fullfile(scratch, sprintf('abl_rt_%d.m', k));
    fid = fopen(fn_out,'w');
    fprintf(fid, '%% round-trip probe, format %s\r\n', fmts{k});
    for i = 1:numel(names)
        v = base.(names{i});
        if isa(v,'function_handle') || isstruct(v); continue; end
        if isnumeric(v) && isscalar(v)
            fprintf(fid, ['params0.%s = ' fmts{k} ';\r\n'], names{i}, v);
        elseif ischar(v) && ~isempty(v)
            fprintf(fid, 'params0.%s = ''%s'';\r\n', names{i}, v);
        elseif iscell(v)
            fprintf(fid, 'params0.%s = {''%s''};\r\n', names{i}, strjoin(v, ''', '''));
        elseif isnumeric(v) && ~isempty(v)
            fprintf(fid, 'params0.%s = [%s];\r\n', names{i}, ...
                    strjoin(arrayfun(@(x) sprintf(fmts{k},x), v(:)', 'uni', 0), ' '));
        end
    end
    fprintf(fid, 'params0.mods = {};\r\nparams0.g = [];\r\n'); fclose(fid);

    rt = getParams(); s2 = loadParams(fn_out);
    f2 = fieldnames(s2); for i = 1:numel(f2); rt.(f2{i}) = s2.(f2{i}); end
    rt.RunForceVelocity=true; rt.RunKtr=false; rt.RunSlack=true; rt.RunStairs=false;
    rt.RunForceVelocityTime=false; rt.RunSlackPassive=true; rt.EvalFeatures=true;
    rt.BreakOnODEUnstable=false; rt.PlotEachSeparately=0; rt.PlotFeatureFitting=0;
    rt.RunSlackSegments='AllPar'; rt.MaxRunTime=150; rt.OptimizeOn='Feats';
    rt.mods = {}; rt.g = ones(1,0);
    cRT(k) = evaluateBakersExp(1, rt, true, false);
    fprintf('  %-6s -> %.10f   (delta vs rep1 %+.3g)\n', fmts{k}, cRT(k), cRT(k)-cA(1));
end

%% ---- which params lose precision at %g ---------------------------------
fprintf('\n=== worst relative precision loss at %%g (6 s.f.) ===\n');
nm = {}; rel = [];
for i = 1:numel(names)
    v = base.(names{i});
    if ~isnumeric(v) || ~isscalar(v) || v == 0 || ~isfinite(v); continue; end
    back = str2double(sprintf('%g', v));
    r = abs(back - v)/abs(v);
    if r > 0; nm{end+1} = names{i}; rel(end+1) = r; end %#ok<SAGROW>
end
[rel, o] = sort(rel,'descend'); nm = nm(o);
for i = 1:min(12,numel(nm))
    fprintf('  %-28s %.10g  ->  rel err %.2e\n', nm{i}, base.(nm{i}), rel(i));
end
fprintf('  %d of %d scalars lose precision; worst %.2e, median %.2e\n', ...
        numel(nm), numel(names), max(rel), median(rel));

save(fullfile(OUTDIR,'repro_diagnostic.mat'), 'cA','cRT','nm','rel');
fprintf('\nwritten: %s\n', fullfile(OUTDIR,'repro_diagnostic.mat'));
