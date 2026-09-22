% RunReproCheck.m -- is the objective reproducible at all?
%
% The refit3 campaign reported 0.4931 for its snapshot #1; re-evaluating that
% same file gave 0.5091 (+3.2 %), which is an order of magnitude larger than a
% round's gain (-0.3 %). If the objective is not reproducible to better than
% the step size, no campaign result means anything.
%
% Three candidate causes, tested here:
%   R1  non-determinism      -- same file, same pool, evaluated twice
%   R2  worker-count effect  -- runSlackExperiment splits the protocol into
%                               parallel chunks; if the SPLIT depends on the
%                               pool the ODE restarts at different points
%   R3  serial vs parallel   -- with no pool at all the whole velocitytable may
%                               run as ONE continuous integration instead of 5
%   R4  %g snapshot loss     -- writeParamsToMFile writes 6 significant digits;
%                               compare against a %.17g rewrite of the same params
%
%   matlab -batch "run('Analyses/RefitStrategy/RunReproCheck.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
TAG = 'refit3_opt';


fprintf('\n== R3: no pool at all ==\n');
delete(gcp('nocreate'));
c_nopool = ev(TAG);
fprintf('  no pool            %.6f\n', c_nopool);

for NW = [3 5]
    delete(gcp('nocreate'));
    cl = parcluster('local');
    jsl = fullfile(tempdir, sprintf('mljob_rc%d_%d', NW, feature('getpid')));
    if ~isfolder(jsl); mkdir(jsl); end
    cl.JobStorageLocation = jsl;
    if cl.NumWorkers < NW; cl.NumWorkers = NW; end
    parpool(cl, NW);
    a = ev(TAG); b = ev(TAG);
    fprintf('\n== %d workers ==\n  eval 1 %.6f\n  eval 2 %.6f   (repeat delta %.2e)\n', NW, a, b, abs(a-b));
end

fprintf('\n== R4: %%g snapshot loss ==\n');
% rewrite the same params at full precision and compare
p = getParams(loadParams(TAG), [], true, false);
tmp = fullfile(tempdir,'refit3_opt_17g.m');
fid = fopen(tmp,'w'); fprintf(fid,'%% full-precision rewrite\r\n');
nm = fieldnames(p);
for i=1:numel(nm)
    v = p.(nm{i});
    if isa(v,'function_handle')||isstruct(v); continue; end
    if isnumeric(v)&&isscalar(v);     fprintf(fid,'params0.%s = %.17g;\r\n',nm{i},v);
    elseif islogical(v)&&isscalar(v); fprintf(fid,'params0.%s = %d;\r\n',nm{i},v);
    elseif ischar(v)&&~isempty(v);    fprintf(fid,'params0.%s = ''%s'';\r\n',nm{i},v);
    elseif iscell(v)&&~isempty(v);    fprintf(fid,'params0.%s = {''%s''};\r\n',nm{i},strjoin(v,''', '''));
    elseif isnumeric(v)&&~isempty(v)
        fprintf(fid,'params0.%s = [%s];\r\n',nm{i}, ...
                strjoin(arrayfun(@(x) sprintf('%.17g',x),v(:)','uni',0),' '));
    end
end
fprintf(fid,'params0.mods = {};\r\nparams0.g = [];\r\n'); fclose(fid);
addpath(tempdir);
fprintf('  %%g file           %.6f\n', ev(TAG));
fprintf('  %%.17g rewrite     %.6f\n', ev('refit3_opt_17g'));
fprintf('  (both loaded from disk; any difference here IS the serialisation)\n');

fprintf('\n== for reference ==\n  refit3_seed        %.6f   (campaign said 0.4945)\n', ev('refit3_seed'));
fprintf('  campaign said snapshot #1 = 0.4931, best after 8 rounds = 0.4214\n');
fprintf('\nDONE\n');

function c = ev(tag_or_p)
    if ischar(tag_or_p); p = getParams(loadParams(tag_or_p), [], true, false);
    else;                p = tag_or_p; end
    p.PlotEachSeparately=0; p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0;
    p.EvalFeatures=1; p.BreakOnODEUnstable=0; p.MaxRunTime=600;
    p = getParams(p, p.g, true);
    c = evaluateBakersExp([], p, false, false);
end
