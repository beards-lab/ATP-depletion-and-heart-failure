% BuildRefit2Seed.m -- second refit attempt, per fjezek 2026-09-21.
%
% RESTART POINT: params/ktrfence_seed.m, NOT ktrfence_opt.
%   Measured on the same objective, opt is better on slack+FV (0.289 vs 0.462).
%   It is still the wrong place to start: that number was earned on an objective
%   with no A/steady term and a ktr fence anchored on the wrong preparation, and
%   opt's advantage sits entirely in peak1_dSL and vall_t, which remain 2.7 and
%   2.5 preparation-widths out and will be re-optimised from scratch anyway.
%   What opt did acquire is a broken A/steady (see below), which the seed does
%   not have. Starting from the seed does not throw away a solved feature.
%
% OBJECTIVE CHANGES
%   1. NEW TERM  A_over_steady|fl:0.01008|1.757
%      The post-slack plateau as a fraction of the restretched plateau -- the
%      thing that reads optically and that scoring A and steady separately
%      leaves free. Earns a tight deadband honestly: preparation strength is
%      common-mode and cancels, so its relative SD across the three preps is
%      0.078 against 0.266 (A) and 0.236 (steady), 3.4x tighter. Floor = s^2/3
%      (pooled, n = 3) and w = lambda/s^2, so w*fl = lambda*nEl like every
%      other chi-square term.
%   2. SLACK FLOORS 1/3 -> 1/4 of the original one-preparation-width value
%      (x 0.75 on the current floors). PRAGMATIC, not derived: s^2/n is the
%      statistically defensible stop, and this goes past it. Justified only as
%      "the level terms should start paying sooner than the measurement's own
%      resolution strictly requires"; label it as a choice when reporting.
%   3. FV FLOOR 1/2 -> 1/4 (x 0.5). Verified safe with margin: F-V sits at 0.46
%      preparation-widths and does not pay until the floor drops below 0.21x its
%      current value, so 0.5x leaves ~2x headroom. It stays inert -- narrowing
%      F-V buys nothing while it is this well fitted, and its weight is 1804,
%      so the only thing narrowing can do is detonate. Going further is unsafe.
%   4. ktr FENCE RETARGETED on the 3-preparation mean, 43.19 (was 03/27's 49.4),
%      with floor s^2/3 = 0.00585 and w = 0.6057. This is no longer inert: the
%      seed is 24 % fast against the pooled mean, so the fence now pulls.
%      Model change in runKtrExperiment.m, gated on UsePooledTarget.
%   5. PS_* floors unchanged at 1 -- single-source targets, n = 1.
%
%   matlab -batch "run('Analyses/RefitStrategy/BuildRefit2Seed.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
if ~exist('SEEDTAG','var'); SEEDTAG = 'refit2_seed'; end

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_r2_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 5); end

p = getParams(); s = loadParams('ktrfence_seed');
f = fieldnames(s); for i = 1:numel(f); p.(f{i}) = s.(f{i}); end
p.RunForceVelocity=true; p.RunSlack=true; p.RunStairs=false; p.RunKtr=true;
p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
p.RunSlackSegments='AllPar'; p.MaxRunTime=150; p.OptimizeOn='Feats';
p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0; p.mods={}; p.g=ones(1,0);

fprintf('=== restart point on the OLD objective ===\n');
c_old = evaluateBakersExp(1, p, true, false);
fprintf('  ktrfence_seed = %.4f\n', c_old);

%% ---- rescale the floors -------------------------------------------------
POOLED = {'ktr','A','t0_crossing','restretchSlopeStart','peak1_y','peak1_dSL', ...
          'vall_y','vall_t','peak2','steady','vall2_dy','rsK80'};
fn = p.fn;
for i = 1:numel(fn)
    tk = split(fn{i},'|'); j = find(startsWith(tk,'fl:'),1);
    if isempty(j); continue; end
    base = regexprep(tk{1}, '\[.*$', '');
    if     strcmp(base,'ktrP_k63');        continue;               % replaced below
    elseif strcmp(base,'FV_normpowerAvg'); sc = 0.50;              % 1/2 -> 1/4
    elseif ismember(base, POOLED);         sc = 0.75;              % 1/3 -> 1/4
    else;                                  sc = 1.00;              % PS_*, n = 1
    end
    if sc == 1; continue; end
    old = str2double(extractAfter(tk{j},'fl:'));
    tk{j} = sprintf('fl:%.6g', old*sc);
    fn{i} = strjoin(tk,'|');
    fprintf('  refloor %-22s x%.2f   %.6g -> %.6g\n', tk{1}, sc, old, old*sc);
end

%% ---- retarget the ktr fence --------------------------------------------
LAMBDA = 0.0106294;
S_KTRP = 0.132;                                   % 3-day spread
FENCE  = sprintf('ktrP_k63|fl:%.4g|%.4g', S_KTRP^2/3, LAMBDA/S_KTRP^2);
for i = 1:numel(fn)
    if startsWith(fn{i},'ktrP_k63'); fn{i} = FENCE; end
end
fprintf('  fence -> %s   (target now the 3-prep mean 43.19)\n', FENCE);

%% ---- the new A/steady term ---------------------------------------------
S_RATIO = 0.078; NEL = 5;
RATIO = sprintf('A_over_steady|fl:%.4g|%.4g', (1/3)*NEL*S_RATIO^2, LAMBDA/S_RATIO^2);
fn{end+1} = RATIO;
fprintf('  new   -> %s   (w*fl = %.6g = lambda*nEl)\n', RATIO, ...
        (LAMBDA/S_RATIO^2)*((1/3)*NEL*S_RATIO^2));
p.fn = fn;

%% ---- pool: unchanged, plus nothing new ---------------------------------
sl = loadParams('lean_seed');
pool = sl.leanPool;
if ~ismember('x_M_slack', pool); pool{end+1} = 'x_M_slack'; end
keep = true(1,numel(pool));
for i = 1:numel(pool)
    nm = pool{i}; us = strfind(nm,'__');
    if isempty(us)
        keep(i) = isfield(p,nm) && isscalar(p.(nm)) && p.(nm) ~= 0;
    else
        b = nm(1:us(1)-1); ix = str2double(nm(us(1)+2:end));
        keep(i) = isfield(p,b) && numel(p.(b)) >= ix && p.(b)(ix) ~= 0;
    end
end
p.leanPool = pool(keep);
fprintf('  pool: %d parameters\n', numel(p.leanPool));

%% ---- score --------------------------------------------------------------
rep = checkFeatureCoverage(p);
if ~rep.ok
    fprintf(2,'UNSCORABLE:\n');
    if ~isempty(rep.missingModel); fprintf(2,'  model: %s\n', strjoin(rep.missingModel,', ')); end
    if ~isempty(rep.missingData);  fprintf(2,'  data : %s\n', strjoin(rep.missingData ,', ')); end
    error('BuildRefit2Seed:coverage','objective has unscorable terms');
end
fprintf('\n=== seed on the NEW objective ===\n');
[c_new,~,ct,nm] = evaluateBakersExp(1, p, true, false);
fprintf('  %s = %.4f   (was %.4f on the old one -- NOT comparable)\n', SEEDTAG, c_new, c_old);
[~,o] = sort(ct,'descend');
for i = o(1:min(12,numel(o)))
    if ct(i) <= 0; continue; end
    fprintf('    %-24s %8.4f  (%4.1f %%)\n', nm{i}, ct(i), 100*ct(i)/c_new);
end
for q = {'A_over_steady','ktrP_k63','FV_normpowerAvg','A','steady'}
    i = find(strcmp(nm,q{1}),1);
    if ~isempty(i); fprintf('  %-22s %8.4f\n', q{1}, ct(i)); end
end

% and what the OLD optimum would score on it, for the record
q2 = getParams(loadParams('ktrfence_opt'), [], true, false);
q2.fn = p.fn; q2.RunKtr = 1; q2.PlotEachSeparately=0; q2.ghostSave=''; q2.ghostLoad='';
q2.EvalFeatures=1; q2.BreakOnODEUnstable=0; q2.SaveBest=0; q2 = getParams(q2,q2.g,true);
[c_opt,~,ct2,nm2] = evaluateBakersExp(1, q2, true, false);
iR = find(strcmp(nm2,'A_over_steady'),1);
fprintf('\n  ktrfence_opt on the SAME new objective = %.4f  (A_over_steady %.4f)\n', c_opt, ct2(iR));

%% ---- write --------------------------------------------------------------
out = fullfile(root,'params',[SEEDTAG '.m']);
iWrite17g(out, p, sprintf(['ktrfence_seed + A_over_steady term + slack floors x0.75 + FV floor x0.5 ' ...
    '+ ktr fence retargeted on the 3-prep mean 43.19. cost %.6f'], c_new));
fprintf('\nwritten: %s\n', out);
q3 = getParams(); s3 = loadParams(out);
f3 = fieldnames(s3); for i=1:numel(f3); q3.(f3{i}) = s3.(f3{i}); end
nmf = fieldnames(p); worst = 0;
for i=1:numel(nmf)
    a = p.(nmf{i});
    if ~isnumeric(a)||~isscalar(a)||a==0||~isfinite(a); continue; end
    if ~isfield(q3,nmf{i}); continue; end
    worst = max(worst, abs(q3.(nmf{i})-a)/abs(a));
end
fprintf('round-trip worst relative error: %.3g\n', worst);
fprintf('\nDONE\n');

function iWrite17g(fname, params, comment)
fid = fopen(fname, 'w');
fprintf(fid, '%% Generated file to manipulate the simulation parameters directly\r\n');
fprintf(fid, '%% %s\r\n', comment);
nm = fieldnames(params);
for i = 1:numel(nm)
    v = params.(nm{i});
    if isa(v,'function_handle') || isstruct(v); continue; end
    if isnumeric(v) && isscalar(v);       fprintf(fid, 'params0.%s = %.17g;\r\n', nm{i}, v);
    elseif islogical(v) && isscalar(v);   fprintf(fid, 'params0.%s = %d;\r\n', nm{i}, v);
    elseif ischar(v) && ~isempty(v);      fprintf(fid, 'params0.%s = ''%s'';\r\n', nm{i}, v);
    elseif iscell(v) && ~isempty(v);      fprintf(fid, 'params0.%s = {''%s''};\r\n', nm{i}, strjoin(v, ''', '''));
    elseif isnumeric(v) && ~isempty(v)
        fprintf(fid, 'params0.%s = [%s];\r\n', nm{i}, ...
                strjoin(arrayfun(@(x) sprintf('%.17g',x), v(:)', 'uni', 0), ' '));
    end
end
fprintf(fid, 'params0.mods = {};\r\nparams0.g = [];\r\n');
fclose(fid);
end
