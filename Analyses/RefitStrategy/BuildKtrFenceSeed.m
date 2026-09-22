% BuildKtrFenceSeed.m -- the refit seed, per the three decisions of 2026-09-18.
%
% BASE  params/lean_opt.m  (the MechanismAblation output, lean objective:
%       floored, POOLED slack target, rsK80, no coolDownLS, FV share 0.25)
%
% 1. MODEL -- both ktr fixes, which are FREE on this objective
%      UseMaxwellTensionOnly 0 -> 1, and x_M_slack 0 -> 0.0108565 WITH it
%        (dPUdT_CombinedTransitions.m:322 divides by x_M_slack; the ablation
%         deleted the mechanism because the ktr protocol was not in the
%         objective that measured it as free. It costs 13.9x the data noise.)
%      tau_reg 26.6779 ms -> 0.2 ms   (Analyses/KtrProtocolVsSlackKtr)
%    Together: Et 0.5311 -> 0.5169 AND held-out ktr 13.94x -> 1.27x the noise.
%
% 2. DEADBAND -- floor = s^2/n, n = how many sources the term's target averages.
%    The floors were set at FLOOR_CHI2 = 1 when the target was ONE preparation.
%    UsePooledTarget = 1 made the slack target the MEAN OF THREE, and the
%    uncertainty on a mean of three is s/sqrt(3). So:
%      pooled slack terms  n = 3  -> x 1/3
%      FV_normpowerAvg     n = 2  -> x 1/2   (it averages two lab curves)
%      passive PS_*        n = 1  -> x 1     (single source, unchanged)
%    Measured safe zone: F-V stays inert down to FLOOR ~ 0.12 and takes 69 % of
%    the objective at 0.10, because its pinned weight is 1749x the median
%    w*floor. 1/3 sits well inside that. Live data terms 6 -> 9 of 16.
%
% 3. ktr FENCE -- 'ktrP_k63|fl:0.0174|0.610'
%    s = 0.132 (3-day spread 49.4/42.1/38.1), floor = 1 prep-width = s^2,
%    weight = lambda/s^2 with lambda = 0.0106294, so w*fl = lambda exactly:
%    commensurate with every other term. The model sits 0.12 prep-widths off,
%    so it pays NOTHING at the seed. It is a fence against tau_reg drifting
%    back, not a fitting target. The ktr TRACE stays fully held out.
%
% Pool comes from lean_seed.m: lean_opt.m's leanPool field was overwritten with
% the FEATURE list by the snapshot writer and is unusable.
%
%   matlab -batch "run('Analyses/RefitStrategy/BuildKtrFenceSeed.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
if ~exist('SEEDTAG','var'); SEEDTAG = 'ktrfence_seed'; end
if ~exist('NWORKERS','var'); NWORKERS = 5; end

c = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_seed_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
c.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(c, NWORKERS); end

p = getParams(); s = loadParams('lean_opt');
f = fieldnames(s); for i = 1:numel(f); p.(f{i}) = s.(f{i}); end
p.RunForceVelocity=true; p.RunSlack=true; p.RunStairs=false;
p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
p.RunSlackSegments='AllPar'; p.MaxRunTime=150; p.OptimizeOn='Feats';
p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0;
p.mods = {}; p.g = ones(1,0);

fnOld = p.fn;
fprintf('=== base lean_opt on its own objective ===\n');
c_base = evaluateBakersExp(1, p, true, false);
fprintf('  lean_opt = %.4f\n', c_base);

%% ---- 1. model fixes -----------------------------------------------------
p.UseMaxwellTensionOnly = 1;
p.x_M_slack             = 0.0108565;
p.tau_reg               = 2e-4;
p.RunKtr                = 1;                 % needed to PRODUCE ktrP_k63
assert(~isempty(p.ktr_velocitytableonfile), 'no ktr protocol file on the paramset');

%% ---- 2. re-floor --------------------------------------------------------
% n = sources averaged into each term's target. pooledTargetFeatures.m states
% it replaces the slack features and deliberately leaves FV / ktr / PS_* alone.
POOLED = {'ktr','A','t0_crossing','restretchSlopeStart','peak1_y','peak1_dSL', ...
          'vall_y','vall_t','peak2','steady','vall2_dy','rsK80'};
fn = p.fn; nScaled = 0;
for i = 1:numel(fn)
    tk = split(fn{i},'|');
    j  = find(startsWith(tk,'fl:'), 1);
    if isempty(j); continue; end
    base = regexprep(tk{1}, '\[.*$', '');           % strip a [1:4] selector
    if     strcmp(base,'FV_normpowerAvg'); nSrc = 2;
    elseif ismember(base, POOLED);         nSrc = 3;
    else;                                  nSrc = 1;
    end
    old = str2double(extractAfter(tk{j},'fl:'));
    new = old / nSrc;
    tk{j} = sprintf('fl:%.6g', new);
    fn{i} = strjoin(tk,'|');
    fprintf('  refloor %-22s n=%d   %.6g -> %.6g\n', tk{1}, nSrc, old, new);
    nScaled = nScaled + 1;
end
fprintf('  %d floors rescaled\n', nScaled);

%% ---- 3. the ktr fence ---------------------------------------------------
LAMBDA  = 0.0106294;              % w*floor for a chi-square-consistent term
S_KTRP  = 0.132;                  % 3-day spread of 49.4 / 42.1 / 38.1
FENCE   = sprintf('ktrP_k63|fl:%.4g|%.4g', S_KTRP^2, LAMBDA/S_KTRP^2);
fn{end+1} = FENCE;
fprintf('  fence: %s   (w*fl = %.6g = lambda)\n', FENCE, LAMBDA);
p.fn = fn;

%% ---- 4. pool ------------------------------------------------------------
sl = loadParams('lean_seed');
pool = sl.leanPool;
if ~ismember('x_M_slack', pool); pool{end+1} = 'x_M_slack'; end   % mechanism is live again
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
if any(~keep)
    fprintf('  dropping %d zero-base params: %s\n', nnz(~keep), strjoin(pool(~keep),', '));
end
p.leanPool = pool(keep);
fprintf('  pool: %d parameters (tau_reg in pool: %d, x_M_slack in pool: %d)\n', ...
        numel(p.leanPool), ismember('tau_reg',p.leanPool), ismember('x_M_slack',p.leanPool));

%% ---- 5. score and verify -----------------------------------------------
rep = checkFeatureCoverage(p);
if ~rep.ok
    fprintf(2,'UNSCORABLE fn entries:\n');
    if ~isempty(rep.missingModel); fprintf(2,'  model: %s\n', strjoin(rep.missingModel,', ')); end
    if ~isempty(rep.missingData);  fprintf(2,'  data : %s\n', strjoin(rep.missingData ,', ')); end
    error('BuildKtrFenceSeed:coverage','objective has unscorable terms');
end
fprintf('\n=== seed on the NEW objective ===\n');
t0 = tic; [c_new,~,ct,nm] = evaluateBakersExp(1, p, true, false); t_eval = toc(t0);
fprintf('  %s = %.4f     (one evaluation: %.1f s)\n', SEEDTAG, c_new, t_eval);
[~,o] = sort(ct,'descend');
for i = o(1:min(10,numel(o)))
    if ct(i) <= 0; continue; end
    fprintf('    %-24s %8.4f  (%4.1f %%)\n', nm{i}, ct(i), 100*ct(i)/c_new);
end
iF = find(strcmp(nm,'ktrP_k63'),1);
fprintf('  fence term ktrP_k63 = %.6f  (inert at the seed if 0)\n', ct(iF));

% cost of having RunKtr on at all
q = p; q.RunKtr = 0; q.fn = fn(1:end-1);
t1 = tic; evaluateBakersExp(1, q, true, false); t_noktr = toc(t1);
fprintf('  evaluation without the ktr protocol: %.1f s  -> the fence costs %+.0f %% per evaluation\n', ...
        t_noktr, 100*(t_eval/t_noktr - 1));

%% ---- 6. write -----------------------------------------------------------
out = fullfile(root,'params',[SEEDTAG '.m']);
iWrite17g(out, p, sprintf(['lean_opt + MaxwellTensionOnly restored (x_M_slack %.6g) + tau_reg %.4g; ' ...
    'floors rescaled to s^2/n (pooled 3, FV 2, PS_* 1); ktr fence %s; RunKtr on. ' ...
    'cost %.6f (lean_opt on the old objective %.6f)'], ...
    p.x_M_slack, p.tau_reg, FENCE, c_new, c_base));
fprintf('\nwritten: %s\n', out);

q2 = getParams(); s2 = loadParams(out);
f2 = fieldnames(s2); for i = 1:numel(f2); q2.(f2{i}) = s2.(f2{i}); end
nmf = fieldnames(p); worst = 0; wname = '';
for i = 1:numel(nmf)
    a = p.(nmf{i});
    if ~isnumeric(a) || ~isscalar(a) || a == 0 || ~isfinite(a); continue; end
    if ~isfield(q2, nmf{i}); continue; end
    r = abs(q2.(nmf{i}) - a)/abs(a);
    if r > worst; worst = r; wname = nmf{i}; end
end
fprintf('round-trip worst relative error: %.3g (%s)\n', worst, wname);
V = validateKtrProtocol(SEEDTAG);
fprintf('\nDONE\n');

function iWrite17g(fname, params, comment)
fid = fopen(fname, 'w');
fprintf(fid, '%% Generated file to manipulate the simulation parameters directly\r\n');
fprintf(fid, '%% %s\r\n', comment);
nm = fieldnames(params);
for i = 1:numel(nm)
    v = params.(nm{i});
    if isa(v,'function_handle') || isstruct(v); continue; end
    if isnumeric(v) && isscalar(v)
        fprintf(fid, 'params0.%s = %.17g;\r\n', nm{i}, v);
    elseif islogical(v) && isscalar(v)
        fprintf(fid, 'params0.%s = %d;\r\n', nm{i}, v);
    elseif ischar(v) && ~isempty(v)
        fprintf(fid, 'params0.%s = ''%s'';\r\n', nm{i}, v);
    elseif iscell(v) && ~isempty(v)
        fprintf(fid, 'params0.%s = {''%s''};\r\n', nm{i}, strjoin(v, ''', '''));
    elseif isnumeric(v) && ~isempty(v)
        fprintf(fid, 'params0.%s = [%s];\r\n', nm{i}, ...
                strjoin(arrayfun(@(x) sprintf('%.17g',x), v(:)', 'uni', 0), ' '));
    end
end
fprintf(fid, 'params0.mods = {};\r\nparams0.g = [];\r\n');
fclose(fid);
end
