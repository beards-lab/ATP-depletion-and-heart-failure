% BuildRefit4Seed.m -- rebuild the objective on the re-curated pooled target.
%
% UPSTREAM (this session):
%   * extractSlackAttributes' peak prominence was 0.1 kPa = 0.07-0.11 sigma of
%     the recordings' noise, so findpeaks returned 11-15 ripples per re-stretch
%     and peak1_t(1) took the earliest. 7 of 15 preparation x segment cells
%     locked onto the rising limb instead of the overshoot. Now
%     max(3*sigma, 5 % of the window's force range).
%   * all three 8 mM files re-curated through ONE code path (04/03 had never
%     been through UpdateSlackFeatures at all).
%   * pooled target rebuilt. peak1_dSL s 52.4 -> 19.5 %, vall_t 57.1 -> 19.8 %.
%   * buildPooledTarget's fitFDist now returns d0 = Inf -- with the two broken
%     features fixed the remaining spreads are statistically indistinguishable
%     at n = 3 -- so d0 is capped at 5 and the weights come out near-uniform.
%     That is a real consequence: the old differential weighting was partly an
%     artefact of the same bug.
%
% THIS SCRIPT: keep refit3_opt's PARAMETERS (they banked A/steady, the ktr
% fence and a halved held-out staircase), and re-derive every pooled slack
% term's floor and weight from the NEW spreads, with the same rules as before:
%     fl = (1/n) * nEl * s^2       n = 3 for the pooled slack block
%     w  = lambda / s^2            lambda = 0.0106294, so w*fl = lambda*nEl
% Non-pooled terms are untouched: FV (own spread, n = 2), ktrP_k63 (3-day ktr
% spread), A_over_steady (own ratio spread), PS_* (floored at seed x1.2).
%
%   matlab -batch "run('Analyses/RefitStrategy/BuildRefit4Seed.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
if ~exist('SEEDTAG','var'); SEEDTAG = 'refit4_seed'; end
LAMBDA = 0.0106294;

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_r4_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 5); end

p = getParams(); s = loadParams('refit3_opt');
f = fieldnames(s); for i = 1:numel(f); p.(f{i}) = s.(f{i}); end
% leanPool comes back from every optimiser snapshot holding the FEATURE list,
% not parameter names -- writeParamsToMFile wrote params.fn for every cell
% field (now fixed). Take the pool from lean_seed, which predates the damage.
sl = loadParams('lean_seed');
pool = sl.leanPool;
if ~ismember('x_M_slack', pool); pool{end+1} = 'x_M_slack'; end
keep = true(1,numel(pool));
for i = 1:numel(pool)
    nmp = pool{i}; us = strfind(nmp,'__');
    if isempty(us)
        keep(i) = isfield(p,nmp) && isscalar(p.(nmp)) && p.(nmp) ~= 0;
    else
        b = nmp(1:us(1)-1); ix = str2double(nmp(us(1)+2:end));
        keep(i) = isfield(p,b) && numel(p.(b)) >= ix && p.(b)(ix) ~= 0;
    end
end
p.leanPool = pool(keep);
fprintf('pool restored: %d parameters (first few: %s)\n', ...
        numel(p.leanPool), strjoin(p.leanPool(1:min(6,end)),', '));

p.RunForceVelocity=true; p.RunSlack=true; p.RunStairs=false; p.RunKtr=true;
p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
p.RunSlackSegments='AllPar'; p.MaxRunTime=150; p.OptimizeOn='Feats';
p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0; p.mods={}; p.g=ones(1,0);

fprintf('=== refit3_opt on the OLD objective + OLD target, for the record ===\n');
c_before = evaluateBakersExp(1, p, true, false);
fprintf('  %.4f  (this now uses the RE-CURATED data, so it is already not the 0.4271 of before)\n', c_before);

%% ---- re-derive the pooled slack terms ----------------------------------
PT = pooledTarget8mM();
POOLED = {'ktr','A','t0_crossing','restretchSlopeStart','peak1_y','peak1_dSL', ...
          'vall_y','vall_t','peak2','steady','vall2_dy','rsK80'};
% t0_crossing is scored against the pooled 't0'
NAMEMAP = containers.Map({'t0_crossing'}, {'t0'});
fn = p.fn;
fprintf('\n%-22s %8s %10s %10s %10s %10s\n','term','nEl','s_new','fl_old','fl_new','w_new');
for i = 1:numel(fn)
    tk = split(fn{i},'|');
    base = regexprep(tk{1}, '\[.*$', '');
    if ~ismember(base, POOLED); continue; end
    key = base; if isKey(NAMEMAP, base); key = NAMEMAP(base); end
    j = find(strcmp(PT.name, key), 1);
    if isempty(j); warning('no pooled entry for %s', base); continue; end
    snew = PT.s(j);
    if ~isfinite(snew) || snew <= 0; warning('bad s for %s', base); continue; end
    % element count: rsK80[1:4] scores 4, everything else 5
    nEl = 5;
    if contains(tk{1},'[')
        inner = extractBetween(tk{1},'[',']');
        nEl = numel(eval(['[' inner{1} ']']));
    end
    flnew = (1/3) * nEl * snew^2;
    wnew  = LAMBDA / snew^2;
    jf = find(startsWith(tk,'fl:'),1);
    flold = NaN; if ~isempty(jf); flold = str2double(extractAfter(tk{jf},'fl:')); end
    if isempty(jf)
        tk = [tk(1:end-1); {sprintf('fl:%.6g',flnew)}; tk(end)];
    else
        tk{jf} = sprintf('fl:%.6g', flnew);
    end
    tk{end} = sprintf('%.6g', wnew);
    fn{i} = strjoin(tk,'|');
    fprintf('%-22s %8d %10.4f %10.4g %10.4g %10.4g\n', base, nEl, snew, flold, flnew, wnew);
end
p.fn = fn;

%% ---- score --------------------------------------------------------------
rep = checkFeatureCoverage(p);
if ~rep.ok
    fprintf(2,'UNSCORABLE:\n');
    if ~isempty(rep.missingModel); fprintf(2,'  model: %s\n', strjoin(rep.missingModel,', ')); end
    if ~isempty(rep.missingData);  fprintf(2,'  data : %s\n', strjoin(rep.missingData ,', ')); end
    error('BuildRefit4Seed:coverage','objective has unscorable terms');
end
fprintf('\n=== seed on the NEW objective ===\n');
[c_new,~,ct,nm] = evaluateBakersExp(1, p, true, false);
fprintf('  %s = %.4f\n', SEEDTAG, c_new);
[~,o] = sort(ct,'descend');
for i = o(1:min(14,numel(o)))
    if ct(i) <= 0; continue; end
    fprintf('    %-24s %8.4f  (%4.1f %%)\n', nm{i}, ct(i), 100*ct(i)/c_new);
end
fprintf('\n  zero-cost terms: %s\n', strjoin(nm(ct<=1e-9), ', '));

%% ---- write --------------------------------------------------------------
out = fullfile(root,'params',[SEEDTAG '.m']);
iWrite17g(out, p, sprintf(['refit3_opt params on the RE-CURATED pooled target (peak prominence fix); ' ...
    'pooled slack floors/weights re-derived from the new spreads. cost %.6f'], c_new));
fprintf('\nwritten: %s\n', out);
q = getParams(); s2 = loadParams(out);
f2 = fieldnames(s2); for i=1:numel(f2); q.(f2{i}) = s2.(f2{i}); end
nmf = fieldnames(p); worst = 0;
for i=1:numel(nmf)
    a = p.(nmf{i});
    if ~isnumeric(a)||~isscalar(a)||a==0||~isfinite(a); continue; end
    if ~isfield(q,nmf{i}); continue; end
    worst = max(worst, abs(q.(nmf{i})-a)/abs(a));
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
