% BuildRefit6Seed.m -- re-price ktr_shape, and give the optimiser the four
% mechanosensing parameters it has never been able to reach.
%
% 1. ktr_shape at s := 1 x s_shape (was 2 x): fl 0.0204 -> 0.0051,
%    w 0.696 -> 2.784, so the slack ladder's slope goes from 4 % of the
%    objective to ~18 %. refit5 ran 6.3 h with three cheap slope levers
%    compulsory and closed 5 % of the gap; the levers exist, so it is a pricing
%    problem, and this is the price.
%
% 2. THE SRX CYCLE HAS TWO DEAD BRANCHES.
%
%        PT -> SR    ksr0 *exp(-F/sigma2)      367.5   live
%        SR -> PT    kmsr *exp(+F/sigma1)      0       DEAD  (sigma1 inert)
%        SRD-> PD    kmsrd*exp(+F/sigma_srd1)  26.8    live
%        PD -> SRD   ksrd *exp(-F/sigma_srd2)  0       DEAD  (sigma_srd2 inert)
%
%    Half the force-gated transitions are switched off, and the two rate
%    constants are ZERO -- which is exactly what optimizeFeatures' zero-base
%    filter drops from the pool, so no campaign of any length can revive them.
%    That is why `sigma1 x0.5` and `sigma1 x2` both scored EXACTLY the baseline
%    in the sensitivity sweep: dead code.
%
%    The structural probe found `kmsr = 50` raises the ktr LEVEL from 47.8 to
%    52.0 (pooled data segment 1 = 51.6) at **zero F-V cost and zero re-stretch
%    cost** -- the only lever in the whole probe that did not sell something
%    else. Every switch that bought slope (UsePassiveForSR, P_bound_max <= 0.6)
%    cost 100+ units of F-V or 6x the re-stretch block.
%
%    This script seeds both dead rates at the SMALLEST value that is clearly
%    nonzero and cheap, so the mechanism is merely *available*, not imposed,
%    and adds kmsr / ksrd / sigma1 / sigma_srd2 to leanPool.
%
% WHY THIS IS THE RIGHT BET FOR THE LOW-ATP CALL. The 2 mM signature is
% ktr x0.54 against force x1.18 -- ktr is the ATP-sensitive observable, and
% SRX turnover is where ATP enters. A model whose force-gated SRX exit is
% identically zero cannot express an ATP effect through that route at all, so
% whatever the 2 mM fit is doing today, it is doing it with the wrong mechanism.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/BuildRefit6Seed.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~exist('SEEDTAG','var'); SEEDTAG = 'refit6_seed'; end
LAMBDA = 0.0106294;

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_r6_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 6); end

p0 = getParams(loadParams('refit5_opt'), [], true, false);
p0.RunForceVelocity=true; p0.RunSlack=true; p0.RunStairs=false; p0.RunKtr=true;
p0.RunForceVelocityTime=false; p0.RunSlackPassive=true; p0.EvalFeatures=true;
p0.BreakOnODEUnstable=false; p0.PlotEachSeparately=0; p0.PlotFeatureFitting=0;
p0.RunSlackSegments='AllPar'; p0.MaxRunTime=300; p0.OptimizeOn='Feats';
p0.ghostSave=''; p0.ghostLoad=''; p0.SaveBest=0; p0.mods={}; p0.g=ones(1,0);

%% ---- 1. re-price ktr_shape ---------------------------------------------
DAYS = {'03_27','04_03','04_10'};
Xk = zeros(3,5);
for d = 1:3
    S = load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d})));
    Xk(d,:) = S.features_data.ktr(:)'/S.features_data.ktr(1);
end
Xk = Xk(:,2:5);
s_shape = mean(std(Xk,0,1,'omitnan')./abs(mean(Xk,1,'omitnan')),'omitnan');
fl_new = (1/3)*4*s_shape^2; w_new = LAMBDA/s_shape^2;
fn = p0.fn;
for i = 1:numel(fn)
    if startsWith(fn{i},'ktr_shape')
        fn{i} = sprintf('ktr_shape[2:5]|fl:%.6g|%.6g', fl_new, w_new);
    end
end
p0.fn = fn;
fprintf('ktr_shape: s %.4f | fl %.6f | w %.4f\n', s_shape, fl_new, w_new);

RESTRETCH = {'peak1_y','peak1_dSL','vall_y','vall_t','vall2_dy','peak2', ...
             'restretchSlopeStart','rsK80'};

%% ---- 2. how cheap is it to switch the dead branches on? ----------------
ACT = { 'kmsr', 0; 'kmsr', 1; 'kmsr', 2; 'kmsr', 5; 'kmsr', 10; ...
        'ksrd', 1; 'ksrd', 2; 'ksrd', 5; 'ksrd', 10 };
nA = size(ACT,1);
AC = nan(nA,1); AR = nan(nA,1); AK = nan(nA,1); AF = nan(nA,1); ASL = nan(nA,1);
parfor i = 1:nA
    p = p0; p.(ACT{i,1}) = ACT{i,2};
    p = getParams(p, p.g, true, false);
    p.RunSlackSegments = 'All';
    try
        [c,~,ct,nm] = evaluateBakersExp(1, p, true, false);
        AC(i) = c;
        AR(i) = sum(ct(ismember(nm, RESTRETCH)));
        j = find(strcmp(nm,'ktr_shape'),1);       if ~isempty(j); AK(i) = ct(j); end
        j = find(strcmp(nm,'FV_normpowerAvg'),1); if ~isempty(j); AF(i) = ct(j); end
        [~,~,fm] = runSlackExperiment(p);
        bb = polyfit(1:5, fm.ktr(:)', 1); ASL(i) = bb(1);
    catch
    end
end
fprintf('\n=== cost of switching a dead SRX branch on ===\n');
fprintf('%-16s %9s %10s %10s %8s %9s\n','setting','TOTAL','restretch','ktr_shape','FV','slope');
for i = 1:nA
    fprintf('%-16s %9.4f %10.4f %10.4f %8.4f %+9.3f\n', ...
            sprintf('%s = %g', ACT{i,1}, ACT{i,2}), AC(i), AR(i), AK(i), AF(i), ASL(i));
end

% pick the cheapest CLEARLY-NONZERO activation for each: we want the mechanism
% available to the search, not imposed by the seed. Accept <= 10 % over the
% kmsr = 0 baseline and no re-stretch or F-V regression worth the name.
base = AC(1);
ik = find(strcmp(ACT(:,1),'kmsr') & cell2mat(ACT(:,2)) > 0 & AC <= base*1.10 & AF < 1e-6, 1, 'last');
if isempty(ik); ik = 2; end
kmsr_seed = ACT{ik,2};
is = find(strcmp(ACT(:,1),'ksrd') & AC <= base*1.10 & AF < 1e-6, 1, 'last');
ksrd_seed = 0; if ~isempty(is); ksrd_seed = ACT{is,2}; end
fprintf('\nseeding kmsr = %g, ksrd = %g\n', kmsr_seed, ksrd_seed);

p0.kmsr = kmsr_seed;
if ksrd_seed > 0; p0.ksrd = ksrd_seed; end
p0 = getParams(p0, p0.g, true, false);

%% ---- 3. extend the pool -------------------------------------------------
pool = p0.leanPool;
NEW = {'kmsr','sigma1'};
if ksrd_seed > 0; NEW = [NEW, {'ksrd','sigma_srd2'}]; end
for i = 1:numel(NEW)
    if ~ismember(NEW{i}, pool) && isfield(p0, NEW{i}) && p0.(NEW{i}) ~= 0
        pool{end+1} = NEW{i}; %#ok<SAGROW>
    end
end
p0.leanPool = pool;
fprintf('pool: %d parameters (added: %s)\n', numel(pool), strjoin(NEW,', '));

%% ---- 4. score and write -------------------------------------------------
rep = checkFeatureCoverage(p0);
if ~rep.ok
    fprintf(2,'UNSCORABLE:\n');
    if ~isempty(rep.missingModel); fprintf(2,'  model: %s\n', strjoin(rep.missingModel,', ')); end
    if ~isempty(rep.missingData);  fprintf(2,'  data : %s\n', strjoin(rep.missingData ,', ')); end
    error('BuildRefit6Seed:coverage','objective has unscorable terms');
end
p0.RunSlackSegments = 'AllPar';
[c_new,~,ct,nm] = evaluateBakersExp(1, p0, true, false);
fprintf('\n=== %s = %.4f ===\n', SEEDTAG, c_new);
[~,o] = sort(ct,'descend');
for i = o(1:min(14,numel(o)))
    if ct(i) <= 0; continue; end
    fprintf('    %-24s %8.4f  (%4.1f %%)\n', nm{i}, ct(i), 100*ct(i)/c_new);
end
fprintf('  re-stretch block: %.4f\n', sum(ct(ismember(nm,RESTRETCH))));

out = fullfile(root,'params',[SEEDTAG '.m']);
iWrite17g(out, p0, sprintf(['refit5_opt params; ktr_shape re-priced at 1x s_shape; ' ...
    'kmsr=%g ksrd=%g seeded and kmsr/sigma1/ksrd/sigma_srd2 added to leanPool. cost %.6f'], ...
    kmsr_seed, ksrd_seed, c_new));
fprintf('\nwritten: %s\n', out);
q = getParams(); s2 = loadParams(out);
f2 = fieldnames(s2); for i=1:numel(f2); q.(f2{i}) = s2.(f2{i}); end
nmf = fieldnames(p0); worst = 0;
for i=1:numel(nmf)
    a = p0.(nmf{i});
    if ~isnumeric(a)||~isscalar(a)||a==0||~isfinite(a); continue; end
    if ~isfield(q,nmf{i}); continue; end
    worst = max(worst, abs(q.(nmf{i})-a)/abs(a));
end
fprintf('round-trip worst relative error: %.3g\n', worst);
fprintf('leanPool round-trips as %d entries (first: %s)\n', numel(q.leanPool), q.leanPool{1});
save(fullfile(resdir,'refit6_seed_build.mat'),'ACT','AC','AR','AK','AF','ASL', ...
     'kmsr_seed','ksrd_seed','s_shape','fl_new','w_new','c_new');
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
