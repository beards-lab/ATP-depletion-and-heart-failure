% BuildRefit7Seed.m -- rebuild the 8 mM objective around the metrics that carry
% the ATP effect.
%
% WHY (Analyses/KtrSegmentSlope/RunAtpRelevance.m, all 12 files re-curated on
% one code path). Paired 2mM/8mM log-ratio across the three preparations:
%
%     rsK80    x0.38   t = -26.2      re-stretch recovery RATE
%     ktr      x0.53   t = -10.5      hold-redevelopment rate (old estimator)
%     ktr80    x0.57   t =  -8.3      the same, audited estimator
%     rsT20    x1.72   t =  +4.2      re-stretch recovery DELAY
%     t0       x1.42   t =  +3.7
%     force SHAPE (peak1_y, vall_y, peak2, A over steady)   x1.01-1.06, flat
%   and the ktr LADDER is itself an ATP effect: -0.2 / -7.4 / -4.9 %/segment at
%   8 mM against -10.0 / -11.1 / -16.2 at 2 mM, steeper in every preparation.
%
% The 2 mM data enter as RATIOS, so the 8 mM fit has to get the BASELINE of
% every ATP-carrying metric right before any ATP mechanism can be fitted. And
% the two strongest carriers sat inside deadbands built from a MODERATED
% spread about twice their measured one: rsK80 0.179 vs 0.090 raw, ktr 0.181
% vs 0.113 -- which is why `ktr` never bit.
%
% CHANGES
%   ATP carriers, scored on their OWN between-preparation SD (raw, no
%   moderation), fl = (1/3)*nEl*s^2, w = lambda/s^2:
%     ktr80[1:5]        replaces `ktr` (threshold-biased free-amplitude fit)
%     ktr80_shape[1:5]  replaces ktr_shape (was /ktr(1), amplified seg-1 noise)
%     rsK80[1:4]        re-floored from 0.179 to its raw spread
%     rsT20[1:4]        NEW -- was not scored at all
%     ktrP_k80          replaces ktrP_k63: the ktr PROTOCOL, pooled over the
%                       three days by pooledKtrTarget, its own raw spread
%   ATP-flat force SHAPE terms become GUARDS (floored at 1.2x their current raw
%   cost): they cost nothing unless they get >20 % worse. Their error is the
%   same at both ATP levels, so it does not bias an ATP ratio.
%     A, peak1_y, vall_y, peak2, vall2_dy, restretchSlopeStart, A_over_steady
%   Kept as they are: steady (the force scale -- the one force ATP effect) and
%   its level guard, t0_crossing, vall_t, peak1_dSL, FV, PS_*, bounds.
%
% Also tested here: kmsr = 5 vs 0. kmsr did nothing for the slope at 1-10,
% cost +0.023 at the seed, and refit6 evaluations ran ~2x slower than refit5's.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/BuildRefit7Seed.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~exist('SEEDTAG','var'); SEEDTAG = 'refit7_seed'; end
LAMBDA = 0.0106294;
GUARD_X = 1.2;

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_r7_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 6); end

p0 = getParams(loadParams('refit6_opt'), [], true, false);
p0.RunForceVelocity=true; p0.RunSlack=true; p0.RunStairs=false; p0.RunKtr=true;
p0.RunForceVelocityTime=false; p0.RunSlackPassive=true; p0.EvalFeatures=true;
p0.BreakOnODEUnstable=false; p0.PlotEachSeparately=0; p0.PlotFeatureFitting=0;
p0.RunSlackSegments='AllPar'; p0.MaxRunTime=150; p0.OptimizeOn='Feats';
p0.ghostSave=''; p0.ghostLoad=''; p0.SaveBest=0; p0.mods={}; p0.g=ones(1,0);

%% ---- raw between-prep spreads, straight from the three curated files ---
DAYS = {'03_27','04_03','04_10'};
X = struct('ktr80',[],'ktr80_shape',[],'rsK80',[],'rsT20',[]);
for d = 1:3
    S = load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d})));
    f = S.features_data;
    X.ktr80(d,:)       = f.ktr80(:)';
    X.ktr80_shape(d,:) = f.ktr80(:)' / mean(f.ktr80);
    X.rsK80(d,:)       = f.rsK80(1:4);
    X.rsT20(d,:)       = f.rsT20(1:4);
end
cv = @(M) mean(std(M,0,1,'omitnan') ./ abs(mean(M,1,'omitnan')), 'omitnan');
SR = struct('ktr80',cv(X.ktr80), 'ktr80_shape',cv(X.ktr80_shape), ...
            'rsK80',cv(X.rsK80), 'rsT20',cv(X.rsT20));
PK = pooledKtrTarget();
jk = find(strcmp(PK.name,'ktrP_k80'));
SR.ktrP_k80 = PK.s(jk);
fprintf('raw between-prep s:  ktr80 %.3f | ktr80_shape %.3f | rsK80 %.3f | rsT20 %.3f | ktrP_k80 %.3f\n', ...
        SR.ktr80, SR.ktr80_shape, SR.rsK80, SR.rsT20, SR.ktrP_k80);
fprintf('ktr protocol, per prep k80: %s -> pooled %.2f\n', mat2str(round(PK.perPrep(:,jk)',2)), PK.val(jk));

term = @(spec, nEl, s) sprintf('%s|fl:%.6g|%.6g', spec, (1/3)*nEl*s^2, LAMBDA/s^2);

%% ---- the new fn ----------------------------------------------------------
fn0  = p0.fn;
DROP = {'ktr','ktr_shape','ktrP_k63','rsK80'};
GUARD = {'A','peak1_y','vall_y','peak2','vall2_dy','restretchSlopeStart','A_over_steady'};
base = @(e) regexprep(regexprep(e,'\|.*$',''), '\[.*$', '');
keep = true(1,numel(fn0));
for i = 1:numel(fn0)
    if contains(fn0{i}, ','); continue; end          % grouped boundary entries
    if ismember(base(fn0{i}), DROP); keep(i) = false; end
end
fn = fn0(keep);
fn{end+1} = term('ktr80',          5, SR.ktr80);
fn{end+1} = term('ktr80_shape',    5, SR.ktr80_shape);
fn{end+1} = term('rsK80[1:4]',     4, SR.rsK80);
fn{end+1} = term('rsT20[1:4]',     4, SR.rsT20);
fn{end+1} = term('ktrP_k80',       1, SR.ktrP_k80);

% guards: pass 1 measures each term's RAW cost at the seed (floor stripped)
iG = find(cellfun(@(e) ~contains(e,',') && ismember(base(e), GUARD), fn));
fnRaw = fn;
for i = iG
    tk = split(fn{i},'|'); tk = tk(~startsWith(tk,'fl:'));
    fnRaw{i} = strjoin(tk,'|');
end
pr = p0; pr.fn = fnRaw;
[~,~,ctR,nmR] = evaluateBakersExp(1, pr, true, false);
fprintf('\n%-22s %10s %10s %10s\n','guard term','raw cost','new floor','weight');
for i = iG
    tk = split(fn{i},'|');
    w  = str2double(tk{end});
    raw = ctR(i) / w;                                % unfloored, unweighted
    tk = tk(~startsWith(tk,'fl:'));
    tk = [tk(1:end-1); {sprintf('fl:%.6g', GUARD_X*raw)}; tk(end)];
    fn{i} = strjoin(tk,'|');
    fprintf('%-22s %10.5f %10.5f %10.4g\n', base(fn{i}), raw, GUARD_X*raw, w);
end
p0.fn = fn;

% Coverage: evalFeatureCost charges MISSING_FEATURE_COST = 100 per absent
% feature, so pass 1 already is the coverage check -- no second (serial,
% slow) evaluation needed.
bad = find(ctR >= 100);
if ~isempty(bad)
    error('BuildRefit7Seed:coverage','unscorable terms: %s', strjoin(nmR(bad),', '));
end

%% ---- kmsr 5 vs 0: cost AND wall time ------------------------------------
VAR = {'kmsr as fitted', p0.kmsr; 'kmsr = 0', 0};
C = nan(2,1); TT = nan(2,1); CT = cell(2,1); NM = {};
for v = 1:2
    pv = p0; pv.kmsr = VAR{v,2}; pv = getParams(pv, pv.g, true, false);
    t0 = tic;
    [C(v),~,CT{v},NM] = evaluateBakersExp(1, pv, true, false);
    TT(v) = toc(t0);
    fprintf('%-16s cost %.4f   wall %.0f s\n', VAR{v,1}, C(v), TT(v));
end
% prefer kmsr = 0 unless it costs more than the objective's ~1 % reproducibility
useZero = C(2) <= C(1) * 1.01;
if useZero
    p0.kmsr = 0; p0 = getParams(p0, p0.g, true, false);
    p0.leanPool = setdiff(p0.leanPool, {'kmsr','sigma1'}, 'stable');
    ct = CT{2}; c_new = C(2);
    fprintf('-> kmsr = 0 (sigma1 dead again, both dropped from the pool)\n');
else
    ct = CT{1}; c_new = C(1);
    fprintf('-> kmsr kept at %g\n', p0.kmsr);
end

%% ---- composition ---------------------------------------------------------
ATPC = {'ktr80','ktr80_shape','rsK80','rsT20','ktrP_k80','steady','t0_crossing'};
lab = cellfun(@(s) base(s), NM, 'UniformOutput', false);
fprintf('\n=== %s = %.4f ===\n', SEEDTAG, c_new);
[~,o] = sort(ct,'descend');
for i = o(:)'
    if ct(i) <= 1e-9; continue; end
    tag = ''; if ismember(lab{i}, ATPC); tag = '   ATP'; end
    if ismember(lab{i}, GUARD); tag = '   guard'; end
    fprintf('    %-22s %8.4f  (%4.1f %%)%s\n', lab{i}, ct(i), 100*ct(i)/c_new, tag);
end
fprintf('  ATP-carrying share: %.1f %%\n', 100*sum(ct(ismember(lab,ATPC)))/c_new);
fprintf('  zero-cost: %s\n', strjoin(unique(lab(ct<=1e-9)),', '));

%% ---- write ---------------------------------------------------------------
out = fullfile(root,'params',[SEEDTAG '.m']);
iWrite17g(out, p0, sprintf(['refit6_opt params; objective rebuilt around the ATP-carrying ' ...
    'metrics (ktr80, ktr80_shape, rsK80, rsT20, pooled ktrP_k80 on raw spreads; force ' ...
    'shape -> guards x%.1f); kmsr %g. cost %.6f'], GUARD_X, p0.kmsr, c_new));
q = getParams(); s2 = loadParams(out);
f2 = fieldnames(s2); for i=1:numel(f2); q.(f2{i}) = s2.(f2{i}); end
fprintf('written %s | pool %d | leanPool round-trips %d | fn %d terms\n', out, ...
        numel(p0.leanPool), numel(q.leanPool), numel(q.fn));
save(fullfile(resdir,'refit7_seed_build.mat'),'SR','PK','C','TT','ct','NM','c_new','fn');
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
