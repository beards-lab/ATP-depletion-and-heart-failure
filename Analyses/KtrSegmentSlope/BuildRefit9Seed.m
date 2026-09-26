% BuildRefit9Seed.m -- refit8's ATP-carrying objective, rebuilt under the fixed
% numerics (Analyses/StrainWindowNoise/conclusions.md).
%
% WHAT CHANGED UNDERNEATH. refit7/refit8's noise (19.5 % of the objective) was
% not the strain window: features were read off ode15s STEP POINTS at a
% hard-coded RelTol 1e-2, so a 1e-6 parameter change moved the step sequence
% and with it every crossing-time feature. Tighter tolerances used to stall on
% a corner in the tension-only titin Maxwell rectifier at isometric rest. The
% fix, all opt-in parameters carried in the snapshot:
%     SmoothRectifierEps = 0.01   odeRelTol = 1e-3   odeAbsTol = 1e-6
%     odeOutputDt = 5e-4          ClampStrainDepExtrapolation = 1
% gives 0.3 % noise on the full objective at refit8_opt.
%
% CONSEQUENCES HANDLED HERE
%   * old costs are not comparable (refit8_opt ~0.14 -> 0.356), and part of
%     refit7/8's progress was sampling artefact -> re-score both candidate
%     seeds under the new numerics and seed from the better one on the
%     ATP-carrying block;
%   * refit8's floors added s_num^2 measured under the OLD noise (1-9 %) ->
%     back to the plain chi-2 floor nEl*s^2/3 (s = own raw between-prep spread);
%   * guards (ATP-flat force shape, vall_t, peak1_dSL) re-derived at 1.2x of
%     the chosen seed's raw cost under the new numerics.
% Objective design is unchanged from refit8 (fjezek, 2026-09-23): ATP-carrying
% metrics scored, force shape as guards, ktr protocol pooled on its own raw
% spread, staircase visual only.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/BuildRefit9Seed.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~exist('SEEDTAG','var'); SEEDTAG = 'refit9_seed'; end
LAMBDA = 0.0106294; GUARD_X = 1.2;
NUM = struct('SmoothRectifierEps',0.01, 'odeRelTol',1e-3, 'odeAbsTol',1e-6, ...
             'odeOutputDt',5e-4, 'ClampStrainDepExtrapolation',1);
CAND = {'refit8_seed','refit8_opt'};

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_r9_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 6); end

%% ---- the objective, floors from the between-prep spread only -----------
DAYS = {'03_27','04_03','04_10'};
X = struct('ktrfit',[],'ktrfit_shape',[],'rsK80',[],'rsT20',[]);
for d = 1:3
    S = load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d})));
    f = S.features_data;
    X.ktrfit(d,:)       = f.ktrfit(:)';
    X.ktrfit_shape(d,:) = f.ktrfit(:)' / mean(f.ktrfit);
    X.rsK80(d,:)        = f.rsK80(1:4);
    X.rsT20(d,:)        = f.rsT20(1:4);
end
cv = @(M) mean(std(M,0,1,'omitnan') ./ abs(mean(M,1,'omitnan')), 'omitnan');
SB = struct('ktrfit',cv(X.ktrfit), 'ktrfit_shape',cv(X.ktrfit_shape), ...
            'rsK80',cv(X.rsK80), 'rsT20',cv(X.rsT20));
PK = pooledKtrTarget(); SB.ktrP_k80 = PK.s(strcmp(PK.name,'ktrP_k80'));
term = @(spec, nEl, s) sprintf('%s|fl:%.6g|%.6g', spec, nEl*s^2/3, LAMBDA/s^2);

S0 = loadParams('refit8_seed'); fn0 = S0.fn;
ATPT  = {'ktrfit','ktrfit_shape','rsK80','rsT20','ktrP_k80'};
GUARD = {'A','peak1_y','vall_y','peak2','vall2_dy','restretchSlopeStart','A_over_steady','vall_t','peak1_dSL'};
ATPC  = [ATPT, {'steady','t0_crossing'}];
base  = @(e) regexprep(regexprep(e,'\|.*$',''), '\[.*$', '');
keep = true(1,numel(fn0));
for i = 1:numel(fn0)
    if contains(fn0{i}, ','); continue; end
    if ismember(base(fn0{i}), ATPT); keep(i) = false; end
end
fn = fn0(keep);
fn{end+1} = term('ktrfit',       5, SB.ktrfit);
fn{end+1} = term('ktrfit_shape', 5, SB.ktrfit_shape);
fn{end+1} = term('rsK80[1:4]',   4, SB.rsK80);
fn{end+1} = term('rsT20[1:4]',   4, SB.rsT20);
fn{end+1} = term('ktrP_k80',     1, SB.ktrP_k80);
iG = find(cellfun(@(e) ~contains(e,',') && ismember(base(e), GUARD), fn));
fnRaw = fn;
for i = iG
    tk = split(fn{i},'|'); tk = tk(~startsWith(tk,'fl:')); fnRaw{i} = strjoin(tk,'|');
end
fprintf('%-14s %9s\n','ATP term','s (raw)');
for q = ATPT; fprintf('%-14s %9.4f\n', q{1}, SB.(q{1})); end

%% ---- re-score both candidates under the new numerics --------------------
R = struct();
for c = 1:numel(CAND)
    p = getParams(loadParams(CAND{c}), [], true, false);
    for q = fieldnames(NUM)'; p.(q{1}) = NUM.(q{1}); end
    assert(p.dS == 0.002 && p.MaxStrainArraySize == 60 && mod(p.MaxStrainArraySize,2) == 0 ...
           && p.WindowsOverflowStepCount == -1, 'grid settings differ from the verified ones');
    p.RunForceVelocity=true; p.RunSlack=true; p.RunStairs=false; p.RunKtr=true;
    p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
    p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
    p.RunSlackSegments='AllPar'; p.MaxRunTime=300; p.OptimizeOn='Feats';
    p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0; p.mods={}; p.g=ones(1,0);
    p.fn = fnRaw; p = getParams(p, p.g, true, false);
    t0 = tic; [c0,~,ct,NM] = evaluateBakersExp(1, p, true, false); wt = toc(t0);
    if any(ct >= 100); error('BuildRefit9Seed:coverage','unscorable term for %s', CAND{c}); end
    lab = cellfun(@(s) base(s), NM, 'UniformOutput', false);
    R(c).p = p; R(c).ct = ct; R(c).lab = lab; R(c).wall = wt;
    R(c).atp = sum(ct(ismember(lab, ATPC)));
    R(c).guardRaw = ct(iG) ./ cellfun(@(e) str2double(extractAfter(e, find(e=='|',1,'last'))), fn(iG));
    fprintf('\n%-12s ATP block %.4f | total with raw guards %.4f | %.0f s\n', CAND{c}, R(c).atp, c0, wt);
    for j = find(ismember(lab, ATPC))
        fprintf('    %-14s %8.4f\n', lab{j}, ct(j));
    end
end
[~, ib] = min([R.atp]);
fprintf('\n-> seeding from %s (ATP block %.4f vs %.4f)\n', CAND{ib}, R(ib).atp, R(3-ib).atp);
p = R(ib).p;

%% ---- guards at 1.2x of the chosen point's raw cost ----------------------
for k = 1:numel(iG)
    i = iG(k);
    tk = split(fn{i},'|'); tk = tk(~startsWith(tk,'fl:'));
    tk = [tk(1:end-1); {sprintf('fl:%.6g', GUARD_X*R(ib).guardRaw(k))}; tk(end)];
    fn{i} = strjoin(tk,'|');
end
p.fn = fn; p = getParams(p, p.g, true, false);

%% ---- the seed, twice: unperturbed and at a 1e-6 perturbation -------------
[c1,~,ct1,NM] = evaluateBakersExp(1, p, true, false);
pp = p;
for i = 1:numel(pp.leanPool)
    nm = pp.leanPool{i};
    if isempty(strfind(nm,'__')) && isfield(pp,nm) && isscalar(pp.(nm)) && pp.(nm) ~= 0 %#ok<STREMP>
        pp.(nm) = pp.(nm)*(1 + 1e-6*randn);
    end
end
pp = getParams(pp, pp.g, true, false);
c2 = evaluateBakersExp(1, pp, true, false);
lab = cellfun(@(s) base(s), NM, 'UniformOutput', false);
fprintf('\n=== %s = %.4f  (1e-6 neighbour %.4f, %+.2f %%) ===\n', SEEDTAG, c1, c2, 100*(c2-c1)/c1);
[~,o] = sort(ct1,'descend');
for i = o(:)'
    if ct1(i) <= 1e-9; continue; end
    tag = ''; if ismember(lab{i}, ATPC); tag = '   ATP'; elseif ismember(lab{i}, GUARD); tag = '   guard'; end
    fprintf('    %-22s %8.4f  (%4.1f %%)%s\n', lab{i}, ct1(i), 100*ct1(i)/c1, tag);
end
fprintf('  ATP-carrying share: %.1f %%\n', 100*sum(ct1(ismember(lab,ATPC)))/c1);

%% ---- write ---------------------------------------------------------------
out = fullfile(root,'params',[SEEDTAG '.m']);
writeParamsToMFile(out, p, [], sprintf(['%s under the fixed numerics (SiLU rectifier, RelTol 1e-3, dense ' ...
    '0.5 ms output, extrapolation clamp); ATP-carrying objective, floors nEl*s^2/3, guards 1.2x. cost %.6f'], ...
    CAND{ib}, c1));
q = getParams(); s2 = loadParams(out);
f2 = fieldnames(s2); for i=1:numel(f2); q.(f2{i}) = s2.(f2{i}); end
for fq = fieldnames(NUM)'
    assert(isfield(q, fq{1}) && q.(fq{1}) == NUM.(fq{1}), 'numerics setting %s did not round-trip', fq{1});
end
fprintf('written %s | seeded from %s | pool %d | fn %d | numerics round-trip OK\n', out, CAND{ib}, ...
        numel(p.leanPool), numel(q.fn));
seed = struct('snap',SEEDTAG,'cost',c1,'featCost',ct1(:)','featNames',{NM(:)'}, ...
              'when',datestr(now,'yyyy-mm-dd HH:MM'));
save(fullfile(resdir,'refit9_seed_breakdown.mat'),'seed');
save(fullfile(resdir,'refit9_seed_build.mat'),'SB','R','c1','c2','fn','CAND','ib');
fprintf('\nDONE\n');
