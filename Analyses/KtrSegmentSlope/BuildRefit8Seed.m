% BuildRefit8Seed.m -- refit7's objective made robust to the model's own
% numerical noise, seeded from refit7's round-1 point (the one that reproduced).
%
% WHAT WENT WRONG IN refit7 (Analyses/KtrSegmentSlope/RunNumericalNoise.m):
%   round 2 accepted 0.3355 and its snapshot re-scored 0.3707. Nine evaluations
%   of that point at random 1e-6 perturbations span 0.3355-0.4160 (SD 7.2 %),
%   and the accepted value was the MINIMUM of the nine -- a winner's curse.
%   All of the scatter sits in the three rate terms: ktr80_shape 23 %, rsK80
%   27 %, rsT20 19 % of their own cost; vall_t / steady / t0 move 0.1-0.3 %.
%   Cause (memory objective-noise-floor): discrete decisions inside the model
%   -- event-driven strain-grid shifts, onset/window choices -- make the
%   recovery traces jump ~1 %; with the model still 5-10 % off the data the
%   cost moves as 2*residual*jitter.
%
% THREE FIXES, each measured:
%   1. RATES FROM THE WHOLE BAND. ktrfit / ktrfit_shape / rsKfit replace
%      ktr80 / ktr80_shape / rsK80: a weighted log-linear regression over the
%      20-80 % band (bandRateFit) instead of two crossing times. Numerical
%      jitter 1.09 -> 0.72 %, 1.03 -> 0.67 %, 1.43 -> 1.00 %. They carry the
%      ATP effect as well or better (ktrfit x0.57, t -9.7; rsKfit x0.39,
%      t -15.8). rsT20 stays on crossings: its regression version was NOISIER
%      (5.1 vs 3.2 %).
%   2. NOISE IN THE DEADBAND. fl = nEl * (s_between^2/3 + s_num^2): a term
%      stops paying inside the model's own reproducibility as well as inside
%      the uncertainty of the 3-heart mean. Weights stay lambda/s_between^2.
%   3. NOISE-AWARE ACCEPTANCE in the campaign (optimizeFeatures CONFIRM_N): a
%      candidate is banked only if the MEAN of it and two 1e-6-perturbed
%      copies beats the incumbent's mean. This is what actually stops the
%      winner's curse -- (1) and (2) only shrink the noise.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/BuildRefit8Seed.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~exist('SEEDTAG','var'); SEEDTAG = 'refit8_seed'; end
LAMBDA = 0.0106294;
GUARD_X = 1.2;

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_r8_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 6); end

p0 = getParams(loadParams(fullfile(root,'params','refit7_iter','refit7_iter_001.m')), [], true, false);
p0.RunForceVelocity=true; p0.RunSlack=true; p0.RunStairs=false; p0.RunKtr=true;
p0.RunForceVelocityTime=false; p0.RunSlackPassive=true; p0.EvalFeatures=true;
p0.BreakOnODEUnstable=false; p0.PlotEachSeparately=0; p0.PlotFeatureFitting=0;
p0.RunSlackSegments='AllPar'; p0.MaxRunTime=150; p0.OptimizeOn='Feats';
p0.ghostSave=''; p0.ghostLoad=''; p0.SaveBest=0; p0.mods={}; p0.g=ones(1,0);

%% ---- spreads: between-prep (raw, from the curated files) and numerical --
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
N = load(fullfile(resdir,'numerical_noise.mat'),'SNUM'); SN = N.SNUM;
fprintf('%-14s %10s %10s\n','term','s_between','s_num');
for q = {'ktrfit','ktrfit_shape','rsK80','rsT20','ktrP_k80'}
    fprintf('%-14s %10.4f %10.4f\n', q{1}, SB.(q{1}), SN.(q{1}));
end
term = @(spec, nEl, sb, sn) sprintf('%s|fl:%.6g|%.6g', spec, nEl*(sb^2/3 + sn^2), LAMBDA/sb^2);

%% ---- the objective -------------------------------------------------------
fn0  = p0.fn;
DROP = {'ktr80','ktr80_shape','rsK80','rsT20','ktrP_k80'};
GUARD = {'A','peak1_y','vall_y','peak2','vall2_dy','restretchSlopeStart','A_over_steady','vall_t','peak1_dSL'};
% vall_t and peak1_dSL too: not ATP-carrying (paired t -2.3 / +2.6, n = 3), and
% vall_t has topped every objective since refit4 without moving (the 2-state
% wall). At full weight they were 55 % of the objective and would let the
% search trade ATP metrics for them.
base = @(e) regexprep(regexprep(e,'\|.*$',''), '\[.*$', '');
keep = true(1,numel(fn0));
for i = 1:numel(fn0)
    if contains(fn0{i}, ','); continue; end
    if ismember(base(fn0{i}), DROP); keep(i) = false; end
end
fn = fn0(keep);
fn{end+1} = term('ktrfit',       5, SB.ktrfit,       SN.ktrfit);
fn{end+1} = term('ktrfit_shape', 5, SB.ktrfit_shape, SN.ktrfit_shape);
fn{end+1} = term('rsK80[1:4]',   4, SB.rsK80,        SN.rsK80);
fn{end+1} = term('rsT20[1:4]',   4, SB.rsT20,        SN.rsT20);
fn{end+1} = term('ktrP_k80',     1, SB.ktrP_k80,     SN.ktrP_k80);

% guards re-derived at THIS seed: floor at 1.2x each term's raw cost
iG = find(cellfun(@(e) ~contains(e,',') && ismember(base(e), GUARD), fn));
fnRaw = fn;
for i = iG
    tk = split(fn{i},'|'); tk = tk(~startsWith(tk,'fl:')); fnRaw{i} = strjoin(tk,'|');
end
pr = p0; pr.fn = fnRaw;
[~,~,ctR,nmR] = evaluateBakersExp(1, pr, true, false);
bad = find(ctR >= 100);
if ~isempty(bad); error('BuildRefit8Seed:coverage','unscorable: %s', strjoin(nmR(bad),', ')); end
for i = iG
    tk = split(fn{i},'|'); w = str2double(tk{end}); raw = ctR(i)/w;
    tk = tk(~startsWith(tk,'fl:'));
    tk = [tk(1:end-1); {sprintf('fl:%.6g', GUARD_X*raw)}; tk(end)];
    fn{i} = strjoin(tk,'|');
end
p0.fn = fn;

%% ---- score, three times at 1e-6 perturbations: what the campaign will see
ATPC = {'ktrfit','ktrfit_shape','rsK80','rsT20','ktrP_k80','steady','t0_crossing'};
rng(11);
C = nan(3,1); CT = [];
for r = 1:3
    pv = p0;
    if r > 1
        for i = 1:numel(pv.leanPool)
            nm = pv.leanPool{i};
            if isempty(strfind(nm,'__')) && isfield(pv,nm) && isscalar(pv.(nm)) && pv.(nm) ~= 0 %#ok<STREMP>
                pv.(nm) = pv.(nm)*(1 + 1e-6*randn);
            end
        end
        pv = getParams(pv, pv.g, true, false);
    end
    [C(r),~,ct,NM] = evaluateBakersExp(1, pv, true, false);
    CT(r,:) = ct(:)'; %#ok<AGROW>
end
lab = cellfun(@(s) base(s), NM, 'UniformOutput', false);
ct = mean(CT,1); c_new = mean(C);
fprintf('\n=== %s: %s -> mean %.4f (SD %.4f, %.1f %%) ===\n', SEEDTAG, mat2str(round(C',4)), ...
        c_new, std(C), 100*std(C)/c_new);
[~,o] = sort(ct,'descend');
for i = o(:)'
    if ct(i) <= 1e-9; continue; end
    tag = ''; if ismember(lab{i}, ATPC); tag = '   ATP'; end
    fprintf('    %-22s %8.4f  (%4.1f %%)  sd %.4f%s\n', lab{i}, ct(i), 100*ct(i)/c_new, std(CT(:,i)), tag);
end
fprintf('  ATP-carrying share: %.1f %%\n', 100*sum(ct(ismember(lab,ATPC)))/c_new);

%% ---- write ---------------------------------------------------------------
out = fullfile(root,'params',[SEEDTAG '.m']);
writeParamsToMFile(out, p0, [], sprintf(['refit7 round-1 point; rates from the whole 20-80 %% band ' ...
    '(ktrfit, ktrfit_shape; rsK80 kept on the audited crossings), floors nEl*(s^2/3 + s_num^2). mean cost %.6f'], c_new));
q = getParams(); s2 = loadParams(out);
f2 = fieldnames(s2); for i=1:numel(f2); q.(f2{i}) = s2.(f2{i}); end
fprintf('written %s | pool %d | fn %d terms (round-trips %d)\n', out, numel(p0.leanPool), numel(fn), numel(q.fn));
save(fullfile(resdir,'refit8_seed_build.mat'),'SB','SN','C','CT','NM','c_new','fn');
seed = struct('snap',SEEDTAG,'cost',c_new,'featCost',ct,'featNames',{NM(:)'}, ...
              'when',datestr(now,'yyyy-mm-dd HH:MM'));
save(fullfile(resdir,'refit8_seed_breakdown.mat'),'seed');
fprintf('\nDONE\n');
