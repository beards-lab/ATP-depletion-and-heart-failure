% RunStructuralProbe.m -- build the refit6 objective (ktr_shape at 1 x s_shape)
% and then ask, on THAT objective, whether any structural switch buys the slack
% ladder's slope more cheaply than the strain-dependence levers did.
%
% WHY THIS BEFORE A LONG RUN. refit5 ran 6.3 h with three cheap slope levers
% compulsory every round and closed 5 % of the slope gap. The levers exist
% (-1.4 to -3.8 per +15 %) so it is a pricing problem -- but raising the price
% only helps if SOME reachable configuration can pay it. Two of the mechanisms
% that should produce a length-dependent RATE are switched off by things the
% optimiser can never touch:
%
%   kmsr = 0            ->  RSR2PT = kmsr*exp(F_SR/sigma1)*P_SR is identically
%                           zero, so the direct force-gated SRX exit is OFF and
%                           sigma1 is dead code. The only way out of SRX is
%                           ksr2srd (force-INdependent) then kmsrd.
%   UsePassiveForSR = 0 ->  the SRX gate reads total force, not passive force.
%                           Passive force is what actually carries the length
%                           signal at a shortened SL; this is the canonical
%                           thick-filament mechanosensing route.
%
% Neither is in leanPool and kmsr is excluded anyway by the zero-base filter, so
% no campaign of any length can find them.
%
% WHY IT MATTERS BEYOND THIS FEATURE. The 2 mM signature is ktr x0.54 against
% force x1.18 (Docs, [[lowatp-analysis]]). ktr is the ATP-sensitive observable.
% A model whose ktr has the wrong dependence on LENGTH is a model whose ktr
% mechanism is wrong, and that dependence is what the low-ATP fit has to lean
% on. Chasing the slope is not cosmetic.
%
% Reported per variant: the full objective, the ktr ladder slope, and the
% RE-STRETCH subtotal -- peak1_y, peak1_dSL, vall_y, vall_t, vall2_dy, peak2,
% restretchSlopeStart, rsK80 -- because the earlier analysis concluded the
% low-ATP fit is spoiled by exactly those high-ATP re-stretch errors, so a
% variant that buys slope by selling re-stretch is worse than useless.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunStructuralProbe.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~isfolder(resdir); mkdir(resdir); end
LAMBDA = 0.0106294;

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_sp_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 6); end

%% ---- the refit6 objective: ktr_shape at 1 x s_shape ---------------------
p0 = getParams(loadParams('refit5_opt'), [], true, false);
p0.RunForceVelocity=true; p0.RunSlack=true; p0.RunStairs=false; p0.RunKtr=true;
p0.RunForceVelocityTime=false; p0.RunSlackPassive=true; p0.EvalFeatures=true;
p0.BreakOnODEUnstable=false; p0.PlotEachSeparately=0; p0.PlotFeatureFitting=0;
p0.RunSlackSegments='AllPar'; p0.MaxRunTime=300; p0.OptimizeOn='Feats';
p0.ghostSave=''; p0.ghostLoad=''; p0.SaveBest=0; p0.mods={}; p0.g=ones(1,0);

DAYS = {'03_27','04_03','04_10'};
Xk = zeros(3,5);
for d = 1:3
    S = load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d})));
    Xk(d,:) = S.features_data.ktr(:)'/S.features_data.ktr(1);
end
Xk = Xk(:,2:5);
s_shape = mean(std(Xk,0,1,'omitnan')./abs(mean(Xk,1,'omitnan')),'omitnan');
s_new = s_shape;                       % 1 x, was 2 x
fl_new = (1/3)*4*s_new^2; w_new = LAMBDA/s_new^2;

fn = p0.fn; found = false;
for i = 1:numel(fn)
    if startsWith(fn{i},'ktr_shape')
        fn{i} = sprintf('ktr_shape[2:5]|fl:%.6g|%.6g', fl_new, w_new); found = true;
    end
end
if ~found; error('RunStructuralProbe:noterm','ktr_shape not in fn'); end
p0.fn = fn;
fprintf('ktr_shape re-priced: s %.4f -> %.4f | fl %.5f -> %.5f | w %.4f -> %.4f\n', ...
        2*s_shape, s_new, (1/3)*4*(2*s_shape)^2, fl_new, LAMBDA/(2*s_shape)^2, w_new);

RESTRETCH = {'peak1_y','peak1_dSL','vall_y','vall_t','vall2_dy','peak2', ...
             'restretchSlopeStart','rsK80'};

%% ---- variants ----------------------------------------------------------
% {label, {field,value; field,value; ...}}
V = { 'BASELINE refit5_opt',        {}
      'UsePassiveForSR=1',          {'UsePassiveForSR',1}
      'kmsr=10',                    {'kmsr',10}
      'kmsr=50',                    {'kmsr',50}
      'kmsr=10 + PassiveForSR',     {'kmsr',10;'UsePassiveForSR',1}
      'kmsr=50 + PassiveForSR',     {'kmsr',50;'UsePassiveForSR',1}
      'useHalfActiveForSR=1',       {'useHalfActiveForSR',1}
      'UseDirectSRXTransition=1',   {'UseDirectSRXTransition',1}
      'P_bound_max=0.8',            {'P_bound_max',0.8}
      'P_bound_max=0.6',            {'P_bound_max',0.6}
      'UseOverlapFactor=0',         {'UseOverlapFactor',0} };
nV = size(V,1);
C = nan(nV,1); SLP = nan(nV,1); RS = nan(nV,1); KS = nan(nV,1); FV = nan(nV,1);
KTR = nan(nV,5); ERR = cell(nV,1);

parfor i = 1:nV
    p = p0; M = V{i,2};
    for r = 1:size(M,1); p.(M{r,1}) = M{r,2}; end
    p = getParams(p, p.g, true, false);
    p.RunSlackSegments = 'All';        % serial inside the parfor
    try
        [c,~,ct,nm] = evaluateBakersExp(1, p, true, false);
        C(i) = c;
        RS(i) = sum(ct(ismember(nm, RESTRETCH)));
        j = find(strcmp(nm,'ktr_shape'),1);   if ~isempty(j); KS(i) = ct(j); end
        j = find(strcmp(nm,'FV_normpowerAvg'),1); if ~isempty(j); FV(i) = ct(j); end
        [~,~,fm] = runSlackExperiment(p);
        KTR(i,:) = fm.ktr(:)';
        bb = polyfit(1:5, KTR(i,:), 1); SLP(i) = bb(1);
    catch ME
        ERR{i} = ME.message;
    end
end

fprintf('\n=== structural probe on the refit6 objective (data slope -1.848) ===\n');
fprintf('%-26s %9s %9s %9s %9s %9s   %s\n', ...
        'variant','TOTAL','restretch','ktr_shape','FV','slope','ktr seg1..seg5');
for i = 1:nV
    if ~isempty(ERR{i}); fprintf('%-26s  FAILED: %s\n', V{i,1}, ERR{i}); continue; end
    fprintf('%-26s %9.4f %9.4f %9.4f %9.4f %+9.3f   %s\n', V{i,1}, C(i), RS(i), KS(i), FV(i), SLP(i), ...
            sprintf('%6.1f', KTR(i,:)));
end
fprintf('%-26s %9s %9s %9s %9s %+9.3f   %s\n','POOLED DATA','-','-','-','-',-1.848, ...
        sprintf('%6.1f',[51.61 49.50 47.87 46.32 43.96]));

%% ---- verdict ------------------------------------------------------------
base = 1;
fprintf('\n=== relative to baseline ===\n');
fprintf('%-26s %10s %10s %10s\n','variant','d TOTAL','d restr','d slope');
for i = 2:nV
    if ~isempty(ERR{i}); continue; end
    fprintf('%-26s %+10.4f %+10.4f %+10.3f\n', V{i,1}, C(i)-C(base), RS(i)-RS(base), SLP(i)-SLP(base));
end
good = find(SLP < SLP(base)-0.15 & C < C(base)*1.5 & RS < RS(base)*1.2);
if isempty(good)
    fprintf('\nNo structural switch buys slope at an acceptable price.\n');
else
    fprintf('\nCandidates (slope down, total within 1.5x, re-stretch within 1.2x): %s\n', ...
            strjoin(V(good,1)', ' | '));
end

save(fullfile(resdir,'structural_probe.mat'),'V','C','SLP','RS','KS','FV','KTR', ...
     's_shape','fl_new','w_new');
fprintf('\nDONE\n');
