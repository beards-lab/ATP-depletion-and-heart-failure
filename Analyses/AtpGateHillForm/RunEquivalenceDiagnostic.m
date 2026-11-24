% RunEquivalenceDiagnostic.m -- why did the 2 mM cost move by 0.4 % when the gate
% factors are identical to 3e-16?
%
% Two candidate explanations, and this script separates them:
%   (H1) the conversion is not quite exact -- RunHillGateConversion drops gates with
%        |p| < PTHRESH (= 1e-3) as optimiser noise, which perturbs those rates by a
%        few ppm. Test: convert EVERY gate (PTHRESH = 0) and re-compare.
%   (H2) the slack pipeline is not bit-reproducible, so a few-ppm perturbation (or
%        none at all) shows up as a visible feature change. Test: run the SAME
%        parameter set twice and compare it with itself.
%
% H2 is the one that matters beyond this analysis: it puts a NOISE FLOOR under every
% cost difference this project quotes, including the power-gate optimisation itself.
%
% Run: matlab -batch "run('Analyses/AtpGateHillForm/RunEquivalenceDiagnostic.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end
if isempty(gcp('nocreate')); parpool('local',6); end

BASE = 'lowATP_powergate2_opt'; ATP2 = 2;
GATED = {'k2','kmsrd','kah','kamh','ka','kd','k1','k_1', ...
         'ksr0','ksr2srd','ksrd2sr','tau_reg','A0','eta_M', ...
         'kstiff2','kSE','kSE_M','mu_neg','dr','dr2','P_bound_max'};

b0 = getParams(loadParams(BASE), [], true, false); Aref = b0.MgATP_ref;

%% ---- 0. exact size of the perturbation the PTHRESH cut introduces -------------
fprintf('\n=== per-rate factor at %g mM: power law vs the converted gate ===\n', ATP2);
fprintf('%-14s %12s %14s %14s %12s\n','rate','p','power factor','hill(kept)','rel. diff');
worstKept = 0; worstAll = 0;
hillAll = struct();
for i = 1:numel(GATED)
    p = b0.(['pAtp_' GATED{i}]);
    if p == 0, continue; end
    fp = (ATP2/Aref)^p;
    [Ka, na] = solveHillK(fp, ATP2, Aref);           % convert EVERY gate, no threshold
    hillAll.(GATED{i}) = [na Ka];
    fa = atpGateFactor(ATP2, Aref, Ka, na);
    if abs(p) >= 1e-3, fk = fa; else, fk = 1; end     % what the PTHRESH cut actually used
    fprintf('%-14s %12.4g %14.10f %14.10f %12.3g\n', GATED{i}, p, fp, fk, abs(fk-fp)/fp);
    worstKept = max(worstKept, abs(fk-fp)/fp);
    worstAll  = max(worstAll,  abs(fa-fp)/fp);
end
fprintf('\nlargest per-rate perturbation:  with the PTHRESH cut %.3g   converting all %.3g\n', ...
        worstKept, worstAll);

%% ---- 1. simulate: self-vs-self, PTHRESH-cut, and full conversion --------------
q = iCfg(b0);
fprintf('\n=== slack runs at %g mM ===\n', ATP2);
t=tic; [~,~,fA,fd] = runSlackExperiment(q);  fprintf('  power (run 1)      %4.0f s\n', toc(t));
t=tic; [~,~,fB,~ ] = runSlackExperiment(q);  fprintf('  power (run 2)      %4.0f s\n', toc(t));

qh = q; qh.AtpGateForm = 'hill';                      % full conversion, no threshold
for i = 1:numel(GATED)
    if isfield(hillAll, GATED{i})
        qh.(['pAtp_' GATED{i}]) = hillAll.(GATED{i})(1);
        qh.(['KAtp_' GATED{i}]) = hillAll.(GATED{i})(2);
    else
        qh.(['pAtp_' GATED{i}]) = 0; qh.(['KAtp_' GATED{i}]) = Inf;
    end
end
t=tic; [~,~,fC,~] = runSlackExperiment(qh); fprintf('  hill  (all gates)  %4.0f s\n', toc(t));

qk = getParams(loadParams('lowATP_hillgate_opt'), [], true, false); qk = iCfg(qk);
t=tic; [~,~,fD,~] = runSlackExperiment(qk); fprintf('  hill  (PTHRESH cut)%4.0f s\n', toc(t));

FN = q.fn;
cA = sum(evalFeatureCost(fd,fA,FN,2)); cB = sum(evalFeatureCost(fd,fB,FN,2));
cC = sum(evalFeatureCost(fd,fC,FN,2)); cD = sum(evalFeatureCost(fd,fD,FN,2));
fprintf('\ncost: power#1 %.6f | power#2 %.6f | hill(all) %.6f | hill(cut) %.6f\n', cA,cB,cC,cD);
fprintf('  H2 reproducibility floor  |power#1 - power#2| = %.3g\n', abs(cA-cB));
fprintf('  full conversion           |hill(all) - power#1| = %.3g\n', abs(cC-cA));
fprintf('  thresholded conversion    |hill(cut) - power#1| = %.3g\n', abs(cD-cA));

iCmp(fA,fB,'power#1 vs power#2   (pure noise floor)');
iCmp(fA,fC,'power#1 vs hill(all) (exact conversion)');
iCmp(fA,fD,'power#1 vs hill(cut) (thresholded)');

save(fullfile(resdir,'equivalence.mat'),'cA','cB','cC','cD','fA','fB','fC','fD','fd','FN','hillAll');
fprintf('\nEQUIVALENCE DIAGNOSTIC DONE\n');

%% ---- helpers -----------------------------------------------------------------
function q = iCfg(q)
q.RunForceVelocity=false; q.RunForceVelocityTime=false; q.RunKtr=false; q.RunStairs=false;
q.RunSlackPassive=false;  q.RunSlack=true; q.EvalFeatures=true; q.BreakOnODEUnstable=false;
q.RunSlackSegments='AllPar'; q.MaxRunTime=300;
q.PlotEachSeparately=0; q.PlotFeatureFitting=0;
q.UseAtpDetachR2=0; q.UseAtpKmsrd=0; q.UseStrainDepAtpGate=0; q.UseR2Ceiling=0; q.UseD0State=0;
q.MgATP=2; q.MgADP=0; q.Pi=0.8;
q.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
fn = q.fn; [~,iu] = unique(fn,'stable'); fn = fn(iu);
fn = cellfun(@(s) regexprep(s,'^rsK(\|.*)?$','rsK|0.1'), fn, 'UniformOutput', false);
q.fn = fn(~(startsWith(fn,'FV') | startsWith(fn,'PS_')));
end

function iCmp(fa, fb, label)
fns = intersect(fieldnames(fa), fieldnames(fb));
worst = 0; wname = '-';
for i = 1:numel(fns)
    x = fa.(fns{i}); y = fb.(fns{i});
    if ~isnumeric(x) || ~isnumeric(y) || ~isequal(size(x),size(y)) || isempty(x), continue; end
    d = max(abs(double(x(:)) - double(y(:))) ./ max(abs(double(x(:))), 1e-12));
    if isfinite(d) && d > worst, worst = d; wname = fns{i}; end
end
fprintf('  %-42s max rel. feature difference %.3g  (worst: %s)\n', label, worst, wname);
end
