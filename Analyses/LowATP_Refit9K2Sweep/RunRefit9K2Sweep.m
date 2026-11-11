% RunRefit9K2Sweep.m -- low ATP on the current 8 mM incumbent, the rough way:
% only A2 detachment feels ATP, and we sweep how much.
%
%   matlab -batch "run('Analyses/LowATP_Refit9K2Sweep/RunRefit9K2Sweep.m')"
%
% DESIGN
%   [MgATP] takes exactly two values: 8 mM (reference) and 2 mM. The single
%   knob is g = R2(2 mM)/R2(8 mM), the factor on the A2-detachment rate. It is
%   implemented with the merged power gate pAtp_k2 (evaluateModel), which is
%   exactly 1 at MgATP_ref = 8 for any exponent, so the 8 mM leg is untouched:
%       g = (2/8)^p   ->   p = log(g)/log(1/4)
%   Hydrolysis (kah/kamh) is NOT gated -- the physiological choice.
%
% SCORING -- fit the 2 mM data, not the ratio
%   A 2/8 ratio doubles every 8 mM miss, so the 2 mM leg is scored in ABSOLUTE
%   terms against the pooled 2 mM target (pooledTarget2mM, pooledKtrTarget(2)),
%   with the incumbent's objective moved onto the 2 mM spread
%   (rescaleFnToCondition). Two numbers per run:
%     score_plain  sum of all 2 mM chi-2 terms
%     score_rel    each term x a = min(z,3)/3, z = |r-1|/sqrt(s^2+sr^2) the
%                  data's own ATP-effect z-score -> terms the data say do not
%                  move with ATP drop out; boundary/guard terms keep weight 1
%   plus FV2 (the in-preparation 2 ML/s force-velocity point, pooled over the
%   two preparations that have it) as an extra chi-2 term, and the 2/8 ratios
%   as a diagnostic.
%
% Outputs -> results/k2sweep_run_<k>.mat (one per run, written as it finishes),
%            results/k2sweep.mat, results/k2sweep.txt, results/k2sweep_summary.json,
%            results/k2sweep_*.png (PlotK2Sweep.m)

if ~exist('CORE','var');  CORE  = 'robust1_seed'; end
if ~exist('GLIST','var'); GLIST = [1 0.75 0.54 0.44 0.35 0.30]; end
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resK2 = fullfile(root,'Analyses','LowATP_Refit9K2Sweep','results');
diary(fullfile(resK2,'k2sweep.log')); diary on;
fprintf('k2 sweep on %s, started %s\n', CORE, datestr(now));

clK2 = parcluster('local');
jslK2 = fullfile(tempdir, sprintf('mljob_k2_%d', feature('getpid')));
if ~isfolder(jslK2); mkdir(jslK2); end
clK2.JobStorageLocation = jslK2;
if isempty(gcp('nocreate')); parpool(clK2, 5); end

%% ---- the two legs -------------------------------------------------------
pK2 = getParams(loadParams(CORE), [], true, false);
pK2.RunForceVelocity=true; pK2.RunSlack=true; pK2.RunStairs=false; pK2.RunKtr=true;
pK2.RunForceVelocityTime=false; pK2.RunSlackPassive=false; pK2.EvalFeatures=true;
pK2.BreakOnODEUnstable=false; pK2.PlotEachSeparately=0; pK2.PlotFeatureFitting=0;
pK2.RunSlackSegments='AllPar'; pK2.MaxRunTime=400; pK2.OptimizeOn='Feats';
pK2.ghostSave=''; pK2.ghostLoad=''; pK2.SaveBest=0; pK2.mods={}; pK2.g=ones(1,0);
pK2.recalculateDataFeats=false; pK2.UsePooledTarget=true;
pK2.UseAtpPowerGates=true; pK2.MgATP_ref=8; pK2.pAtp_k2=0;
pK2.UseAtpDetachR2=0; pK2.UseAtpKmsrd=0;

primK2 = @(c) cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), c, 'uni', 0);
fn8K2 = pK2.fn; [~,iuK2] = unique(fn8K2,'stable'); fn8K2 = fn8K2(iuK2);
fn8K2 = fn8K2(~startsWith(primK2(fn8K2),'PS_'));          % passive (PNB) leg not run here
fn2K2 = rescaleFnToCondition(fn8K2, 2);
% Two 8 mM-only terms leave the 2 mM leg:
%   FV_normpowerAvg -- runFVExperiment always simulates ATP_c(1) (8 mM) against
%                      Baker's 8 mM curve, so it is not a 2 mM measurement; the
%                      FV2 point below is the 2 mM force-velocity term instead.
%   steady[1] guard -- its 52.5-87.7 kPa window is an 8 mM window; the pooled
%                      2 mM steady is ~90 kPa, so it would punish a correct fit.
%                      The steady chi-2 term still scores the 2 mM level.
drop2 = startsWith(fn2K2, 'FV_normpowerAvg') | startsWith(fn2K2, 'steady[1]');
fn2K2 = fn2K2(~drop2);

p8K2 = pK2; p8K2.MgATP = 8; p8K2.fn = fn8K2;
p8K2.velocitytableonfile     = 'protocol_03_27_2026_8mM_slack.mat';
p8K2.ktr_velocitytableonfile = 'protocol_03_27_2026_velocitytable_ktr.mat';
p2K2 = p8K2; p2K2.MgATP = 2; p2K2.fn = fn2K2; p2K2.RunForceVelocity = false;
p2K2.velocitytableonfile     = 'protocol_03_27_2026_2mM_slack.mat';
p2K2.ktr_velocitytableonfile = 'protocol_03_27_2026_velocitytable_ktr_2mM.mat';

%% ---- ATP relevance per term: a = min(z,3)/3 -------------------------------
PT8K2 = pooledTarget8mM(); PT2K2 = pooledTarget2mM();
KT8K2 = pooledKtrTarget(8); KT2K2 = pooledKtrTarget(2);
nmK2  = primK2(fn2K2);
aK2   = ones(1, numel(fn2K2)); zK2 = nan(1, numel(fn2K2));
for iK2 = 1:numel(fn2K2)
    isBnd = contains(fn2K2{iK2}, ',') || ~isempty(regexp(fn2K2{iK2}, '\|-?[\d.]+-{1}-?[\d.]+(\||$)', 'once'));
    if isBnd; continue; end                               % guards / bounds keep weight 1
    nmZ = nmK2{iK2}; if strcmp(nmZ, 't0_crossing'); nmZ = 't0'; end   % same measurement
    j = find(strcmp(PT8K2.name, nmZ), 1);
    if ~isempty(j)
        zK2(iK2) = PT8K2.z(j);
    else
        j = find(strcmp(KT8K2.name, nmK2{iK2}), 1);
        if ~isempty(j)
            r  = KT2K2.val(j)/KT8K2.val(j);
            sr = r*sqrt(KT2K2.s(j)^2 + KT8K2.s(j)^2);
            zK2(iK2) = abs(r-1)/sqrt(KT8K2.s(j)^2 + sr^2);
        elseif any(strcmp(nmK2{iK2}, {'ktrfit_shape','ktr80_shape','A_over_steady'}))
            zK2(iK2) = NaN;                               % no measured ratio -> keep 1
        end
    end
    if isfinite(zK2(iK2)); aK2(iK2) = min(zK2(iK2),3)/3; end
end

%% ---- FV2: the in-preparation force-velocity point, 2 preparations ------------
FV2K2 = fv2Target();                      % pooled P/P0 at 2 ML/s, 8 and 2 mM
fprintf('FV2 target P/P0: 8 mM %.3f +- %.3f | 2 mM %.3f +- %.3f  (n = %d)\n', ...
        FV2K2.m8, FV2K2.s8, FV2K2.m2, FV2K2.s2, FV2K2.n);

%% ---- run ------------------------------------------------------------------------
t0K2 = tic;
[fm8K2, fd8K2, tr8K2] = iRun(p8K2);
c8K2 = evalFeatureCost(fd8K2, fm8K2, fn8K2);
fv8K2 = fv2Model(p8K2);
fprintf('8 mM: cost %.4f  FV2 %.3f  (%.0f s)\n', sum(c8K2), fv8K2, toc(t0K2));
save(fullfile(resK2,'k2sweep_run_8mM.mat'), 'fm8K2','fd8K2','tr8K2','c8K2','fv8K2','fn8K2');

RK2 = struct('g',{},'p',{},'fm',{},'fd',{},'tr',{},'cost',{},'fv2',{},'score_plain',{},'score_rel',{},'sec',{});
for kK2 = 1:numel(GLIST)
    tK2 = tic;
    gK2 = GLIST(kK2);
    q   = p2K2; q.pAtp_k2 = log(gK2)/log(1/4);
    [fmK2, fdK2, trK2] = iRun(q);
    cK2 = evalFeatureCost(fdK2, fmK2, fn2K2);
    fvK2 = fv2Model(q);
    cFV  = ((fvK2 - FV2K2.m2)/FV2K2.s2)^2;                 % chi-2 of the FV2 point
    zFV  = abs(FV2K2.r-1)/sqrt(FV2K2.s8r^2 + FV2K2.sr^2);
    RK2(kK2).g = gK2; RK2(kK2).p = q.pAtp_k2; RK2(kK2).fm = fmK2; RK2(kK2).fd = fdK2; RK2(kK2).tr = trK2;
    RK2(kK2).cost = cK2; RK2(kK2).fv2 = fvK2;
    RK2(kK2).score_plain = sum(cK2) + cFV;
    RK2(kK2).score_rel   = sum(cK2(:)'.*aK2) + min(zFV,3)/3*cFV;
    RK2(kK2).sec = toc(tK2);
    fprintf('g = %.2f (p = %.3f): plain %.4f  rel %.4f  FV2 %.3f  (%.0f s)\n', gK2, q.pAtp_k2, ...
            RK2(kK2).score_plain, RK2(kK2).score_rel, fvK2, RK2(kK2).sec);
    R1 = RK2(kK2); save(fullfile(resK2, sprintf('k2sweep_run_%d.mat', kK2)), 'R1'); %#ok<NASGU>
end
save(fullfile(resK2,'k2sweep.mat'), 'RK2','fm8K2','fd8K2','tr8K2','c8K2','fv8K2','fn8K2','fn2K2', ...
     'aK2','zK2','FV2K2','CORE','GLIST');
fprintf('sweep done in %.0f s\n', toc(t0K2));
SummarizeK2Sweep;           % text table + JSON + figures
fprintf('DONE\n'); diary off;

%% ---- local functions ------------------------------------------------------------
function [fm, fd, tr] = iRun(p)
    % RunBakersExp is a SCRIPT; a function workspace keeps it from clobbering ours.
    params0 = p;                                               %#ok<NASGU>
    LoadData; RunBakersExp;                                    %#ok<*ASGLU>
    fm = features_model; fd = features_data;
    % thinned slack trace (model + data) for the figures; never let it kill a run
    tr = struct();
    try
    k  = makeMonotonous(out_slack.t(:));
    tt = out_slack.t(k); tq = (tt(1):5e-4:tt(end))';
    tr.t = tq; tr.F = interp1(tt, out_slack.Force(k), tq, 'linear');
    tr.SL = interp1(tt, out_slack.SL(k), tq, 'linear');
    if isfield(out_slack,'datatable'); tr.data = out_slack.datatable(1:4:end,:); end
    catch eTr
        warning('trace capture failed: %s', eTr.message);
    end
end
