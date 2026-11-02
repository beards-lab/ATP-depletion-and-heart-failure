% RunSrxFamily.m -- is the SRX fraction during contraction identifiable from the
% 8 mM battery (slack-restretch, ktr protocol, F-V, passive)?
%
%   matlab -batch "run('Analyses/SRXDuringContraction/RunSrxFamily.m')"
%
% Three one-parameter families move the SRX on top of robust1_seed:
%   L  (size)          ksr0 and ksrd x lam       -- more recruitment at every load
%   S  (mechanosensing) sigma2 and sigma_srd1 x s -- weaker load dependence: active
%                                                   SRX rises, rest SRX barely moves
%   M  (speed, partial) ksr0, ksrd, kmsrd x mu    -- entry/escape only; the pool SHRINKS
%                                                   (0.030 -> 0.008, see srxApplyVariant)
%   MA (speed)         every SRX rate x mu        -- pool size kept, slower exchange
%                                                   (mu -> 0 is the mant-like frozen pool)
% Every variant is COMPENSATED: ka is re-solved so the isometric force at SL 2.2 equals
% the incumbent's (secant, 3 s warm start), and SRXT_0/SRXD_0 are set to the variant's
% own isometric steady state, so a slow pool starts where it would sit anyway. The
% reference (lam = s = mu = 1) is run through the same path as a regression check.
% This is a ONE-KNOB compensation (ka only), not a refit: the cost differences are an
% UPPER bound on how much the data resist a given SRX level.
%
% Scoring: the incumbent objective (8 mM, pooled target), reported (a) as is and (b) as
% the DATA cost = without the SRX_ss literature bound, which is what identifiability
% means here.
% Outputs -> results/family_run_<k>.mat, results/family.mat (SummarizeSrxFamily.m)
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsF = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
diary(fullfile(rsF, 'family.log')); diary on;
fprintf('SRX family on robust1_seed, started %s\n', datestr(now));

clF = parcluster('local');
jslF = fullfile(tempdir, sprintf('mljob_srx_%d', feature('getpid')));
if ~isfolder(jslF); mkdir(jslF); end
clF.JobStorageLocation = jslF;
if isempty(gcp('nocreate')); parpool(clF, 5); end

pF = getParams(loadParams('robust1_seed'), [], true, false);
pF.RunForceVelocity=true; pF.RunSlack=true; pF.RunStairs=false; pF.RunKtr=true;
pF.RunForceVelocityTime=false; pF.RunSlackPassive=true; pF.EvalFeatures=true;
pF.BreakOnODEUnstable=false; pF.PlotEachSeparately=0; pF.PlotFeatureFitting=0;
pF.RunSlackSegments='AllPar'; pF.MaxRunTime=400; pF.OptimizeOn='Feats';
pF.ghostSave=''; pF.ghostLoad=''; pF.SaveBest=0; pF.mods={}; pF.g=ones(1,0);
fnF = pF.fn;
fnNoSRX = regexprep(fnF, 'SRX_ss\[1\]\|[^,]*,?', '');       % data cost: no SRX literature bound

if ~exist('VARS', 'var')
    VARS = struct('fam', {'ref', 'L','L','L','L', 'S','S','S','S','S', 'M','M','M', 'MA','MA','MA'}, ...
                  'x',   { 1,     0.25, 3, 10, 30, 1.3, 1.7, 2.2, 3.0, 4.0, 0.1, 0.01, 0.001, 0.1, 0.01, 0.001});
end

%% ---- compensation (parallel over variants, isometric only) ----------------------
[F_target, ss_ref] = srxIsoForce(pF, pF.ka, [], 3);
fprintf('incumbent isometric F = %.3f kPa, SRX %.4f\n', F_target, ss_ref(1) + ss_ref(2));
nV = numel(VARS); kaV = zeros(1, nV); icV = zeros(nV, 2); FV_ = zeros(1, nV);
parfor iV = 1:nV
    q = srxApplyVariant(pF, VARS(iV));
    ic = [];
    if startsWith(VARS(iV).fam, 'M'); ic = ss_ref; end        % slow pool: start at the reference state
    ka = q.ka; [F1, s1] = srxIsoForce(q, ka, ic, 3); ic = s1;
    ka2 = ka * F_target / F1; [F2, s2] = srxIsoForce(q, ka2, ic, 3); ic = s2;
    for it = 1:4                                          % secant on ka
        if abs(F2 - F_target) < 0.05; break; end
        ka3 = ka2 + (F_target - F2) * (ka2 - ka) / (F2 - F1);
        ka = ka2; F1 = F2; ka2 = max(ka3, 1); [F2, s2] = srxIsoForce(q, ka2, ic, 3); ic = s2;
    end
    kaV(iV) = ka2; icV(iV, :) = s2; FV_(iV) = F2;
end
for iV = 1:nV
    fprintf('%-3s x=%6.3g  ka %.2f (x%.3f)  F %.3f  SRX_iso %.4f\n', VARS(iV).fam, VARS(iV).x, ...
            kaV(iV), kaV(iV)/pF.ka, FV_(iV), sum(icV(iV, :)));
end

%% ---- full 8 mM battery per variant ----------------------------------------------
RF = struct('fam',{},'x',{},'ka',{},'ic',{},'fm',{},'fd',{},'tr',{},'cost',{},'costNoSRX',{}, ...
            'total',{},'dataTotal',{},'sec',{});
for iV = 1:nV
    fR = fullfile(rsF, sprintf('family_run_%d.mat', iV));
    if isfile(fR)                                         % resume: keep a finished variant
        L = load(fR);
        if strcmp(L.R1.fam, VARS(iV).fam) && abs(L.R1.x - VARS(iV).x) < 1e-12
            RF(iV) = L.R1; fprintf('%-3s x=%6.3g: kept from %s\n', VARS(iV).fam, VARS(iV).x, fR); continue
        end
    end
    if isempty(gcp('nocreate')); parpool(clF, 5); end     % a lost pool would silently run serially
    tV = tic;
    q = srxApplyVariant(pF, VARS(iV)); q.ka = kaV(iV); q.SRXT_0 = icV(iV, 1); q.SRXD_0 = icV(iV, 2);
    try
        [fm, fd, tr] = iRun(q);
        c = evalFeatureCost(fd, fm, fnF); cN = evalFeatureCost(fd, fm, fnNoSRX);
    catch eV
        warning('variant %d failed: %s', iV, eV.message); fm = struct(); fd = struct(); tr = struct();
        c = NaN; cN = NaN;
    end
    RF(iV).fam = VARS(iV).fam; RF(iV).x = VARS(iV).x; RF(iV).ka = kaV(iV); RF(iV).ic = icV(iV, :);
    RF(iV).fm = fm; RF(iV).fd = fd; RF(iV).tr = tr; RF(iV).cost = c; RF(iV).costNoSRX = cN;
    RF(iV).total = sum(c); RF(iV).dataTotal = sum(cN); RF(iV).sec = toc(tV);
    srx = NaN; if isfield(fm, 'SRX_ss'); srx = fm.SRX_ss(1); end
    fprintf('%-3s x=%6.3g: cost %.4f  data cost %.4f  SRX_ss %.4f  (%.0f s)\n', RF(iV).fam, RF(iV).x, ...
            RF(iV).total, RF(iV).dataTotal, srx, RF(iV).sec);
    R1 = RF(iV); save(fullfile(rsF, sprintf('family_run_%d.mat', iV)), 'R1'); %#ok<NASGU>
end
save(fullfile(rsF, 'family.mat'), 'RF', 'VARS', 'fnF', 'fnNoSRX', 'F_target', 'ss_ref');
fprintf('DONE %s\n', datestr(now)); diary off;

%% ---- local functions ------------------------------------------------------------
function [fm, fd, tr] = iRun(p)
    % RunBakersExp is a SCRIPT; a function workspace keeps it from clobbering ours.
    params0 = p;                                               %#ok<NASGU>
    LoadData; RunBakersExp;                                    %#ok<*ASGLU>
    fm = features_model; fd = features_data;
    % NB: in 'AllPar' mode out_slack is not the slack cycles (it comes back as one flat
    % isometric segment), so tr is not a usable SRX trace; RunSrxProtocols has the traces.
    tr = struct();
    try
        k  = makeMonotonous(out_slack.t(:));
        tt = out_slack.t(k); tq = (tt(1):5e-4:tt(end))';
        tr.t = tq; tr.F = interp1(tt, out_slack.Force(k), tq, 'linear');
        tr.SL = interp1(tt, out_slack.SL(k), tq, 'linear');
        tr.SRX = interp1(tt, out_slack.SR(k) + out_slack.SRD(k), tq, 'linear');
        tr.PT = interp1(tt, out_slack.PuATP(k), tq, 'linear');
        tr.att = interp1(tt, out_slack.p1_0(k) + out_slack.p2_0(k), tq, 'linear');
        if isfield(out_slack, 'datatable'); tr.data = out_slack.datatable(1:4:end, :); end
    catch eTr
        warning('trace capture failed: %s', eTr.message);
    end
end
