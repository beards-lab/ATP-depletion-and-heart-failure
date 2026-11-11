% RunSrxCompensation.m -- can other knobs hide a different SRX pool?
%
%   matlab -batch "run('Analyses/SRXDuringContraction/RunSrxCompensation.m')"
%
% The total SRX pool (isometric SR+SRD at SL 2.2) is PINNED at 50 % and 200 % of the
% incumbent (0.030 -> 0.015 / 0.060) by the inflow (ksr0 & ksrd together). Other knobs are
% then refitted to bring the 8 mM features back to the data. The incumbent pool (x1) is put
% through the same refit as a CONTROL, so the verdict compares like with like. Two sets:
%   A  cross-bridge knobs only: ka kd k1 k_1 k2 kah kstiff2 A0 k_cr   (+ inflow, holds the pool)
%   B  A + the SRX KINETICS: outflow kmsrd, its load sensitivity sigma_srd1, inflow load
%      sensitivity sigma2                                           (pool still pinned)
% Method: knobs are log multipliers theta. Jacobian of the residuals (17 features, data-SD
% units, srxResiduals) by forward differences at the incumbent (step x1.15). Then two damped
% Gauss-Newton steps per (target, set), each verified by a full battery run. The pool
% constraint is a residual with SD 5 % of the target. The fit is judged by
% chi2 = sum z^2 over the 17 features: if a compensated fibre returns to the incumbent's
% chi2, the pool size is not identifiable from these data.
% Every evaluation is cached (results/comp/<tag>.mat), so a rerun resumes.
% Output -> results/compensation.mat (PlotSrxCompensation.m)
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsC = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results'); cdir = fullfile(rsC, 'comp');
if ~isfolder(cdir); mkdir(cdir); end
diary(fullfile(rsC, 'compensation.log')); diary on;
fprintf('SRX compensation on robust1_seed, started %s\n', datestr(now));
clC = parcluster('local'); jslC = fullfile(tempdir, sprintf('mljob_srxc_%d', feature('getpid')));
if ~isfolder(jslC); mkdir(jslC); end
clC.JobStorageLocation = jslC;
if isempty(gcp('nocreate')); parpool(clC, 5); end

pC = getParams(loadParams('robust1_seed'), [], true, false);
pC.RunForceVelocity=true; pC.RunSlack=true; pC.RunStairs=false; pC.RunKtr=true;
pC.RunForceVelocityTime=false; pC.RunSlackPassive=true; pC.EvalFeatures=true;
pC.BreakOnODEUnstable=false; pC.PlotEachSeparately=0; pC.PlotFeatureFitting=0;
pC.RunSlackSegments='AllPar'; pC.MaxRunTime=400; pC.OptimizeOn='Feats';
pC.ghostSave=''; pC.ghostLoad=''; pC.SaveBest=0; pC.mods={}; pC.g=ones(1,0);

KN = {'ka', {'ka'}; 'kd', {'kd'}; 'k1', {'k1'}; 'k_1', {'k_1'}; 'k2', {'k2'}; 'kah', {'kah'}; ...
      'kstiff2', {'kstiff2'}; 'A0', {'A0'}; 'k_cr', {'k_cr'}; 'inflow', {'ksr0', 'ksrd'}; ...
      'kmsrd', {'kmsrd'}; 'sigma_srd1', {'sigma_srd1'}; 'sigma2', {'sigma2'}};
nK = size(KN, 1); iIn = find(strcmp(KN(:, 1), 'inflow'));
SETS = struct('name', {'A', 'B'}, 'idx', {[1:9, iIn], 1:nK});
TGT = [0.5 1 2];   % x1 = CONTROL: the incumbent pool refitted the same way
H = log(1.15); LAMBDA = 0.5; STEPMAX = log(4); WPOOL = 1/0.05;
ev = @(th, tag) iEval(pC, KN, th, tag, cdir);

%% ---- incumbent and Jacobian ------------------------------------------------------
E0 = ev(zeros(nK, 1), 'base');
fprintf('base: chi2 %.3f  cost %.4f  SRX %.4f  F %.2f\n', E0.chi2, E0.cost, E0.srx, E0.F);
nz = numel(E0.z); Jz = zeros(nz, nK); Js = zeros(1, nK);
for j = 1:nK
    th = zeros(nK, 1); th(j) = H; Ej = ev(th, sprintf('jac_%s', KN{j, 1}));
    dz = (Ej.z - E0.z) / H; dz(~isfinite(dz)) = 0; Jz(:, j) = dz;
    Js(j) = (log(Ej.srx) - log(E0.srx)) / H;
    fprintf('J %-10s  dlogSRX/dlog %+6.2f  |dz| %6.2f\n', KN{j, 1}, Js(j), norm(dz));
end

%% ---- pin the pool, then compensate -------------------------------------------------
R = struct('target', {}, 'set', {}, 'theta', {}, 'E', {}, 'hist', {});
U = cell(1, numel(TGT));
for it = 1:numel(TGT)
    T = TGT(it) * E0.srx;
    % inflow factor for the target pool, isometric runs only (cheap)
    th = zeros(nK, 1); th(iIn) = log(TGT(it)) / Js(iIn);
    for k = 1:4
        s = srxIsoState(srxApplyTheta(pC, KN, th)); e = log(s.SRX) - log(T);
        if abs(e) < 0.01; break; end
        th(iIn) = th(iIn) - e / Js(iIn);
    end
    U{it} = ev(th, sprintf('pin_%g', TGT(it)));
    fprintf('pool x%g: inflow x%.3f  SRX %.4f  chi2 %.2f (uncompensated)\n', TGT(it), exp(th(iIn)), U{it}.srx, U{it}.chi2);
    for iS = 1:numel(SETS)
        S = SETS(iS).idx; thS = th; E = U{it}; hist = {E};
        for step = 1:2
            r = [E.z; WPOOL * (log(E.srx) - log(T))]; ok = isfinite(r);
            J = [Jz(:, S); WPOOL * Js(S)]; J = J(ok, :); r = r(ok);
            d = -(J' * J + LAMBDA * eye(numel(S))) \ (J' * r);
            d = max(-STEPMAX, min(STEPMAX, d));
            thN = thS; thN(S) = thN(S) + d;
            EN = ev(thN, sprintf('comp_%g_%s_%d', TGT(it), SETS(iS).name, step)); hist{end+1} = EN; %#ok<SAGROW>
            fprintf('  set %s step %d: chi2 %.2f  SRX %.4f  F %.2f\n', SETS(iS).name, step, EN.chi2, EN.srx, EN.F);
            if EN.chi2 < E.chi2; thS = thN; E = EN; else; break; end   % keep the better point
        end
        R(end+1) = struct('target', TGT(it), 'set', SETS(iS).name, 'theta', thS, 'E', E, 'hist', {hist}); %#ok<SAGROW>
    end
end

%% ---- traces of the final points ------------------------------------------------------
TH = [{zeros(nK, 1)}, cellfun(@(u) u.theta, U, 'uni', 0), {R.theta}];
TL = [{'incumbent'}, arrayfun(@(t) sprintf('pool x%g, uncompensated', t), TGT, 'uni', 0), ...
      arrayfun(@(r) sprintf('pool x%g, set %s', r.target, r.set), R, 'uni', 0)];
TR = cell(numel(TH), 3);
parfor i = 1:numel(TH)
    [s, a, b] = srxTraces(srxApplyTheta(pC, KN, TH{i}), root); TR(i, :) = {s, a, b};
end
save(fullfile(rsC, 'compensation.mat'), 'KN', 'SETS', 'TGT', 'E0', 'Jz', 'Js', 'U', 'R', 'TH', 'TL', 'TR', ...
     'H', 'LAMBDA', 'WPOOL');
fprintf('DONE %s\n', datestr(now)); diary off;

%% ---- local functions ------------------------------------------------------------
function E = iEval(p, KN, th, tag, cdir)
    f = fullfile(cdir, [tag '.mat']);
    if isfile(f); L = load(f); if isequal(L.E.theta, th); E = L.E; return; end; end
    if isempty(gcp('nocreate')); parpool(5); end
    t0 = tic; q = srxApplyTheta(p, KN, th);
    s = srxIsoState(q); q.SRXT_0 = s.SR; q.SRXD_0 = s.SRD;
    [fm, fd] = iRun(q);
    [z, names] = srxResiduals(fm, fd);
    E = struct('tag', tag, 'theta', th, 'z', z, 'names', {names}, 'chi2', sum(z(isfinite(z)).^2), ...
               'cost', sum(evalFeatureCost(fd, fm, p.fn)), 'srx', s.SRX, 'F', s.F, 'R2T', s.R2T, ...
               'fm', fm, 'fd', fd, 'sec', toc(t0));
    save(f, 'E');
    fprintf('    [%s] chi2 %.2f  SRX %.4f  (%.0f s)\n', tag, E.chi2, E.srx, E.sec);
end

function [fm, fd] = iRun(p)
    % RunBakersExp is a SCRIPT; a function workspace keeps it from clobbering ours.
    params0 = p; %#ok<NASGU>
    LoadData; RunBakersExp; %#ok<*ASGLU>
    fm = features_model; fd = features_data;
end
