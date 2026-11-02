% RunSrxTwin.m -- would the proposed identification protocol pin the SRX pool to +-20 %?
%
%   matlab -batch "run('Analyses/SRXDuringContraction/RunSrxTwin.m')"
%
% TWIN (synthetic-data) test. Truth = robust1_seed (pool 0.030); synthetic "data" = the
% model's own features, so the truth sits at chi2 = 0 and every residual is pure
% identifiability, never model misfit. Noise = one preparation:
%   battery features   the pooled between-prep SD (srxResiduals)
%   protocol B k63(T)  8 %;  F20(T) 5 %          (release-hold-restretch, T = 5..80 ms)
%   protocol A t50/t90 10 %                     (activation step from rest)
% Three analyses:
%   1  Fisher information (linearised): SD of log(pool) with all 13 knobs free, for data
%      sets battery / +A / +B / +A+B, and how it scales with preparations.
%   2  Twin profile: pool pinned at x0.8 and x1.2 (inflow), all 13 knobs refitted by
%      damped Gauss-Newton (2 verified steps) to (i) the battery only, (ii) battery + B.
%      Every point is scored on battery, B and A: does B expose the battery-only refit?
%   3  The proposed "noised seed" test in linearised form (Monte Carlo): measurement noise
%      + a 5 % noised seed on the 12 nuisance knobs, with and without the seed acting as a
%      prior, to show what that test measures.
% Battery runs reuse the RunSrxCompensation cache (results/comp/); protocol runs are cached
% in results/twin/. Output -> results/twin.mat (PlotSrxTwin.m)
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsT = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
cdirC = fullfile(rsT, 'comp'); cdirT = fullfile(rsT, 'twin'); if ~isfolder(cdirT); mkdir(cdirT); end
diary(fullfile(rsT, 'twin.log')); diary on;
fprintf('SRX twin identifiability, started %s\n', datestr(now));
clT = parcluster('local'); jslT = fullfile(tempdir, sprintf('mljob_srxt_%d', feature('getpid')));
if ~isfolder(jslT); mkdir(jslT); end
clT.JobStorageLocation = jslT;
if isempty(gcp('nocreate')); parpool(clT, 5); end

pT = getParams(loadParams('robust1_seed'), [], true, false);
pT.RunForceVelocity=true; pT.RunSlack=true; pT.RunStairs=false; pT.RunKtr=true;
pT.RunForceVelocityTime=false; pT.RunSlackPassive=true; pT.EvalFeatures=true;
pT.BreakOnODEUnstable=false; pT.PlotEachSeparately=0; pT.PlotFeatureFitting=0;
pT.RunSlackSegments='AllPar'; pT.MaxRunTime=400; pT.OptimizeOn='Feats';
pT.ghostSave=''; pT.ghostLoad=''; pT.SaveBest=0; pT.mods={}; pT.g=ones(1,0);

KN = {'ka', {'ka'}; 'kd', {'kd'}; 'k1', {'k1'}; 'k_1', {'k_1'}; 'k2', {'k2'}; 'kah', {'kah'}; ...
      'kstiff2', {'kstiff2'}; 'A0', {'A0'}; 'k_cr', {'k_cr'}; 'inflow', {'ksr0', 'ksrd'}; ...
      'kmsrd', {'kmsrd'}; 'sigma_srd1', {'sigma_srd1'}; 'sigma2', {'sigma2'}};   % = RunSrxCompensation
nK = size(KN, 1); iIn = find(strcmp(KN(:, 1), 'inflow')); H = log(1.15);
TB = [5 10 20 40 80] * 1e-3; nT = numel(TB);
TGT = [0.8 1.2]; LAMBDA = 0.5; STEPMAX = log(4); WPOOL = 1/0.02;     % pool pinned to 2 %
ev = @(th, tag) iEval(pT, KN, th, tag, cdirC, cdirT, TB);

%% ---- truth and Jacobian ------------------------------------------------------------
E0 = ev(zeros(nK, 1), 'base');
[zb0, ~, ~, mb0, ~, sdb] = srxResiduals(E0.fm, E0.fd);  %#ok<ASGLU>
p0 = E0.P; sdP = [0.08 * p0.k63, 0.05 * p0.F20]; sdA = 0.10 * [p0.t50, p0.t90];
vecB = @(E) [E.P.k63, E.P.F20]; vecA = @(E) [E.P.t50, E.P.t90];
res = @(E) struct('bat', (iM(E) - mb0) ./ sdb, 'B', ((vecB(E) - vecB(E0)) ./ sdP)', ...
                  'A', ((vecA(E) - vecA(E0)) ./ sdA)');
fprintf('truth: SRX %.4f  F %.2f  k63(T) %s  F20 %s  t50 %.1f ms\n', E0.srx, E0.F, mat2str(p0.k63, 3), mat2str(p0.F20, 2), 1e3*p0.t50);
Jb = zeros(numel(mb0), nK); JB = zeros(2*nT, nK); JA = zeros(2, nK); Js = zeros(1, nK);
for j = 1:nK
    th = zeros(nK, 1); th(j) = H; Ej = ev(th, sprintf('jac_%s', KN{j, 1})); r = res(Ej);
    Jb(:, j) = r.bat / H; JB(:, j) = r.B / H; JA(:, j) = r.A / H; Js(j) = (log(Ej.srx) - log(E0.srx)) / H;
    fprintf('J %-10s |dz| battery %6.2f  B %6.2f  A %5.2f  dlogSRX %+5.2f\n', KN{j, 1}, sqrt(sum(Jb(:, j).^2, 'omitnan')), norm(JB(:, j)), norm(JA(:, j)), Js(j));
end
Jb(~isfinite(Jb)) = 0; JB(~isfinite(JB)) = 0; JA(~isfinite(JA)) = 0;   % a failed run carries no information

%% ---- 1 Fisher information: SD of log(pool) --------------------------------------------
CFG = struct('name', {'battery', 'battery + A', 'battery + B', 'battery + A + B'}, ...
             'J', {Jb, [Jb; JA], [Jb; JB], [Jb; JA; JB]});
for c = 1:numel(CFG)
    [CFG(c).sd, CFG(c).corr, CFG(c).cond] = iPoolSD(CFG(c).J, Js, 0);
    CFG(c).sd_prior = iPoolSD(CFG(c).J, Js, 1/log(2)^2);    % weak prior: nuisance known within x2
    fprintf('FIM %-16s SD(log pool) %.3f (x%.2f)   with x2 prior %.3f   cond %.1e\n', CFG(c).name, ...
            CFG(c).sd, exp(CFG(c).sd), CFG(c).sd_prior, CFG(c).cond);
end

%% ---- 2 twin profile at x0.8 / x1.2 ----------------------------------------------------
PR = struct('target', {}, 'fit', {}, 'theta', {}, 'E', {}, 'r', {}, 'chi2', {}, 'hist', {});
PIN = cell(1, numel(TGT));
for it = 1:numel(TGT)
    T = TGT(it) * E0.srx; th = zeros(nK, 1); th(iIn) = log(TGT(it)) / Js(iIn);
    for k = 1:4
        s = srxIsoState(srxApplyTheta(pT, KN, th)); e = log(s.SRX) - log(T);
        if abs(e) < 0.005; break; end
        th(iIn) = th(iIn) - e / Js(iIn);
    end
    Eu = ev(th, sprintf('twin_pin_%g', TGT(it))); PIN{it} = struct('theta', th, 'E', Eu, 'r', res(Eu));
    fprintf('pool x%g pinned (SRX %.4f): chi2 battery %.2f  B %.2f  A %.2f\n', TGT(it), Eu.srx, ...
            iC(res(Eu).bat), iC(res(Eu).B), iC(res(Eu).A));
    for fit = {'bat', 'batB'}
        thS = th; E = Eu; r = res(E); hist = {E};
        J = Jb; if strcmp(fit{1}, 'batB'); J = [Jb; JB]; end
        rv = @(r_) iStack(r_, fit{1});
        for step = 1:2
            rr = [rv(r); WPOOL * (log(E.srx) - log(T))]; ok = isfinite(rr);
            JJ = [J; WPOOL * Js]; JJ = JJ(ok, :); rr = rr(ok);
            d = -(JJ' * JJ + LAMBDA * eye(nK)) \ (JJ' * rr); d = max(-STEPMAX, min(STEPMAX, d));
            thN = thS + d; EN = ev(thN, sprintf('twin_%g_%s_%d', TGT(it), fit{1}, step)); rN = res(EN); hist{end+1} = EN; %#ok<SAGROW>
            fprintf('  fit %-4s step %d: chi2 battery %.2f  B %.2f  A %.2f  SRX %.4f\n', fit{1}, step, ...
                    iC(rN.bat), iC(rN.B), iC(rN.A), EN.srx);
            if iC(rv(rN)) < iC(rv(r)); thS = thN; E = EN; r = rN; else; break; end
        end
        PR(end+1) = struct('target', TGT(it), 'fit', fit{1}, 'theta', thS, 'E', E, 'r', r, ...
                           'chi2', [iC(r.bat), iC(r.B), iC(r.A)], 'hist', {hist}); %#ok<SAGROW>
    end
end

%% ---- 3 the noised-seed test, linearised Monte Carlo ------------------------------------
rng(1); NMC = 4000; nui = setdiff(1:nK, iIn); SEED = 0.05;
MC = struct('cfg', {}, 'prior', {}, 'truth', {}, 'est', {});
for c = [1 3]
    J = CFG(c).J; nr = size(J, 1);
    for tr = log(TGT)
        th0 = zeros(nK, 1); th0(iIn) = tr / Js(iIn);             % truth: pool moved by the inflow
        for pr = {'none', 'seed 5 %'}
            est = zeros(1, NMC);
            for n = 1:NMC
                y = J * th0 + randn(nr, 1);                        % synthetic data + one-prep noise
                seed = th0; seed(nui) = seed(nui) + SEED * randn(numel(nui), 1);
                P = 1e-8 * eye(nK); if strcmp(pr{1}, 'seed 5 %'); P(nui, nui) = eye(numel(nui)) / SEED^2; end
                thh = (J' * J + P) \ (J' * y + P * seed);
                est(n) = Js * thh;
            end
            MC(end+1) = struct('cfg', CFG(c).name, 'prior', pr{1}, 'truth', tr, 'est', est); %#ok<SAGROW>
            fprintf('MC %-12s prior %-8s truth x%.2f: estimate x%.2f, SD %.3f (x%.2f)\n', CFG(c).name, pr{1}, ...
                    exp(tr), exp(mean(est)), std(est), exp(std(est)));
        end
    end
end
save(fullfile(rsT, 'twin.mat'), 'KN', 'TB', 'TGT', 'E0', 'Jb', 'JB', 'JA', 'Js', 'CFG', 'PIN', 'PR', 'MC', ...
     'sdb', 'sdP', 'sdA', 'mb0', 'WPOOL', 'LAMBDA');
fprintf('DONE %s\n', datestr(now)); diary off;

%% ---- local functions ------------------------------------------------------------
function c = iC(z); c = sum(z(isfinite(z)).^2); end
function v = iStack(r, fit); if strcmp(fit, 'bat'); v = r.bat; else; v = [r.bat; r.B]; end; end
function m = iM(E); [~, ~, ~, m] = srxResiduals(E.fm, E.fd); end

function [sd, corr, cn] = iPoolSD(J, c, prior)
    % SD of the pool functional c*theta with all knobs free (Cramer-Rao, unit noise)
    F = J' * J + prior * eye(size(J, 2)) + 1e-10 * eye(size(J, 2));
    S = inv(F); sd = sqrt(c * S * c');
    corr = (S * c') ./ sqrt(diag(S) * (c * S * c'));   % correlation of each knob with the pool
    cn = cond(J' * J);
end

function E = iEval(p, KN, th, tag, cdirC, cdirT, TB)
    fB = fullfile(cdirC, [tag '.mat']); fP = fullfile(cdirT, [tag '_proto.mat']);
    E = [];
    if isfile(fB); L = load(fB); if isequal(L.E.theta, th); E = L.E; end; end
    q = srxApplyTheta(p, KN, th);
    if isempty(E)
        if isempty(gcp('nocreate')); parpool(5); end
        t0 = tic; s = srxIsoState(q); q2 = q; q2.SRXT_0 = s.SR; q2.SRXD_0 = s.SRD;
        [fm, fd] = iRun(q2);
        E = struct('tag', tag, 'theta', th, 'fm', fm, 'fd', fd, 'srx', s.SRX, 'F', s.F, 'sec', toc(t0));
        save(fB, 'E');
    end
    if isfile(fP); L = load(fP); P = L.P; else; P = srxProtoFeatures(q, TB); save(fP, 'P'); end
    E.P = P;
    fprintf('    [%s] SRX %.4f  k63(5ms) %.1f  t50 %.1f ms\n', tag, E.srx, P.k63(1), 1e3*P.t50);
end

function [fm, fd] = iRun(p)
    % RunBakersExp is a SCRIPT; a function workspace keeps it from clobbering ours.
    params0 = p; %#ok<NASGU>
    LoadData; RunBakersExp; %#ok<*ASGLU>
    fm = features_model; fd = features_data;
end
