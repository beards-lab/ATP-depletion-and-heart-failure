% RunCondTestVerify.m — is the 50 ms stiffness dip real?
%
% RunConditioningTest found the second ramp's onset stiffness is ENHANCED at
% short recovery (ratio 1.99 at 5 ms) but DEPRESSED at 50 ms (0.878) before
% returning to 1.000 by 350 ms. A non-monotonic minimum on a 10-point sweep is
% exactly the shape a single bad point makes, so:
%
%   [1] refine the recovery interval around the dip;
%   [2] re-run the dip point with 64x denser solver output;
%   [3] report the attached fraction at each onset, to see whether the dip has
%       a state-variable correlate or is purely a force-trace artefact.
%
% Writes results/condtest_verify.mat

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
resdir = fullfile(here, 'results');
C = load(fullfile(resdir, 'characterize.mat'));

params0 = getParams(loadParams('params/estiff125_opt.m'), [], true, false);
params0.ghostLoad = ''; params0.ghostSave = ''; params0.SaveBest = false;
params0.PlotEachSeparately = false; params0.PlotProbsOnFig = false;
params0.ValuesInTime = true; params0.MaxRunTime = 320;

L0 = 1.00; FRAC = 0.04; T_WARM = 2.00; T_POST = 1.00; SRS_EPS = 0.0045;
V_VERB = FRAC / 0.100;

%% [1] Refined sweep through the dip --------------------------------------
fine = [0.020 0.025 0.030 0.035 0.040 0.045 0.050 0.055 0.060 0.070 0.085 0.100];
fprintf('== [1] refined recovery sweep, verbatim (v = %.3f ML/s) ==\n', V_VERB);
fprintf('%8s | %8s %7s | %8s %8s | %8s %8s\n', ...
        'tRec/ms', 'SRSeps2', 'ratio', 'F@on2', 'ratio/F', 'att@on1', 'att@on2');
Fi = struct();
for j = 1:numel(fine)
    [vt, meta] = buildConditioningTestProtocol(L0, FRAC, V_VERB, T_WARM, ...
                                               fine(j), T_POST, params0.ML);
    [s, Fon, aOn, ~] = iMeasure(params0, vt, meta, SRS_EPS, ones(1, size(vt,1)-1));
    Fi(j).tRec = fine(j); Fi(j).s = s; Fi(j).Fon = Fon; Fi(j).aOn = aOn;
    fprintf('%8.0f | %8.2f %7.3f | %8.3f %8.3f | %8.4f %8.4f\n', ...
            1000*fine(j), s(2), s(2)/s(1), Fon(2), (s(2)/s(1))/(Fon(2)/Fon(1)), ...
            aOn(1), aOn(2));
end

%% [2] Sampling convergence at the dip ------------------------------------
fprintf('\n== [2] solver-output convergence at the dip (tRec = 50 ms) ==\n');
nsubs = [1 4 16 64];
Cv = struct();
for q = 1:numel(nsubs)
    [vt, meta] = buildConditioningTestProtocol(L0, FRAC, V_VERB, T_WARM, ...
                                               0.050, T_POST, params0.ML);
    [s, Fon, aOn, npts] = iMeasure(params0, vt, meta, SRS_EPS, ...
                                   [1 nsubs(q) nsubs(q) 1 nsubs(q) nsubs(q) 1]);
    Cv(q).nsub = nsubs(q); Cv(q).s = s; Cv(q).ratio = s(2)/s(1);
    fprintf('  nsub=%3d: pts in win = [%d %d]  SRS1 = %8.3f  SRS2 = %8.3f  ratio = %.4f\n', ...
            nsubs(q), npts(1), npts(2), s(1), s(2), s(2)/s(1));
end
d = 100*(Cv(1).ratio - Cv(end).ratio)/Cv(end).ratio;
fprintf('  --> ratio changes %+.3f%% from coarsest to 64x denser output\n', d);

%% [3] Same check at an enhancement point ---------------------------------
fprintf('\n== [3] convergence at an enhancement point (tRec = 5 ms) ==\n');
for q = [1 numel(nsubs)]
    [vt, meta] = buildConditioningTestProtocol(L0, FRAC, V_VERB, T_WARM, ...
                                               0.005, T_POST, params0.ML);
    [s, ~, ~, npts] = iMeasure(params0, vt, meta, SRS_EPS, ...
                               [1 nsubs(q) nsubs(q) 1 nsubs(q) nsubs(q) 1]);
    fprintf('  nsub=%3d: pts in win = [%d %d]  ratio = %.4f\n', ...
            nsubs(q), npts(1), npts(2), s(2)/s(1));
end

save(fullfile(resdir, 'condtest_verify.mat'), 'Fi', 'Cv', 'fine', 'nsubs');
fprintf('\nSaved results/condtest_verify.mat\n');


% =========================================================================
function [s, Fon, aOn, npts] = iMeasure(params0, vt, meta, srsEps, nsub)
    [Td, Vd] = iDensify(vt, nsub);
    p = params0;
    p.SL0 = 2 * vt(1, 4); p.Slim_l = 1.4; p.Slim_r = 2.3;
    if isfield(p, 'PU0'); p = rmfield(p, 'PU0'); end
    p = getParams(p, p.g, true);
    p.Velocity = Vd(:);
    lastwarn('');
    [~, out] = evaluateModel(str2func(p.modelFcn), Td(:), p);
    if ~isempty(lastwarn); fprintf('    [warn] %s\n', lastwarn); end

    t = out.t(:); F = out.Force(:);
    F0 = interp1(t, F, meta.t0, 'linear', 'extrap');
    Fn = F / F0;
    dur = min(srsEps * meta.L0 / meta.v, meta.dt);
    att = out.p1_0(:) + out.p2_0(:);
    s = zeros(1,2); Fon = zeros(1,2); aOn = zeros(1,2); npts = zeros(1,2);
    for e = 1:2
        ton = meta.onsets(e);
        s(e)    = iSlope(t, Fn, vt, ton, dur);
        Fon(e)  = interp1(t, Fn,  ton, 'linear', 'extrap');
        aOn(e)  = interp1(t, att, ton, 'linear', 'extrap');
        npts(e) = nnz(t >= ton & t <= ton + dur);
    end
end

function b = iSlope(t, Fn, vt, ton, dur)
    tq = linspace(ton, ton + dur, 201)';
    Fq = interp1(t, Fn, tq, 'linear', 'extrap');
    Lq = iIntegrateVT(vt, tq);
    p  = polyfit(Lq, Fq, 1);
    b  = p(1);
end

function Lq = iIntegrateVT(vt, tq)
    tk = vt(:,1); vk = vt(:,2);
    Lk = zeros(numel(tk),1); Lk(1) = vt(1,4);
    for k = 2:numel(tk)
        Lk(k) = Lk(k-1) + vk(k-1)*(tk(k) - tk(k-1));
    end
    Lq = interp1(tk, Lk, tq(:), 'linear', 'extrap');
end

function [Td, Vd] = iDensify(vt, nsub)
    Td = []; Vd = [];
    for k = 1:size(vt, 1) - 1
        tt = linspace(vt(k,1), vt(k+1,1), nsub(k) + 1);
        Td = [Td, tt(1:end-1)];                  %#ok<AGROW>
        Vd = [Vd, repmat(vt(k,2), 1, nsub(k))];  %#ok<AGROW>
    end
    Td(end+1) = vt(end,1); Vd(end+1) = vt(end,2);
end
