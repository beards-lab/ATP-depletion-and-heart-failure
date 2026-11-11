% RunVerify.m — do not trust the numbers in RunSawtoothCardiac until these pass.
%
% Three things could make the headline results artefacts rather than physics:
%
%   [1] REPEATABILITY. The 0.4 ML/s protocol is run twice by RunSawtoothCardiac
%       (once as case A, once as a sweep point) and reported peak_t of 19.52 ms
%       and 19.0 ms. ode15s is deterministic, so identical inputs must give
%       identical output; a difference means the two calls are NOT identical and
%       something is leaking between runs.
%
%   [2] OUTPUT SAMPLING. The SRS regression windows held only ~21-23 raw solver
%       points. ode15s places steps where the solution is stiff, not where the
%       measurement window is, so the slope may be reading the interpolant.
%       Fix: subdivide the ramp into sub-segments of the SAME velocity. That
%       forces solver output at each boundary without changing the physics
%       (each restart continues from PU(end,:)), giving a convergence test.
%
%   [3] WARM-START CONVERGENCE. Every force here is normalised by F0 at t = 0.
%       If the 2 s isometric hold has not settled, every ratio is wrong.
%
% Writes results/verify.mat

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
resdir = fullfile(here, 'results');
C = load(fullfile(resdir, 'characterize.mat'));

params0 = getParams(loadParams('params/estiff125_opt.m'), [], true, false);
params0.ghostLoad = ''; params0.ghostSave = ''; params0.SaveBest = false;
params0.PlotEachSeparately = false; params0.PlotProbsOnFig = false;
params0.ValuesInTime = true; params0.MaxRunTime = 320;

L0 = 1.00; FRAC = 0.04; T_WARM = 2.00; T_POST = 1.00; SRS_EPS = 0.0045;
V = [FRAC/0.100, 0.35*C.Vmax];              % verbatim, kinetically matched
tags = {'verbatim', 'kinmatched'};

%% [1] Repeatability -------------------------------------------------------
fprintf('== [1] repeatability, v = %.3f ML/s ==\n', V(1));
[vt, meta] = buildSawtoothProtocol(L0, FRAC, V(1), T_WARM, T_POST, params0.ML);
r = struct();
for k = 1:2
    out = iRun(params0, vt(:,1), vt(:,2));
    a   = extractSawtoothAttributes(out, meta, SRS_EPS);
    r(k).n = numel(out.t); r(k).F0 = a.F0_kPa; r(k).peak_t = a.peak_t;
    r(k).SRS_t10 = a.SRS_t10; r(k).SRS_eps = a.SRS_eps; r(k).t = out.t;
    fprintf('  run %d: n_t = %5d  F0 = %.6f kPa  peak_t = %.4f ms  SRS10 = %.4f\n', ...
            k, r(k).n, r(k).F0, 1000*r(k).peak_t, r(k).SRS_t10);
end
same_t = isequal(size(r(1).t), size(r(2).t)) && max(abs(r(1).t - r(2).t)) == 0;
fprintf('  identical time vectors: %d   |dF0| = %.3g   |dpeak_t| = %.3g ms\n', ...
        same_t, abs(r(1).F0 - r(2).F0), 1000*abs(r(1).peak_t - r(2).peak_t));

%% [2] Output-sampling convergence ----------------------------------------
% Densify ONLY the two ramp segments; the long isometric holds do not need it.
nsubs = [1 4 16 64];
conv = struct();
for i = 1:2
    [vt, meta] = buildSawtoothProtocol(L0, FRAC, V(i), T_WARM, T_POST, params0.ML);
    fprintf('\n== [2] sampling convergence, %s (v = %.3f ML/s, ramp %.2f ms) ==\n', ...
            tags{i}, V(i), 1000*meta.dt);
    for q = 1:numel(nsubs)
        [Td, Vd] = iDensify(vt, [1 nsubs(q) nsubs(q) 1]);
        out = iRun(params0, Td, Vd);
        a   = extractSawtoothAttributes(out, meta, SRS_EPS);
        conv(i,q).nsub = nsubs(q); conv(i,q).a = a; conv(i,q).nt = numel(out.t);
        fprintf('  nsub=%3d: n_t=%5d  pts(10ms)=%4d pts(eps)=%4d | SRS10=%8.3f  SRSeps=%8.3f  peak=%.4f at %.3f ms\n', ...
                nsubs(q), numel(out.t), a.n_srs_t10, a.n_srs_eps, ...
                a.SRS_t10, a.SRS_eps, a.peak_Fn, 1000*a.peak_t);
    end
    b = conv(i,end).a;
    fprintf('  --> vs finest: dSRS10 = %+.2f%%  dSRSeps = %+.2f%%  dpeak = %+.3f%%\n', ...
            100*(conv(i,1).a.SRS_t10 - b.SRS_t10)/b.SRS_t10, ...
            100*(conv(i,1).a.SRS_eps - b.SRS_eps)/b.SRS_eps, ...
            100*(conv(i,1).a.peak_Fn - b.peak_Fn)/b.peak_Fn);
end

%% [3] Warm-start convergence ---------------------------------------------
fprintf('\n== [3] warm-start convergence (2 s isometric hold) ==\n');
[vt, meta] = buildSawtoothProtocol(L0, FRAC, V(1), T_WARM, T_POST, params0.ML);
out = iRun(params0, vt(:,1), vt(:,2));
tp = [-1.5 -1.0 -0.5 -0.2 -0.05 0];
Fp = interp1(out.t, out.Force, tp, 'linear', 'extrap');
for k = 1:numel(tp)
    fprintf('  F(t=%+6.2f s) = %.4f kPa   (%.4f%% from t=0)\n', ...
            tp(k), Fp(k), 100*(Fp(k) - Fp(end))/Fp(end));
end
attp = interp1(out.t, out.p1_0 + out.p2_0, tp, 'linear', 'extrap');
SRp  = interp1(out.t, out.SR, tp, 'linear', 'extrap');
fprintf('  attached: %.5f -> %.5f   SRX: %.5f -> %.5f  (t=-1.5 s -> 0)\n', ...
        attp(1), attp(end), SRp(1), SRp(end));
warm_drift = abs(Fp(3) - Fp(end)) / Fp(end);           % last 0.5 s of the hold
fprintf('  drift over the final 0.5 s of the hold: %.4f%%\n', 100*warm_drift);

save(fullfile(resdir, 'verify.mat'), 'r', 'conv', 'nsubs', 'tp', 'Fp', ...
     'attp', 'SRp', 'same_t', 'warm_drift', 'V', 'tags');
fprintf('\nSaved results/verify.mat\n');


% =========================================================================
function out = iRun(params0, T, Vel)
    p = params0;
    p.SL0 = 2.0; p.Slim_l = 1.4; p.Slim_r = 2.3;
    if isfield(p, 'PU0'); p = rmfield(p, 'PU0'); end
    p = getParams(p, p.g, true);
    p.Velocity = Vel(:);
    lastwarn('');
    [~, out] = evaluateModel(str2func(p.modelFcn), T(:), p);
    if ~isempty(lastwarn); fprintf('    [warn] %s\n', lastwarn); end
end

function [Td, Vd] = iDensify(vt, nsub)
%IDENSIFY Split each velocity segment into nsub(k) sub-segments of equal
%   duration and IDENTICAL velocity. Physics unchanged; solver output denser.
    Td = []; Vd = [];
    for k = 1:size(vt, 1) - 1
        tt = linspace(vt(k,1), vt(k+1,1), nsub(k) + 1);
        Td = [Td, tt(1:end-1)];               %#ok<AGROW>
        Vd = [Vd, repmat(vt(k,2), 1, nsub(k))]; %#ok<AGROW>
    end
    Td(end+1) = vt(end,1); Vd(end+1) = vt(end,2);
end
