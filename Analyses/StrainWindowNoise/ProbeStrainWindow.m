%% ProbeStrainWindow.m -- what do strain-window shifts delete, and why is it not empty?
%
%   P1  per-head transition rates vs strain, evaluated exactly as the RHS does
%       (fast-path pchip coefficients, including extrapolation past the knots)
%   P2  one serial slack run with the shift logger on: mass and active force
%       carried by every block of bins a shift zeroes
%
% Set PARAMSET (default 'refit8_opt') and RUN (default [1 2]) before calling.
if ~exist('RUN', 'var'), RUN = [1 2]; end
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit8_opt'; end
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
addpath(genpath(root));
res  = fullfile(here, 'results');

P = getParams(loadParams(PARAMSET), [], true, false);
P.RunForceVelocity = false; P.RunSlack = true; P.RunStairs = false; P.RunKtr = false;
P.RunForceVelocityTime = false; P.RunSlackPassive = true; P.EvalFeatures = true;
P.BreakOnODEUnstable = false; P.PlotEachSeparately = 0; P.PlotFeatureFitting = 0;
P.RunSlackSegments = 'All'; P.MaxRunTime = 600; P.OptimizeOn = 'Feats';
P.ghostSave = ''; P.ghostLoad = ''; P.SaveBest = 0; P.mods = {}; P.g = ones(1,0);

%% P1 -- rates vs strain ------------------------------------------------------
if any(RUN == 1)
    s = (-0.15:0.0005:0.15)';
    ev = @(brk, cof, x) evalPP(brk, cof, x);
    C = P.cached;
    R.s    = s;
    R.R1D  = P.kd  * ev(C.pwsdR1D_brk, C.pwsdR1D_cof, s);
    R.R12  = P.k1  * ev(C.pwsd_brk,    C.pwsd_cof,    s);
    R.R21  = P.k_1 * ev(C.pwsdR21_brk, C.pwsdR21_cof, s);
    R.R2   = P.k2  * ev(C.pwsd2_brk,   C.pwsd2_cof,   s + P.dr2 - P.dr);
    R.rip  = P.k2rip * min(1e9, max(0, exp(P.alphaRip*(s + P.dr3))));
    R.knots = struct('R1D', C.pwsdR1D_brk([1 end])', 'R12', C.pwsd_brk([1 end])', ...
                     'R21', C.pwsdR21_brk([1 end])', 'R2', C.pwsd2_brk([1 end])' - (P.dr2 - P.dr));
    R.grid = [min(P.s) max(P.s)];
    fprintf('grid s range %.4f..%.4f um (%d bins of %.4f)\n', R.grid, P.ss, P.dS);
    nm = {'R1D','R12','R21','R2','rip'};
    at = [-0.15 -0.10 -0.06 -0.03 -0.015 0 0.015 0.03 0.06 0.10 0.15];
    fprintf('%-6s', 'rate'); fprintf('%9.0f', 1e3*at); fprintf('   <- strain, nm\n');
    for k = 1:numel(nm)
        v = interp1(s, R.(nm{k}), at);
        fprintf('%-6s', nm{k}); fprintf('%9.3g', v); fprintf('\n');
    end
    fprintf('knot ranges (strain, um): R1D %s  R12 %s  R21 %s  R2 %s\n', mat2str(R.knots.R1D,3), ...
        mat2str(R.knots.R12,3), mat2str(R.knots.R21,3), mat2str(R.knots.R2,3));
    save(fullfile(res, ['rates_' PARAMSET '.mat']), 'R');
end

%% P2 -- what the shifts delete in one slack run --------------------------------
if any(RUN == 2)
    global WINDOW_SHIFT_LOG %#ok<GVMIS>
    WINDOW_SHIFT_LOG = [];
    Pl = P; Pl.LogWindowShifts = true;
    tic; [~, out, fm] = runSlackExperiment(Pl); wall = toc;
    L = WINDOW_SHIFT_LOG;
    fprintf('\n%d shifts in one slack run (%.0f s)\n', numel(L), wall);
    fprintf('%8s %3s %8s %9s %9s %9s %18s %18s\n', 't (s)', 'ie', 'bound', 'del mass', 'del frac', 'dF (kPa)', 'deleted s (nm)', 'occupied s (nm)');
    for k = 1:numel(L)
        fprintf('%8.4f %3d %8.4f %9.2e %8.1f%% %9.3f %8.0f..%-8.0f %8.0f..%-8.0f\n', L(k).t, L(k).ie, L(k).bound, L(k).mass, ...
            100*L(k).massFrac, L(k).dFact, 1e3*L(k).sDel, 1e3*L(k).sOccupied);
    end
    save(fullfile(res, ['shiftlog_' PARAMSET '.mat']), 'L', 'fm', 'wall');
end

function y = evalPP(brk, cof, x)
% Same segment rule as the RHS fast path: interior breaks pick the segment,
% the end segments are extrapolated as cubics.
b = brk(:);
si = max(1, min(size(cof,1), sum(x >= b(1:end-1)', 2)));
dx = x - b(si);
y = ((cof(si,1).*dx + cof(si,2)).*dx + cof(si,3)).*dx + cof(si,4);
end
