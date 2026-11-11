% RunConditioningTest.m
% =========================================================================
% Stiffness slope at the onset of the FIRST and SECOND ramp-up.
%
% This is the measurement van der Zee et al. 2026 actually report: short-range
% stiffness on a conditioning stretch, then on a test stretch after a recovery
% hold. Their result is that the second is DEPRESSED (muscle thixotropy).
%
% Analyses/SawtoothSRS/conclusions.md §6 predicted, from the single-sawtooth
% recovery, that this model has no memory beyond ~350 ms and so should show an
% effect ONLY at short recovery intervals. That is what this script tests, so
% the recovery interval is swept rather than fixed -- densely between 20 and
% 100 ms, where the sign of the effect turns out to change.
%
% Both velocity cases from the main analysis are run:
%   verbatim    0.400 ML/s (100 ms ramp)  - the published input dynamics
%   kin-matched 2.625 ML/s (15.2 ms ramp) - the same v/Vmax the paper used
%
% Slopes are reported on two windows (see §4 of conclusions.md for why the
% literal 10 ms recipe stops being "short range" once the ramp speed changes):
%   SRS_eps  strain-matched to the skeletal window (0.45 % L0)   <- primary
%   SRS_t10  the literal first-10-ms regression                  <- secondary
%
% Convergence and the reality of the 40-100 ms depression are checked
% separately in RunCondTestVerify.m.
%
% Writes results/condtest.mat, results/condtest.png, results/condtest_table.csv
%
% Run: matlab -batch "run('Analyses/SawtoothSRS/RunConditioningTest.m')"
% =========================================================================

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
resdir = fullfile(here, 'results');
C = load(fullfile(resdir, 'characterize.mat'));      % Vmax, ktr

params0 = getParams(loadParams('params/estiff125_opt.m'), [], true, false);
params0.ghostLoad = ''; params0.ghostSave = ''; params0.SaveBest = false;
params0.PlotEachSeparately = false; params0.PlotProbsOnFig = false;
params0.ValuesInTime = true; params0.MaxRunTime = 320;

L0 = 1.00; FRAC = 0.04; T_WARM = 2.00; T_POST = 1.00; SRS_EPS = 0.0045;
V_VERB = FRAC / 0.100;
V      = [V_VERB, 0.35*C.Vmax];
tags   = {'verbatim', 'kinmatched'};

% Dense through 20-100 ms: that is where the effect crosses from enhancement to
% depression, and a coarse grid there reads as a single stray point.
tRecs = [0.005 0.010 0.020 0.025 0.030 0.035 0.040 0.045 0.050 0.055 ...
         0.060 0.070 0.085 0.100 0.150 0.200 0.350 0.500 1.000 2.000];

fprintf('Vmax = %.3f ML/s, ktr = %.2f 1/s (tau = %.1f ms)\n\n', C.Vmax, C.ktr, 1000/C.ktr);

rows = {};
R = struct();
for i = 1:numel(V)
    fprintf('=== %s: v = %.3f ML/s, ramp = %.2f ms, amplitude +4%% L0 (40 nm/HS) ===\n', ...
            tags{i}, V(i), 1000*FRAC*L0/V(i));
    fprintf('%8s | %8s %8s %7s | %8s %8s %7s | %7s %8s %8s\n', ...
            'tRec/ms', 'SRSeps1', 'SRSeps2', 'ratio', 'SRS10_1', 'SRS10_2', 'ratio', ...
            'F@on2', 'att@on2', 'ratio/F');
    for j = 1:numel(tRecs)
        [vt, meta] = buildConditioningTestProtocol(L0, FRAC, V(i), T_WARM, ...
                                                   tRecs(j), T_POST, params0.ML);
        out = iRun(params0, vt);

        t  = out.t(:); F = out.Force(:);
        F0 = interp1(t, F, meta.t0, 'linear', 'extrap');    % isometric, 1st onset
        Fn = F / F0;

        % Window durations, clipped to the ramp so a slope never straddles the apex.
        d_eps = min(SRS_EPS * meta.L0 / meta.v, meta.dt);
        d_t10 = min(0.010, meta.dt);

        s_eps = [iSlope(t, Fn, vt, meta.onsets(1), d_eps), ...
                 iSlope(t, Fn, vt, meta.onsets(2), d_eps)];
        s_t10 = [iSlope(t, Fn, vt, meta.onsets(1), d_t10), ...
                 iSlope(t, Fn, vt, meta.onsets(2), d_t10)];

        % Force the muscle actually starts the test stretch from. A changed
        % slope could just be a changed starting force rather than a genuine
        % stiffness change, so divide that out as well.
        Fon = [interp1(t, Fn, meta.onsets(1), 'linear', 'extrap'), ...
               interp1(t, Fn, meta.onsets(2), 'linear', 'extrap')];

        att = out.p1_0(:) + out.p2_0(:);
        aOn = [interp1(t, att, meta.onsets(1), 'linear', 'extrap'), ...
               interp1(t, att, meta.onsets(2), 'linear', 'extrap')];

        R(i,j).tag = tags{i}; R(i,j).v = V(i); R(i,j).tRec = tRecs(j);
        R(i,j).meta = meta; R(i,j).F0 = F0;
        R(i,j).SRS_eps = s_eps; R(i,j).SRS_t10 = s_t10;
        R(i,j).ratio_eps = s_eps(2)/s_eps(1);
        R(i,j).ratio_t10 = s_t10(2)/s_t10(1);
        R(i,j).Fon = Fon; R(i,j).attOn = aOn;
        R(i,j).ratio_perF = (s_eps(2)/s_eps(1)) / (Fon(2)/Fon(1));
        R(i,j).t = t; R(i,j).Fn = Fn;

        fprintf('%8.0f | %8.2f %8.2f %7.3f | %8.2f %8.2f %7.3f | %7.3f %8.4f %8.3f\n', ...
                1000*tRecs(j), s_eps(1), s_eps(2), R(i,j).ratio_eps, ...
                s_t10(1), s_t10(2), R(i,j).ratio_t10, Fon(2), aOn(2), ...
                R(i,j).ratio_perF);

        rows(end+1, :) = {tags{i}, V(i), 1000*tRecs(j), s_eps(1), s_eps(2), ...
                          R(i,j).ratio_eps, s_t10(1), s_t10(2), R(i,j).ratio_t10, ...
                          Fon(2), aOn(2), R(i,j).ratio_perF};   %#ok<SAGROW>
    end
    fprintf('\n');
end

T = cell2table(rows, 'VariableNames', {'case','v_MLs','tRec_ms','SRSeps_ramp1', ...
     'SRSeps_ramp2','ratio_eps','SRSt10_ramp1','SRSt10_ramp2','ratio_t10', ...
     'F_at_onset2','att_at_onset2','ratio_per_onsetForce'});
writetable(T, fullfile(resdir, 'condtest_table.csv'));
fprintf('Wrote results/condtest_table.csv\n');

%% Figure -----------------------------------------------------------------
f = figure('Visible', 'off', 'Position', [30 30 1320 780]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
co = [0 0.447 0.741; 0.85 0.325 0.098];
co3 = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125];
gx = @() set(gca, 'XScale', 'log');

% (1) absolute slopes
nexttile; hold on; box on;
for i = 1:numel(V)
    e1 = arrayfun(@(k) R(i,k).SRS_eps(1), 1:numel(tRecs));
    e2 = arrayfun(@(k) R(i,k).SRS_eps(2), 1:numel(tRecs));
    plot(1000*tRecs, e2, 'o-', 'Color', co(i,:), 'LineWidth', 1.7, ...
         'MarkerSize', 4, 'DisplayName', sprintf('%s, 2nd ramp', tags{i}));
    yline(e1(1), '--', 'Color', co(i,:), 'LineWidth', 1.2, 'HandleVisibility', 'off');
end
gx(); xlabel('recovery interval (ms)'); ylabel('SRS, strain-matched (F_0/L_0)');
legend('Location', 'best', 'FontSize', 7);
title('onset slope, 2nd ramp (dashed = 1st)', 'FontSize', 9);

% (2) ratio
nexttile; hold on; box on;
for i = 1:numel(V)
    rr = arrayfun(@(k) R(i,k).ratio_eps,  1:numel(tRecs));
    rf = arrayfun(@(k) R(i,k).ratio_perF, 1:numel(tRecs));
    plot(1000*tRecs, rr, 'o-', 'Color', co(i,:), 'LineWidth', 1.7, ...
         'MarkerSize', 4, 'DisplayName', sprintf('%s', tags{i}));
    plot(1000*tRecs, rf, 's:', 'Color', co(i,:), 'LineWidth', 1.2, ...
         'MarkerSize', 3, 'DisplayName', sprintf('%s, per onset F', tags{i}));
end
yline(1, ':k', 'HandleVisibility', 'off');
gx(); xlabel('recovery interval (ms)'); ylabel('SRS_2 / SRS_1');
legend('Location', 'best', 'FontSize', 7);
title('enhancement above 1, depression below', 'FontSize', 9);

% (3) attached fraction at the second onset
nexttile; hold on; box on;
for i = 1:numel(V)
    aa = arrayfun(@(k) R(i,k).attOn(2), 1:numel(tRecs));
    plot(1000*tRecs, aa, 'o-', 'Color', co(i,:), 'LineWidth', 1.7, ...
         'MarkerSize', 4, 'DisplayName', tags{i});
end
yline(R(1,1).attOn(1), '--k', 'isometric', 'FontSize', 7, 'HandleVisibility', 'off');
gx(); xlabel('recovery interval (ms)'); ylabel('attached fraction at 2nd onset');
legend('Location', 'best', 'FontSize', 7);
title('bridge count at the 2nd onset', 'FontSize', 9);

% (4) traces, verbatim
nexttile; hold on; box on;
pick = [find(tRecs == 0.005), find(tRecs == 0.070), find(tRecs == 0.500)];
% Longest recovery first: the CONDITIONING event is identical in all three, so
% whichever is drawn last hides the others over 0-200 ms. Drawing longest-first
% leaves the shortest on top, where the interesting second event is.
for k = numel(pick):-1:1
    r = R(1, pick(k));
    plot(1000*r.t, r.Fn, '-', 'Color', co3(k,:), 'LineWidth', 1.4, ...
         'DisplayName', sprintf('t_{rec} = %.0f ms', 1000*r.tRec));
end
yline(1, ':k', 'HandleVisibility', 'off');
xlabel('t (ms)'); ylabel('F / F_0'); xlim([-20 900]);
legend('Location', 'northeast', 'FontSize', 7);
title(sprintf('verbatim traces, v = %.2f ML/s (1st event common)', V(1)), 'FontSize', 9);

% (5)-(6) force-length of the two ramps, at enhancement and at the dip
for k = 1:2
    nexttile; hold on; box on;
    r = R(1, pick(k));
    for e = 1:2
        ton = r.meta.onsets(e);
        d   = min(SRS_EPS*r.meta.L0/r.meta.v, r.meta.dt);
        tz  = linspace(ton, ton + r.meta.dt, 400)';
        Lz  = iCommandedL(r.meta, tz);
        plot(100*(Lz - r.meta.L0)/r.meta.L0, interp1(r.t, r.Fn, tz), '-', ...
             'Color', co(e,:), 'LineWidth', 1.7, ...
             'DisplayName', sprintf('ramp %d: %.1f F_0/L_0', e, r.SRS_eps(e)));
        tw = [ton, ton + d]; Lw = iCommandedL(r.meta, tw);
        Fa = interp1(r.t, r.Fn, ton);
        plot(100*(Lw - r.meta.L0)/r.meta.L0, Fa + r.SRS_eps(e)*(Lw - r.meta.L0), ...
             '--k', 'LineWidth', 1.4, 'HandleVisibility', 'off');
    end
    xlabel('strain from L_0 (% L_0)'); ylabel('F / F_0');
    legend('Location', 'southeast', 'FontSize', 7);
    title(sprintf('t_{rec} = %.0f ms  (ratio %.2f)', 1000*r.tRec, r.ratio_eps), ...
          'FontSize', 9);
end

sgtitle(sprintf('Conditioning and test ramp-up: onset stiffness   (+4%% L_0 from SL 2.0 \\mum, k_{tr} = %.0f s^{-1}, \\tau = %.0f ms)', ...
        C.ktr, 1000/C.ktr), 'FontSize', 10);
try
    exportgraphics(f, fullfile(resdir, 'condtest.png'), 'Resolution', 150);
catch
    print(f, fullfile(resdir, 'condtest.png'), '-dpng', '-r150');
end
close(f);

save(fullfile(resdir, 'condtest.mat'), 'R', 'tRecs', 'V', 'tags', ...
     'FRAC', 'L0', 'SRS_EPS', 'T', '-v7.3');
fprintf('Saved results/condtest.mat and results/condtest.png\n');


% =========================================================================
function out = iRun(params0, vt)
    p = params0;
    p.SL0 = 2 * vt(1, 4);
    p.Slim_l = 1.4; p.Slim_r = 2.3;      % same grid as the stairs runner
    if isfield(p, 'PU0'); p = rmfield(p, 'PU0'); end
    p = getParams(p, p.g, true);
    p.Velocity = vt(:, 2);
    lastwarn('');
    [~, out] = evaluateModel(str2func(p.modelFcn), vt(:, 1), p);
    if ~isempty(lastwarn); fprintf('    [warn] %s\n', lastwarn); end
end

function b = iSlope(t, Fn, vt, ton, dur)
%ISLOPE Regression slope dFn/dL over [ton, ton+dur], uniform resample.
%   Length comes from integrating velocitytable column 2 - the trajectory the
%   model actually follows - not from column 4, which is bookkeeping.
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

function Lq = iCommandedL(meta, tq)
%ICOMMANDEDL Length (ML) for the two-sawtooth protocol, from meta alone.
    tq = tq(:);
    up1 = max(0, min(tq, meta.t1) - meta.t0);
    dn1 = max(0, min(tq, meta.t2) - meta.t1);
    up2 = max(0, min(tq, meta.t4) - meta.t3);
    dn2 = max(0, min(tq, meta.t5) - meta.t4);
    Lq  = meta.L0 + meta.v * (up1 - dn1 + up2 - dn2);
end
