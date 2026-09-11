% RunSawtoothCardiac.m
% =========================================================================
% Port the skeletal short-range-stiffness sawtooth to cardiac, and evaluate it.
%
% SOURCE PROTOCOL  van der Zee TJ, Simha SN, Milburn GN, Campbell KS, Ting LH,
%   De Groote F. "Cross-bridge model for predicting muscle short-range stiffness
%   during movement." PLOS Comput Biol 2026, 10.1371/journal.pcbi.1014748.
%   Permeabilized rat soleus, L0 <-> SL 2.60 um, isokinetic stretch of ~3.83% L0
%   at 0.45 L0/s (stated as ~35% of unloaded shortening velocity, so Vmax ~1.3
%   L0/s), released isokinetically at the same speed, 22 C. Outcome measure:
%   short-range stiffness = regression slope of normalised force on length over
%   the first 0.010 s after stretch onset.
%
% CARDIAC PORT     Idealised to the round numbers: +4% L0 in 100 ms, released at
%   the same speed. Reference length is SL 2.0 um total (L0 = 1.0 ML, since
%   params0.ML = 2), which is also the length every other runner in this repo
%   hardcodes for its isometric reference.
%
%   Case A  UNMATCHED / verbatim input dynamics - the published numbers applied
%           to cardiac as-is: 4% L0, 100 ms, v = 0.400 ML/s.
%   Case B  KINETICALLY MATCHED - same amplitude, but the ramp speed set to the
%           same fraction of Vmax the paper used (0.35), using this model's own
%           Vmax from RunCharacterize.m.
%   Sweep   Same amplitude across a range of v, so the transition between the
%           two regimes is visible rather than inferred from two points.
%
% 8 mM ATP only, on params/estiff125_opt.m (the current best high-ATP set).
%
% Writes results/sawtooth.mat, results/*.png, and the protocol files
% data/protocol_sawtooth_cardiac_{verbatim,kinmatched}.mat
%
% Run: matlab -batch "run('Analyses/SawtoothSRS/RunSawtoothCardiac.m')"
% =========================================================================

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
resdir = fullfile(here, 'results');

C = load(fullfile(resdir, 'characterize.mat'));   % Vmax, ktr, F0, derived rates

params0 = getParams(loadParams('params/estiff125_opt.m'), [], true, false);
params0.ghostLoad = ''; params0.ghostSave = ''; params0.SaveBest = false;
params0.PlotEachSeparately = false; params0.PlotProbsOnFig = false;
params0.ValuesInTime  = true;
params0.MaxRunTime    = 320;

% Protocol geometry, shared by every case.
L0      = 1.00;     % ML  -> SL 2.0 um total
FRAC    = 0.04;     % +4% L0 -> +0.08 um SL -> 40 nm per half-sarcomere
T_WARM  = 2.00;     % s isometric warm-start (>> 1/ktr = 25 ms)
T_POST  = 1.00;     % s recovery hold

% The strain the paper's own 10 ms window covers, used as the strain-matched
% SRS window so cardiac and skeletal slopes are taken over the same excursion:
%   0.45 L0/s * 0.010 s = 0.0045 L0.
SRS_EPS = 0.0045;

fprintf('Vmax = %.3f ML/s, ktr = %.2f 1/s (tau = %.1f ms)\n', C.Vmax, C.ktr, 1000/C.ktr);

%% [1] The two headline cases ---------------------------------------------
% Named once and reused in the sweep below. Writing the literal 0.4 there
% instead differs from FRAC/0.100 by one ulp, which is enough to shift the
% solver's step placement - harmless, but it makes the two runs of the "same"
% protocol print different peak times.
V_VERB = FRAC / 0.100;

cases = struct('name', {}, 'tag', {}, 'v', {});
cases(1) = struct('name', 'A unmatched (verbatim 4%/100ms)', 'tag', 'verbatim',   'v', V_VERB);
cases(2) = struct('name', 'B kinetically matched (35% Vmax)', 'tag', 'kinmatched', 'v', 0.35*C.Vmax);

for i = 1:numel(cases)
    [vt, meta] = buildSawtoothProtocol(L0, FRAC, cases(i).v, T_WARM, T_POST, params0.ML);
    [out, p]   = iRunSawtooth(params0, vt);
    a          = extractSawtoothAttributes(out, meta, SRS_EPS);

    cases(i).vt = vt; cases(i).meta = meta; cases(i).out = out; cases(i).a = a;
    cases(i).vrel = cases(i).v / C.Vmax;

    iSaveProtocol(vt, meta, out, a, ['protocol_sawtooth_cardiac_' cases(i).tag '.mat']);

    fprintf('\n=== %s ===\n', cases(i).name);
    fprintf('  v = %.4f ML/s  (%.1f%% Vmax)   ramp = %.2f ms   amplitude = %.3f um SL (%.0f nm/HS)\n', ...
            meta.v, 100*cases(i).vrel, 1000*meta.dt, meta.dSL, meta.nm_per_HS);
    fprintf('  isometric F0        = %.2f kPa\n', a.F0_kPa);
    fprintf('  SRS (10 ms window)  = %.3f F0/L0   [window = %.2f ms, %.1f nm/HS, %.2f ktr tau]\n', ...
            a.SRS_t10, 1000*a.srs_t10_dur, a.srs_t10_nmHS, a.srs_t10_dur*C.ktr);
    fprintf('  SRS (%.2f%% L0 window)= %.3f F0/L0   [strain-matched to skeletal]\n', ...
            100*a.srs_eps_strain, a.SRS_eps);
    fprintf('  peak F/F0           = %.3f at %.2f ms (%.1f nm/HS)   yield before apex: %d\n', ...
            a.peak_Fn, 1000*a.peak_t, a.peak_nmHS, a.yields);
    fprintf('  F/F0 at apex        = %.3f    min during release = %.3f (slack: %d)\n', ...
            a.Fn_apex, a.Fn_minRel, a.slackedRel);
    fprintf('  F/F0 back at L0     = %.3f    after %.2g s hold = %.3f  (RFE = %+.3f)\n', ...
            a.Fn_return, T_POST, a.Fn_end, a.RFE);
    fprintf('  post-release OVERSHOOT = %.3f F0 at %.1f ms after return to L0\n', ...
            a.Fn_postMax, 1000*a.Fn_postMax_t);
    fprintf('  work in/out         = %.3f / %.3f kPa*ML   hysteresis = %.3f (elastic return %.0f%%)\n', ...
            a.W_stretch, a.W_release, a.W_hyst, 100*a.elasticReturn);
    fprintf('  attached fraction   = %.3f iso -> %.3f apex -> %.3f at L0 (min %.3f)\n', ...
            a.att_t0, a.att_t1, a.att_t2, a.att_min);
    fprintf('  passive share of apex force = %.1f%%\n', 100*a.pasFrac_apex);
    fprintf('  solver points: %d in ramp, %d in 10 ms window, %d in strain window\n', ...
            a.n_ramp, a.n_srs_t10, a.n_srs_eps);
end

%% [2] Velocity sweep at fixed amplitude ----------------------------------
% Fixed 4% L0, ramp speed spanning both headline cases and beyond, so the
% regime change is a curve rather than two isolated claims.
vsw = [0.1 0.2 V_VERB 0.8 1.5 0.35*C.Vmax 4.0 6.0];
S = struct('v', num2cell(vsw));
for j = 1:numel(vsw)
    [vt, meta] = buildSawtoothProtocol(L0, FRAC, vsw(j), T_WARM, T_POST, params0.ML);
    [out, ~]   = iRunSawtooth(params0, vt);
    a          = extractSawtoothAttributes(out, meta, SRS_EPS);
    S(j).meta = meta; S(j).a = a; S(j).vrel = vsw(j)/C.Vmax;
    S(j).t = out.t; S(j).Fn = out.Force / a.F0_kPa;
    fprintf('sweep v=%5.2f ML/s (%4.1f%% Vmax): SRS10=%7.2f  SRSeps=%7.2f  peak=%.3f at %.1f nm/HS  hyst=%+.3f  postMax=%.3f\n', ...
            vsw(j), 100*S(j).vrel, a.SRS_t10, a.SRS_eps, a.peak_Fn, a.peak_nmHS, ...
            a.W_hyst, a.Fn_postMax);
end

%% [3] Figures -------------------------------------------------------------
iFigTraces(cases, C, resdir);
iFigWorkLoops(cases, S, resdir);
iFigSweep(S, cases, C, resdir);

%% [4] Save ----------------------------------------------------------------
save(fullfile(resdir, 'sawtooth.mat'), 'cases', 'S', 'C', 'L0', 'FRAC', ...
     'T_WARM', 'T_POST', 'SRS_EPS', 'vsw', '-v7.3');
fprintf('\nSaved results/sawtooth.mat\n');


% =========================================================================
function [out, p] = iRunSawtooth(params0, vt)
%IRUNSAWTOOTH One sawtooth simulation on a grid wide enough for the excursion.
    p = params0;
    p.SL0 = 2 * vt(1, 4);        % um; ML = 2 so SL = 2 * L
    % Same grid the file-driven stairs runner uses. It was sized for stretches
    % up to SL 2.2 um, so it comfortably spans this 2.0 -> 2.08 excursion, and
    % UseSpaceExtension catches anything the attached distribution does beyond.
    p.Slim_l = 1.4;
    p.Slim_r = 2.3;
    if isfield(p, 'PU0'); p = rmfield(p, 'PU0'); end
    p = getParams(p, p.g, true);
    p.Velocity = vt(:, 2);       % after getParams: it does not set Velocity

    lastwarn('');
    [~, out] = evaluateModel(str2func(p.modelFcn), vt(:, 1), p);
    w = lastwarn;
    if ~isempty(w)
        fprintf('  [warn] %s\n', w);
    end
end

function iSaveProtocol(vt, meta, out, a, fname)
%ISAVEPROTOCOL Write the protocol file the standard runners consume.
%   datatable = [t(s), SL(um), F(normalised by isometric force at stretch
%   onset)] - the same normalisation runStairsExperiment applies, so this file
%   can be replayed through that machinery unchanged. Force here is the model's
%   own output (there is no cardiac sawtooth measurement yet), which makes the
%   file a round-trip self-test rather than data.
    Lof = @(tt) meta.L0 + meta.v * max(0, min(tt, meta.t1) - meta.t0) ...
                        - meta.v * max(0, min(tt, meta.t2) - meta.t1);
    tq = (vt(1,1) : 1e-4 : vt(end,1))';        % 10 kHz: resolves a 15 ms ramp
    datatable = [ tq, ...
                  2 * Lof(tq), ...
                  interp1(out.t, out.Force / a.F0_kPa, tq, 'linear', 'extrap') ];
    velocitytable = vt;
    features_data = struct();
    protocolMeta  = meta;
    if ~exist('data', 'dir'); mkdir('data'); end
    save(fullfile('data', fname), 'velocitytable', 'datatable', ...
         'features_data', 'protocolMeta');
    fprintf('  wrote data/%s (%d samples)\n', fname, size(datatable, 1));
end

function iFigTraces(cases, C, resdir)
    f = figure('Visible', 'off', 'Position', [50 50 1180 760]);
    tiledlayout(2, 2, 'Padding', 'compact');
    for i = 1:numel(cases)
        m = cases(i).meta; a = cases(i).a; o = cases(i).out;
        Lof = @(tt) m.L0 + m.v*max(0, min(tt, m.t1) - m.t0) - m.v*max(0, min(tt, m.t2) - m.t1);

        nexttile; hold on; box on;
        yyaxis right
        tq = linspace(-0.02, m.t2 + 0.2, 4000)';
        plot(1000*tq, 2*Lof(tq), '-', 'LineWidth', 1.2);
        ylabel('SL (\mum)');
        yyaxis left
        plot(1000*o.t, o.Force / a.F0_kPa, '-', 'LineWidth', 1.8);
        yline(1, ':k');
        xlim(1000*[-0.02, m.t2 + 0.2]);
        xlabel('t (ms)'); ylabel('F / F_0');
        title(sprintf('%s\nv=%.3f ML/s (%.0f%% V_{max}), ramp %.1f ms', ...
                      cases(i).name, m.v, 100*cases(i).vrel, 1000*m.dt), 'FontSize', 9);

        nexttile; hold on; box on;
        % Zoom on the SRS window with the fitted slope drawn on the F-L plane.
        tz = linspace(m.t0, min(m.t0 + max(0.010, m.dt), m.t1), 500)';
        Lz = Lof(tz);
        plot(100*(Lz - m.L0)/m.L0, interp1(o.t, o.Force/a.F0_kPa, tz), 'k-', ...
             'LineWidth', 1.8, 'DisplayName', 'model F-L');
        for which = 1:2
            if which == 1
                dur = a.srs_t10_dur; b = a.SRS_t10; c = [0.85 0.33 0.10]; nm = '10 ms';
            else
                dur = min(a.srs_eps_strain*m.L0/m.v, m.dt); b = a.SRS_eps; c = [0 0.45 0.74]; nm = 'strain-matched';
            end
            tw = [m.t0, m.t0 + dur];
            Lw = Lof(tw);
            Fa = interp1(o.t, o.Force/a.F0_kPa, m.t0);
            plot(100*(Lw - m.L0)/m.L0, Fa + b*(Lw - m.L0), '--', 'Color', c, ...
                 'LineWidth', 1.6, 'DisplayName', sprintf('%s: %.1f F_0/L_0', nm, b));
        end
        legend('Location', 'southeast', 'FontSize', 8);
        xlabel('strain (% L_0)'); ylabel('F / F_0');
        title('short-range stiffness fits', 'FontSize', 9);
    end
    sgtitle(sprintf('Cardiac sawtooth: +4%% L_0 from SL 2.0 \\mum, released at the same speed   (k_{tr} = %.0f s^{-1}, V_{max} = %.1f ML/s)', C.ktr, C.Vmax), 'FontSize', 10);
    exportgraphics_compat(f, fullfile(resdir, 'traces.png'));
end

function iFigWorkLoops(cases, S, resdir)
    f = figure('Visible', 'off', 'Position', [50 50 1080 420]);
    tiledlayout(1, 2, 'Padding', 'compact');

    nexttile; hold on; box on;
    for i = 1:numel(cases)
        m = cases(i).meta; a = cases(i).a; o = cases(i).out;
        Lof = @(tt) m.L0 + m.v*max(0, min(tt, m.t1) - m.t0) - m.v*max(0, min(tt, m.t2) - m.t1);
        tq = linspace(m.t0, m.t2, 3000)';
        plot(2*Lof(tq), interp1(o.t, o.Force/a.F0_kPa, tq), '-', 'LineWidth', 1.8, ...
             'DisplayName', sprintf('%s (hyst %.3f)', cases(i).tag, a.W_hyst));
    end
    xlabel('SL (\mum)'); ylabel('F / F_0'); legend('Location', 'northwest', 'FontSize', 8);
    title('work loops, headline cases');

    nexttile; hold on; box on;
    cm = parula(numel(S) + 1);
    for j = 1:numel(S)
        m = S(j).meta;
        Lof = @(tt) m.L0 + m.v*max(0, min(tt, m.t1) - m.t0) - m.v*max(0, min(tt, m.t2) - m.t1);
        tq = linspace(m.t0, m.t2, 3000)';
        plot(2*Lof(tq), interp1(S(j).t, S(j).Fn, tq), '-', 'Color', cm(j,:), ...
             'LineWidth', 1.4, 'DisplayName', sprintf('%.2f ML/s', m.v));
    end
    xlabel('SL (\mum)'); ylabel('F / F_0'); legend('Location', 'northwest', 'FontSize', 7);
    title('work loops vs ramp speed');
    exportgraphics_compat(f, fullfile(resdir, 'workloops.png'));
end

function iFigSweep(S, cases, C, resdir)
    vrel = [S.vrel];
    A = [S.a];
    f = figure('Visible', 'off', 'Position', [50 50 1180 700]);
    tiledlayout(2, 2, 'Padding', 'compact');

    nexttile; hold on; box on;
    plot(100*vrel, [A.SRS_t10], 'o-', 'LineWidth', 1.6, 'DisplayName', 'SRS, 10 ms window');
    plot(100*vrel, [A.SRS_eps], 's-', 'LineWidth', 1.6, 'DisplayName', 'SRS, strain-matched');
    iMarkCases(cases, C);
    set(gca, 'XScale', 'log'); xlabel('ramp speed (% V_{max})'); ylabel('SRS (F_0 / L_0)');
    legend('Location', 'best', 'FontSize', 8); title('short-range stiffness');

    nexttile; hold on; box on;
    plot(100*vrel, [A.peak_Fn], 'o-', 'LineWidth', 1.6, 'DisplayName', 'peak F/F_0');
    plot(100*vrel, [A.Fn_apex], 's-', 'LineWidth', 1.6, 'DisplayName', 'F/F_0 at apex');
    plot(100*vrel, [A.Fn_end],  '^-', 'LineWidth', 1.6, 'DisplayName', 'F/F_0 after recovery');
    iMarkCases(cases, C);
    set(gca, 'XScale', 'log'); xlabel('ramp speed (% V_{max})'); ylabel('F / F_0');
    legend('Location', 'best', 'FontSize', 8); title('force landmarks');

    nexttile; hold on; box on;
    plot(100*vrel, [A.peak_nmHS], 'o-', 'LineWidth', 1.6);
    iMarkCases(cases, C);
    set(gca, 'XScale', 'log'); xlabel('ramp speed (% V_{max})'); ylabel('strain at peak (nm / HS)');
    title('yield point');

    nexttile; hold on; box on;
    yyaxis left
    plot(100*vrel, [A.W_hyst], 'o-', 'LineWidth', 1.6); ylabel('hysteresis (kPa\cdotML)');
    yyaxis right
    plot(100*vrel, [A.att_min], 's-', 'LineWidth', 1.6); ylabel('min attached fraction');
    set(gca, 'XScale', 'log'); xlabel('ramp speed (% V_{max})');
    title('dissipation and cross-bridge detachment');

    sgtitle('Fixed +4% L_0 amplitude, ramp speed swept', 'FontSize', 10);
    exportgraphics_compat(f, fullfile(resdir, 'sweep.png'));
end

function iMarkCases(cases, C)
    for i = 1:numel(cases)
        xline(100*cases(i).v/C.Vmax, ':', cases(i).tag, 'FontSize', 7, ...
              'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    end
    xline(35, '-.', 'skeletal 35%', 'FontSize', 7, 'HandleVisibility', 'off');
end

function exportgraphics_compat(f, fname)
    try
        exportgraphics(f, fname, 'Resolution', 150);
    catch
        print(f, fname, '-dpng', '-r150');
    end
    close(f);
    fprintf('  wrote %s\n', fname);
end
