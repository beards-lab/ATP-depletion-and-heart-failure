% CreateFVRampProtocolVelocityTable.m
% Generates the FV2 protocol files from the raw <cond>_slack.txt bursts.
%
% WHAT FV2 IS
% The 04-series slack recordings open with an event that is NOT a slack. It is
% a four-phase load-clamp manoeuvre (times from the 04/10 8mM file):
%
%   1. fast step release   50-52 ms   1.1004 -> ~1.056 ML  (~-13 ML/s)
%   2. ISOVELOCITY RAMP    52-86 ms   1.056  -> ~0.995 ML  (-2.01 ML/s, flat)
%   3. slack release       86-92 ms   0.995  -> 0.799 ML   (up to -58 ML/s)
%   4. hold, then restretch 110-131 ms 0.799 -> 1.1004 ML  (+15.5 ML/s)
%
% Only after this do the five true slack releases (-27 to -52 ML/s) begin.
%
% Phase 2 is the measurement: the step in phase 1 drops force onto the load the
% fibre can bear at 2 ML/s, and force then holds a clean PLATEAU (~16 kPa, i.e.
% 0.31*Fss, at 8mM) for the whole 34 ms ramp. That is a genuine sub-Vmax point
% on the force-velocity curve at a single velocity -- hence FV2, the
% force-velocity ramp at 2 ML/s (-4.03 um/s in SL). Phase 3 then unloads the
% fibre completely so phase 4 measures redevelopment from zero.
%
% The commanded velocity is identical to three decimals across 2mM, 8mM and
% PNB+Mava, so the plateau forces are directly comparable between conditions.
%
% It is carved out here as its own experiment so the slack protocol files keep
% the 22-row / 5-slack layout that runSlackExperiment indexes literally.
%
% CONVENTIONS (match CreateKtrProtocolVelocityTable / CreateStairs...):
%   - datatable force is NORMALISED by the pre-ramp isometric Fss, so the file
%     is drop-in for the existing runKtrExperiment / runStairsExperiment
%     scoring. The absolute Fss (kPa) is saved alongside, so absolute force is
%     recovered as datatable(:,3) * Fss.
%   - datatable length is in um (SL ~ L x 2).
%
% Produces:
%   data/protocol_04_10_2026_velocitytable_FV2.mat          (8mM ATP)
%   data/protocol_04_10_2026_velocitytable_FV2_2mM.mat      (2mM ATP)
%   data/protocol_04_10_2026_velocitytable_FV2_PNBMava.mat  (8mM + PNB/Mava)

configs = {
    '../data/04 10 2026 Male 2-8/8mM_slack.txt',          '../data/protocol_04_10_2026_velocitytable_FV2.mat';
    '../data/04 10 2026 Male 2-8/2mM_slack.txt',          '../data/protocol_04_10_2026_velocitytable_FV2_2mM.mat';
    '../data/04 10 2026 Male 2-8/8mM_slack_PNB_Mava.txt', '../data/protocol_04_10_2026_velocitytable_FV2_PNBMava.mat';
};

% Slack .txt Time column is in MILLISECONDS; divide by 1000 for a seconds axis.
fvSegment  = [0, 0.45];   % clip window (s): ends before the first true slack (0.500)
preWin     = [0, 0.045];  % pre-ramp isometric window, for L_base and Fss
warmupLead = 5;           % isometric warm-up before the ramp (s), so the model
                          % reaches steady state before the length change
thresh     = 0.005;       % 0.5% ML transition-detection threshold

if ~exist('fvRunIdx', 'var'); fvRunIdx = 1:size(configs, 1); end

for ci = fvRunIdx(:)'
    fpath   = configs{ci, 1};
    outPath = configs{ci, 2};

    M = readmatrix(fpath);
    idx0 = find(M(:,1) >= -0.1, 1);
    t = M(idx0:end, 1) / 1000;   % ms -> s
    L = M(idx0:end, 2);
    F = M(idx0:end, 3);

    clip = t >= fvSegment(1) & t <= fvSegment(2);

    % ── Detect the phase boundaries ─────────────────────────────────────────
    L_base = median(L(t >= preWin(1) & t < preWin(2)));
    L_min  = min(L(clip));

    % The isovelocity ramp is the one sustained stretch of MODERATE shortening
    % speed: the step and the slack release are an order of magnitude faster,
    % the holds an order slower. Take the longest such run inside the clip.
    dt_s   = median(diff(t));
    v_s    = gradient(movmean(L, max(1, round(5e-4/dt_s))), t);   % 0.5 ms smoothing
    isRamp = clip & v_s < -0.5 & v_s > -5;
    d  = diff([0; isRamp; 0]);
    rs = find(d == 1); re = find(d == -1) - 1;
    [~, iLong] = max(t(re) - t(rs));
    t_r1 = t(rs(iLong));  L_r1 = L(rs(iLong));   % ramp start (after the step)
    t_r2 = t(re(iLong));  L_r2 = L(re(iLong));   % ramp end (before the slack)

    % Step release: last sample still at baseline before the ramp
    i_step = find(clip & t < t_r1 & L >= L_base - thresh, 1, 'last');
    t_step = t(i_step);

    % Slack bottom: first sample within thresh of the minimum
    i_bot = find(t > t_r2 & L <= L_min + thresh, 1, 'first');
    t_bot = t(i_bot);

    % Restretch start: last sample near the minimum before L rises
    i_rise = find(t > t_bot & L > L_min + thresh, 1, 'first');
    t_rs   = t(max(1, i_rise - 1));

    % Top: first sample back at the baseline plateau
    i_top = find(t > t_rs & L >= L_base - thresh, 1, 'first');
    t_top = t(i_top);

    % ── Segment velocities ──────────────────────────────────────────────────
    v_step      = (L_r1   - L_base) / (t_r1  - t_step);
    v_ramp      = (L_r2   - L_r1)   / (t_r2  - t_r1);
    v_slack     = (L_min  - L_r2)   / (t_bot - t_r2);
    v_restretch = (L_base - L_min)  / (t_top - t_rs);

    % ── 8-row velocitytable: [t, vel(ML/s), vel_um(um/s), L(ML)] ────────────
    velocitytable = [
        t_step - warmupLead, 0,            0,              L_base;
        t_step,              v_step,       v_step*2,       L_base;
        t_r1,                v_ramp,       v_ramp*2,       L_r1;
        t_r2,                v_slack,      v_slack*2,      L_r2;
        t_bot,               0,            0,              L_min;
        t_rs,                v_restretch,  v_restretch*2,  L_min;
        t_top,               0,            0,              L_base;
        fvSegment(2),        0,            0,              L_base;
    ];

    % ── Datatable (downsampled, force normalised by pre-ramp Fss) ───────────
    dsf       = 5;
    datatable = [downsample(t, dsf), downsample(L, dsf), downsample(F, dsf)];

    Fss = mean(datatable(datatable(:,1) >= preWin(1) & datatable(:,1) < preWin(2), 3));
    if Fss == 0; Fss = 1; end
    datatable(:, 3) = datatable(:, 3) / Fss;

    % ML fraction -> um (SL ~ L x 2)
    if mean(datatable(:,2)) < 1.5
        datatable(:,2) = datatable(:,2) * 2;
    end
    datatable = datatable(datatable(:,1) >= fvSegment(1) & datatable(:,1) <= fvSegment(2), :);

    % ── Features (absolute kPa; the physiology the protocol measures) ───────
    % The plateau is read over the middle 60% of the ramp, skipping the
    % transients where the step settles and where the slack release begins.
    ramp_span = t_r2 - t_r1;
    wPlateau  = t >= t_r1 + 0.2*ramp_span & t <= t_r2 - 0.2*ramp_span;
    wRestr    = t >= t_rs & t <= fvSegment(2);
    F_plateau = mean(F(wPlateau));

    % Redevelopment half-time. The restretch ends on a large ELASTIC peak that
    % already exceeds Fss, so "time to x% of Fss" measured from t_top is
    % degenerate. Measure instead from the minimum force reached after that
    % peak has decayed: half-time to recover from there back to Fss.
    wPost   = find(t > t_top & t <= fvSegment(2));
    [Fdip, iRel] = min(movmean(F(wPost), max(1, round(1e-3/dt_s))));
    iDip    = wPost(iRel);
    Fhalf   = Fdip + 0.5*(Fss - Fdip);
    iHalf   = find(t > t(iDip) & F >= Fhalf, 1, 'first');
    if isempty(iHalf); trec50 = NaN; else; trec50 = t(iHalf) - t(iDip); end

    features_data = struct( ...
        'FV2_v',             v_ramp, ...            % commanded ramp velocity (ML/s)
        'FV2_Fss',           Fss, ...               % pre-ramp isometric force (kPa)
        'FV2_Fplateau',      F_plateau, ...         % force borne during the ramp (kPa)
        'FV2_Fnorm',         F_plateau / Fss, ...   % ... relative to isometric
        'FV2_Fslack',        mean(F(t >= t_bot & t <= t_rs)), ... % force while slack (kPa)
        'FV2_restretchPeak', max(F(wRestr)), ...    % peak force after the restretch (kPa)
        'FV2_postDip',       Fdip, ...              % force minimum after the elastic peak (kPa)
        'FV2_trec50',        trec50);               % half-time, dip -> Fss (s)

    save(outPath, 'velocitytable', 'datatable', 'features_data', 'Fss');

    fprintf('=== %s ===\n', fpath);
    fprintf('  Phases: step %.1f-%.1f ms | RAMP %.1f-%.1f ms | slack %.1f-%.1f ms | restretch %.1f-%.1f ms\n', ...
            1000*[t_step t_r1 t_r1 t_r2 t_r2 t_bot t_rs t_top]);
    fprintf('  L %.4f -> %.4f -(ramp)-> %.4f -> %.4f -> %.4f ML\n', L_base, L_r1, L_r2, L_min, L_base);
    fprintf('  v: step=%.2f  RAMP=%.3f ML/s (%.2f um/s SL)  slack=%.2f  restretch=%.2f ML/s\n', ...
            v_step, v_ramp, v_ramp*2, v_slack, v_restretch);
    fprintf('  Fss=%.2f kPa | RAMP PLATEAU=%.2f kPa (%.3f of Fss) | slack=%.2f kPa\n', ...
            Fss, F_plateau, F_plateau/Fss, features_data.FV2_Fslack);
    fprintf('  restretchPeak=%.2f kPa | postDip=%.2f kPa | trec50=%.1f ms\n', ...
            features_data.FV2_restretchPeak, Fdip, 1000*trec50);
    fprintf('  Saved %dx4 velocitytable, %dx3 datatable to %s\n\n', ...
            size(velocitytable,1), size(datatable,1), outPath);

    % ── Diagnostic plot ─────────────────────────────────────────────────────
    L_ideal = zeros(size(velocitytable,1), 1);
    L_ideal(1) = velocitytable(1, 4);
    for k = 2:size(velocitytable, 1)
        L_ideal(k) = L_ideal(k-1) + velocitytable(k-1,2) * (velocitytable(k,1) - velocitytable(k-1,1));
    end

    figure(70+ci); clf; tiledlayout(2,1);
    ax1 = nexttile;
    plot(t, L, 'Color',[0.7 0.7 0.7]); hold on;
    plot(velocitytable(:,1), L_ideal, 'r-', 'LineWidth', 1.5);
    plot(velocitytable(:,1), velocitytable(:,4), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    title(sprintf('FV2 VT (%d pts): %s', size(velocitytable,1), fpath), 'Interpreter', 'none');
    ylabel('Length (Lo)'); xlabel('Time (s)');
    legend('Raw L', 'Ideal L', 'VT points', 'Location', 'best');

    ax2 = nexttile;
    plot(t, F/Fss, 'Color',[0.7 0.7 0.7]); hold on;
    plot(datatable(:,1), datatable(:,3), 'r-', 'LineWidth', 1.5);
    ylabel('Force (norm.)'); xlabel('Time (s)');
    legend('Raw F/Fss', 'datatable', 'Location', 'best');
    linkaxes([ax1 ax2], 'x');
    xlim(ax1, fvSegment + [-0.02 0.02]);
end
