% CreateKtrProtocolVelocityTable.m
% Generates ktr velocity-table protocol files from raw stiff_ktr.txt data.
%
% The ktr protocol: slow ramp to 80% ML → hold → slow restretch to 105% ML
% → brief hold → return to 100% ML → isometric force redevelopment.
%
% Produces exactly 8 velocitytable rows:
%   [warm-start, shortening, hold-80%, restretch, peak-105%, return, recovery, end]
% Force in datatable is normalized by pre-shortening isometric Fss.
%
% Note: 04/03 files have 4 columns (t, L, F, SL); only cols 1-3 are used.
%
% Produces:
%   data/protocol_03_27_2026_velocitytable_ktr.mat      (8mM ATP)
%   data/protocol_03_27_2026_velocitytable_ktr_2mM.mat  (2mM ATP)
%   data/protocol_04_03_2026_velocitytable_ktr.mat      (8mM ATP)
%   data/protocol_04_03_2026_velocitytable_ktr_2mM.mat  (2mM ATP)

configs = {
    '../data/03 27 2026 M/8mM_stiff_ktr.txt', '../data/protocol_03_27_2026_velocitytable_ktr.mat';
    '../data/03 27 2026 M/2mM_stiff_ktr.txt', '../data/protocol_03_27_2026_velocitytable_ktr_2mM.mat';
    '../data/04 03 2026 F/8mM_stiff_ktr.txt', '../data/protocol_04_03_2026_velocitytable_ktr.mat';
    '../data/04 03 2026 F/2mM_stiff_ktr.txt', '../data/protocol_04_03_2026_velocitytable_ktr_2mM.mat';
};

% stiff_ktr Time column is in MILLISECONDS (header "(ms)"); the merged/stairs
% files are in seconds. Divide by 1000 so every protocol shares a seconds
% time axis and the ktr event runs at its true (fast) timescale.
ktrSegment  = [1.583, 1.660]; % clip window around the ktr event (s)
warmupLead  = 5;              % isometric warm-up before shortening (s) — lets
                              % the model reach steady state (the old 15 s
                              % "warm-up" was an artifact of ms-as-seconds)

for ci = 1:size(configs, 1)
    fpath   = configs{ci, 1};
    outPath = configs{ci, 2};

    M = readmatrix(fpath);
    idx0 = find(M(:,1) >= -0.1, 1);
    t = M(idx0:end, 1) / 1000;   % ms → s
    L = M(idx0:end, 2);
    F = M(idx0:end, 3);

    % ── Detect 6 key transition times ───────────────────────────────────────
    L_base = median(L(t >= 1.570 & t < ktrSegment(1)));
    thresh = 0.005;   % 0.5% ML — threshold for detecting transitions

    % 1. Shortening start: last sample at baseline before L drops
    i_drop = find(L < L_base - thresh, 1, 'first');
    t_ss   = t(max(1, i_drop - 1));

    % 2. Hold at 80%: first sample within thresh of the global minimum
    [L_min, i_min] = min(L);
    i_hold = find(L(1:i_min) <= L_min + thresh, 1, 'first');
    t_hs   = t(i_hold);

    % 3. Restretch start: last sample near minimum before L rises
    i_rise = find(t > t_hs & L > L_min + thresh, 1, 'first');
    t_rs   = t(max(1, i_rise - 1));

    % 4. Peak: global maximum
    [L_max, i_max] = max(L);
    t_pk   = t(i_max);

    % 5. Return start: last sample near peak before L falls
    i_fall = find(t > t_pk & L < L_max - thresh, 1, 'first');
    t_ret  = t(max(1, i_fall - 1));

    % 6. Recovery: first sample where L has fallen back below L_base + thresh
    %    (L must drop below this level after coming down from the overshoot)
    i_recov = find(t > t_pk & L < L_base + thresh, 1, 'first');
    t_recov = t(i_recov);

    % ── Segment velocities ────────────────────────────────────────────────────
    v_short     = (L_min  - L_base) / (t_hs   - t_ss);
    v_restretch = (L_max  - L_min)  / (t_pk   - t_rs);
    v_return    = (L_base - L_max)  / (t_recov - t_ret);

    % ── 8-row velocitytable: [t, vel(ML/s), vel_um(µm/s), L(ML)] ────────────
    %    Warm-start row leads shortening by warmupLead seconds of isometric hold.
    velocitytable = [
        t_ss - warmupLead, 0,        0,              L_base;
        t_ss,          v_short,      v_short*2,      L_base;
        t_hs,          0,            0,              L_min;
        t_rs,          v_restretch,  v_restretch*2,  L_min;
        t_pk,          0,            0,              L_max;
        t_ret,         v_return,     v_return*2,     L_max;
        t_recov,       0,            0,              L_base;
        ktrSegment(2), 0,            0,              L_base;
    ];

    % ── Datatable (downsampled, normalized force) ─────────────────────────────
    dsf       = 5;
    datatable = [downsample(t, dsf), downsample(L, dsf), downsample(F, dsf)];

    % Normalize force by pre-shortening Fss (isometric baseline before t_ss)
    i_pre = datatable(:,1) >= ktrSegment(1) & datatable(:,1) < t_ss;
    Fss   = mean(datatable(i_pre, 3));
    if Fss == 0; Fss = 1; end
    datatable(:, 3) = datatable(:, 3) / Fss;

    % Convert ML fraction to µm (SL ≈ L × 2 µm)
    if mean(datatable(:,2)) < 1.5
        datatable(:,2) = datatable(:,2) * 2;
    end
    datatable = datatable(datatable(:,1) >= ktrSegment(1) & datatable(:,1) <= ktrSegment(2), :);

    features_data = struct();
    save(outPath, 'velocitytable', 'datatable', 'features_data');

    fprintf('=== %s ===\n', fpath);
    fprintf('  Keypoints: t_ss=%.2f  t_hs=%.2f  t_rs=%.2f  t_pk=%.2f  t_ret=%.2f  t_recov=%.2f\n', ...
            t_ss, t_hs, t_rs, t_pk, t_ret, t_recov);
    fprintf('  Velocities: v_short=%.5f  v_restretch=%.5f  v_return=%.5f ML/s\n', ...
            v_short, v_restretch, v_return);
    fprintf('  Saved %dx4 velocitytable, %dx3 datatable to %s\n\n', ...
            size(velocitytable,1), size(datatable,1), outPath);

    % ── Diagnostic plot ───────────────────────────────────────────────────────
    L_ideal = zeros(size(velocitytable,1), 1);
    L_ideal(1) = velocitytable(1, 4);
    for k = 2:size(velocitytable, 1)
        dt_k    = velocitytable(k,1) - velocitytable(k-1,1);
        L_ideal(k) = L_ideal(k-1) + velocitytable(k-1,2) * dt_k;
    end

    figure(50+ci); clf; tiledlayout(2,1);
    ax1 = nexttile;
    plot(t, L, 'Color',[0.7 0.7 0.7]); hold on;
    plot(velocitytable(:,1), L_ideal, 'r-', 'LineWidth', 1.5);
    plot(velocitytable(:,1), velocitytable(:,4), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    xlim([ktrSegment(1)-5, ktrSegment(2)+5]);
    title(['Ktr VT (8 pts): ' fpath], 'Interpreter', 'none');
    ylabel('Length (Lo)'); xlabel('Time (s)');
    legend('Raw L', 'Ideal L', 'VT points', 'Location', 'best');

    ax2 = nexttile;
    plot(t, F/Fss, 'Color',[0.7 0.7 0.7]); hold on;
    plot(datatable(:,1), datatable(:,3), 'r-', 'LineWidth', 1.5);
    xlim([ktrSegment(1)-5, ktrSegment(2)+5]);
    ylabel('Force (norm.)'); xlabel('Time (s)');
    legend('Raw F/Fss', 'datatable', 'Location', 'best');
    linkaxes([ax1 ax2], 'x');
end
