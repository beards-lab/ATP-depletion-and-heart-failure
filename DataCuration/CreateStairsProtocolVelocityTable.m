% CreateStairsProtocolVelocityTable.m
% Generates stairs (ramp-up) velocity-table protocol files from
% 02_Merged_*mM_Active.txt data.
%
% The stairs protocol is an incremental staircase stretch (analogous to the
% legacy bakers_rampup8): the fiber is held isometric, then stretched in a
% series of small steps (each a fast ramp + hold) up to ~110% ML, and the
% force response to each step is recorded.
%
%   03/27 (t≈72.47–73.18): 8 steps of +0.01 ML  → levels 1.01 … 1.07, 1.10
%   04/03 (t≈72.53–73.15): 4 steps (+0.01,+0.03,+0.03,+0.03) → 1.01,1.04,1.07,1.10
%
% Steps are detected generically from the (smoothed) length velocity, so the
% same code adapts to either step count. Force in datatable is normalized by
% pre-stretch isometric Fss; length is stored in µm (SL ≈ L × 2).
%
% Note: 04/03 files have 4 columns (t, L, F, SL); only cols 1-3 are used.
%
% Produces:
%   data/protocol_03_27_2026_velocitytable_stairs.mat
%   data/protocol_04_03_2026_velocitytable_stairs.mat

% {file, clip window [s], warm-start t [s], pre-stretch Fss window [s], output}
configs = {
  '../data/03 27 2026 M/02_Merged_8mM_Active.txt', [72.3 73.6], 67.0, [72.30 72.45], '../data/protocol_03_27_2026_velocitytable_stairs.mat';
  '../data/04 03 2026 F/02_Merged_8mM_Active.txt', [72.3 73.6], 67.0, [72.30 72.50], '../data/protocol_04_03_2026_velocitytable_stairs.mat';
};

for c = 1:size(configs, 1)
    fpath = configs{c,1}; seg = configs{c,2}; t_ws = configs{c,3};
    prewin = configs{c,4}; outp = configs{c,5};

    M = readmatrix(fpath); t = M(:,1); L = M(:,2); F = M(:,3);
    mask = t >= seg(1) & t <= seg(2);
    tz = t(mask); Lz = L(mask); Fz = F(mask);

    L_base = median(Lz(tz >= prewin(1) & tz < prewin(2)));
    Fss    = mean(Fz(tz  >= prewin(1) & tz < prewin(2)));

    % ── Detect staircase ramps from smoothed velocity ───────────────────────────
    dt  = median(diff(tz));
    win = round(0.001 / dt);            % 1 ms smoothing (holds noisy at 50 µs sampling)
    Ls  = movmean(Lz, win);
    v   = gradient(Ls, tz);             % ML/s
    isRamp = abs(v) > 0.3;              % ramps ~2.5 ML/s; hold noise < 0.1 ML/s

    d  = diff([0; isRamp; 0]);
    rs = find(d == 1); re = find(d == -1) - 1;
    keep = (tz(re) - tz(rs)) >= 0.0008; % drop < 0.8 ms noise spikes
    rs = rs(keep); re = re(keep);

    grpS = []; grpE = []; ms = rs(1); me = re(1);
    for i = 2:numel(rs)
        if tz(rs(i)) - tz(me) < 0.02    % merge sub-ramps within one step (< 20 ms)
            me = re(i);
        else
            grpS(end+1) = ms; grpE(end+1) = me; ms = rs(i); me = re(i);
        end
    end
    grpS(end+1) = ms; grpE(end+1) = me;

    % ── Build velocitytable: warm-start hold, then [ramp, hold] per step ─────────
    VT = [t_ws, 0, 0, L_base];
    Lprev = L_base;
    for i = 1:numel(grpS)
        t0 = tz(grpS(i)); t1 = tz(grpE(i));
        if i < numel(grpS); iHE = grpS(i+1) - 1; else; iHE = numel(Lz); end
        L_lvl = median(Lz(grpE(i):iHE));      % plateau level after this ramp
        vseg  = (L_lvl - Lprev) / (t1 - t0);
        VT = [VT; t0, vseg, vseg*2, Lprev];   % ramp
        VT = [VT; t1, 0,    0,      L_lvl];   % hold
        Lprev = L_lvl;
    end
    VT = [VT; seg(2), 0, 0, Lprev];           % final end row
    velocitytable = VT;

    % ── Datatable (downsampled, normalized force, µm length) ─────────────────────
    dsf = 5; mask_dt = t >= prewin(1) & t <= seg(2);
    datatable = [downsample(t(mask_dt), dsf), ...
                 downsample(L(mask_dt), dsf) * 2, ...
                 downsample(F(mask_dt), dsf) / Fss];

    features_data = struct();
    save(outp, 'velocitytable', 'datatable', 'features_data');

    fprintf('%s\n  Fss=%.2f mN, %d steps, L: %.3f → %.3f, VT=%dx4, DT=%dx3\n', ...
        outp, Fss, numel(grpS), L_base, Lprev, size(VT,1), size(datatable,1));
    fprintf('  levels: '); fprintf('%.3f ', VT(3:2:end,4)); fprintf('\n');
end
