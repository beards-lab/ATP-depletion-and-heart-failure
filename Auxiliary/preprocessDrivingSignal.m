%% preprocessDrivingSignal.m
%
% Loads 06_Merged_8mM_Active_PNB_Mava.txt, resamples to a uniform 500 Hz
% grid, then compares four smoothing strategies on the burst region (69–78s):
%
%   1. Butterworth LP  — classic zero-phase IIR
%   2. Savitzky-Golay  — local polynomial (wide window = 200 ms)
%   3. Envelope midline — (upper + lower pchip envelope) / 2
%   4. LOESS           — locally weighted regression (5% window)
%
% The script ends with a figure showing all four on the full signal and three
% burst zoom panels so you can pick the one that matches "pencil midline".
%
% After choosing a method, export to .mat via the save() line at the bottom
% and use sigLookup (local function) for O(1) interpolation inside dxdt.
%
% Output struct 'sig' fields:
%   t0     — start time (s)
%   Fs     — output sample rate (Hz)
%   N      — number of samples
%   y      — smoothed signal (N×1)  ← swap in your chosen method
%   m      — Catmull-Rom slopes (N×1) for O(1) cubic Hermite interp
%
% See also: sigLookup (local function at bottom)

clear; clc;
cd(fileparts(mfilename('fullpath')));
addpath(genpath('..'));

%% ── 1. Load ──────────────────────────────────────────────────────────────
fname = '..\data\03 27 2026 M\06_Merged_8mM_Active_PNB_Mava.txt';
fid   = fopen(fname);
C     = textscan(fid, '%f%f%f');
fclose(fid);

t_raw = C{1};   % time (s), non-uniform
y_raw = C{3};
fprintf('Loaded: %d pts  t=[%.1f .. %.1f]s\n', numel(t_raw), t_raw(1), t_raw(end));

%% ── 2. Resample to uniform 500 Hz grid ───────────────────────────────────
Fs_out = 1000;
t0     = t_raw(1);
t_end  = t_raw(end);
t_uni  = (t0 : 1/Fs_out : t_end)';
N      = numel(t_uni);
y_uni  = interp1(t_raw, y_raw, t_uni, 'pchip');
fprintf('Uniform %g Hz grid: %d pts\n', Fs_out, N);

%% ── 3. Four smoothing methods ────────────────────────────────────────────

% --- 3a. Butterworth LP (200 Hz zero-phase) ---
lp_hz      = 400;
[b, a]     = butter(4, lp_hz / (Fs_out/2), 'low');
y_butter   = filtfilt(b, a, y_uni);

% --- 3b. Savitzky-Golay, wide window (200 ms = 101 samples at 500 Hz) ---
sg_win     = 101;   % must be odd; increase for more smoothing
sg_order   = 3;
y_sg       = sgolayfilt(y_uni, sg_order, sg_win);

% --- 3c. LOESS (locally weighted regression, 5% bandwidth) ---
loess_span = 0.05;
y_loess    = smooth(t_uni, y_uni, loess_span, 'loess');

% --- 3d. Envelope midline (upper+lower pchip / 2) ---
% Peaks at least 10 ms apart to avoid noise peaks
min_dist = round(Fs_out * 0.01);
[~, idx_max] = findpeaks( y_uni, 'MinPeakDistance', min_dist);
[~, idx_min] = findpeaks(-y_uni, 'MinPeakDistance', min_dist);

if numel(idx_max) > 3 && numel(idx_min) > 3
    env_upper   = interp1(t_uni(idx_max), y_uni(idx_max), t_uni, 'pchip', 'extrap');
    env_lower   = interp1(t_uni(idx_min), y_uni(idx_min), t_uni, 'pchip', 'extrap');
    y_envelope  = (env_upper + env_lower) / 2;
else
    warning('Not enough peaks for envelope — falling back to SG');
    y_envelope = y_sg;
end
fprintf('Envelope peaks: upper=%d  lower=%d\n', numel(idx_max), numel(idx_min));

%% ── 4. Choose method for export ──────────────────────────────────────────
% Swap y_filt to whichever method you prefer:
%   y_butter  — Butterworth (fastest transitions, may ring slightly)
%   y_sg      — Savitzky-Golay wide (smooth midline, no ringing)
%   y_loess   — LOESS (similar to SG, slightly heavier compute)
%   y_envelope — envelope midline (best "pencil" in active burst, poor in noise-only regions)
y_filt = y_sg;   % ← CHANGE THIS

%% ── 5. Catmull-Rom slopes for O(1) cubic Hermite interp ─────────────────
m      = zeros(N, 1);
m(2:end-1) = (y_filt(3:end) - y_filt(1:end-2)) * (Fs_out/2);
m(1)       = (y_filt(2) - y_filt(1))       * Fs_out;
m(end)     = (y_filt(end) - y_filt(end-1)) * Fs_out;

%% ── 6. Package and save ──────────────────────────────────────────────────
sig.t0  = t0;
sig.Fs  = Fs_out;
sig.N   = N;
sig.y   = y_filt;
sig.m   = m;

out_mat = '..\data\03 27 2026 M\06_Merged_8mM_Active_PNB_Mava_proc.mat';
% save(out_mat, 'sig');
% fprintf('Saved: %s\n', out_mat);

%% ── 7. Comparison figure ─────────────────────────────────────────────────
colors = struct( ...
    'butter',   [0.18 0.45 0.70], ...
    'sg',       [0.85 0.33 0.10], ...
    'loess',    [0.77 0.13 0.77], ...
    'envelope', [0.17 0.63 0.17]);

figure('Name','Smoothing method comparison', 'Position',[40 40 1300 900]);

% ── Panel 1: full signal ──
subplot(3,2,[1 2]);
plot(t_raw, y_raw, 'Color',[0.82 0.82 0.82], 'DisplayName','raw'); hold on;
plot(t_uni, y_butter,  '-',  'Color',colors.butter,   'LineWidth',1.4, 'DisplayName',sprintf('Butterworth %gHz LP',lp_hz));
plot(t_uni, y_sg,      '-',  'Color',colors.sg,       'LineWidth',1.4, 'DisplayName',sprintf('SG win=%d (%.0fms)', sg_win, sg_win/Fs_out*1e3));
plot(t_uni, y_loess,   '-',  'Color',colors.loess,    'LineWidth',1.4, 'DisplayName',sprintf('LOESS %.0f%%', loess_span*100));
plot(t_uni, y_envelope,'-',  'Color',colors.envelope, 'LineWidth',1.8, 'DisplayName','Envelope midline');
legend('Location','northwest'); grid on;
title('Full signal — all methods'); ylabel('signal'); xlabel('time (s)');
xline(69, 'k--', 'burst start'); xline(78, 'k--', 'burst end');

% ── Panel 2: burst overview (69–78 s) ──
subplot(3,2,[3 4]);
bm = t_uni >= 69 & t_uni <= 78;
plot(t_uni(bm), y_uni(bm),      'Color',[0.82 0.82 0.82]); hold on;
plot(t_uni(bm), y_butter(bm),   '-', 'Color',colors.butter,   'LineWidth',1.4);
plot(t_uni(bm), y_sg(bm),       '-', 'Color',colors.sg,       'LineWidth',1.4);
plot(t_uni(bm), y_loess(bm),    '-', 'Color',colors.loess,    'LineWidth',1.4);
plot(t_uni(bm), y_envelope(bm), '-', 'Color',colors.envelope, 'LineWidth',1.8);
grid on; title('Burst region (69–78 s)'); ylabel('signal'); xlabel('time (s)');

% ── Panel 3: early burst / onset zoom (69–70.5 s) ──
subplot(3,2,5);
zm1 = t_uni >= 69 & t_uni <= 70.5;
plot(t_uni(zm1), y_uni(zm1),      'Color',[0.82 0.82 0.82]); hold on;
plot(t_uni(zm1), y_butter(zm1),   '-', 'Color',colors.butter,   'LineWidth',1.4);
plot(t_uni(zm1), y_sg(zm1),       '-', 'Color',colors.sg,       'LineWidth',1.4);
plot(t_uni(zm1), y_loess(zm1),    '-', 'Color',colors.loess,    'LineWidth',1.4);
plot(t_uni(zm1), y_envelope(zm1), '-', 'Color',colors.envelope, 'LineWidth',1.8);
grid on; title('Onset zoom (69–70.5 s)'); ylabel('signal'); xlabel('time (s)');

% ── Panel 4: active peaks zoom (74–76 s) ──
subplot(3,2,6);
zm2 = t_uni >= 74 & t_uni <= 76;
plot(t_uni(zm2), y_uni(zm2),      'Color',[0.82 0.82 0.82]); hold on;
plot(t_uni(zm2), y_butter(zm2),   '-', 'Color',colors.butter,   'LineWidth',1.4);
plot(t_uni(zm2), y_sg(zm2),       '-', 'Color',colors.sg,       'LineWidth',1.4);
plot(t_uni(zm2), y_loess(zm2),    '-', 'Color',colors.loess,    'LineWidth',1.4);
plot(t_uni(zm2), y_envelope(zm2), '-', 'Color',colors.envelope, 'LineWidth',1.8);
grid on; title('Active peaks zoom (74–76 s)'); ylabel('signal'); xlabel('time (s)');

% shared legend via invisible axes
legend_labels = {'raw @500Hz', ...
    sprintf('Butterworth %gHz', lp_hz), ...
    sprintf('SG wide (%.0fms window)', sg_win/Fs_out*1e3), ...
    sprintf('LOESS %.0f%%', loess_span*100), ...
    'Envelope midline'};
sgtitle(sprintf('Smoothing comparison  |  active method → %s  (change y_filt in section 4)', ...
    'y_sg'));

%% ── 8. Envelope components detail ───────────────────────────────────────
figure('Name','Envelope detail', 'Position',[40 500 900 400]);
bm = t_uni >= 69 & t_uni <= 78;
plot(t_uni(bm), y_uni(bm), 'Color',[0.82 0.82 0.82]); hold on;
plot(t_uni(bm), env_upper(bm), 'r-', 'LineWidth',1.2, 'DisplayName','upper envelope');
plot(t_uni(bm), env_lower(bm), 'b-', 'LineWidth',1.2, 'DisplayName','lower envelope');
plot(t_uni(bm), y_envelope(bm), 'g-', 'LineWidth',2, 'DisplayName','midline');
% mark the peak points used
in_burst_max = idx_max(t_uni(idx_max) >= 69 & t_uni(idx_max) <= 78);
in_burst_min = idx_min(t_uni(idx_min) >= 69 & t_uni(idx_min) <= 78);
plot(t_uni(in_burst_max), y_uni(in_burst_max), 'r^', 'MarkerSize',4, 'DisplayName','maxima used');
plot(t_uni(in_burst_min), y_uni(in_burst_min), 'bv', 'MarkerSize',4, 'DisplayName','minima used');
legend('Location','best'); grid on;
title('Envelope midline — components (burst 69–78 s)');
xlabel('time (s)'); ylabel('signal');

%% ── 9. How to use inside dxdt ────────────────────────────────────────────
fprintf('\n── dxdt usage ──────────────────────────────────────────\n');
fprintf('  load(''06_Merged_8mM_Active_PNB_Mava_proc.mat'', ''sig'');\n');
fprintf('  params.driveSig = sig;\n');
fprintf('  %% inside dPUdT_CombinedTransitions:\n');
fprintf('  y_t = sigLookup(params.driveSig, t);   %% O(1), no search\n');
fprintf('────────────────────────────────────────────────────────\n');

%% ── Local function: O(1) cubic Hermite lookup ────────────────────────────
function y = sigLookup(sig, t)
% SIGLOOKUP  O(1) cubic Hermite interpolation on a uniform grid.
%
%   Clamps t to [t0, t0+(N-1)/Fs] — safe to call outside the signal range.
%   Uses Catmull-Rom slopes precomputed in sig.m.

    k   = max(1, min(sig.N - 1, 1 + floor((t - sig.t0) * sig.Fs)));
    tau = (t - (sig.t0 + (k-1)/sig.Fs)) * sig.Fs;   % in [0, 1)

    t2  = tau*tau;  t3 = t2*tau;
    h00 =  2*t3 - 3*t2 + 1;
    h10 =    t3 - 2*t2 + tau;
    h01 = -2*t3 + 3*t2;
    h11 =    t3 -   t2;

    inv_Fs = 1 / sig.Fs;
    y = h00*sig.y(k) + h10*sig.m(k)*inv_Fs ...
      + h01*sig.y(k+1) + h11*sig.m(k+1)*inv_Fs;
end
