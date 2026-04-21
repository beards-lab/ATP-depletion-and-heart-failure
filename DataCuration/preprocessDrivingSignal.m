%% preprocessDrivingSignal.m
%
% Loads a raw merged txt file, resamples to a uniform grid, applies three
% smoothing methods (Butterworth LP, Savitzky-Golay, LOESS), and exports
% one .mat per method × subsignal combination.
%
% Output modes (controlled by sigOnly):
%   sigOnly = true  → <base>_<method>_<subname>_sig.mat  (sig struct only)
%   sigOnly = false → <base>_<method>_<subname>_all.mat  (sig + velocitytable
%                     + datatable with filtered force column)
%
% sig struct fields (sigLookup-compatible):
%   t0  — clip start time (s)
%   Fs  — sample rate (Hz)
%   N   — number of samples
%   y   — filtered signal (N×1)
%   m   — Catmull-Rom slopes (N×1)
%   t   — time vector (optional, only if exportTime = true)

clear; clc;
cd(fileparts(mfilename('fullpath')));
addpath(genpath('..'));

%% ── 0. Configuration ────────────────────────────────────────────────────
dataDir  = '..\data\03 27 2026 M\';
fname    = '06_Merged_8mM_Active_PNB_Mava.txt';
% fname = '01_Merged_Relax.txt';
outDir   = dataDir;

Fs_out     = 1000;    % resample rate (Hz)
lp_hz      = 200;    % Butterworth LP cutoff (Hz)
sg_win     = 15;     % SG window (samples, must be odd)
sg_order   = 3;
loess_span = 0.0001; % LOESS fraction of data width

% Subsignals: each entry defines one clipped output window.
% vtPath: path to matching velocitytable mat (used only when sigOnly=false).
subsignals(1).name   = 'slack';
subsignals(1).tspan  = [74, 77.5];
subsignals(1).vtPath = '..\data\protocol_03_27_2026_velocitytable_slack.mat';
subsignals(2).name   = 'stairs';
subsignals(2).tspan  = [72, 74];
subsignals(2).vtPath = '';

exportMethods = {'butter', 'sg', 'loess'};

% sigOnly = true  → save sig struct only, suffix _sig
% sigOnly = false → save sig + velocitytable + datatable, suffix _all
sigOnly    = true;
exportTime = false;  % include sig.t (time vector) in output
dryRun = false; % run without exporting files

%% ── 1. Load ─────────────────────────────────────────────────────────────
fid   = fopen(fullfile(dataDir, fname));
C     = textscan(fid, '%f%f%f');
fclose(fid);


t_raw = C{1};   % time (s), non-uniform
[t_raw i_tuniq] = unique(t_raw);
y_raw = C{3}(i_tuniq);   % signal (force / length)
fprintf('Loaded: %d pts  t=[%.1f .. %.1f]s\n', numel(t_raw), t_raw(1), t_raw(end));

%% ── 2. Resample to uniform grid ─────────────────────────────────────────
t0    = t_raw(1);
t_end = t_raw(end);
t_uni = (t0 : 1/Fs_out : t_end)';
N     = numel(t_uni);
y_uni = interp1(t_raw, y_raw, t_uni, 'pchip');
fprintf('Uniform %g Hz grid: %d pts\n', Fs_out, N);

%% ── 3. Smoothing methods ────────────────────────────────────────────────

% 3a. Butterworth LP (zero-phase)
[b, a]   = butter(4, lp_hz / (Fs_out/2), 'low');
y_butter = filtfilt(b, a, y_uni);

% 3b. Savitzky-Golay
y_sg = sgolayfilt(y_uni, sg_order, sg_win);

% 3c. LOESS
y_loess = smooth(t_uni, y_uni, loess_span, 'loess');

% 3d. Envelope midline (figures only — not exported)
min_dist = round(Fs_out * 0.005);
[~, idx_max] = findpeaks( y_uni, 'MinPeakDistance', min_dist);
[~, idx_min] = findpeaks(-y_uni, 'MinPeakDistance', min_dist);
if numel(idx_max) > 3 && numel(idx_min) > 3
    env_upper  = interp1(t_uni(idx_max), y_uni(idx_max), t_uni, 'pchip', 'extrap');
    env_lower  = interp1(t_uni(idx_min), y_uni(idx_min), t_uni, 'pchip', 'extrap');
    y_envelope = (env_upper + env_lower) / 2;
else
    warning('Not enough peaks for envelope — falling back to SG');
    y_envelope = y_sg;
    env_upper  = y_sg;
    env_lower  = y_sg;
end
fprintf('Envelope peaks: upper=%d  lower=%d\n', numel(idx_max), numel(idx_min));

%% ── 4. Export: all methods × all subsignals ─────────────────────────────
[~, base, ~] = fileparts(fname);
methodMap    = struct('butter', y_butter, 'sg', y_sg, 'loess', y_loess);

for mi = 1:numel(exportMethods)
    mname  = exportMethods{mi};
    y_filt = methodMap.(mname);

    for si = 1:numel(subsignals)
        sname = subsignals(si).name;
        ts    = subsignals(si).tspan;

        k1 = max(1, round((ts(1) - t0) * Fs_out) + 1);
        k2 = min(N, round((ts(2) - t0) * Fs_out) + 1);
        y_clip = y_filt(k1:k2);
        Nc     = numel(y_clip);

        % Catmull-Rom slopes — computed fresh on the clipped signal
        mc = zeros(Nc, 1);
        mc(2:end-1) = (y_clip(3:end) - y_clip(1:end-2)) * (Fs_out/2);
        mc(1)       = (y_clip(2)     - y_clip(1))       * Fs_out;
        mc(end)     = (y_clip(end)   - y_clip(end-1))   * Fs_out;

        sig.t0 = t0 + (k1-1) / Fs_out;
        sig.Fs = Fs_out;
        sig.N  = Nc;
        sig.y  = y_clip;
        sig.m  = mc;
        if exportTime
            sig.t = (sig.t0 : 1/Fs_out : sig.t0 + (Nc-1)/Fs_out)';
        end

        if dryRun
            fprintf('Skipping exporting [t=%.2f..%.2f s  N=%d] \n', sig.t0, sig.t0+(Nc-1)/Fs_out, Nc);
            continue;
        end
        if sigOnly
            outpath = fullfile(outDir, sprintf('%s_%s_%s_sig.mat', base, mname, sname));
            save(outpath, 'sig');
        else
            outpath = fullfile(outDir, sprintf('%s_%s_%s_all.mat', base, mname, sname));
            vtp = subsignals(si).vtPath;
            if ~isempty(vtp) && isfile(vtp)
                tmp           = load(vtp, 'velocitytable', 'datatable');
                velocitytable = tmp.velocitytable;
                datatable     = tmp.datatable;
                % Replace raw force (col 3) with filtered signal at datatable times
                datatable(:,3) = interp1(t_uni, y_filt, datatable(:,1), 'pchip', 'extrap');
                save(outpath, 'sig', 'velocitytable', 'datatable');
            else
                warning('vtPath missing for subsignal "%s" — saving sig only', sname);
                save(outpath, 'sig');
            end
        end
        fprintf('Saved  [t=%.2f..%.2f s  N=%d]  %s\n', sig.t0, sig.t0+(Nc-1)/Fs_out, Nc, outpath);
    end
end

%% ── 5. Comparison figure ─────────────────────────────────────────────────
colors = struct( ...
    'butter',   [0.18 0.45 0.70], ...
    'sg',       [0.85 0.33 0.10], ...
    'loess',    [0.77 0.13 0.77], ...
    'envelope', [0.17 0.63 0.17]);

figure('Name','Smoothing method comparison', 'Position',[40 40 1300 900]);

subplot(2,2,[1 2]);
plot(t_raw, y_raw, 'Color',[0.82 0.82 0.82], 'DisplayName','raw'); hold on;
plot(t_uni, y_butter,  '-', 'Color',colors.butter,   'LineWidth',1.4, 'DisplayName',sprintf('Butterworth %gHz LP',lp_hz));
plot(t_uni, y_sg,      '-', 'Color',colors.sg,       'LineWidth',1.4, 'DisplayName',sprintf('SG win=%d (%.0fms)',sg_win,sg_win/Fs_out*1e3));
plot(t_uni, y_loess,   '-', 'Color',colors.loess,    'LineWidth',1.4, 'DisplayName',sprintf('LOESS %.2f%%',loess_span*100));
plot(t_uni, y_envelope,'-', 'Color',colors.envelope, 'LineWidth',1.8, 'DisplayName','Envelope midline');
legend('Location','northwest'); grid on;
title('Full signal — all methods'); ylabel('signal'); xlabel('time (s)');
for ssi = 1:numel(subsignals)
    xline(subsignals(ssi).tspan(1), 'k--');
    xline(subsignals(ssi).tspan(2), 'k--', subsignals(ssi).name);
end

subplot(2,2,[3 4]);
bm = t_uni >= 69 & t_uni <= 78;
plot(t_uni(bm), y_uni(bm),      'Color',[0.82 0.82 0.82]); hold on;
plot(t_uni(bm), y_butter(bm),   '-', 'Color',colors.butter,   'LineWidth',1.4);
plot(t_uni(bm), y_sg(bm),       '-', 'Color',colors.sg,       'LineWidth',1.4);
plot(t_uni(bm), y_loess(bm),    '-', 'Color',colors.loess,    'LineWidth',1.4);
plot(t_uni(bm), y_envelope(bm), '-', 'Color',colors.envelope, 'LineWidth',1.8);
grid on; title('Burst region (69–78 s)'); ylabel('signal'); xlabel('time (s)');
sgtitle(sprintf('Smoothing comparison  |  methods: %s', strjoin(exportMethods, ', ')));

return
%% ── 6. Envelope components detail ───────────────────────────────────────
figure('Name','Envelope detail', 'Position',[40 500 900 400]);
bm = t_uni >= 69 & t_uni <= 78;
plot(t_uni(bm), y_uni(bm),      'Color',[0.82 0.82 0.82]); hold on;
plot(t_uni(bm), env_upper(bm),  'r-', 'LineWidth',1.2, 'DisplayName','upper envelope');
plot(t_uni(bm), env_lower(bm),  'b-', 'LineWidth',1.2, 'DisplayName','lower envelope');
plot(t_uni(bm), y_envelope(bm), 'g-', 'LineWidth',2,   'DisplayName','midline');
in_burst_max = idx_max(t_uni(idx_max) >= 69 & t_uni(idx_max) <= 78);
in_burst_min = idx_min(t_uni(idx_min) >= 69 & t_uni(idx_min) <= 78);
plot(t_uni(in_burst_max), y_uni(in_burst_max), 'r^', 'MarkerSize',4, 'DisplayName','maxima');
plot(t_uni(in_burst_min), y_uni(in_burst_min), 'bv', 'MarkerSize',4, 'DisplayName','minima');
legend('Location','best'); grid on;
title('Envelope midline — components (burst 69–78 s)');
xlabel('time (s)'); ylabel('signal');
