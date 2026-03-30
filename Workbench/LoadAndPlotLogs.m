%% Load and plot all Log_*.txt files from a data directory
close all; clear; clc;

dataDir = '../data/03 27 2026 M/';
S = dir([dataDir '*Log_*.txt']);
[~, idx] = sort({S.name});
S = S(idx);

n = length(S);
colors = lines(n);
labels = cell(n, 1);

%% Load all files and detect sync anchor points
% Anchor events (in t-space after initial drop sync):
%   t1: big L drop in 68-72s
%   t2: start of gradual L step-up in 72-73s
%   t3: big L drop in 74-75s (may not exist in all files)

data = struct('t', {}, 'F', {}, 'L', {}, 'anchors', {}, 't_shifted', {}, 'F_shifted', {}, 'L_shifted', {});
for i = 1:n
    fpath = fullfile(S(i).folder, S(i).name);
    M = readmatrix(fpath, 'NumHeaderLines', 4);

    % Initial sync: t=0 at first big L drop
    idx0 = find(diff(M(:,2)) < -0.01, 1);
    t = (M(:,1) - M(idx0,1)) / 1000;

    % t1: big L drop in 68-72s window
    w1 = find(t >= 68 & t <= 72);
    i1 = find(diff(M(w1,2)) < -0.05, 1);
    t1 = t(w1(i1));

    % t2: first step-up (start of ramp) in 72-73s window
    w2 = find(t >= 72 & t <= 73);
    i2 = find(diff(M(w2,2)) > 0.003, 1);
    t2 = t(w2(i2));

    % t3: big L drop in 74-75s window (optional)
    w3 = find(t >= 73.5 & t <= 75.5);
    if ~isempty(w3)
        [min_diff, i3] = min(diff(M(w3,2)));
        if min_diff < -0.01
            t3 = t(w3(i3));
        else
            t3 = NaN;
        end
    else
        t3 = NaN;
    end

    data(i).t       = t;
    data(i).F       = M(:,3);
    data(i).L       = M(:,2);
    data(i).anchors = [t1, t2, t3];
    labels{i} = strrep(S(i).name(1:end-4), '_', ' ');
end

%% Use file 06 (last file) as reference time axis for piecewise shift
% figure(106);hold on;
refIdx = n; % 06_Log_8mM_Active_PNB_Mava is last after sorting
a_ref  = data(refIdx).anchors; % [t1_ref, t2_ref, t3_ref]

for i = 1:n
    t_i = data(i).t;
    a_i = data(i).anchors;

    % Piecewise shift: shift each segment so its anchor aligns with ref
    t_shifted = t_i; % before 65s: no extra shift (initial sync is fine)
    
    mid1 = 71.5; % Midpoint between t1 (~70) and t2 (~72.5)
    mid2 = 73.8; % Midpoint between t2 (~72.5) and t3 (~74.5)

    % Segment 1 (65 to 71.5): shift by (t1_ref - t1)
    seg1 = t_i >= 65 & t_i < mid1;
    t_shifted(seg1) = t_i(seg1) + (a_ref(1) - a_i(1));

    % Segment 2 (71.5 to 73.8): shift by (t2_ref - t2)
    seg2 = t_i >= mid1 & t_i < mid2;
    t_shifted(seg2) = t_i(seg2) + (a_ref(2) - a_i(2));

    % Segment 3 (73.8 to end): shift by (t3_ref - t3) or dt2 if missing
    seg3 = t_i >= mid2;
    if ~isnan(a_i(3)) && ~isnan(a_ref(3))
        t_shifted(seg3) = t_i(seg3) + (a_ref(3) - a_i(3));
    else
        t_shifted(seg3) = t_i(seg3) + (a_ref(2) - a_i(2));
    end

    % Ensure strict monotonicity
    [t_shifted, iu] = unique(t_shifted, 'stable');
    F_i = data(i).F(iu);
    L_i = data(i).L(iu);
    ok = [true; diff(t_shifted) > 0];
    
    data(i).t_shifted_orig = t_shifted(ok);
    data(i).F_shifted_orig = F_i(ok);
    data(i).L_shifted_orig = L_i(ok);
    
    data(i).t_shifted = t_shifted(ok);
    data(i).F_shifted = F_i(ok);
    data(i).L_shifted = L_i(ok);

    % ax1 = nexttile(1);hold on;plot(data(i).t_shifted, data(i).F_shifted);
    % ax2 = nexttile(2);hold on;plot(data(i).t_shifted, data(i).L_shifted);
    % linkaxes([ax1 ax2], 'x');
end

%% Process High-Sampled Files and Merge
allFiles = dir([dataDir '*.txt']);
isLog = contains({allFiles.name}, 'Log_');
highResFiles = allFiles(~isLog);
isMerged = contains({highResFiles.name}, 'Merged_');
highResFiles = highResFiles(~isMerged);

fprintf('Aligning high-resolution burst files...\n');
hi_res_ranges = cell(n, 1);  % track injection time ranges per log
for ii = 1:n; hi_res_ranges{ii} = {}; end
for k = 1:length(highResFiles)
    fpath = fullfile(highResFiles(k).folder, highResFiles(k).name);
    try
        M_hi = readmatrix(fpath, 'NumHeaderLines', 4);

    catch
        continue;
    end
    t_hi = M_hi(:,1)/1000; % ms to s
    L_hi = M_hi(:,2);
    F_hi = M_hi(:,3);
    
    % Ensure t_hi sample points are unique and strictly monotonic for interp1
    [t_hi, ui] = unique(t_hi, 'stable');
    L_hi = L_hi(ui);
    F_hi = F_hi(ui);
    ok = [true; diff(t_hi) > 0];
    t_hi = t_hi(ok);
    L_hi = L_hi(ok);
    F_hi = F_hi(ok);
    
    % figure(k);clf; nexttile(1);ax1 = plot(t_hi, F_hi);ax2 = nexttile(2);plot(t_hi, L_hi);title(highResFiles(k).name);linkaxes([ax1 ax2], 'x');
    
    
    % Find correct Log file by explicit name pairing
    name_hi = highResFiles(k).name;
    best_log_idx = -1;
    
    for i = 1:n
        name_log = S(i).name;
        if contains(name_hi, '2mM') && contains(name_log, '2mM')
            best_log_idx = i; break;
        elseif contains(name_hi, '8mM') && contains(name_hi, 'PNB_Mava') && contains(name_log, 'PNB_Mava')
            best_log_idx = i; break;
        elseif contains(name_hi, '8mM') && ~contains(name_hi, 'PNB_Mava') && contains(name_log, '8mM') && ~contains(name_log, 'PNB_Mava') && ~contains(name_log, 'repeat')
            best_log_idx = i; break;
        end
    end
    
    if best_log_idx > 0
        i = best_log_idx;
        
        % Determine burst type and corresponding log time window
        if contains(name_hi, 'slack')
            t_win = [73, 78];  % slack L-drop happens ~74s (high-force)
            burst_type = 'slack';
        elseif contains(name_hi, 'stiff') || contains(name_hi, 'ktr')
            t_win = [70, 76];  % ktr L-drop happens ~72s
            burst_type = 'ktr';
        elseif contains(name_hi, 'stretch')
            t_win = [72, 78];  % stretch L-drop happens ~74s
            burst_type = 'stretch';
        else
            fprintf('  Unknown burst type for %s, skipping.\n', name_hi);
            continue;
        end
        
        % Extract the log L-signal in the constrained window
        t_sh = data(i).t_shifted_orig;
        L_sh = data(i).L_shifted_orig;
        win_idx = t_sh >= t_win(1) & t_sh <= t_win(2);
        t_log_win = t_sh(win_idx);
        L_log_win = L_sh(win_idx);
        
        if isempty(t_log_win)
            fprintf('  No log data in window [%.0f, %.0f]s for %s, skipping.\n', t_win(1), t_win(2), name_hi);
            continue;
        end
        
        % Resample both to 1kHz for xcorr
        dt = 0.001;
        t_grid_log = (min(t_log_win):dt:max(t_log_win))';
        L_log_1k = interp1(t_log_win, L_log_win, t_grid_log, 'linear', 'extrap');
        
        t_grid_hi = (min(t_hi):dt:max(t_hi))';
        L_hi_1k = interp1(t_hi, L_hi, t_grid_hi, 'linear', 'extrap');
        
        % xcorr on mean-subtracted L
        [R, lags] = xcorr(L_log_1k - mean(L_log_1k), L_hi_1k - mean(L_hi_1k));
        [~, i_max] = max(R);
        lag_samples = lags(i_max);
        time_offset = lag_samples * dt + min(t_log_win) - min(t_hi);
        t_hi_aligned = t_hi + time_offset;
        
        fprintf('  Mapped %-30s -> %-30s | type: %-7s | offset: %.2fs | range: [%.1f, %.1f]s\n', ...
            name_hi, S(i).name, burst_type, time_offset, min(t_hi_aligned), max(t_hi_aligned));
        
        % % --- Diagnostic plot: show alignment quality ---
        % figure(100+k); clf;
        % subplot(2,1,1); hold on;
        % plot(data(i).t_shifted_orig, data(i).L_shifted_orig, 'b', 'LineWidth', 2);
        % plot(t_hi_aligned, L_hi, 'r', 'LineWidth', 1);
        % legend('Log (original)', ['Hi-res: ' name_hi], 'Interpreter', 'none');
        % ylabel('L'); title(['Alignment: ' name_hi ' -> ' S(i).name ' (' burst_type ')'], 'Interpreter', 'none');
        % xline(min(t_hi_aligned), '--r'); xline(max(t_hi_aligned), '--r');
        % xlim(t_win + [-2, 2]);
        % 
        % subplot(2,1,2); hold on;
        % plot(data(i).t_shifted_orig, data(i).F_shifted_orig, 'b', 'LineWidth', 2);
        % plot(t_hi_aligned, F_hi, 'r', 'LineWidth', 1);
        % ylabel('F'); xlabel('Time (s)');
        % xline(min(t_hi_aligned), '--r'); xline(max(t_hi_aligned), '--r');
        % xlim(t_win + [-2, 2]);
        
        
        % Record the injection range and name
        hi_res_ranges{i}{end+1} = struct('t_start', min(t_hi_aligned), ...
            't_end', max(t_hi_aligned), 'name', name_hi);
        
        % Inject high-res data into data(i).t_shifted, L_shifted, F_shifted
        keep = data(i).t_shifted < min(t_hi_aligned) | data(i).t_shifted > max(t_hi_aligned);
        
        new_t = [data(i).t_shifted(keep); t_hi_aligned];
        new_L = [data(i).L_shifted(keep); L_hi];
        new_F = [data(i).F_shifted(keep); F_hi];
        
        % Sort by time
        [new_t, sortIdx] = sort(new_t);
        new_L = new_L(sortIdx);
        new_F = new_F(sortIdx);
        
        % Remove exact duplicates
        [new_t, uniqueIdx] = unique(new_t, 'stable');
        new_L = new_L(uniqueIdx);
        new_F = new_F(uniqueIdx);
        
        data(i).t_shifted = new_t;
        data(i).L_shifted = new_L;
        data(i).F_shifted = new_F;
    else
        fprintf('  Could not find explicit mapping for %s\n', name_hi);
    end
end

%% Save Merged Output Files
fprintf('\nSaving merged outputs...\n');
for i = 1:n
    outNameRaw = strrep(S(i).name, 'Log_', 'Merged_');
    outName = fullfile(dataDir, outNameRaw);
    
    % Format: [Time(s), L, F]
    outData = [data(i).t_shifted, data(i).L_shifted, data(i).F_shifted];
    
    try
        writematrix(outData, outName, 'Delimiter', '\t');
        fprintf('  Saved %s\n', outNameRaw);
    catch ME
        fprintf('  Failed to save %s: %s\n', outNameRaw, ME.message);
    end
end

%% Plot Merged Traces vs Original Logs
figure(101);clf;
for i = 1:n
    figure(101);
    ax1_1 = nexttile(1);hold on;
    plot(data(i).t_shifted, data(i).F_shifted, 'LineWidth', 1);
    
    ax1_2 = nexttile(2);hold on;
    plot(data(i).t_shifted, data(i).L_shifted, 'LineWidth', 1);
    
    
    figure(i+10); clf;
    
    ax1 = subplot(2,1,1); hold on;
    % Plot merged trace (thin gray)
    plot(ax1, data(i).t_shifted, data(i).F_shifted, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    % Plot FULL original log on top (thick colored) — always visible
    plot(ax1, data(i).t_shifted_orig, data(i).F_shifted_orig, 'Color', colors(i,:), 'LineWidth', 2);
    ylabel(ax1, 'Force (kPa)');
    
    outNameRaw = strrep(S(i).name, 'Log_', 'Merged_');
    title_str = ['Dataset: ', strrep(strrep(outNameRaw, '_', ' '), '.txt', '')];
    title(ax1, title_str);
    
    ax2 = subplot(2,1,2); hold on;
    % Plot merged trace (thin gray)
    plot(ax2, data(i).t_shifted, data(i).L_shifted, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    % Plot FULL original log on top (thick colored)
    plot(ax2, data(i).t_shifted_orig, data(i).L_shifted_orig, 'Color', colors(i,:), 'LineWidth', 2);
    ylabel(ax2, 'Length (Lo)');
    xlabel(ax2, 'Time (s)');
    
    % Mark hi-res injection boundaries with vertical lines + label
    for jj = 1:length(hi_res_ranges{i})
        hr = hi_res_ranges{i}{jj};
        xline(ax1, hr.t_start, '--g', 'LineWidth', 1.5);
        xline(ax1, hr.t_end,   '--g', 'LineWidth', 1.5);
        xline(ax2, hr.t_start, '--g', 'LineWidth', 1.5);
        xline(ax2, hr.t_end,   '--g', 'LineWidth', 1.5);
        % Annotate the burst name at the top of F axis
        yl = ylim(ax1);
        text(ax1, (hr.t_start + hr.t_end)/2, yl(2), ...
            strrep(hr.name, '_', ' '), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 7, 'Color', [0 0.5 0], 'Interpreter', 'none');
    end
    
    legend(ax2, {'Merged (with Hi-Res)', 'Original Log'}, 'Location', 'best');
    linkaxes([ax1 ax2], 'x');
end
linkaxes([ax1_1 ax1_2], 'x');
%% compensate for run-down
figure(297);clf;
r828 = 80.09/71;
ax1 = nexttile(1);hold on;
plot(data(2).t_shifted, data(2).F_shifted, 'LineWidth', 1);
plot(data(3).t_shifted, data(3).F_shifted, 'LineWidth', 1);
plot(data(4).t_shifted, data(4).F_shifted*r828, 'LineWidth', 1);

ax2 = nexttile(2);hold on;
plot(data(2).t_shifted, data(2).L_shifted, 'LineWidth', 1);
plot(data(3).t_shifted, data(3).L_shifted, 'LineWidth', 1);
linkaxes([ax1 ax2], 'x');