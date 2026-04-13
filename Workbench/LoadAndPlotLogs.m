%% Load and plot all Log_*.txt files from a data directory
close all; clear; clc;
rewriteData = false;
plotAll = false;

dataDir = '../data/03 27 2026 M/';
fpath = [dataDir '02_Merged_8mM_Active.txt'];
fpath = [dataDir '06_Merged_8mM_Active_PNB_Mava.txt'];
data_ref = readmatrix(fpath, 'NumHeaderLines', 4);


% dataset 1
dataDir = '../data/03 27 2026 M/';

% dataset 2
% dataDir = '../data/04 03 2026 F/';
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

data = struct('t', {}, 'F', {}, 'L', {}, 'anchors', {}, 't_shifted', {}, 'F_shifted', {}, 'L_shifted', {}, 'HiResOffset', {});
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
    if plotAll
        ax1 = nexttile(1);hold on;plot(data(i).t_shifted, data(i).F_shifted);
        ax2 = nexttile(2);hold on;plot(data(i).t_shifted, data(i).L_shifted);
        linkaxes([ax1 ax2], 'x');
    end

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

    if plotAll
        % figure(1);clf; nexttile(1);ax1 = plot(t_hi, F_hi);ax2 = nexttile(2);plot(t_hi, L_hi);title(highResFiles(k).name);linkaxes([ax1 ax2], 'x');
    end
    
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
            t_win = [74.5, 78.5];  % only the HIGH-force slack event (second one)
            burst_type = 'slack';
        elseif contains(name_hi, 'stiff') || contains(name_hi, 'ktr')
            t_win = [69, 73];  % ktr L-drop happens ~72s
            burst_type = 'ktr';
        elseif contains(name_hi, 'stretch')
            t_win = [72, 74];  % stretch L-drop happens ~74s
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
        
        if plotAll
        % --- Diagnostic plot: show alignment quality ---
            figure(100+k); clf;
            subplot(2,1,1); hold on;
            plot(data(i).t_shifted_orig, data(i).L_shifted_orig, 'b', 'LineWidth', 2);
            plot(t_hi_aligned, L_hi, 'r', 'LineWidth', 1);
            legend('Log (original)', ['Hi-res: ' name_hi], 'Interpreter', 'none');
            ylabel('L'); title(['Alignment: ' name_hi ' -> ' S(i).name ' (' burst_type ')'], 'Interpreter', 'none');
            xline(min(t_hi_aligned), '--r'); xline(max(t_hi_aligned), '--r');
            xlim(t_win + [-2, 2]);
    
            subplot(2,1,2); hold on;
            plot(data(i).t_shifted_orig, data(i).F_shifted_orig, 'b', 'LineWidth', 2);
            plot(t_hi_aligned, F_hi, 'r', 'LineWidth', 1);
            ylabel('F'); xlabel('Time (s)');
            xline(min(t_hi_aligned), '--r'); xline(max(t_hi_aligned), '--r');
            xlim(t_win + [-2, 2]);
        end
        
        
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

        data(i).HiResOffset = [data(i).HiResOffset; time_offset];

    else
        fprintf('  Could not find explicit mapping for %s\n', name_hi);
    end
end

%% Repeat Alignment pass with High-Res data (Fine Alignment)
fprintf('\nPerforming fine-alignment on merged data using xcorr boundaries...\n');

% last one is the reference
refIdx = n;
dt_grid   = 0.0005; % 0.5ms interpolation grid
max_shift = 0.020;  % ±20ms search range
shifts = (-max_shift : dt_grid : max_shift)';

plotAll = true;
if plotAll
    figure(1001); clf;
end
for i = 1:n
    data(i).t_fineShifted = data(i).t_shifted;
    data(i).L_fineShifted = data(i).L_shifted;
    data(i).F_fineShifted = data(i).F_shifted;
    % if i == refIdx; continue; end

    dts = [0,0,0];
    for jj = 1:length(hi_res_ranges{i})
        hr = hi_res_ranges{i}{jj};

        i_ref_win = data(refIdx).t_shifted >= hr.t_start & data(refIdx).t_shifted <= hr.t_end;
        i_hr_win  = data(i).t_shifted     >= hr.t_start & data(i).t_shifted     <= hr.t_end;

        if ~any(i_ref_win) || ~any(i_hr_win)
            fprintf('  i=%d jj=%d: no data in window [%.2f, %.2f], skipping\n', i, jj, hr.t_start, hr.t_end);
            continue;
        end

        t_ref = data(refIdx).t_shifted(i_ref_win);
        L_ref = data(refIdx).L_shifted(i_ref_win);
        t_cur = data(i).t_shifted(i_hr_win);
        L_cur = data(i).L_shifted(i_hr_win);

        % Interpolate reference onto uniform grid
        t_grid = (min(t_ref) : dt_grid : max(t_ref))';
        L_ref_g = interp1(t_ref, L_ref, t_grid, 'linear', 'extrap');

        % Grid search: shift cur by each candidate, measure MSE on overlap
        mse = nan(size(shifts));
        for k = 1:length(shifts)
            L_cur_sh = interp1(t_cur, L_cur, t_grid + shifts(k), 'linear', NaN);
            valid = ~isnan(L_cur_sh);
            if sum(valid) < 10; continue; end
            mse(k) = mean((L_ref_g(valid) - L_cur_sh(valid)).^2);
        end

        [~, best_k] = min(mse);
        dts(jj) = shifts(best_k);
        fprintf('  i=%d jj=%d [%.2f-%.2f s]: best shift = %.4f s (MSE=%.3e)\n', ...
            i, jj, hr.t_start, hr.t_end, dts(jj), mse(best_k));

        
        data(i).t_fineShifted(i_hr_win) = data(i).t_fineShifted(i_hr_win) - dts(jj);
        
        % Diagnostic: before/after overlay
        if plotAll
        figure(200+i*10+jj); clf;
        subplot(2,1,1); hold on;
        plot(t_ref, L_ref, 'k', 'DisplayName', 'ref');
        plot(t_cur, L_cur, 'r', 'DisplayName', sprintf('cur (before, i=%d)', i));
        % plot(data(i).t_fineShifted(i_hr_win), data(i).L_fineShifted(i_hr_win), 'b--', 'DisplayName', 'cur (after)');
        plot(data(i).t_fineShifted, data(i).L_fineShifted, 'b--', 'DisplayName', 'cur (after)');
        legend; ylabel('L (Lo)'); title(sprintf('Fine alignment i=%d jj=%d', i, jj));
        subplot(2,1,2);
        plot(shifts, mse); xlabel('shift (s)'); ylabel('MSE'); xline(dts(jj), 'r--');
        end

       
        % i_data_win = data(i).t_shifted >= hr.t_start & data(i).t_shifted <= hr.t_end;
        % i_shift_win = data(i).t_shifted >= hr.t_start + dts(jj) & data(i).t_shifted <= hr.t_end + dts(jj);
        % datanew_t(i_shift_win) = data(i).t_shifted(i_data_win);

    end

    % Ensure strict monotonicity
    [t_fineShifted, iu] = unique(data(i).t_fineShifted, 'stable');
    ok = [true; diff(t_fineShifted) > 0];      
    F_fineShifted = data(i).F_shifted(iu);
    L_fineShifted = data(i).L_shifted(iu);
    data(i).t_fineShifted = t_fineShifted(ok);
    data(i).F_fineShifted = F_fineShifted(ok);
    data(i).L_fineShifted = L_fineShifted(ok);

    fprintf('  %s fine shifts: ktr=%.3fs, stairs=%.3fs, slack=%.3fs\n', S(i).name, dts(1), dts(2), dts(3));
    if plotAll        
        figure(1001); 
        ax1 = subplot(2,1,1); hold on;
        plot(data(i).t_fineShifted, data(i).F_fineShifted, DisplayName=sprintf('Fine alignment i=%d', i));
        legend; ylabel('F'); title('Fine alignment ALL');
        ax2 = subplot(2,1,2);hold on;
        plot(data(i).t_fineShifted, data(i).L_fineShifted);
        xlabel('shift (s)'); ylabel('L');
        linkaxes([ax1 ax2], 'x');
    end
    
    % 
    % % Apply standard piecewise shift with these high-res offsets
    % t_new = data(i).t_shifted;
    % 
    % mid1 = 71.5;
    % mid2 = 73.8;
    % 
    % seg1 = t_new >= 65 & t_new < mid1;
    % t_new(seg1) = t_new(seg1) + dts(1);
    % 
    % seg2 = t_new >= mid1 & t_new < mid2;
    % t_new(seg2) = t_new(seg2) + dts(2);
    % 
    % seg3 = t_new >= mid2;
    % t_new(seg3) = t_new(seg3) + dts(3);
    % 
    % % Same safety cleanup as initial coarse piecewise shift
    % [t_new, iu] = unique(t_new, 'stable');
    % F_new = data(i).F_shifted(iu);
    % L_new = data(i).L_shifted(iu);
    % 
    % ok = [true; diff(t_new) > 0];
    % 
    % data(i).t_shifted = t_new(ok);
    % data(i).F_shifted = F_new(ok);
    % data(i).L_shifted = L_new(ok);
end

%% Save Merged Output Files
fprintf('\nSaving merged outputs...\n');
for i = 1:n
    outNameRaw = strrep(S(i).name, 'Log_', 'Merged_');
    outName = fullfile(dataDir, outNameRaw);
    
    % Format: [Time(s), L, F]
    outData = [data(i).t_fineShifted, data(i).L_fineShifted, data(i).F_fineShifted];
    
    try
        if rewriteData
            writematrix(outData, outName, 'Delimiter', '\t');
            fprintf('  Saved %s\n', outNameRaw);
        else
            fprintf('  Skipping writing %s due to rewriteData = false \n', outNameRaw);
        end
    catch ME
        fprintf('  Failed to save %s: %s\n', outNameRaw, ME.message);
    end
end

save(fullfile(dataDir, 'AllDataMerged'), 'data');

%% Plot Merged Traces vs Original Logs
figure(101);clf;
for i = 1:n
    figure(101);
    ax1_1 = nexttile(1);hold on;
    plot(data(i).t_fineShifted, data(i).F_fineShifted, 'LineWidth', 1);
    
    ax1_2 = nexttile(2);hold on;
    plot(data(i).t_fineShifted, data(i).L_fineShifted, 'LineWidth', 1);
    linkaxes([ax1_1 ax1_2], 'x');
    
    
    figure(i+10); clf;
    
    ax1 = subplot(2,1,1); hold on;
    % Plot merged trace (thin gray)
    plot(ax1, data(i).t_fineShifted, data(i).F_fineShifted, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    % Plot FULL original log on top (thick colored) — always visible
    plot(ax1, data(i).t_shifted_orig, data(i).F_shifted_orig, 'Color', colors(i,:), 'LineWidth', 2);
    ylabel(ax1, 'Force (kPa)');
    
    outNameRaw = strrep(S(i).name, 'Log_', 'Merged_');
    title_str = ['Dataset: ', strrep(strrep(outNameRaw, '_', ' '), '.txt', '')];
    title(ax1, title_str);
    
    ax2 = subplot(2,1,2); hold on;
    % Plot merged trace (thin gray)
    plot(ax2, data(i).t_fineShifted, data(i).L_fineShifted, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
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

return
%% graveyeard below

%% compensate for run-down


refData = data(2);
t_ref = refData.t_fineShifted;
F_ref = refData.F_fineShifted;
L_ref = refData.L_fineShifted;
t_ref = refData.t_shifted_orig;
F_ref = refData.F_shifted_orig;
L_ref = refData.L_shifted_orig;


% run-down data
rundownData = data(4);
F_rd = interp1(rundownData.t_fineShifted, rundownData.F_fineShifted, t_ref);
L_rd = interp1(rundownData.t_fineShifted, rundownData.L_fineShifted, t_ref);

% passive reference
F_pas = 0; %interp1(data(5).t_fineShifted, data(5).F_fineShifted, t_ref);
F_rdact = F_rd - F_pas;
F_refact = F_ref - F_pas;

L0 = 1.0;
k  = -0.6;   % SL slope: increase until high-SL portion of corrected trace aligns with F_refact
r0 = (68)/(56);
r0=1.214286;
% r828 = 80.09/71;
F_rdact_comp = scaleRundown(F_rdact, L_ref, r0, L0, k);
t_rd_comp = rundownData.t_fineShifted;
F_rd_comp = interp1(t_ref, F_rdact_comp + F_pas, t_rd_comp);
%%
figure(297);clf;
% plot(t_ref, F_refact, t_ref, F_rdact*r0, t_ref, F_rdact_comp,'LineWidth', 1);
ax1 = nexttile(1);plot(t_ref, F_refact + F_pas, t_rd_comp, F_rd_comp,'LineWidth', 1);
% ax2 = nexttile(2);plot(t_ref, L_ref, t_ref, L_rd,'LineWidth', 1);
% linkaxes([ax1, ax2], 'x')
xlim([66.6300   78.5624])
fprintf('With k = %g;r0=%f; cost = %g\n', k, r0, nansum(F_rdact_comp - F_refact));
legend('First', 'Third rescaled with f(L)');
xy_F = linspace(1, 100, 10);%[1:100];
xy_L = linspace(0.6, 1.2, 10);%:0.1:1.2];
sr = scaleRundown(xy_F, xy_L', r0, L0, k)';
% nexttile(2);cla;surf(xy_F, xy_L, sr);
% nexttile(2); plot(xy_F', scaleRundown(xy_F, L_ref, r0, L0, k)');
% plot(t_ref, F_refact, t_ref, F_rdact)
% F_rdact
%% rundown corrig
t = refData.t_shifted;
F = refData.F_shifted;
L = refData.L_shifted;
zones = [71.5, 72.4;77.4, 78.4];
% at 2.0
mask = any(t >= zones(:,1)' & t <= zones(:,2)', 2);
p1 = polyfit(t(mask), F(mask), 1);

% at 2.2
L_smooth = movmean(L, 20, 'omitnan');
falling = find(diff(L_smooth >= 1.1) == -1) + 1;  % first sample below threshold
t_drops = t(falling);
mask = false(size(t));
for i = 1:numel(t_drops)
    mask = mask | (t >= t_drops(i) - 0.15 & t < t_drops(i));
end
p2 = polyfit(t(mask), F(mask), 1);
clf;
plot(t, F, t(mask), F(mask), t, polyval(p1, t), t, polyval(p2, t), linewidth = 1)

%%

%%
% plot(rundownData.t_shifted, rundownData.F_shifted*r828, 'LineWidth', 1);

ax2 = nexttile(2);hold on;
plot(refData.t_shifted, refData.L_shifted, 'LineWidth', 1);
plot(data(3).t_shifted, data(3).L_shifted, 'LineWidth', 1);

%% Report Figure A — Comparison 2, upscaled-4, 3, 5
clr2 = [0.00, 0.45, 0.70];  % blue:   8mM Active
clr3 = [0.84, 0.10, 0.11];  % red:    2mM Active
clr4 = [0.49, 0.18, 0.56];  % purple: 8mM repeat upscaled
clr5 = [0.20, 0.63, 0.17];  % green:  PNB+Mava reference

% Upscale i=4 to match i=2 using scaleRundown
t_ref2 = refData.t_fineShifted;
F_ref2 = refData.F_fineShifted;
L_ref2 = refData.L_fineShifted;
F_rd4  = interp1(rundownData.t_fineShifted, rundownData.F_fineShifted, t_ref2, 'linear', NaN);
F_pas5 = interp1(data(5).t_fineShifted, data(5).F_fineShifted, t_ref2, 'linear', NaN);
F_rdact4      = F_rd4 - F_pas5;
r0_up = 68/56;  L0_up = 1.0;  k_up = -0.8;
F_rdact4_comp = scaleRundown(F_rdact4, L_ref2, r0_up, L0_up, k_up);
F_4_upscaled  = F_rdact4_comp + F_pas5;  % total force, rundown-compensated

zoom_xlims  = {[66, 79], [69.5, 74], [74, 77.5], [66, 79]};
zoom_titles = {'Overview (66–79 s)', 'ktr (69.5–74 s)', 'Slack (74–77.5 s)', 'Length protocol'};

fA = figure('Name', 'Report_Comparison_235', 'Units', 'centimeters', 'Position', [2 2 20 24]);
tlA = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
axA = gobjects(4, 1);

for ii = 1:4
    axA(ii) = nexttile(ii); hold on;
    if ii < 4
        plot(refData.t_fineShifted, refData.F_fineShifted, 'Color', clr2, 'LineWidth', 1.2, 'DisplayName', '8mM Active (i=2)');
        plot(t_ref2, F_4_upscaled,                         '--', 'Color', clr4, 'LineWidth', 1.0, 'DisplayName', '8mM repeat upscaled (i=4)');
        plot(data(3).t_fineShifted, data(3).F_fineShifted, 'Color', clr3, 'LineWidth', 1.2, 'DisplayName', '2mM Active (i=3)');
        plot(data(5).t_fineShifted, data(5).F_fineShifted, 'Color', clr5, 'LineWidth', 1.2, 'DisplayName', 'PNB+Mava (i=5)');
        ylabel('Force (kPa)');
        if ii == 1; legend('Location', 'northwest', 'FontSize', 7); end
    else
        plot(data(5).t_fineShifted, data(5).L_fineShifted, 'Color', clr5, 'LineWidth', 1.2);
        ylabel('Length (Lo)'); xlabel('Time (s)');
    end
    xlim(zoom_xlims{ii});
    title(zoom_titles{ii});
    xline(72.0, '--k', 'Alpha', 0.3, 'HandleVisibility', 'off');
    xline(74.5, '--k', 'Alpha', 0.3, 'HandleVisibility', 'off');
    box on;
end

title(tlA, '03/27/2026 — Force comparison: 8mM vs 2mM vs PNB+Mava');
exportgraphics(fA, fullfile(dataDir, 'Report_Comparison_235.png'), 'Resolution', 300);
fprintf('Saved Report_Comparison_235.png\n');

%% Report Figure B — Rundown linear fits i=2 and i=3
ds_order = [3, 2];
ds_labels = {'2mM Active (i=3)', '8mM Active (i=2)'};

fB = figure('Name', 'Report_Rundown_23', 'Units', 'centimeters', 'Position', [2 2 18 16]);

for ii_loop = 1:2
    ds_idx = ds_order(ii_loop);
    t = data(ds_idx).t_fineShifted;
    F = data(ds_idx).F_fineShifted;
    L = data(ds_idx).L_fineShifted;

    % Mask 1: SL=2.0 plateau zones (stable force before each L-drop)
    zones = [71.5, 72.4; 77.4, 78.4];
    mask1 = any(t >= zones(:,1)' & t <= zones(:,2)', 2);
    p1 = polyfit(t(mask1), F(mask1), 1);

    % Mask 2: pre-drop points at SL=2.2 (L >= 1.1, 150ms window before each drop)
    L_smooth = movmean(L, 20, 'omitnan');
    falling  = find(diff(L_smooth >= 1.1) == -1) + 1;
    t_drops  = t(falling);
    mask2 = false(size(t));
    for jj = 1:numel(t_drops)
        mask2 = mask2 | (t >= t_drops(jj) - 0.15 & t < t_drops(jj));
    end
    p2 = polyfit(t(mask2), F(mask2), 1);

    t_line = linspace(71, 79, 500);

    subplot(2, 1, ii_loop); hold on;
    win = t >= 70 & t <= 79.5;
    plot(t(win), F(win), 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'DisplayName', 'Raw trace');
    plot(t(mask1), F(mask1), 'o', 'Color', [0.2 0.2 0.8], 'MarkerSize', 3, 'LineStyle', 'none', 'DisplayName', 'SL=2.0 pts');
    plot(t_line, polyval(p1, t_line), '--', 'Color', [0.2 0.2 0.8], 'LineWidth', 1.5, 'DisplayName', sprintf('SL=2.0: %.3f kPa/s', p1(1)));
    plot(t(mask2), F(mask2), 's', 'Color', [0.8 0.2 0.2], 'MarkerSize', 3, 'LineStyle', 'none', 'DisplayName', 'SL=2.2 pts');
    plot(t_line, polyval(p2, t_line), ':', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5, 'DisplayName', sprintf('SL=2.2: %.3f kPa/s', p2(1)));
    xlim([70 79.5]); ylabel('Force (kPa)'); xlabel('Time (s)');
    title(sprintf('%s — Rundown', ds_labels{ii_loop}));
    legend('Location', 'northeast', 'FontSize', 7); box on;
    yl = ylim;
    text(70.2, yl(2) - 0.05*(yl(2)-yl(1)), sprintf('slope_{2.0} = %.3f kPa/s', p1(1)), ...
        'Color', [0.2 0.2 0.8], 'FontSize', 8, 'VerticalAlignment', 'top');
    text(70.2, yl(2) - 0.15*(yl(2)-yl(1)), sprintf('slope_{2.2} = %.3f kPa/s', p2(1)), ...
        'Color', [0.8 0.2 0.2], 'FontSize', 8, 'VerticalAlignment', 'top');
end

sgtitle('Rundown fits: 2mM Active vs 8mM Active');
exportgraphics(fB, fullfile(dataDir, 'Report_Rundown_23.png'), 'Resolution', 300);
fprintf('Saved Report_Rundown_23.png\n');

% -------------------------------------------------------------------------
function F_corrected = scaleRundown(F_rdact, L_ref, r0, L0, k, f0)
% scaleRundown  Scale run-down active force with SL-dependent compensation.
%
%   F_corrected = scaleRundown(F_rdact, L_ref, r0, L0, k)
%
%   Inputs:
%     F_rdact - run-down active force vector (same length as L_ref)
%     L_ref   - sarcomere length vector (same time axis as F_rdact)
%     r0      - base scale factor at pivot length L0 (e.g. 75/66)
%     L0      - pivot sarcomere length (e.g. mean(L_ref))
%     k       - SL slope [1/length_unit]: positive = more compensation at higher SL
%
%   The scale factor varies linearly with SL:
%     scale(L) = r0 + k*(L - L0)
%
%   Start with k=0 to recover flat-ratio behaviour, then increase k until
%   the corrected trace aligns with the reference active force.

    scale = r0 + k .* (L_ref - L0);
    F_corrected = (F_rdact) .* scale;
end