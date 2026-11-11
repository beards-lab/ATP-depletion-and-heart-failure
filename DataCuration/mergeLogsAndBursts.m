function [data, hi_res_ranges, S] = mergeLogsAndBursts(dataDir, opts)
%MERGELOGSANDBURSTS Align the slow Log_*.txt traces and inject the hi-res bursts.
%
%   [data, hi_res_ranges, S] = mergeLogsAndBursts(dataDir, opts)
%
%   One protocol day produces two kinds of file:
%     * `NN_Log_<condition>.txt`  - the whole 100 s protocol at 10 ms (slow, low-res)
%     * `<atp>_<burst>.txt`       - a 1-3.5 s window at 50 us (fast, hi-res)
%   The two are recorded by separate acquisitions with independent clocks, so the
%   hi-res bursts have to be located inside the log time base before anything can
%   be extracted from them. This function does that, and hands back one merged
%   [t, L, F] trace per condition on a common time axis.
%
%   Alignment happens in three passes:
%     1. Coarse : t=0 at the first big L drop, then a piecewise shift onto the
%                 reference log using three protocol anchors (ktr drop, start of
%                 the staircase, first slack drop).
%     2. Injection: each burst file is located by xcorr of L against the log
%                 inside the burst's expected window, then spliced in, replacing
%                 the log samples it covers.
%     3. Fine  : a +/-20 ms grid search per injected range, minimising L MSE
%                 against the reference condition, so all conditions share one
%                 time base to sub-millisecond accuracy.
%
%   Inputs:
%     dataDir - Folder holding the Log_*.txt and burst .txt files (trailing
%               separator optional).
%     opts    - Optional struct. Fields (defaults in brackets):
%       rewriteData    [false] Write NN_Merged_*.txt (Log_ -> Merged_) to dataDir.
%                              Left off by default: these files are inputs to the
%                              protocol builders, so overwriting them silently
%                              re-bases every downstream .mat.
%       saveAllData    [false] Write AllDataMerged.mat (the returned `data`).
%       plotAll        [false] Per-burst alignment diagnostics.
%       plotFineAlign  [plotAll] Fine-alignment before/after overlays + MSE curves.
%       winSlack   [74.5 78.5] Log-time window used to locate the slack burst.
%                              Widen the left edge for protocols whose slack file
%                              opens with an extra event before the first slack
%                              (e.g. the 04-series FV2 ramp at ~74.4 s), otherwise
%                              that event falls outside the xcorr window.
%       winKtr       [69   73] Log-time window for the ktr / stiff burst.
%       winStretch   [72   74] Log-time window for the staircase stretch burst.
%       refIdx           [end] Index (into the name-sorted log list) of the
%                              condition every other condition is aligned to.
%
%   Outputs:
%     data          - 1xN struct array, one entry per log, name-sorted. Merged
%                     trace in t_fineShifted / L_fineShifted / F_fineShifted;
%                     the log-only trace is kept in *_shifted_orig for QC.
%     hi_res_ranges - 1xN cell; per log, a cell of structs {t_start,t_end,name}
%                     recording where each burst was spliced in.
%     S             - The name-sorted dir() listing of the log files.
%
%   Extracted from LoadAndPlotLogs.m so the same merge can be run for any
%   protocol day without editing a hardcoded path.
%
%   See also: LoadAndPlotLogs, CreateProtocolVelocityTable

    if nargin < 2; opts = struct(); end
    opts = iDefaults(opts, struct( ...
        'rewriteData',   false, ...
        'saveAllData',   false, ...
        'plotAll',       false, ...
        'winSlack',      [74.5, 78.5], ...
        'winKtr',        [69, 73], ...
        'winStretch',    [72, 74], ...
        'refIdx',        []));
    if ~isfield(opts, 'plotFineAlign'); opts.plotFineAlign = opts.plotAll; end

    if ~isempty(dataDir) && ~ismember(dataDir(end), '\/'); dataDir = [dataDir filesep]; end

    S = dir([dataDir '*Log_*.txt']);
    if isempty(S)
        error('mergeLogsAndBursts:noLogs', 'No *Log_*.txt files in %s', dataDir);
    end
    [~, idx] = sort({S.name});
    S = S(idx);
    n = length(S);
    if isempty(opts.refIdx); opts.refIdx = n; end

    %% 1. Load logs and detect the protocol anchors
    % Anchor events (in t-space after the initial-drop sync):
    %   t1: big L drop in 68-72 s        (ktr release)
    %   t2: start of the L step-up in 72-73 s (staircase)
    %   t3: big L drop in 73.5-75.5 s    (first slack; may be absent)
    data = struct('t', {}, 'F', {}, 'L', {}, 'anchors', {}, 't_shifted', {}, ...
                  'F_shifted', {}, 'L_shifted', {}, 'HiResOffset', {});
    for i = 1:n
        M = readmatrix(fullfile(S(i).folder, S(i).name), 'NumHeaderLines', 4);

        % Initial sync: t=0 at the first big L drop
        idx0 = find(diff(M(:,2)) < -0.01, 1);
        t = (M(:,1) - M(idx0,1)) / 1000;

        w1 = find(t >= 68 & t <= 72);
        i1 = find(diff(M(w1,2)) < -0.05, 1);
        t1 = t(w1(i1));

        w2 = find(t >= 72 & t <= 73);
        i2 = find(diff(M(w2,2)) > 0.003, 1);
        t2 = t(w2(i2));

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
    end

    %% 2. Coarse piecewise shift onto the reference log
    a_ref = data(opts.refIdx).anchors;
    mid1  = 71.5;   % between t1 (~70) and t2 (~72.5)
    mid2  = 73.8;   % between t2 (~72.5) and t3 (~74.5)

    for i = 1:n
        t_i = data(i).t;
        a_i = data(i).anchors;

        t_shifted = t_i;   % before 65 s the initial sync is good enough

        seg1 = t_i >= 65 & t_i < mid1;
        t_shifted(seg1) = t_i(seg1) + (a_ref(1) - a_i(1));

        seg2 = t_i >= mid1 & t_i < mid2;
        t_shifted(seg2) = t_i(seg2) + (a_ref(2) - a_i(2));

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
        ok  = [true; diff(t_shifted) > 0];

        data(i).t_shifted_orig = t_shifted(ok);
        data(i).F_shifted_orig = F_i(ok);
        data(i).L_shifted_orig = L_i(ok);

        data(i).t_shifted = t_shifted(ok);
        data(i).F_shifted = F_i(ok);
        data(i).L_shifted = L_i(ok);
    end

    %% 3. Locate and inject the hi-res burst files
    allFiles     = dir([dataDir '*.txt']);
    isLog        = contains({allFiles.name}, 'Log_');
    highResFiles = allFiles(~isLog);
    isMerged     = contains({highResFiles.name}, 'Merged_');
    highResFiles = highResFiles(~isMerged);

    fprintf('Aligning high-resolution burst files...\n');
    hi_res_ranges = cell(n, 1);
    for ii = 1:n; hi_res_ranges{ii} = {}; end

    for k = 1:length(highResFiles)
        fpath = fullfile(highResFiles(k).folder, highResFiles(k).name);
        try
            M_hi = readmatrix(fpath, 'NumHeaderLines', 4);
        catch
            continue;
        end
        t_hi = M_hi(:,1)/1000;   % ms to s
        L_hi = M_hi(:,2);
        F_hi = M_hi(:,3);

        % interp1 needs a strictly increasing sample grid
        [t_hi, ui] = unique(t_hi, 'stable');
        L_hi = L_hi(ui);
        F_hi = F_hi(ui);
        ok   = [true; diff(t_hi) > 0];
        t_hi = t_hi(ok);  L_hi = L_hi(ok);  F_hi = F_hi(ok);

        % Pair the burst with its log by explicit name matching
        name_hi      = highResFiles(k).name;
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
        if best_log_idx < 0
            fprintf('  Could not find explicit mapping for %s\n', name_hi);
            continue;
        end
        i = best_log_idx;

        % Burst type sets the log window the xcorr searches in
        if contains(name_hi, 'slack')
            t_win = opts.winSlack;      burst_type = 'slack';
        elseif contains(name_hi, 'stiff') || contains(name_hi, 'ktr')
            t_win = opts.winKtr;        burst_type = 'ktr';
        elseif contains(name_hi, 'stretch')
            t_win = opts.winStretch;    burst_type = 'stretch';
        else
            fprintf('  Unknown burst type for %s, skipping.\n', name_hi);
            continue;
        end

        t_sh    = data(i).t_shifted_orig;
        L_sh    = data(i).L_shifted_orig;
        win_idx = t_sh >= t_win(1) & t_sh <= t_win(2);
        t_log_win = t_sh(win_idx);
        L_log_win = L_sh(win_idx);
        if isempty(t_log_win)
            fprintf('  No log data in window [%.0f, %.0f]s for %s, skipping.\n', t_win(1), t_win(2), name_hi);
            continue;
        end

        % xcorr on mean-subtracted L, both resampled to 1 kHz
        dt = 0.001;
        t_grid_log = (min(t_log_win):dt:max(t_log_win))';
        L_log_1k   = interp1(t_log_win, L_log_win, t_grid_log, 'linear', 'extrap');
        t_grid_hi  = (min(t_hi):dt:max(t_hi))';
        L_hi_1k    = interp1(t_hi, L_hi, t_grid_hi, 'linear', 'extrap');

        [R, lags]    = xcorr(L_log_1k - mean(L_log_1k), L_hi_1k - mean(L_hi_1k));
        [~, i_max]   = max(R);
        time_offset  = lags(i_max) * dt + min(t_log_win) - min(t_hi);
        t_hi_aligned = t_hi + time_offset;

        fprintf('  Mapped %-30s -> %-30s | type: %-7s | offset: %.2fs | range: [%.1f, %.1f]s\n', ...
            name_hi, S(i).name, burst_type, time_offset, min(t_hi_aligned), max(t_hi_aligned));

        if opts.plotAll
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

        hi_res_ranges{i}{end+1} = struct('t_start', min(t_hi_aligned), ...
            't_end', max(t_hi_aligned), 'name', name_hi);

        % Splice: drop the log samples the burst covers, insert the burst
        keep  = data(i).t_shifted < min(t_hi_aligned) | data(i).t_shifted > max(t_hi_aligned);
        new_t = [data(i).t_shifted(keep); t_hi_aligned];
        new_L = [data(i).L_shifted(keep); L_hi];
        new_F = [data(i).F_shifted(keep); F_hi];

        [new_t, sortIdx] = sort(new_t);
        new_L = new_L(sortIdx);
        new_F = new_F(sortIdx);
        [new_t, uniqueIdx] = unique(new_t, 'stable');
        new_L = new_L(uniqueIdx);
        new_F = new_F(uniqueIdx);

        data(i).t_shifted   = new_t;
        data(i).L_shifted   = new_L;
        data(i).F_shifted   = new_F;
        data(i).HiResOffset = [data(i).HiResOffset; time_offset];
    end

    %% 4. Fine alignment of each injected range against the reference condition
    fprintf('\nPerforming fine-alignment on merged data using xcorr boundaries...\n');
    dt_grid   = 0.0005;   % 0.5 ms interpolation grid
    max_shift = 0.020;    % +/-20 ms search range
    shifts    = (-max_shift : dt_grid : max_shift)';

    if opts.plotFineAlign; figure(1001); clf; end

    for i = 1:n
        data(i).t_fineShifted = data(i).t_shifted;
        data(i).L_fineShifted = data(i).L_shifted;
        data(i).F_fineShifted = data(i).F_shifted;

        dts = [0, 0, 0];
        for jj = 1:length(hi_res_ranges{i})
            hr = hi_res_ranges{i}{jj};

            i_ref_win = data(opts.refIdx).t_shifted >= hr.t_start & data(opts.refIdx).t_shifted <= hr.t_end;
            i_hr_win  = data(i).t_shifted           >= hr.t_start & data(i).t_shifted           <= hr.t_end;
            if ~any(i_ref_win) || ~any(i_hr_win)
                fprintf('  i=%d jj=%d: no data in window [%.2f, %.2f], skipping\n', i, jj, hr.t_start, hr.t_end);
                continue;
            end

            t_ref = data(opts.refIdx).t_shifted(i_ref_win);
            L_ref = data(opts.refIdx).L_shifted(i_ref_win);
            t_cur = data(i).t_shifted(i_hr_win);
            L_cur = data(i).L_shifted(i_hr_win);

            t_grid  = (min(t_ref) : dt_grid : max(t_ref))';
            L_ref_g = interp1(t_ref, L_ref, t_grid, 'linear', 'extrap');

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

            if opts.plotFineAlign
                figure(200+i*10+jj); clf;
                subplot(2,1,1); hold on;
                plot(t_ref, L_ref, 'k', 'DisplayName', 'ref');
                plot(t_cur, L_cur, 'r', 'DisplayName', sprintf('cur (before, i=%d)', i));
                plot(data(i).t_fineShifted, data(i).L_fineShifted, 'b--', 'DisplayName', 'cur (after)');
                legend; ylabel('L (Lo)'); title(sprintf('Fine alignment i=%d jj=%d', i, jj));
                subplot(2,1,2);
                plot(shifts, mse); xlabel('shift (s)'); ylabel('MSE'); xline(dts(jj), 'r--');
            end
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

        if opts.plotFineAlign
            figure(1001);
            ax1 = subplot(2,1,1); hold on;
            plot(data(i).t_fineShifted, data(i).F_fineShifted, 'DisplayName', sprintf('Fine alignment i=%d', i));
            legend; ylabel('F'); title('Fine alignment ALL');
            ax2 = subplot(2,1,2); hold on;
            plot(data(i).t_fineShifted, data(i).L_fineShifted);
            xlabel('Time (s)'); ylabel('L');
            linkaxes([ax1 ax2], 'x');
        end
    end

    %% 5. Save
    fprintf('\nSaving merged outputs...\n');
    for i = 1:n
        outNameRaw = strrep(S(i).name, 'Log_', 'Merged_');
        outName    = fullfile(dataDir, outNameRaw);
        outData    = [data(i).t_fineShifted, data(i).L_fineShifted, data(i).F_fineShifted];  % [t(s), L, F]
        try
            if opts.rewriteData
                writematrix(outData, outName, 'Delimiter', '\t');
                fprintf('  Saved %s\n', outNameRaw);
            else
                fprintf('  Skipping writing %s due to rewriteData = false \n', outNameRaw);
            end
        catch ME
            fprintf('  Failed to save %s: %s\n', outNameRaw, ME.message);
        end
    end

    if opts.saveAllData
        save(fullfile(dataDir, 'AllDataMerged'), 'data');
        fprintf('  Saved AllDataMerged.mat\n');
    end
end

% ----------------------------------------------------------------------------
function o = iDefaults(o, d)
%IDEFAULTS Fill every field of d that o does not already define.
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(o, fn{i}) || isempty(o.(fn{i}))
            o.(fn{i}) = d.(fn{i});
        end
    end
end
