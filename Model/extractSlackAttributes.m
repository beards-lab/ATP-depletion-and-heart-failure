function features = extractSlackAttributes(data_t, data_y, data_SL, velocitytable, features, out, plotResults, useSmoothing, datatable)
% EXTRACTSLACKATTRIBUTES  Extract force/SL features from a slack-release protocol.
%
%   FEATURES = EXTRACTSLACKATTRIBUTES(DATA_T, DATA_Y, DATA_SL, VELOCITYTABLE,
%                                      FEATURES, OUT, PLOTRESULTS, USESMOOTHING)
%   Segments the slack-release time series by velocity table entries, fits
%   exponential recovery curves, and extracts timing/amplitude features.
%
%   Inputs:
%     DATA_T        - time vector (s)
%     DATA_Y        - force vector (kPa)
%     DATA_SL       - sarcomere length vector (um)
%     VELOCITYTABLE - Nx2 matrix [time, velocity]; rows with velocity<0 mark
%                     the start of each slack-release segment
%     FEATURES      - struct to accumulate results into (optional, default [])
%     OUT           - evaluateModel output struct (optional)
%     PLOTRESULTS   - if true, plot intermediate fits (default false)
%     USESMOOTHING  - true to apply Gaussian smoothing before peak/valley
%                     detection (default false — appropriate for simulation
%                     output which is already smooth; set true for noisy
%                     experimental data).
%
%   Outputs:
%     FEATURES - struct with extracted features (peak force, timing, ktr, SL, etc.)
%
%   See also: fitSlackForceOnset, fitRecovery, evaluateModel

% Named constants
FORCE_THRESHOLD_KPA   = 10;    % minimum force (kPa) to include in active-phase window
PEAK_MIN_DISTANCE_S   = 0.01;  % minimum time between detected peaks (s)
PEAK_MIN_PROMINENCE   = 0.5;   % minimum peak prominence for findpeaks (kPa)
MARKER_SIZE           = 12;    % plot marker size
% Smoothing windows (data path only, useSmoothing=true).
% FAST:  Savitzky-Golay for the re-stretch peak section.  11 ms preserves
%   peaks wider than ~5 ms (the physiological rapid overshoot) while
%   removing sub-5 ms noise.  SG is preferred over Gaussian here because
%   it fits a local polynomial and does not round off peak heights.
% LIGHT: narrow movmedian used only to locate the undershoot minimum.
%   The dip is several kPa deep and tens of ms wide, well above the noise
%   floor, so a short window suffices to reject point noise without
%   smearing the dip into the much longer recovery tail.
% SLOW:  robust LOESS for the recovery/overshoot tail (from the located
%   minimum onward).  Wider window is safe there because the dynamics are
%   slow (50-200 ms); rloess handles occasional spikes without artefact.
SMOOTH_WIN_FAST_S     = 0.011;  % SG window for rapid-peak region   (≈11 ms)
SMOOTH_WIN_LIGHT_S    = 0.005;  % movmedian window for undershoot min (≈5 ms)
SMOOTH_WIN_SLOW_S     = 0.025;  % rloess window for overshoot tail    (≈25 ms)

    if nargin < 5
        features = [];
    end
    if nargin < 6        
        out = [];
    end
    if nargin < 7
        plotResults = false;
    end
    if nargin < 8 || isempty(useSmoothing)
        useSmoothing = false;
    end
    if nargin < 9
        datatable = [];
    end
    % (smoothing parameters are now named constants at the top of the function)
    if plotResults
        if ~isempty(get(gcf, 'Children'))
            figure(406); clf;
        end
        hold on;
        plot(data_t, data_y, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
        xlabel('Time (s)'); ylabel('Force (kPa)'); box on;
        title('Slack feature extraction');
    end


    feats = struct(); 
    if size(data_t, 1) < size(data_t, 2)
        data_t = data_t';
    end
    if size(data_y, 1) < size(data_y, 2)
        data_y = data_y';
    end
    if size(data_SL, 1) < size(data_SL, 2)
        data_SL = data_SL';
    end

    % Resample to 1 kHz if input is undersampled — ensures consistent
    % smoothing window sizes and peak detection behaviour regardless of
    % acquisition rate.
    TARGET_FS = 1000;  % Hz
    actual_fs = 1 / median(diff(data_t));
    if actual_fs < TARGET_FS * 0.99
        t_rs    = (data_t(1) : 1/TARGET_FS : data_t(end))';
        data_y  = interp1(data_t, data_y,  t_rs, 'linear');
        data_SL = interp1(data_t, data_SL, t_rs, 'linear');
        data_t  = t_rs;
    end

    % each set starts with negative velocity
    segments = find(velocitytable(:, 2) < -1)';

    ms = MARKER_SIZE;

    for i_seg =  1:length(segments)
        if segments(i_seg) + 4 > length(velocitytable)
            warning('Out of segments, out of determination..');
            break;
        end
        
        
        % set 1
        velocity_segment = velocitytable(segments(i_seg):segments(i_seg)+4, :);
        t_seg = velocity_segment(1); 
        feats.t_seg = t_seg;

        % win = data_t > velocity_segment(2) & data_t < velocity_segment(3);
        % t = data_t(win); y = data_y(win);

        % cut off below FORCE_THRESHOLD_KPA (excludes near-zero force data)
        win = data_t > velocity_segment(2) & data_t < velocity_segment(3) & data_y > FORCE_THRESHOLD_KPA; % & data_t < t_seg + 0.048;
        t = data_t(win); y = data_y(win); SL = data_SL(win);
        t = t - t_seg;
    
        % F0 = min(y);
        F0 = 0;
        init_tail = 0:0.001:0.15;

        % ── Fit 1: single exponential (original) ─────────────────────────
        % F(τ) = F0 + A·(1 − exp(−k·τ)),  τ = max(t−t0, 0)
        y_exp = @(A, k, t0, a, x) F0 + A .* (1 - exp(-k .* max(x - t0, 0))) + 0*a*(x-t0);
        try
            [ae, gof_e] = fit(t, y, y_exp, ...
                'StartPoint', [100, 50, 0.01, 0], ...
                'Lower',      [10,  0.01, 0.0, -100], ...
                'Upper',      [200, 500,  0.1 , 100]);

            if plotResults
                plot(init_tail + t_seg, ae(init_tail), '-', t + t_seg, ae(t), 'LineWidth', 1.5);
                plot(ae.t0 + t_seg, F0, '*', 'LineWidth', 2, 'MarkerSize', ms);
            end

            feats.ktr     = ae.k;
            feats.A       = ae.A + F0;
            feats.t0      = ae.t0;
            feats.ktr_rmse = gof_e.rmse;
        catch
            feats.ktr      = NaN;
            feats.A        = NaN;
            feats.t0       = NaN;
            feats.ktr_rmse = NaN;
        end

        % ── Fit 2: additive underdamped sinusoid ─────────────────────────
        % F(τ) = F0 + A·(1−exp(−k·τ)) + B·exp(−gam·τ)·sin(ω·τ)
        %   τ = max(t−t0, 0);  F(0)=F0, F(∞)=F0+A ✓
        %   B≥0: overshoot first (sin>0), then diminishing oscillations
        %   gam independent of k → oscillations can persist many cycles (gam ≪ k)
        %   To disable oscillation: set B lower bound to upper bound (both 0).
        if ~isnan(feats.ktr)
            sp_A = feats.A - F0;  sp_k = feats.ktr;  sp_t0 = feats.t0;
        else
            sp_A = max(y) - F0;   sp_k = 50;          sp_t0 = 0.01;
        end
        y_osc = @(A, k, B, gam, omega, t0, x) ...
            F0 + A .* (1 - exp(-k .* max(x - t0, 0))) ...
            + B .* exp(-gam .* max(x - t0, 0)) .* sin(omega .* max(x - t0, 0));
        try
            [ao, gof_o] = fit(t, y, y_osc, ...
                'StartPoint', [sp_A,  sp_k,  10,  sp_k/5, 20,  sp_t0], ...
                'Lower',      [10,    0.01,  0,   0.1,    1,   0.0  ], ...
                'Upper',      [200,   500,   100, 500,    300, 0.1  ]);

            if plotResults
                plot(init_tail + t_seg, ao(init_tail), '--', t + t_seg, ao(t), 'LineWidth', 1.5);
            end

            feats.ktr2_k     = ao.k;
            feats.ktr2_A     = ao.A + F0;
            feats.ktr2_B     = ao.B;
            feats.ktr2_gam   = ao.gam;
            feats.ktr2_omega = ao.omega;
            feats.ktr2_rmse  = gof_o.rmse;
            % Amplitude of first oscillation peak above recovery curve
            if ao.B > 0.5 && ao.omega > 0.5
                tau_peak = atan(ao.omega / ao.gam) / ao.omega;
                feats.ktr2_overshoot = ao.B * exp(-ao.gam * tau_peak) * ao.omega / sqrt(ao.gam^2 + ao.omega^2);
            else
                feats.ktr2_overshoot = 0;
            end
        catch
            feats.ktr2_k        = NaN;
            feats.ktr2_A        = NaN;
            feats.ktr2_B        = NaN;
            feats.ktr2_gam      = NaN;
            feats.ktr2_omega    = NaN;
            feats.ktr2_rmse     = NaN;
            feats.ktr2_overshoot = NaN;
        end

        % comparing data slack to model slack by leasts squares
        % this does not work
        if false && ~isempty(datatable)
            %%
            mdatay = datatable(:, 3); mdatat = datatable(:, 1);
            win = mdatat > velocity_segment(2) & mdatat < velocity_segment(3);
            mt = mdatat(win); my = mdatay(win); 
            mt = mt - t_seg;
            
            mys = smooth(my, length(my)/30);
            [mds]  = downsample([mt, mys], 1);
            mtds = mds(:, 1);myds = mds(:, 2);
            yds = interp1(t, y, mtds);
            plot(mt, my, mt, mys, mtds, myds, '|-', ...
                t, y, '-', mtds-tau_opt, yds, '|-', LineWidth=2);

% tic
            maxShiftTime = 0.10;  % 20 ms physiological bound — tune this
            shiftFun = @(tau) mean((interp1(mtds, yds, mtds + tau, 'linear', NaN) - myds).^2, 'omitnan');
            tau_opt = fminbnd(shiftFun, -maxShiftTime, maxShiftTime);
            cost = shiftFun(tau_opt);
            % toc
            % feats.ktrSqrDiff      = cost;
            % feats.ktrTimeShift    = tau_opt;            
            %%
            clf;
            st = (-10:100)/1000;
            ost = [];
            for ist = 1:length(st)
                ost(ist) = shiftFun(st(ist));
            end
            plot(st, ost)

            %%
            dsN = length(yds);
            sh = 1;
            shift = 1:dsN + sh;
            costs = (yds(max(1, min(dsN, shift))) - myds).^2;
            cost = nansum(costs)
            plot(mt, my, mt, mys, mtds, myds, '|-', ...
                t, y, '-', mtds, yds, '|-', mtds, costs, LineWidth=2);



        else
            feats.ktrSqrDiff = NaN;
        end


        if isempty(y)
            y = NaN;
            SL = NaN;
        end
        feats.Am      = y(end);
        feats.SLslack = SL(end);
        feats.SLdiff  = max(data_SL) - SL(end);

        % peaks and valley
        win = data_t >= velocity_segment(3) & data_t <= velocity_segment(4);
        t = data_t(win); y = data_y(win); SL = data_SL(win);
        t = t - t_seg;
        % plot(t, y, '-|b', 'LineWidth', 1);
        
        
        %% detect init slope stiffness
        
        rngStart = t > t(1) & t < t(1) + 1e-3;
        
        % assert there are at least 3 datapoitns
        rngStart(1:3) = true;
        slpStart = polyfit(SL(rngStart), y(rngStart), 1);
        feats.restretchSlopeStart = slpStart(1);

        rngEnd = t > t(end) - 1e-3 & t < t(end);
        rngEnd(end-5:end-1) = true;
        slpEnd = polyfit(SL(rngEnd), y(rngEnd), 1);
        feats.restretchSlopeEnd = slpEnd(1);

        if plotResults
            % Fitted slope lines evaluated at the SL points used, shown in time domain
            plot(t(rngStart) + t_seg, slpStart(1)*SL(rngStart) + slpStart(2), 'g-', 'LineWidth', 2.5);
            plot(t(rngEnd)   + t_seg, slpEnd(1)  *SL(rngEnd)   + slpEnd(2),   'm-', 'LineWidth', 2.5);
        end

        % feats.Sl_V_restretch = (SL(end)-SL(1))/(t(end)-t(1));
        

        % 1. Detect Peaks — always use raw signal for the restretch section.
        %    The rapid overshoot is only 2-5 ms wide; any smoothing window
        %    wide enough to suppress noise also smears the peak into the
        %    valley, making them indistinguishable.  Raw findpeaks with a
        %    prominence threshold is already reliable here.
        feats.v_restretch = velocity_segment(3, 2);
        y_sm = y;   % raw — smoothing not applied to this fast section
        [peak1_y, peak1_t] = findpeaks(y_sm, t, 'MinPeakDistance', PEAK_MIN_DISTANCE_S, 'MinPeakProminence', PEAK_MIN_PROMINENCE);

        if ~isempty(peak1_y)
            p1_time_abs = peak1_t(1);

            feats.peak1_y   = peak1_y(1);
            feats.peak1_t   = p1_time_abs - t(1);
            feats.peak1_SL  = SL(find(t >= p1_time_abs, 1));
            feats.peak1_dSL = feats.peak1_SL - feats.SLslack;

            if plotResults                
                plot(t + t_seg, y_sm, 'k-');
                plot(p1_time_abs + t_seg, peak1_y(1), 'b^', 'MarkerSize', ms, 'LineWidth', 2);
            end

            % Valley: min after peak on (possibly smoothed) signal
            idx_after = t > p1_time_abs;
            if any(idx_after)
                [min_val, min_idx_rel] = min(y_sm(idx_after));
                t_after = t(idx_after);
                feats.vall_y = min_val;
                feats.vall_t = t_after(min_idx_rel) - t(1);
                if plotResults
                    plot(t_after(min_idx_rel) + t_seg, min_val, 'bv', 'MarkerSize', ms, 'LineWidth', 2);
                end
            else
                feats.vall_y = NaN;
                feats.vall_t = NaN;
            end

        else
            feats.peak1_y   = NaN;
            feats.peak1_t   = NaN;
            feats.peak1_SL  = NaN;
            feats.peak1_dSL = NaN;
            feats.vall_y    = NaN;
            feats.vall_t    = NaN;
        end

        feats.peak2 = y(end);
        if plotResults
            plot(t(end) + t_seg, y(end), 'b>', 'MarkerSize', ms, 'LineWidth', 2);
        end

        % steady state
        win = data_t >= velocity_segment(5) - 0.02 & data_t <= velocity_segment(5);
        t = data_t(win); y = data_y(win);
        t = t - t_seg;
        if isempty(y)
            feats.steady = NaN;
        else
            feats.steady = median(y);
            if plotResults
                plot([t(1) + t_seg, t(end) + t_seg], [feats.steady, feats.steady], 'k--', 'LineWidth', 2);
            end
        end

        % undershoot, overshoot - relative to steady state
        win = data_t >= velocity_segment(4) & data_t <= velocity_segment(5);
        t = data_t(win); y = data_y(win);
        t = t - t_seg;

        % 1. Locate the undershoot minimum with light, narrow smoothing only.
        %    The dip is several kPa deep and tens of ms wide — well above the
        %    noise floor — so it doesn't need the wide trendline window; using
        %    the full SMOOTH_WIN_SLOW_S window here would blend the sharp early
        %    dip with the much longer, flatter recovery tail and bias/shift the
        %    detected minimum. A short movmedian only rejects single-sample
        %    spikes without smearing the dip's shape or location.
        if useSmoothing
            light_win = max(3, round(SMOOTH_WIN_LIGHT_S / median(diff(t))));
            y_light = smoothdata(y, 'movmedian', light_win);
        else
            y_light = y;
        end
        [u_v, u_i] = min(y_light);
        feats.vall2_dy = u_v - feats.steady;
        feats.vall2_t  = t(u_i) - t(1);

        % 2. Trendline-smooth only the tail FROM the minimum onward to find
        %    the overshoot — this is the slow recovery phase, where the wide
        %    robust-LOESS window belongs.  Cutting the smoothing here (rather
        %    than across the whole undershoot+recovery window) keeps the
        %    trendline from being dragged down by the sharp dip located above.
        tail = u_i:numel(y);
        if useSmoothing
            rl_win = max(5, round(SMOOTH_WIN_SLOW_S / median(diff(t(tail)))));
            y_sm2_tail = smoothdata(y(tail), 'rloess', rl_win);
        else
            y_sm2_tail = y(tail);
        end
        [o_v, o_i_rel] = max(y_sm2_tail);
        o_i = tail(1) + o_i_rel - 1;
        feats.ovrsht_dy = o_v - feats.steady;
        feats.ovrsht_t  = t(o_i);

        if plotResults
            plot(t(1:u_i) + t_seg, y_light(1:u_i), 'k:');
            plot(t(tail) + t_seg, y_sm2_tail, 'k-');
            plot(t(u_i) + t_seg, u_v, 'rv', 'MarkerSize', ms, 'LineWidth', 2);
            plot(t(o_i) + t_seg, o_v, 'r^', 'MarkerSize', ms, 'LineWidth', 2);
        end

        % ── Output-state features (from evaluateModel out struct) ─────────────
        if ~isempty(out) && isfield(out, 't')
            t_out  = out.t(:);
            t_on   = velocity_segment(2);    % slack onset time (absolute)
            t_rstr = velocity_segment(3);    % restretch start time

            % XTOR: RTD at final time point (steady state of last segment)
            feats.XTOR = out.RTD(end);

            % Pre-slack steady state window (last 20 ms before slack onset)
            mask_ss = t_out >= (t_on - 0.02) & t_out < t_on;
            if any(mask_ss)
                SRD_ss = zeros(size(out.SR(mask_ss)));
                if isfield(out, 'SRD'); SRD_ss = out.SRD(mask_ss); end
                feats.SRX_ss      = mean(out.SR(mask_ss) + SRD_ss, 'omitnan');
                p3_ss = 0;  % 3-state: include the post-stroke p3 population in attached fraction
                if isfield(out, 'p3_0') && ~isempty(out.p3_0); p3_ss = out.p3_0(mask_ss); end
                feats.attached_ss = mean(out.p1_0(mask_ss) + out.p2_0(mask_ss) + p3_ss, 'omitnan');
                feats.PT_ss       = mean(out.PuATP(mask_ss), 'omitnan');
            else
                feats.SRX_ss      = NaN;
                feats.attached_ss = NaN;
                feats.PT_ss       = NaN;
            end

            % Slack-phase features (between slack onset and restretch)
            % Use LXB = SL - LSE (contractile element length).
            mask_slack = t_out >= t_on & t_out <= t_rstr;
            if any(mask_slack) && isfield(out, 'LXB')
                LXB_sl = out.LXB(mask_slack);
                t_sl   = t_out(mask_slack);
                RTD_sl = out.RTD(mask_slack);

                % XTOR_vmax: RTD at time of maximum shortening speed during slack
                dLXB_dt = gradient(LXB_sl, t_sl);
                [~, idx_vmax]   = min(dLXB_dt);
                feats.XTOR_vmax = RTD_sl(idx_vmax);
            else
                feats.XTOR_vmax   = NaN;
            end

            % t0_crossing: time from shortening start to when LXB first crosses SL
            % (i.e., LSE → 0, force → 0). Extends search back to shortening onset.
            % we start from end of slack instead of start - that gives it enough time to buckle
            t_short_start = velocity_segment(2);  % row 2 col 1 of 5×4 segment matrix
            mask_cross = t_out >= t_short_start & t_out <= t_rstr;
            if any(mask_cross) && isfield(out, 'LXB') && isfield(out, 'SL')
                LXB_c = out.LXB(mask_cross);
                SL_c  = out.SL(mask_cross);
                t_c   = t_out(mask_cross);
                % First point where LXB ≥ SL (crossed from below: LSE ≤ 0, force ≤ 0)
                cross_idx = find(LXB_c <= SL_c, 1, 'first');
                if ~isempty(cross_idx)
                    feats.t0_crossing = t_c(cross_idx) - t_seg;
                else
                    feats.t0_crossing = NaN;
                end
            else
                feats.t0_crossing = NaN;
            end
        else
            feats.XTOR        = NaN;
            feats.SRX_ss      = NaN;
            feats.attached_ss = NaN;
            feats.PT_ss       = NaN;
            feats.XTOR_vmax   = NaN;
            feats.t0_crossing = NaN;
        end

        
    
        % get rid of empty, substitute with NaN instead    
        fn = fieldnames(feats);
        for i_elem = 1:length(fn)
            current_value = feats.(fn{i_elem});
            if isempty(current_value)
                feats.(fn{i_elem}) = NaN; % Replace empty with NaN
            end
        end
        
    
        if isempty(features)
            features = struct();   % plain struct; fields accumulate as row vectors below
        end

        % features(i_seg) = feats;
        % i_seg is implicit via (end+1)
        fieldNames = fieldnames(feats);
        for k = 1:numel(fieldNames)
            fname = fieldNames{k};
            if ~isfield(features, fname)
                features.(fname) = feats.(fname);
            elseif ~isempty(feats.(fname))
                features.(fname)(end+1) = feats.(fname)(1);
            else
                features.(fname)(end+1) = NaN;
            end
        end

    end

end
