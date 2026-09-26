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
PEAK_MIN_DISTANCE_S   = 0.001;  % minimum time between detected peaks (s)
% PEAK PROMINENCE. This used to be a fixed 0.1 kPa, which is 0.07-0.11 SIGMA of
% the recordings' own noise (1.36 / 1.32 / 0.91 kPa on 03/27 / 04/03 / 04/10).
% findpeaks then returned 11-15 "peaks" per re-stretch, almost all noise
% ripples, and the code below takes peak1_t(1) -- the FIRST in time, not the
% largest -- so on 7 of 15 preparation x segment cells it locked onto a ripple
% on the rising limb instead of the overshoot. The real overshoot has
% prominence 12-19 kPa; the ripples taken instead had 0.21-2.01 kPa. That put
% peak1_dSL at 0.5-1.8 SAMPLES of ramp travel on those cells and manufactured a
% non-monotonic "arc" in the pooled target that the model was then penalised for
% not reproducing. See Analyses/RestretchPeakAudit.
%
% The rule is now noise-relative, with an absolute floor as a fraction of the
% window's own force range so it still behaves on a SMOOTH model trace where
% sigma ~ 0. A threshold sweep gives identical results at 3, 4 and 6 sigma
% (between-preparation spread 0.641 -> 0.171), i.e. a plateau, so k = 3 sits
% safely inside it.
PEAK_PROM_K_SIGMA     = 3;      % prominence floor in units of the trace's noise sigma
PEAK_PROM_MIN_FRAC    = 0.05;   % ... and at least this fraction of the window's force range
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
% Guardrail scale factors (calibrated at optfull3_opt baseline to sum ~0.5 and ~1)
DOUBLEPEAK_SCALE      = 0.1;    % grading factor for spurious 2nd peak prominence penalty (×1: no 2nd peaks at baseline)
COOLDOWN_LS_SCALE     = 0.15; % scaling factor for cool-down overshoot LS error (calibrated to sum ≈1)
RESTRETCH_VALL_WIN_S  = 0.080;  % search span for the post-restretch minimum (vall2) that anchors the recovery fit

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

    % Noise level of THIS trace, for the peak-prominence rule above. High-pass
    % by subtracting a 9-point moving median, then take a MAD-based sigma so the
    % handful of large residuals at the ramp edges do not inflate it. A smooth
    % model trace gives sigma ~ 0, and the PEAK_PROM_MIN_FRAC floor takes over.
    hp = data_y(:) - movmedian(data_y(:), 9);
    noise_sigma = 1.4826 * median(abs(hp - median(hp)));

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

        % ── ktr80: the AUDITED hold-redevelopment rate ───────────────────
        % `ktr` above is a free-amplitude fit on the part of the hold where
        % force exceeds FORCE_THRESHOLD_KPA. That cut removes the smooth
        % model's fast early foot but hardly any of the noisy data's, so it is
        % biased 9-11 % on model traces against 1.6 % on data
        % (Analyses/KtrProtocolVsSlackKtr), and the free amplitude trades off
        % against k, so a tiny change in the trace moves k: refit6's best
        % scored ktr_shape 0.0827 in the campaign and 0.0978 on re-evaluation.
        % ktr80 is rsK80's construction applied to the hold: amplitude from
        % the hold's own plateau (13 time constants in, measured outside the
        % rise), rate from the 20-80 % band only, no threshold, no fit.
        winH = data_t > velocity_segment(2) & data_t < velocity_segment(3);
        [feats.ktr80, feats.ktrfit] = holdRate80(data_t(winH), data_y(winH));

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
        peak_prom = max(PEAK_PROM_K_SIGMA * noise_sigma, ...
                        PEAK_PROM_MIN_FRAC * (max(y_sm) - min(y_sm)));
        [peak1_y, peak1_t, ~, peak1_prom] = findpeaks(y_sm, t, 'MinPeakDistance', PEAK_MIN_DISTANCE_S, 'MinPeakProminence', peak_prom);

        if ~isempty(peak1_y)
            p1_time_abs = peak1_t(1);

            feats.peak1_y   = peak1_y(1);
            feats.peak1_t   = p1_time_abs - t(1);
            feats.peak1_SL  = SL(find(t >= p1_time_abs, 1));
            feats.peak1_dSL = feats.peak1_SL - feats.SLslack;

            % Guardrail: penalize spurious 2nd peak if present (graded by prominence)
            if numel(peak1_prom) >= 2
                feats.doublePeak = DOUBLEPEAK_SCALE * sum(peak1_prom(2:end));
            else
                feats.doublePeak = 0;
            end

            if plotResults
                plot(t + t_seg, y_sm, 'k-');
                plot(p1_time_abs + t_seg, peak1_y(1), 'b^', 'MarkerSize', ms, 'LineWidth', 2);
                if numel(peak1_prom) >= 2
                    plot(peak1_t(2:end) + t_seg, peak1_y(2:end), 'rx', 'MarkerSize', ms*2, 'LineWidth', 4);
                end
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
            feats.doublePeak = 0;
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

        % ── Post-restretch force redevelopment: first-order fit ───────────
        % The rate at which force recovers AFTER the re-stretch ramp has ended.
        % Scored as a RATE rather than by least squares on the trace: the
        % window is ~280 ms of which only the first ~30 ms carries the rate
        % information, so a trace-wise LSQE is dominated by the (already
        % correct) plateau and is nearly blind to a 2-3x rate error.
        %
        % Conventions match Analyses/RestretchVsKtrRecovery/recoveryWindows.m
        % so the values are directly comparable to that analysis' tables:
        %   baseline B = the vall2 minimum within RESTRETCH_VALL_WIN_S of the
        %                ramp end (the recovery starts there, NOT at the ramp
        %                end - the first ~7 ms are still falling),
        %   plateau  P = median over the last 15% of the window,
        %   amplitude A = P - B held FIXED, so the fitted k is a pure rate and
        %                does not absorb a force-scale error (that is already
        %                scored by vall2_dy / steady).
        % Only k and t0 are free -> the estimate is span-insensitive
        % (verified: data 48.5 s^-1 at a 30 ms span vs 48.4 at full span).
        % y_light (light movmedian on the data path, raw on the sim path) is used
        % only to LOCATE the anchor minimum; the fit itself is on raw y.
        %
        % rsK80/rsT20 are the AUDITED replacements (Analyses/RestretchMetricAudit).
        % Two things were wrong with rsK as a rate, and neither is the
        % exponential:
        %   (1) the amplitude came from the tail of the fit window, so anything
        %       that moved the window's END moved the rate -- truncating to
        %       60 ms shifted rsK by +77 %, and a known time-dilation was
        %       recovered to only 24 %. Sourcing the amplitude from `steady`
        %       (measured separately, at the reference length) instead drops
        %       those to 9 % and 8 %.
        %   (2) the whole-window least squares spends most of its leverage on
        %       the plateau. ln(4)/(t80-t20) uses only the part of the rise
        %       that carries rate information, and is toe- and plateau-free.
        % Together: between-preparation SD 17.2 % -> 6.1 %, 2 mM/8 mM
        % discriminability z 3.2 -> 6.8, and the value is invariant to the
        % crossing smoother (47.7/47.5/47.5 at 2/5/10 ms).
        [feats.rsK, feats.rsA, feats.rsT0, feats.rsR2, feats.rsK63, ...
         feats.rsK80, feats.rsT20, feats.rsKfit, feats.rsTfit] = ...
            fitRestretchRecovery(t, y, y_light, feats.steady, RESTRETCH_VALL_WIN_S);

        % Guardrail: least-squares error of cool-down window vs data
        if ~isempty(datatable)
            dwin = datatable(:,1) >= velocity_segment(4) & datatable(:,1) <= velocity_segment(5) & datatable(:,1) <= velocity_segment(4) + 0.15;
            td = datatable(dwin,1);  yd = datatable(dwin,3);
            if ~isempty(td)
                % Downsample data to 1/10th sampling, interpolate sim onto data grid
                td10 = td(1:10:end);  yd10 = yd(1:10:end);
                % Use only strictly increasing times (matching sim interpolation approach)
                nonrep = makeMonotonous(data_t);
                ysim = interp1(data_t(nonrep), data_y(nonrep), td10, 'linear', 'extrap');
                resid = yd10 - ysim;
                feats.coolDownLS = COOLDOWN_LS_SCALE * mean(resid.^2, 'omitnan');
            else
                feats.coolDownLS = 0;
            end
        else
            feats.coolDownLS = 0;
        end

        % Time-based LSQ error (RMSE) of the WHOLE slack cycle, model vs data,
        % over [slack onset, steady end]. RMSE (not MSE) so RMSE^2 in the
        % feature cost recovers the window MSE; data-vs-data self-error is 0.
        if ~isempty(datatable)
            sw = datatable(:, 1) >= velocity_segment(1) & datatable(:, 1) <= velocity_segment(5);
            tsd = datatable(sw, 1); ysd = datatable(sw, 3);
            if ~isempty(tsd)
                nonrepS = makeMonotonous(data_t);
                ysimS   = interp1(data_t(nonrepS), data_y(nonrepS), tsd, 'linear', 'extrap');
                feats.slackLSQE = sqrt(mean((ysd - ysimS).^2, 'omitnan'));
            else
                feats.slackLSQE = NaN;
            end
        else
            feats.slackLSQE = NaN;
        end

        if plotResults
            plot(t(1:u_i) + t_seg, y_light(1:u_i), 'k:');
            plot(t(tail) + t_seg, y_sm2_tail, 'k-');
            plot(t(u_i) + t_seg, u_v, 'rv', 'MarkerSize', ms, 'LineWidth', 2);
            plot(t(o_i) + t_seg, o_v, 'r^', 'MarkerSize', ms, 'LineWidth', 2);
            if exist('ysim', 'var')
                plot(td10, ysim, 'r-|', td10, yd10, 'b-|')
                plot(td10, COOLDOWN_LS_SCALE*cumsum(resid.^2), 'r-');
            end
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

        % A/steady -- the post-slack plateau as a FRACTION of the restretched
        % plateau. Scoring A and steady separately leaves this unconstrained,
        % and it is what reads optically: a model can sit inside both deadbands
        % and still have a visibly wrong slack plateau. It also deserves a much
        % tighter deadband than either parent, because the preparation's own
        % strength is common-mode and cancels in the ratio -- measured across
        % the three 8 mM preparations, rel SD 0.078 against 0.266 for A and
        % 0.236 for steady, i.e. 3.4x tighter
        % (Analyses/RefitStrategy/PlanRefit2.m).
        if isfield(feats,'A') && isfield(feats,'steady') && ...
                isscalar(feats.A) && isscalar(feats.steady) && feats.steady ~= 0
            feats.A_over_steady = feats.A / feats.steady;
        else
            feats.A_over_steady = NaN;
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

    % ktr_shape -- ktr across the ladder, divided by its own first segment.
    %
    % The five slack segments are not repeats: they are a ladder of increasing
    % release depth (SL held at 2.04 / 2.00 / 1.96 / 1.92 / 1.88 um), so `ktr`
    % is really ktr(SL). All three 8 mM preparations fall across it (-4.9 %,
    % -23.3 %, -15.2 %); the model rises. Scoring `ktr` alone cannot see that,
    % because its between-preparation spread is dominated by each heart's own
    % overall rate (rel SD 0.104 on the level against 0.041 on the shape), and
    % the resulting deadband is ~2.6x wider than the shape it is meant to
    % police. Dividing by segment 1 cancels the per-heart rate, exactly as
    % A_over_steady cancels the per-heart force.
    %
    % Element 1 is identically 1 by construction and carries no information --
    % score with a selector, e.g. 'ktr_shape[2:5]'.
    if isfield(features,'ktr') && numel(features.ktr) > 1 && ...
            isfinite(features.ktr(1)) && features.ktr(1) ~= 0
        features.ktr_shape = features.ktr(:)' / features.ktr(1);
    elseif isfield(features,'ktr')
        features.ktr_shape = nan(1, numel(features.ktr));
    end

    % ktr80_shape -- the ladder from the robust rate, divided by its MEAN over
    % the segments rather than by segment 1. Dividing by one segment puts that
    % segment's noise into every element at full weight (all four ktr_shape
    % elements moved together when ktr(1) moved); the mean spreads it at 1/n.
    % No element is identically 1, so all of them are scored.
    if isfield(features,'ktr80') && numel(features.ktr80) > 1
        mk = mean(features.ktr80, 'omitnan');
        if isfinite(mk) && mk > 0
            features.ktr80_shape = features.ktr80(:)' / mk;
        else
            features.ktr80_shape = nan(1, numel(features.ktr80));
        end
    end
    % the same ladder from the band-regression rate (bandRateFit)
    if isfield(features,'ktrfit') && numel(features.ktrfit) > 1
        mk = mean(features.ktrfit, 'omitnan');
        if isfinite(mk) && mk > 0
            features.ktrfit_shape = features.ktrfit(:)' / mk;
        else
            features.ktrfit_shape = nan(1, numel(features.ktrfit));
        end
    end

end


% =========================================================================
function [k, A, t0, r2, k63, k80, t20, kfit, tfit] = fitRestretchRecovery(t, y, y_anchor, steady, vallWin)
%FITRESTRETCHRECOVERY First-order fit of the force rise after a re-stretch.
%
%   [K, A, T0, R2, K63, K80, T20] = fitRestretchRecovery(T, Y, Y_ANCHOR, STEADY, VALLWIN)
%
%   T, Y     - the post-restretch window (t relative to the segment origin;
%              t(1) is the end of the re-stretch ramp).
%   Y_ANCHOR - same window, lightly smoothed on the data path (identical to Y
%              on the sim path). Used ONLY to locate the starting minimum, so
%              a single noise spike cannot set the baseline. The fit itself
%              always runs on the raw Y.
%   STEADY   - plateau force for this cycle (feats.steady), used only as a
%              fallback when the window is too short for its own plateau.
%   VALLWIN  - span (s) after the ramp end within which the recovery's
%              starting minimum (vall2) is located.
%
%   Returns the rate K (1/s) of  F = B + A*(1-exp(-k*max(t-t0,0)))  with the
%   baseline B (the vall2 minimum) and the amplitude A = P-B both FIXED, so
%   only K and T0 are fitted. Fixing A is what makes K span-insensitive and
%   keeps a force-scale error out of the rate.
%
%   K63 is the model-free crossing rate 1/(t63 - t_on), reported for
%   cross-checking against Analyses/RestretchVsKtrRecovery; it is a
%   piecewise-constant function of the parameters and so is a diagnostic,
%   not a fit target.
%
%   K80 and T20 are the AUDITED rate and delay (Analyses/RestretchMetricAudit):
%
%       A_ext = STEADY - B                  amplitude from OUTSIDE the window
%       z(t)  = (y - B) / A_ext             the shape, 0 -> 1
%       K80   = ln(4) / (t80 - t20)         the rate  (t_x = z first reaches x)
%       T20   = t20                         the delay
%
%   Two changes from K, both measured rather than argued:
%     * the amplitude no longer comes from the tail of the fit window, so the
%       rate no longer moves when the window does (60 ms truncation: +77 % on
%       K, exactly 0 on K80; recovery of a known time-dilation: 24 % error on
%       K, 6 % on K80);
%     * only the 20-80 % band of the rise is used, which is where the rate
%       information is -- the toe and the plateau are excluded by construction
%       instead of being fitted around with a free t0.
%   Result on the three 2026 preparations, segments 1-4: between-preparation
%   SD 17.2 % -> 6.1 %, 2 mM/8 mM discriminability z 3.2 -> 6.8.

    k = NaN; A = NaN; t0 = NaN; r2 = NaN; k63 = NaN; k80 = NaN; t20 = NaN;
    kfit = NaN; tfit = NaN;

    t = t(:); y = y(:); y_anchor = y_anchor(:);
    if numel(y_anchor) ~= numel(y); y_anchor = y; end
    ok = isfinite(t) & isfinite(y) & isfinite(y_anchor);
    t = t(ok); y = y(ok); y_anchor = y_anchor(ok);
    if numel(t) < 20 || ~(t(end) > t(1)); return; end

    % Baseline: the minimum within the first VALLWIN of the window. The
    % recovery starts there — the first few ms after the ramp are still falling.
    nWin = nnz(t <= t(1) + vallWin);      % contiguous prefix: t is ascending
    if nWin < 3; nWin = numel(t); end
    [B, iB] = min(y_anchor(1:nWin));      % smoothed value AND location, so a
                                          % single downward spike cannot set B

    tt = t(iB:end) - t(iB);
    yy = y(iB:end) - B;
    if numel(tt) < 20 || ~(tt(end) > 0); return; end

    % Plateau: median of the last 15% of the remaining window.
    P = median(yy(tt >= tt(end) - 0.15*(tt(end) - tt(1))));
    if ~isfinite(P); P = steady - B; end
    A = P;
    if ~isfinite(A) || A <= 0
        % No net rise (force fell and stayed down) — the rate is undefined.
        A = NaN; return;
    end

    obj = @(q) sum((yy - A .* (1 - exp(-abs(q(1)) .* max(tt - max(q(2), 0), 0)))).^2);
    q   = fminsearch(obj, [45, 0.004], ...
                     optimset('Display', 'off', 'MaxFunEvals', 4e3, 'MaxIter', 4e3));
    k  = abs(q(1));
    t0 = max(q(2), 0);

    pred = A .* (1 - exp(-k .* max(tt - t0, 0)));
    r2   = 1 - sum((yy - pred).^2) / max(sum((yy - mean(yy)).^2), eps);

    % Model-free crossing estimator (diagnostic).
    n    = max(3, round(0.002 / median(diff(tt))));
    yn   = movmedian(yy, n) / A;
    below = find(yn < 0.05);
    if isempty(below); i_on = 1; else; i_on = below(end); end
    i63 = find(yn >= 1 - exp(-1) & (1:numel(yn))' > i_on, 1);
    if ~isempty(i63) && tt(i63) > tt(i_on)
        k63 = 1 / (tt(i63) - tt(i_on));
    end

    % ── The audited rate: 20-80 % band, amplitude from `steady` ──────────
    % The crossings are read off a 5 ms movmedian (the same width the caller
    % uses to locate the minimum) so a single sample cannot set a crossing
    % time; the width does not matter much -- the value is 47.7 / 47.5 / 47.5
    % at 2 / 5 / 10 ms across the three 2026 preparations.
    Aext = steady - B;
    if isfinite(Aext) && Aext > 0
        nx = max(3, round(0.005 / median(diff(tt))));
        zx = movmedian(yy, nx) / Aext;
        t20 = firstCrossing(tt, zx, 0.20);
        t80 = firstCrossing(tt, zx, 0.80);
        if isfinite(t20) && isfinite(t80) && t80 > t20
            k80 = log(4) / (t80 - t20);
        end
        % The same rate and delay from EVERY sample in the 20-80 % band, not
        % two crossings: see bandRateFit. Same amplitude, same anchor.
        [kfit, tfit] = bandRateFit(tt, zx, 0.20, 0.80);
        % diagnostics only: inert unless a caller sets the global to a struct
        global RSK80_DIAG %#ok<GVMIS,TLEV>
        if isstruct(RSK80_DIAG)
            RSK80_DIAG(end+1).tB = t(iB); %#ok<NASGU>
            RSK80_DIAG(end).B = B; RSK80_DIAG(end).Aext = Aext;
            RSK80_DIAG(end).t20 = t20; RSK80_DIAG(end).t80 = t80; RSK80_DIAG(end).k80 = k80;
            RSK80_DIAG(end).tt = tt; RSK80_DIAG(end).zx = zx;
        end
    end
end

% =========================================================================
function [k80, kfit] = holdRate80(t, y)
%HOLDRATE80 Force-redevelopment rate over an isometric hold, 20-80 % band.
%
%   K80 = HOLDRATE80(T, Y) for the hold that follows a slack release:
%       B   = minimum of the smoothed trace in the first half of the hold
%       P   = median of the smoothed trace over its last HOLD_TAIL_S
%       z   = (y - B) / (P - B)
%       K80 = ln(4) / (t80 - t20),   t_x = first upward crossing of z = x
%             after the minimum
%   Same estimator for model and data. A 5 ms moving MEDIAN is used for the
%   crossings: on a monotone rise a median returns the centre sample, so it
%   removes noise without shifting the crossing times.
    HOLD_TAIL_S = 0.040;
    SMOOTH_S    = 0.005;
    k80 = NaN; kfit = NaN;
    t = t(:); y = y(:);
    ok = isfinite(t) & isfinite(y); t = t(ok); y = y(ok);
    if numel(t) < 40 || ~(t(end) > t(1)); return; end
    [t, iu] = unique(t); y = y(iu);                 % model output repeats stamps
    nx = max(3, round(SMOOTH_S / median(diff(t))));
    ys = movmedian(y, nx);
    half = find(t <= t(1) + 0.5*(t(end) - t(1)));
    [B, iB] = min(ys(half));
    P = median(ys(t >= t(end) - HOLD_TAIL_S));
    if ~isfinite(B) || ~isfinite(P) || P - B <= 0; return; end
    z = (ys(iB:end) - B) / (P - B);
    tt = t(iB:end);
    t20 = firstCrossing(tt, z, 0.20);
    t80 = firstCrossing(tt, z, 0.80);
    if isfinite(t20) && isfinite(t80) && t80 > t20
        k80 = log(4) / (t80 - t20);
    end
    kfit = bandRateFit(tt, z, 0.20, 0.80);
end

function [k, tlo] = bandRateFit(tt, z, lo, hi)
%BANDRATEFIT First-order rate (and delay) from every sample in a z-band.
%
%   For z = 1 - exp(-k (t - t0)),  ln(1 - z) = k*t0 - k*t  is EXACTLY linear
%   in t, so a regression of ln(1-z) on t over the band returns the same k as
%   ln((1-lo)/(1-hi)) / (t_hi - t_lo) on a clean exponential -- but it uses
%   every sample in the band instead of two crossing times.
%
%   Why that matters (Analyses/KtrSegmentSlope/RunNumericalNoise.m): the
%   model's recovery traces jitter ~1 % under a 1e-6 parameter change, and a
%   two-point rate inherits the jitter of both points in full. With the model
%   still 5-10 % from the data the cost noise is ~2*residual*jitter, which came
%   to 19-27 % of the rsK80 / rsT20 / ktr80_shape terms -- enough for the
%   search to bank a lucky evaluation (refit7 round 2: 0.3355 in memory,
%   0.3707 on reload). A regression averages the jitter away.
%
%   The band is the FIRST PASSAGE: from the first sample at z >= lo to the
%   first at z >= hi, so a later oscillation cannot re-enter the fit.
%   Weights (1-z)^2 undo the 1/(1-z) noise amplification of the log, i.e. the
%   fit is close to least squares in z itself.
%   TLO is the time the fitted curve crosses LO -- the delay, likewise from
%   the whole band.
    k = NaN; tlo = NaN;
    i1 = find(z >= lo, 1);
    i2 = find(z >= hi, 1);
    if isempty(i1) || isempty(i2) || i2 - i1 < 3; return; end
    idx = i1:i2;
    zz = min(z(idx), 1 - 1e-6);
    x = tt(idx); y = log(1 - zz); w = (1 - zz).^2;
    W = sum(w); xm = sum(w.*x)/W; ym = sum(w.*y)/W;
    sxx = sum(w.*(x - xm).^2);
    if ~(sxx > 0); return; end
    slope = sum(w.*(x - xm).*(y - ym)) / sxx;
    if ~(slope < 0); return; end
    k = -slope;
    a = ym - slope*xm;                            % y = a + slope*t
    tlo = (log(1 - lo) - a) / slope;
end

function tx = firstCrossing(tt, z, lev)
%FIRSTCROSSING first upward crossing of LEV, linearly interpolated. NaN if none.
    tx = NaN;
    i = find(z >= lev, 1);
    if isempty(i); return; end
    if i == 1; tx = tt(1); return; end
    z1 = z(i-1); z2 = z(i);
    if z2 == z1
        tx = tt(i);
    else
        tx = tt(i-1) + (lev - z1) / (z2 - z1) * (tt(i) - tt(i-1));
    end
end
