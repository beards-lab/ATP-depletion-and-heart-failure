function features = extractSlackAttributes(data_t, data_y, data_SL, velocitytable, features, out, plotResults)
% EXTRACTSLACKATTRIBUTES  Extract force/SL features from a slack-release protocol.
%
%   FEATURES = EXTRACTSLACKATTRIBUTES(DATA_T, DATA_Y, DATA_SL, VELOCITYTABLE,
%                                      FEATURES, OUT, PLOTRESULTS)
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

    if nargin < 5
        features = [];
    end
    if nargin < 6        
        out = [];
    end
    if nargin < 7
        plotResults = false;
    end
    if plotResults
        if ~isempty(get(gcf, 'Children'))
            figure(406);clf;
        end
        hold on;
        plot(data_t, data_y, '-b', 'LineWidth', 1);                
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
        t_seg = velocity_segment(2); 

        win = data_t > velocity_segment(2) & data_t < velocity_segment(3);
        t = data_t(win); y = data_y(win);

        % cut off below FORCE_THRESHOLD_KPA (excludes near-zero force data)
        win = data_t > velocity_segment(2) & data_t < velocity_segment(3) & data_y > FORCE_THRESHOLD_KPA; % & data_t < t_seg + 0.048;
        t = data_t(win); y = data_y(win); SL = data_SL(win);
        t = t - t_seg;
    
        F0 = min(y);
        init_tail = 0:0.001:0.15;

        % ── Single-exponential fit ────────────────────────────────────────
        y_exp = @(A, k, t0, A0, x) max(0, A0*0 + F0 + A*(1-exp(-(x-t0)*k)));
        try
            [ae, ~] = fit(t, y, y_exp, 'StartPoint', [100, 50, 0.01, 0], 'Lower', [10, 0.01, 0.0, 0], 'Upper', [200 200 0.1, 100]);

            if plotResults
                plot(init_tail + t_seg, ae(init_tail), '--', t + t_seg, ae(t), 'LineWidth', 2);
                plot(ae.t0 + t_seg, 0, '*', 'LineWidth', 2, 'MarkerSize', ms);
            end

            feats.ktr    = ae.k;
            feats.A      = ae.A + F0;
            feats.t0     = ae.t0;
            feats.Am     = y(end);
            feats.SLslack = SL(end);
            feats.SLdiff  = max(data_SL) - SL(end);
        catch
            feats.ktr     = NaN;
            feats.A       = NaN;
            feats.t0      = NaN;
            feats.Am      = NaN;
            feats.SLslack = NaN;
            feats.SLdiff  = NaN;
        end

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

        % if plotResults
        %     nSlp = 100;
        %     plot(SL, y, SL(rngStart), y(rngStart),'x', SL(rngEnd), y(rngEnd), 'x', ...
        %         head(SL, nSlp), slpStart(1)*head(SL, nSlp)+slpStart(2), '-|', ...
        %         tail(SL, nSlp), slpEnd(1)*tail(SL, nSlp)+slpEnd(2), '-|', 'LineWidth', 1);
        %     xlabel('SL (um)');ylabel('Force (kPa)');
        % end

        % feats.Sl_V_restretch = (SL(end)-SL(1))/(t(end)-t(1));
        

        %%
        % 1. Detect Peaks
        [peak1_y, peak1_t] = findpeaks(y, t, 'MinPeakDistance', PEAK_MIN_DISTANCE_S, 'MinPeakProminence', PEAK_MIN_PROMINENCE);
        feats.v_restretch = velocity_segment(3, 2); 
        
        if ~isempty(peak1_y)
            % --- Peak Logic ---
            % Use the first detected peak
            p1_time_abs = peak1_t(1);
            
            feats.peak1_y = peak1_y(1);
            feats.peak1_t = p1_time_abs - t(1);
            
            % Find SL at peak time
            feats.peak1_SL = SL(find(t >= p1_time_abs, 1));
            feats.peak1_dSL = feats.peak1_SL - feats.SLslack;
            
            if plotResults
                plot(peak1_t(1) + t_seg, peak1_y(1), '*', MarkerSize=ms);
            end
        
            % --- Valley Logic (Min After Peak) ---
            % Find indices strictly after the first peak
            idx_after = t > p1_time_abs;
            
            if any(idx_after)
                y_after = y(idx_after);
                t_after = t(idx_after);
                
                % Calculate the minimum value in the segment after the peak
                [min_val, min_idx] = min(y_after);
                
                feats.vall_y = min_val;
                feats.vall_t = t_after(min_idx) - t(1); % Relative to start time
            else
                % Peak is at the very end of the signal; no valley possible
                feats.vall_y = NaN;
                feats.vall_t = NaN;
            end
        
        else
            % --- No Peak Found ---
            feats.peak1_y = NaN;
            feats.peak1_t = NaN;
            feats.peak1_SL = NaN;
            feats.peak1_dSL = NaN;
            
            % If no peak, there is no "valley after peak"
            feats.vall_y = NaN;
            feats.vall_t = NaN;
        end

        feats.peak2 = y(end);

        % steady state
        win = data_t >= velocity_segment(5) - 0.02 & data_t <= velocity_segment(5);
        t = data_t(win); y = data_y(win);
        t = t - t_seg;
        feats.steady = median(y);
        % plot(t, y, LineWidth=1.5);
        % plot([t(1) t(end)], [feats.steady feats.steady], LineWidth=3);

        % undershoot, overshoot and set - realtive to steady state
        win = data_t >= velocity_segment(4) & data_t <= velocity_segment(5);
        t = data_t(win); y = data_y(win);
        t = t - t_seg;
        % plot(t, y);
        [u_v, u_i] = min(y);
        feats.vall2_dy = u_v - feats.steady;
        feats.vall2_t = t(u_i) - t(1);
        
        
        % smooth data from 
        [o_v, o_i] = max(smoothdata(y(u_i:end), 1, "gaussian", 10));
        feats.ovrsht_dy = o_v - feats.steady;
        feats.ovrsht_t = t(o_i+u_i-1);
        % plot(t(u_i), u_v, '*', MarkerSize=ms)
        % % plot(t, smoothdata(y, 1, "gaussian", 10))
        % plot(feats.ovrsht_t, o_v, '*', MarkerSize=ms)

        if ~isempty(out)
            feats.XTOR = out.RTD(end);
        else
            % default "goal"
            feats.XTOR = 10;
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
            fn = fieldnames(feats);
            S = cell2struct(cell(size(fn)), fn, 1);  % struct with same fields, empty values
            features = repmat(S, 0, 1);              % 0x1 struct array
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
