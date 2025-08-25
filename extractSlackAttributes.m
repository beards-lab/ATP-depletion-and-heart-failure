function features = extractSlackAttributes(data_t, data_y, data_SL, velocitytable)    
    features = [];
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
    segments = find(velocitytable(:, 2) < 0)';

    for i_seg =  1:length(segments)
            
        
        
        % set 1
        velocity_segment = velocitytable(segments(i_seg):segments(i_seg)+4, :);
        t_seg = velocity_segment(2); 

        win = data_t > velocity_segment(2) & data_t < velocity_segment(3);
        t = data_t(win); y = data_y(win);
        plot(t - t_seg, y, '-|b', 'LineWidth', 1.2);

        % cut off at 10kPa
        win = data_t > velocity_segment(2) & data_t < velocity_segment(3) & data_y > 10;
        t = data_t(win); y = data_y(win); SL = data_SL(win);
        t = t - t_seg;
    
        y_exp = @(A, k, t0, x) max(0, A*(1-exp(-(x-t0)*k)));
        [ae be] = fit(t, y, y_exp, 'StartPoint', [100, 50, 0.01], 'Lower', [10, 0.01, 0.0], 'Upper', [200 200 0.1]);
        init_tail = [0:0.001:0.05];
        plot(init_tail, ae(init_tail),'--', t, ae(t), LineWidth=2)
        plot(ae.t0, 0, '*', LineWidth=2, MarkerSize=8)
        
        feats.ktr = ae.k;
        feats.A = ae.A;
        feats.t0 = ae.t0;
        feats.Am = y(end);
        feats.SLslack = SL(end);
    
        % peaks and valley
        win = data_t >= velocity_segment(3) & data_t <= velocity_segment(4);
        t = data_t(win); y = data_y(win); SL = data_SL(win);
        t = t - t_seg;
        plot(t, y, '-|b', 'LineWidth', 1);
        [peak1_y, peak1_t] = findpeaks(y, t, MinPeakDistance=0.01, MinPeakProminence=2);
        feats.v_restretch = velocity_segment(3, 2);        
        feats.peak1_y = peak1_y(1);
        feats.peak1_t = peak1_t(1) - t(1);
        feats.peak1_SL = SL(find(t >= peak1_t(1), 1));
        feats.peak1_dSL = feats.peak1_SL - feats.SLslack;
        plot(peak1_t, peak1_y, '*');
        feats.peak2 = y(end);
        plot(t(end), feats.peak2, '*');
        [vall_y, vall_t] = findpeaks(-y, t, MinPeakDistance=0.01, MinPeakProminence=2);
        feats.vall_t = vall_t - t(1);
        feats.vall_y = -vall_y;
        plot(vall_t, -vall_y, '*');
                
        % steady state
        win = data_t >= velocity_segment(5) - 0.02 & data_t <= velocity_segment(5);
        t = data_t(win); y = data_y(win);
        t = t - t_seg;
        feats.steady = median(y);
        % plot(t, y, LineWidth=1.5);
        plot([t(1) t(end)], [feats.steady feats.steady], LineWidth=3);

        % undershoot, overshoot and set - realtive to steady state
        win = data_t >= velocity_segment(4) & data_t <= velocity_segment(5);
        t = data_t(win); y = data_y(win);
        t = t - t_seg;
        plot(t, y);
        [u_v, u_i] = min(y);
        feats.vall2_dy = u_v - feats.steady;
        feats.vall2_t = t(u_i) - t(1);
        
        plot(t(u_i), u_v, '*')
        % plot(t, smoothdata(y, 1, "gaussian", 10))
        
        % smooth data from 
        [o_v, o_i] = max(smoothdata(y(u_i:end), 1, "gaussian", 10));
        feats.ovrsht_dy = o_v - feats.steady;
        feats.ovrsht_t = t(o_i+u_i-1);
        plot(feats.ovrsht_t, o_v, '*')
        
    
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

        features(i_seg) = feats;
    end

end
