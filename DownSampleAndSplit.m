function [datatable, velocitytable] = DownSampleAndSplit(data_table, ts_d, ts_s, ML, dsf, scaleF, saveAs, offset)
% ts_d - time segment data for cost function
% ts_s - time segment simulation - broke by constant velocity segments
% dsf - downsample factor

% offset in ms
if nargin < 8
    offset = 0;
end

    data_table.Properties.VariableNames = {'Time', 'L', 'F'};
    data_table.Properties.VariableUnits = {'ms', 'Lo', 'kPa'};

    if isempty(ts_d)
        ts_d = [data_table.Time(1) data_table.Time(end)];
    end
    
    if isempty(ts_s)
        ts_s = [data_table.Time(1) data_table.Time(end)];
    end
    
    % Relabel and downsample 
    imin_d = find(data_table.Time >= ts_d(1), 1);
    imax_d = find(data_table.Time >= ts_d(end), 1);

    imin_s = find(data_table.Time >= ts_s(1), 1);
    imax_s = find(data_table.Time >= ts_s(end), 1);    

    t = data_table.Time(imin_s:imax_s);
    tf = data_table.Time(imin_d:imax_d);
    td = downsample(tf, dsf);
    l = data_table.L(imin_s:imax_s);
    lf = movmean(data_table.L(imin_d:imax_d),[dsf/2 dsf/2]); % l filtered
    % round to limit the oscillations
    lf = round(lf, 3);
    
    ld = downsample(lf, dsf);


    f = data_table.F(imin_s:imax_s)*scaleF;
    ff = movmean(data_table.F(imin_d:imax_d)*scaleF,[dsf/2 dsf/2]); % force filtered
    fd = downsample(ff, dsf);
    
    datatable = [td/1000  + offset/1000, ld*ML, fd];

    % Split it into segments with const velocities
    vs = [];
    pos = [data_table.L(find(t >= ts_s(1), 1))]*ML;
    for it = 1:(length(ts_s)-1)
        t1 = find(t >= ts_s(it), 1);
        t2 = find(t >= ts_s(it + 1), 1);
        vs(it) = round((l(t1) - l(t2))/(t(t1) - t(t2))*1000, 1);
        pos(it+1) = pos(it) + vs(it)*ML*(t(t2)-t(t1))/1000;
    end
    vsum = vs*ML; % 
    % time (s), velocity ML/s, velocity um/s
    velocitytable = [(ts_s + offset)/1000;[vs 0];[vsum 0];pos(1:end)]'; 
    

    if ~isempty(saveAs)
        % Export the data into modelica-readable format and for identificatoin
        fn = ['data/' saveAs '.mat'];
        save(fn,  'datatable', 'velocitytable');
        disp(['Saved as ' fn])
    end

%     figure();clf;
    subplot(211);
    hold on;title(saveAs, 'Interpreter', 'none');
    % yyaxis right; 
    plot(td + offset, ld, '-');
    
    
%     plot(t + offset, l, tf + offset, lf, td + offset, ld, '|-');   
    plot(velocitytable(:, 1)*1000, velocitytable(:, 4)/ML, 'x-', 'Linewidth', 1, 'MarkerSize', 10)
%     plot([ts_d;ts_d], repmat([min(data_table.L);max(data_table.L)], 1, length(ts_d)))
    xlabel('time (ms)');
    ylabel('Length (ML)')

    
    subplot(212);hold on;
    % yyaxis left;
    plot(td + offset, fd, '-', 'Linewidth', 2, 'MarkerSize', 10);
    
%     plot(t + offset, f, tf + offset, ff, td + offset, fd, '|-', 'Linewidth', 2, 'MarkerSize', 10);
%     plot([ts_d;ts_d], repmat([min(data_table.F);max(data_table.F)], 1, length(ts_d)))
    xlabel('time (ms)');
    ylabel('Force (kPa)')
    
end