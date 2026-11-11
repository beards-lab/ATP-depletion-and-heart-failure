% Process fmaxdata small step-up to determine detachment rates

%% % 1. Configuration - detachment
clf;    nexttile; hold on;
% 1. Create a logical mask for t >= 0

indices = [11, 7, 3];
span_width = 15e-3;
num_cells = length(fmaxdataCurrent);
num_ranges = length(indices);

% 2. Pre-allocate a table or matrix to store b values
% Rows = fmaxdata cell, Columns = velocitytable index
b_results = zeros(num_cells, num_ranges);

% 3. Loop through fmaxdata cells
for c = 1:num_cells
    % Extract and scale data for current cell
    % Multiplication [1 2.0 1] is applied here
    % raw_data = table2array(fmaxdataCurrent{c}(:, ["x_s_", "x_Lo_", "x_kPa_"]));
    raw_data = fmaxdataCurrent{c};
    raw_data = raw_data .* [1, 2.0, 1];
    
    data_t_full = raw_data(:, 1);
    data_f_full = raw_data(:, 3) / raw_data(1, 3);
    data_f_full = raw_data(:, 3);

    if c == 7
        plot(data_t_full, data_f_full, 'b', LineWidth=3);
    else
        plot(data_t_full, data_f_full, LineWidth=1);
    end
    % plot(out.t, out.SL, out.t, out.Force/out.Force(end), LineWidth=2);
    
    
    % Loop through the specific time ranges
    for r = 1:num_ranges
        base_time = velocitytable(indices(r), 1);
        curDropSpan = base_time + [0, span_width];
        
        % Selection logic
        sel = data_t_full >= curDropSpan(1) & data_t_full < curDropSpan(2);
        t_segment = data_t_full(sel);
        f_segment = data_f_full(sel);
        
        if isempty(t_segment), continue; end % Skip if no data in range
        
        % Fit setup
        t_shifted = t_segment - t_segment(1);
        ft = fittype('(a-c)*exp(-b*x) + c');
        opts = fitoptions(ft);
        opts.StartPoint = [max(f_segment), 1/span_width, min(f_segment)];
        
        % Fit and store b
        fr = fit(t_shifted, f_segment, ft, opts);
        plot(t_segment, fr(t_shifted), 'r-', 'LineWidth', 2);

        b_results(c, r) = fr.b;
    end
end



% Convert to table for a clear summary
b_table_detach = array2table(b_results, 'VariableNames', strcat('Index_', string(indices)));
% disp(b_table);

%% --- Visualization ---

% figure('Position', [100, 100, 1000, 400]);

% Subplot 1: Trend of Decay Rate (b) across ranges
nexttile(2);cla;
plot([1 2 3], b_results', '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
% set(gca, 'XDir','reverse'); % If indices represent a countdown/timeline
xlabel('Velocity Table Index');
ylabel('Decay Rate (b)');
title('Decay Rate Trend per Cell');
grid on;
legend(strcat('Cell ', string(1:num_cells)), 'Location', 'best');

% % Subplot 2: Heatmap of Decay Rates
% subplot(1, 2, 2);
% imagesc(b_results);
% colorbar;
% xticks(1:num_ranges);
% xticklabels(string(indices));
% yticks(1:num_cells);
% xlabel('Range Index');
% ylabel('fmaxdata Cell #');
% title('Heatmap of Decay Constants (b)');


%% Attachment

% c = 2;
% clf;    hold on;
nexttile(1); hold on;

% 1. Create a logical mask for t >= 0

indices = ([19, 13, 9, 5]);
span_width = [4*30e-3, 30e-3, 30e-3, 30e-3];
offset_width = [1e-3, 1e-3, 1e-3, 1e-3]*2;
num_cells = length(fmaxdataCurrent);
num_ranges = length(indices);

% 2. Pre-allocate a table or matrix to store b values
% Rows = fmaxdata cell, Columns = velocitytable index
b_results = zeros(num_cells, num_ranges);

% 3. Loop through fmaxdata cells
for c = 1:num_cells
% for c = [5 6]

    % raw_data = table2array(fmaxdataCurrent{c}(:, ["x_s_", "x_Lo_", "x_kPa_"]));
    raw_data = fmaxdataCurrent{c};
    raw_data = raw_data .* [1, 2.0, 1];
    
    data_t_full = raw_data(:, 1);
    % data_f_full = datatable(:, 3) / datatable(1, 3);
    data_f_full = raw_data(:, 3);
    
    if c == 2
        % plot(data_t_full, data_f_full, 'b', LineWidth=3);
    else
        % plot(data_t_full, data_f_full, LineWidth=1);
    end
    % plot(out.t, out.SL, out.t, out.Force/out.Force(end), LineWidth=2);
    
    
    % Loop through the specific time ranges
    for r = 1:num_ranges
        base_time = velocitytable(indices(r), 1);
        curDropSpan = base_time + [offset_width(r), span_width(r)];
        
        % Selection logic
        sel = data_t_full > curDropSpan(1) & data_t_full < curDropSpan(2);
        t_segment = data_t_full(sel);
        f_segment = data_f_full(sel);
        
        if isempty(t_segment), continue; end % Skip if no data in range
        
        % Fit setup
        t_shifted = t_segment - t_segment(1);
        ft = fittype('(a-c)*exp(-b*x) + c');
        opts = fitoptions(ft);
        opts.StartPoint = [min(f_segment), 1/span_width(r), max(f_segment)];
        
        % Fit and store b
        [fr frgof] = fit(t_shifted, f_segment, ft, opts);
        plot(t_segment, fr(t_shifted), 'r-', 'LineWidth', 2);
    
        b_results(c, r) = fr.b;
    end
end

% Convert to table for a clear summary
b_table_attach = array2table(b_results, 'VariableNames', strcat('Index_', string(indices)));
% disp(b_table);
% StatesInTime

%% --- Visualization ---

% figure('Position', [100, 100, 1000, 400]);
nexttile(2);hold on;

% Subplot 1: Trend of Decay Rate (b) across ranges
plot([1 2 3 4], b_results', '--x', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Velocity Table Index');
ylabel('Decay Rate (b)');
title('Decay Rate Trend per Cell');
grid on;
legend(strcat('Cell ', string(1:num_cells)), 'Location', 'best');
