% Escape underscores for LaTeX/TeX interpreters
fn_latex = strrep(params0.fn, '_', '\_');
%% 2. Calculate Sensitivity
% Difference matrix (Gradient proxy)
sensitivityMatrix = featureMatrixPlusDiff - featureMatrixMinusDiff;

% Handle "dead" features (Zero variance check)
std_devs = std(sensitivityMatrix, 'omitnan');
is_dead = std_devs == 0;
if any(is_dead)
    warning('Features %s have zero sensitivity and were excluded.', strjoin(fn_latex(is_dead), ', '));
    sensitivityMatrix(:, is_dead) = [];
    fn_latex(is_dead) = [];
end

% drop insensitive params - rows
parIsDead = sum(sensitivityMatrix, 2, 'omitnan') == 0 | isnan(sum(sensitivityMatrix, 2));

sensitivityMatrix = sensitivityMatrix(not(parIsDead), :);
%%
% heatmap(sensitivityMatrix)
heatmap(sensitivityMatrix)
%% 3. Calculate Correlation
% Correlation of the gradients
featureCorrelations = corr(sensitivityMatrix);
%% test corr
cx = sensitivityMatrix(:, 9);
cy = sensitivityMatrix(:, 5); 
corrcoef(cx, cy)
scatter(cx, cy);


%% 4. Visualization: Standard Heatmap (Figure 1)
figure(1); clf;
set(gcf, 'Color', 'w', 'Name', 'Correlation Matrix');

h = heatmap(fn_latex, fn_latex, featureCorrelations);
h.Title = 'Feature Trade-off Analysis (Correlation)';
h.XLabel = 'Features';
h.YLabel = 'Features';
h.Colormap = jet;
h.ColorLimits = [-1 1]; % Fixed scale from -1 to +1

%% 5. Visualization: Clustergram (Figure 301)
% Clustergram creates its own UI, so we target it to Fig 301 by handle
    
    % Create the clustergram
    cg = clustergram(featureCorrelations, ...
        'RowLabels', fn_latex, ...
        'ColumnLabels', fn_latex, ...
        'Symmetric', true, ...
        'Colormap', jet, 'Standardize', 'None', 'Annotate', true);
    %, ...
     %   'DisplayRange', 1);
    
    addTitle(cg, 'Hierarchical Clustering of Model Costs');
    
    

%% 6. Automated Report
% Thresholds
anti_thresh = -0.7; % Strong negative correlation
corr_thresh = 0.8;  % Strong positive correlation

fprintf('\n========================================\n');
fprintf('       SENSITIVITY ANALYSIS REPORT      \n');
fprintf('========================================\n');

% 1. Detect Trade-offs (Anticorrelated)
fprintf('\n--- [CRITICAL] TRADE-OFFS (Anticorrelated) ---\n');
fprintf('Action: You likely cannot minimize both simultaneously.\n');
[row, col] = find(tril(featureCorrelations, -1) < anti_thresh);

if isempty(row)
    disp('  No critical trade-offs found.');
else
    for i = 1:length(row)
        fprintf('  [Trade-off] %-15s vs %-15s (Corr: %5.2f)\n', ...
            fn_latex{row(i)}, fn_latex{col(i)}, featureCorrelations(row(i), col(i)));
    end
end

% 2. Detect Redundancies (Correlated)
fprintf('\n--- REDUNDANCIES (Highly Correlated) ---\n');
fprintf('Action: These features are physically linked.\n');
[row_p, col_p] = find(tril(featureCorrelations, -1) > corr_thresh);

if isempty(row_p)
    disp('  No strong redundancies found.');
else
    for i = 1:length(row_p)
        fprintf('  [Linked]    %-15s vs %-15s (Corr: %5.2f)\n', ...
            fn_latex{row_p(i)}, fn_latex{col_p(i)}, featureCorrelations(row_p(i), col_p(i)));
    end
end
fprintf('\n');


%%

%% 1. Setup Data
% Assume optim_history is 5000 x N
% fn is your cell array of feature names
fn_latex = strrep(fn, '_', '\_'); 

% Standardize the data
% Since different costs have different units/scales, we must use 
% correlation (which is scale-invariant) or Z-score the data.
optim_history_total = sum(optim_history, 2);
history_corr = corr(optim_history(optim_history_total < 10, :), 'Rows', 'pairwise');

% 2. Visualization (Figure 1: Heatmap)
figure(1); clf;
set(gcf, 'Color', 'w');

h = heatmap(fn_latex, fn_latex, history_corr);
h.Title = 'Optimization History: Feature Inter-Dependencies';
h.Colormap = jet;
h.ColorLimits = [-1 1];

% Fix for heatmap labels
set(struct(h).NodeChildren(3), 'TickLabelInterpreter', 'tex');

    % Create the clustergram
    cg = clustergram(history_corr, ...
        'RowLabels', fn_latex, ...
        'ColumnLabels', fn_latex, ...
        'Symmetric', true, ...
        'Colormap', jet, ...
        'DisplayRange', 1);

% 3. Conflict Analysis (Who is the "troublemaker"?)
% A "Conflict Score" measures how much a feature fights others.
% We sum the absolute values of negative correlations.
conflicts = history_corr;
conflicts(conflicts > 0) = 0; % Keep only negative correlations
conflict_scores = sum(abs(conflicts), 1);

[sorted_scores, idx] = sort(conflict_scores, 'descend');

fprintf('\n==================================================\n');
fprintf('   OPTIMIZATION CONFLICT REPORT (Ranked by Trade-offs)\n');
fprintf('==================================================\n');
for i = 1:length(idx)
    fprintf('%02d. %-20s | Conflict Score: %.2f\n', ...
        i, fn{idx(i)}, sorted_scores(i));
end


%% ordered heatmap from sluctergram
% 1. Create the Correlation Matrix

% 3. Extract the new orders
% These are the indices of the features as they appear from top-to-bottom/left-to-right
cgRows = get(cg, 'RowLabels');
cgCols = get(cg, 'ColumnLabels');

% 3. Map the names back to original indices to find the order
[~, rowIdx] = ismember(cgRows, fn_latex);
[~, colIdx] = ismember(cgCols, fn_latex);

% 4. Now recreate the matrix with values in boxes
clustered_C = featureCorrelations(rowIdx, colIdx);

figure('Color', 'w');clf;
h_map = heatmap(cgCols, cgRows, clustered_C);
h_map.Colormap = jet;
% h_map.ColorLimits = [-1 1];
% h_map.LabelDisplayData = 'all'; % This puts the numbers in the boxes

%% 4. Create the Reordered Correlation Matrix
clusteredCorr = sensitivityMatrix(rowIdx, colIdx);
clusteredFeatNames = featNames(rowIdx);

% 5. Display the grouped values in a readable Table format
GroupedValues = array2table(clusteredCorr, ...
    'VariableNames', clusteredFeatNames, ...
    'RowNames', clusteredFeatNames);

disp('--- Grouped Feature Correlations (Clustergram Order) ---');
disp(GroupedValues);

% 6. Optional: Export to CSV to examine in Excel
% writetable(GroupedValues, 'Grouped_Muscle_Features.csv', 'WriteRowNames', true);

%% eval the 5 most negatively correlated pairs 
% 1. Settings
N_pairs = 5; % Number of top negative pairs to show
% h = evalin('base', 'optim_history');
% featNames = "Feat_" + (1:size(h,2));

% 2. Calculate Correlation and find the most negative ones
h = sensitivityMatrix;
C = corr(h);
% We only need the lower triangle to avoid duplicates (A-B vs B-A)
C_tri = triu(ones(size(C)), 1); 
C_flat = C;
C_flat(~C_tri) = NaN; % Hide diagonal and lower half

% Sort by correlation value
[sorted_vals, sorted_idx] = sort(C_flat(:), 'ascend');

% 3. Plotting
figure('Color', 'w', 'Name', 'Top Negative Tradeoffs');
tiledlayout('flow', 'Padding', 'compact');

for i = 1:N_pairs
    % Get linear index and convert to row/column
    linear_idx = sorted_idx(i);
    [row, col] = ind2sub(size(C), linear_idx);
    
    val = sorted_vals(i);
    if isnan(val), break; end
    
    % Create the scatter plot
    nexttile;
    scatter(h(:, col), h(:, row), 15, 'k', 'filled', 'MarkerFaceAlpha', 0.15);
    hold on;
    
    % Target point (0,0)
    plot(0, 0, 'rp', 'MarkerSize', 15, 'LineWidth', 2);
    
    % Formatting
    xlabel(featNames(col), 'FontWeight', 'bold');
    ylabel(featNames(row), 'FontWeight', 'bold');
    title(sprintf('Corr: %.3f', val));
    grid on;
    axis tight;
end

sgtitle(sprintf('Top %d Structural Tradeoffs (Pareto Fronts)', N_pairs));


%% 
%% Test: Spearmans and Pearsons corr

% 1. Setup
S = sensitivityMatrix; 


featNames = "Feat_" + (1:size(S,2));
featNames = arrayfun(@(i) sprintf('%d: %s', i, fn_latex{i}), 1:numel(fn), ...
'UniformOutput', false)

N_top = 9; % Number of strongest tradeoffs to show

% 2. Calculate Correlations
R_pearson = corr(S, 'Type', 'Pearson');
R_spearman = corr(S, 'Type', 'Spearman');

% 3. Composite Score (Average of both)
% We only care about pairs that are strongly negative in BOTH metrics
R_comp = (R_pearson + R_spearman) / 2;

% Isolate Upper Triangle and filter out positive correlations (the "friends")
R_tri = triu(R_comp, 1); 
R_tri(R_tri >= 0) = NaN; 

% 4. Find the most significant negative indices
[sortedVals, idx] = sort(R_tri(:), 'ascend');
validIdx = idx(~isnan(sortedVals));
N_to_plot = min(N_top, length(validIdx));

% 5. Visualize
figure('Color', 'w', 'Name', 'Robust Structural Tradeoff Analysis');
t = tiledlayout('flow', 'Padding', 'compact');

for k = 1:N_to_plot
    linear_idx = validIdx(k);
    [row, col] = ind2sub(size(R_tri), linear_idx);
    
    nexttile;
    x = S(:, col);
    y = S(:, row);
    
    % Define boundaries for the "Win-Win" Green Zone
    xl = [min([x; -0.1]) max([x; 0.1])];
    yl = [min([y; -0.1]) max([y; 0.1])];
    
    % Draw the Green Zone (Quadrant III: Both costs decreasing)
    patch([xl(1) 0 0 xl(1)], [yl(1) yl(1) 0 0], [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hold on;
    
    % Plot parameter sensitivities
    scatter(x, y, 60, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    
    % Axes for reference
    line([0 0], yl, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    line(xl, [0 0], 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    
    % Annotate the quadrants
    text(xl(1)*0.9, yl(2)*0.9, 'Tradeoff', 'Color', 'r', 'FontSize', 8);
    text(xl(2)*0.5, yl(1)*0.9, 'Tradeoff', 'Color', 'r', 'FontSize', 8);
    text(xl(1)*0.9, yl(1)*0.9, 'WIN-WIN', 'Color', [0 0.5 0], 'FontWeight', 'bold');

    xlabel(featNames(col) + " (\Delta Cost)");
    ylabel(featNames(row) + " (\Delta Cost)");
    title({['Pearson: ' num2str(R_pearson(row,col),2)], ...
           ['Spearman: ' num2str(R_spearman(row,col),2)]});
    grid on;
end

title(t, 'Structural Tradeoff Detection (Sensitivity Gradients)');
subtitle(t, 'Dots in Green Zone = Progress possible | Empty Green Zone = Structural Deadlock');

%% Attempt 2: with ranking by abs cost per feature

% 1. Setup Data
S = sensitivityMatrix; % Sensitivity (delta cost)
% ASSUMPTION: You have a vector of the current absolute residuals
% If not, replace this with: baseCosts = mean(abs(history(end,:)));
baseCosts = featureMatrixPlus(1, :); 

featNames = string(arrayfun(@(i) sprintf('%d: %s', i, fn_latex{i}), 1:numel(fn), ...
'UniformOutput', false));

N_top = 9; % Number of strongest tradeoffs to show

% 2. Calculate Robust Correlation
R_pearson = corr(S, 'Type', 'Pearson');
R_spearman = corr(S, 'Type', 'Spearman');
R_comp = (R_pearson + R_spearman) / 2;
R_comp  = R_spearman;

% 3. Calculate Cost Impact
[C1, C2] = meshgrid(baseCosts, baseCosts);
CostWeight = C1 .* C2; 

% 4. Create Priority Matrix
% We only care about the upper triangle and negative correlations
PriorityMatrix = R_comp;
PriorityMatrix(triu(ones(size(R_comp)), 0) == 0) = NaN; % Keep only upper triangle
PriorityMatrix(PriorityMatrix >= 0) = NaN;             % Keep only negative (conflicts)

% Final Priority: Magnitude of conflict * Importance of features
PriorityMatrix = abs(PriorityMatrix) .* CostWeight;

% --- THE FIX: Robust Indexing ---
% Find linear indices of all non-NaN entries
validLinearIdx = find(~isnan(PriorityMatrix));

if isempty(validLinearIdx)
    error('No negative correlations found. All features are "friends" or independent.');
end

% Get the priority values for these valid entries and sort THEM
validPriorities = PriorityMatrix(validLinearIdx);
[sortedPriorityValues, sortOrder] = sort(validPriorities, 'descend');

% Map back to the original linear indices
sortedLinearIdx = validLinearIdx(sortOrder);

% 5. Visualize
N_to_plot = min(N_top, length(sortedLinearIdx));
figure('Color', 'w', 'Name', 'High-Priority Structural Tradeoffs');
t = tiledlayout('flow', 'Padding', 'compact');

for k = 1:N_to_plot
    linIdx = sortedLinearIdx(k);
    [row, col] = ind2sub(size(PriorityMatrix), linIdx);
    
    nexttile;
    x = S(:, col);
    y = S(:, row);
    
    % Draw the Green Zone
    xl = [min([x; -0.1]) max([x; 0.1])];
    yl = [min([y; -0.1]) max([y; 0.1])];
    patch([xl(1) 0 0 xl(1)], [yl(1) yl(1) 0 0], [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hold on;
    
    % Scatter points
    scatter(x, y, 60, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    
    % Origin lines
    line([0 0], yl, 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    line(xl, [0 0], 'Color', [0.4 0.4 0.4], 'LineWidth', 1);

    xlabel(sprintf('%s\n(Cost: %.2e)', featNames(col), baseCosts(col)));
    ylabel(sprintf('%s\n(Cost: %.2e)', featNames(row), baseCosts(row)));
    
    % Now sortedPriorityValues(k) is guaranteed to be a number
    title(sprintf('Priority: %.2e\nCorr: %.2f', sortedPriorityValues(k), R_comp(row,col)));
    grid on;
end

title(t, 'High-Priority Structural Deadlocks (Cleaned Indexing)');