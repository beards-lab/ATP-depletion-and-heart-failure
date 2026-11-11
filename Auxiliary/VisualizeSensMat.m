%% INPUTS:
% featureMatrixPlus  (nParam x nFeat)
% featureMatrixMinus (nParam x nFeat)
% row 1 is baseline

%% Compute total costs
baseline = sum(featureMatrixPlus(1,:));

costPlus  = sum(featureMatrixPlus,  2) - baseline;
costMinus = sum(featureMatrixMinus, 2) - baseline;
if ~RunDeltaMinus
    costMinus = costPlus;
end
if ~RunDeltaPlus
    costPlus = costMinus;
end

nParam = length(costPlus);

%% Determine improvement & direction
bestChange = min(costPlus, costMinus);

direction = strings(nParam,1);
direction(costPlus<0 & costMinus>=0) = "+delta";
direction(costMinus<0 & costPlus>=0) = "-delta";
direction(costPlus<0 & costMinus<0)  = "both";
direction(costPlus>=0 & costMinus>=0) = "none";

improvingParams = find(bestChange < 0);

%% Interaction potential (nonlinearity indicator)
interactionPotential = abs(costPlus - costMinus);

%% ---- VISUALIZATION ----------------------------------------------------
figure('Name','Parameter Sensitivity Summary');
tiledlayout(2,1);

% 1) Grouped bar chart: full +δ and −δ effect
nexttile;
bar(1:nParam, [costMinus, costPlus], 'grouped');
yline(0,'k--','LineWidth',0.8);
xlabel('Parameter index');
ylabel('\Delta total cost');
legend('-\delta', '+\delta');
title('Cost change for each parameter (+/- delta)');
xticks(1:nParam);
xticklabels(paramNames);
xtickangle(45);   % or 30 / 90 depending on length

% 2) Scatter: best improvement & interaction potential
nexttile;
scatter(1:nParam, bestChange, 60, interactionPotential, 'filled');
yline(0,'k--','LineWidth',0.8);
xlabel('Parameter index');
ylabel('Best improvement');
title('Best improvement colored by interaction potential');
colormap turbo;
colorbar; ylabel('Interaction potential');

%% ---- SUMMARY PRINT ----------------------------------------------------
fprintf('\nParameter Sensitivity Summary:\n');
for p = 1:nParam
    fprintf('Param %2d: best %+7.4f   direction %-6s   interaction %.4f\n', ...
        p, bestChange(p), direction(p), interactionPotential(p));
end

fprintf('\nImproving parameters: ');
disp(improvingParams');

%% Worth reoptimizing

%% INPUTS:
% featureMatrixPlus  (nParam x nFeat)
% featureMatrixMinus (nParam x nFeat)
% Row 1 is baseline.

%% Compute total cost changes
baseline = sum(featureMatrixPlus(1,:));
costPlus  = sum(featureMatrixPlus,  2) - baseline;
costMinus = sum(featureMatrixMinus, 2) - baseline;

nParam = length(costPlus);

%% 1) Direct improvement
bestChange = min(costPlus, costMinus);      % negative = improvement
improvingParams = find(bestChange < 0);

%% 2) Nonlinear / opposite-sign behavior
% (indicates potentially valuable when other params vary)
nonlinearParams = find(sign(costPlus) ~= sign(costMinus));

%% 3) Interaction potential: large difference between + and - effects
interactionPotential = abs(costPlus - costMinus);

% Threshold: relative to median effect range
ipThresh = 0.25 * median(interactionPotential);   % adjustable
highInteractionParams = find(interactionPotential > ipThresh);

%% Combine all:
shouldReoptimize = unique([improvingParams; nonlinearParams; highInteractionParams]);

%% ---- SUMMARY PRINT ----------------------------------------------------
fprintf('\n=== Parameters Worth Re-Optimizing ===\n');
disp(shouldReoptimize');

fprintf('\nImproving parameters: ');
disp(improvingParams');

fprintf('Nonlinear/opposite-sign parameters: ');
disp(nonlinearParams');

fprintf('High-interaction parameters: ');
disp(highInteractionParams');

%% sorted names
%% INPUT:
% paramNames = string or cell array like {'p1','p2','p3',...}
% featureMatrixPlus, featureMatrixMinus
% row 1 is baseline

baseline = sum(featureMatrixPlus(1,:));

costPlus  = sum(featureMatrixPlus, 2) - baseline;
costMinus = sum(featureMatrixMinus, 2) - baseline;

% Best improvement for each parameter
bestChange = min(costPlus, costMinus);   % negative is better

%% Sort by improvement magnitude (ascending: most negative first)
[sortedValues, idx] = sort(bestChange, 'ascend');

sortedNames = paramNames(idx);

%% Display as a table
T = table(sortedNames', sortedValues, ...
    'VariableNames', {'Parameter','BestImprovement'});

disp(T);
