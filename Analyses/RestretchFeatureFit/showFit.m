% showFit.m — compact feature-cost + restretch/FV readout for the current params0.
% Assumes params0 (with .fn) already in the base workspace. Runs FV+slack, prints
% the total feature cost, the top cost drivers, and a restretch+FV data-vs-model
% table. Kept terse on purpose (token-efficient iterative loop).
params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0; params0.BreakOnODEUnstable = false;
RunBakersExp;
[cv, ~, ~] = evalFeatureCost(features_data, features_model, params0.fn, 2);
fprintf('\n==== TOTAL feature cost = %.4f ====\n', sum(cv,'omitnan'));
[~,ord] = sort(cv,'descend','MissingPlacement','last');
fprintf('top drivers (w*cost):\n');
for k = ord(:)'
    if cv(k) > 0.03
        lbl = params0.fn{k}; if numel(lbl)>34; lbl=[lbl(1:31) '...']; end
        fprintf('   %-34s %6.3f\n', lbl, cv(k));
    end
end
restr = {'Am','peak1_y','peak1_dSL','vall_y','vall_t','peak2','restretchSlopeStart','steady','vall2_dy','ovrsht_dy','ktr'};
fprintf('\n%-20s %-30s %-30s\n','feature','data','model');
for i=1:numel(restr); f=restr{i};
    d=''; m='';
    if isfield(features_data,f);  d=num2str(round(double(features_data.(f)),3,'significant')); end
    if isfield(features_model,f); m=num2str(round(double(features_model.(f)),3,'significant')); end
    fprintf('%-20s %-30s %-30s\n', f, d, m);
end
if isfield(features_model,'FV_fnorm')
    fprintf('FV_fnorm  data =%s\n          model=%s\n', num2str(round(features_data.FV_fnorm,3)), num2str(round(features_model.FV_fnorm,3)));
end
