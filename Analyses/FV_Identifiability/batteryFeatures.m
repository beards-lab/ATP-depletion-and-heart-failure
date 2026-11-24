function [features_model, features_data] = batteryFeatures(params0)
% batteryFeatures  Run the slack battery and return the RAW feature structs.
%
% evaluateBakersExp returns only the weighted cost per feature; the Jacobian
% needs the feature values themselves. Same two calls as evaluateBakersExp:
% LoadData and RunBakersExp are scripts that read params0 (with params0.mods
% and params0.g applied inside) and leave features_data / features_model here.
params0.PlotEachSeparately = false; params0.PlotFeatureFitting = false;
features_model = []; features_data = []; %#ok<NASGU>
LoadData;
RunBakersExp;
end
