function [Et E] = evaluateBakersExp(g, params0, limitNegative)
if nargin < 3
    limitNegative = true;
end
% Evaluate Bakers' problem
if limitNegative && any(g<0)% || ...
        %params0.kstiff1 < params0.kstiff1_n || params0.kstiff2 < params0.kstiff2_n
    Et = NaN;
    warning(' g < 0! Skipping ...');
    return;
end

params0.PlotEachSeparately = false;
% important to start with the g!!
% params0 = getParams(params0, g, true, true);
params0.g = g;

% params0.SL0 = 2.2;
% % params.Slim = 0.18;
% params0.Slim = 0.3;
% params0.N = 30;
% params0.MgATP = 8;
% 
% params0.UseTORNegShift = false;
% params0.UseMutualPairingAttachment = false;
% params0.UseOverlap = true;
% params0.UsePassive = true;
% params0.UseSerialStiffness = true;
% % params = getParams
% params0.PlotEachSeparately = true;
% 
try
%     Et = 1;
    LoadData;
    lastwarn('', ''); 
    RunBakersExp;
    % Et = sum(E);
    % use feature costs instead
    [cost, ~, cost_raw] = evalFeatureCost(features_data, features_model, params0.fn);
    Et = sum(cost);
    
    if Et == 0
        % sum fin wong
        Et = 1e3;
    end
catch e
    try
        if ~isempty(e.cause)
            % Extract the cause message if present
            causeMessage = e.cause{end}.message;
            warning(causeMessage);
            
            % Parse the error value from the cause message
            Et = 1e3 + str2double(causeMessage);
        else
            Et = 1e9;
        end
    catch ee
        Et = 1e12;
    end
end

Et = abs(Et);