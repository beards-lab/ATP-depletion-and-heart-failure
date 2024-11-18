function [Et E] = evaluateBakersExp(g, params0)
% Evaluate Bakers' problem
if any(g<0) 
    Et = NaN;
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
    Et = sum(E);
catch e
    try
        if ~isempty(e.cause)
            % Extract the cause message if present
            causeMessage = e.cause{end}.message;
            
            % Parse the error value from the cause message
            Et = str2double(causeMessage)*1e6;
        else
            Et = 1e9;
        end
    catch ee
        Et = 1e12;
    end
end