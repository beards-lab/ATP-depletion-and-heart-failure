function [Et E] = evaluateBakersExp(g, params0)
% Evaluate Bakers' problem
if any(g<0) 
    Et = NaN;
    return;
end

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
    t_ss = [0 1];
    t_sl0 = [0 0.1];
    lastwarn('', ''); 
    RunBakersExp;
    Et = sum(E);
catch e
    Et = NaN;
end