function [Et E] = evaluateBakersExp(g, params)
% Evaluate Bakers' problem

% important to start with the g!!
params = getParams([], g);

params.SL0 = 2.2;
% params.Slim = 0.18;
params.Slim = 0.3;
params.N = 30;
params.MgATP = 8;

params.UseTORNegShift = false;
params.UseMutualPairingAttachment = false;
params.UseOverlap = true;
params.UsePassive = true;
params.UseSerialStiffness = true;
% params = getParams
params.PlotEachSeparately = true;


try
%     Et = 1;
    LoadData;
    t_ss = [0 1];
    t_sl0 = [0 0.1];
    RunBakersExp;
    Et = sum(E);
catch e
    Et = NaN;
end