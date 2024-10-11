% Attempt on Bahador's parametrization from the paper 
% change in model states etc..
% Beard et al. 2022 "Reduced cardiac muscle power with low ATP simulating heart failure" https://doi.org/10.1016/j.bpj.2022.07.029

clear;
figure(80085);clf;
params0 = getParams();
ModelParamsInitNiceSlack;

params0.ka = 373;
params0.kd = 103;
params0.k1 = 1/(1/40 + 1/419); % merged k1 and k2
params0.k2 = 44;
params0.ksr0 = 250;
params0.kmsr = 250;
params0.dr = 0.01;
params0.kstiff1 = 1400;
params0.kstiff2 = 1400;
params0.alpha0 
params0.alpha1 = 15;
params0.alpha2_L = 3*2.77; % exp((0.03)*3*2.77) ~= exp((0.03^2)*277)
params0.dr2_L = 0;
params0.dr3 = 0;
params0.dr1 = 0;
params0.sigma1 = 33;
params0.sigma2 = 1e9;
params0.k2_R = 1;
params0.k2_L = 44;
params0.dr2_R = 0;
params0.alpha2_R = 1;
params0.justPlotStateTransitionsFlag = false;
params0.RunForceVelocity = false;
params0.RunSlack = true;
params0.ShowStatePlots = true;

LoadData; 
tic
RunBakersExp;
toc
xl = xlim();
figure(3); clf;
rates = [out.RTD; out.RD1; out.R1D; out.R12; out.R21; out.XB_Ripped; out.RSR2PT; out.RPT2SR; out.XB_TORs];
lgs = {'RTD', 'RD1', 'R1D', 'R12', 'R21', 'XB_{Ripped}', 'SR2PT', 'PT2SR', 'A2 dett'}
ints = [1 2 7 8 9]

plot(out.t, rates(ints, :), 'LineWidth',2);
xlim(xl)
legend(lgs(ints))