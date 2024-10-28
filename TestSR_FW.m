%% test SR within the framework
figure(2);clf;

% from TestSR
g = [2.5670    1.1465    2.2885    1.2476];

params0 = getParams();
ModelParamsInitNiceSlack


params0.WindowsOverflowStepCount = 5;
params0.UseSuperRelaxed = 1;
params0.UseSpaceExtension = true;

params0.justPlotStateTransitionsFlag = false;

% params0.Slim = 0.3;
% params0.LXBpivot = 1.7;
% params0.dS = 0.004;
% params0.Slim_l = 1.7;
% params0.Slim_r = 2.21;


params0.kmsr = 5*g(1);
params0.sigma2 = 15*g(2);
params0.ksr = 1000*g(3);
params0.sigma1 = 1e6;


params0.alpha0 = 0;
params0.alpha1 = 0;
params0.alpha2_L = 0;
params0.alpha2_R = 0;
params0.alpha3 = 0;

params0.ka = 100;
params0.kd = 0;

params0.k1 = 25;
params0.k_1 = 0;

params0.k2 = 10;
params0.k_2 = 0;

params0.kah = 100;
params0.kadh = 0;

params0.dr = 0.01;

params0.kSE = 5e3; 

params0.kstiff1 = 1e4*g(4);
params0.kstiff2 = 1e4*g(4);

params0.UseOverlapFactor = true;
params0.ShowStatePlots = true;
params0.UseTitinInterpolation = false;
params0.UsePassive = true;
 
params0.justPlotStateTransitionsFlag = false;
params0.dS = 0.005;
clf;
simulateThat([], params0, true);

% writeParamsToMFile('ModelParams_SR.m', params0);

%%

params0.vmax1 = 20;
params0.FudgeA = 0;
params0.FudgeB = 121.922;
params0.FudgeC = -209.18;

SL = 1.85:0.05:2.05;
velHS = -(params0.FudgeA*SL.^2 +params0.FudgeB*SL + params0.FudgeC);
clf;plot(SL, -velHS, 'x-', LineWidth=2);
    title('Fudged velocity');xlabel('SL');ylabel('Velocity (um/s)');
    



function Es = simulateThat(g, params, plotThat)

if nargin < 2
    plotThat = false;
end

SL = [2.2 2.0400    2.0000    1.9600    1.9200  1.8800];
df = [76.5 68.3878   65.5438   59.2362   52.5507   43.9655];

% SL = 1.8:0.05:2.1; 
% df = ones(size(SL));

for i = 1:length(SL)
    params.SL0 = SL(i);
    params.Velocity = 0;
    if isfield(params, 'PU0')
        params = rmfield(params, 'PU0');
    end

    % reset the PU0
    params = getParams(params, params.g, true);
    
    % [F out] = evaluateModel(modelFcn, velocitytable(:, 1), params);
    [F out] = evaluateModel(@dPUdTCaSimpleAlternative2State, [0 100], params);
F_total(i) = F;
F_passive(i) = out.FXBPassive(end);
end
 

if plotThat
    %%
    % figure(2); clf;
    hold on;
    plot(SL, df, 'o-', LineWidth=2)
    plot(SL,F_total,'x-', SL,F_passive,'x--', LineWidth=2);
    title('Steady-state tension');xlabel('SL');ylabel('Tension (kPa)');
    legend('Data (estimated)', 'Simulation', 'Passive tension component (model)')

end

E = (F_total - df).^2;
Es = sum(E);

end
