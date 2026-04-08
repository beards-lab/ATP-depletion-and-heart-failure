%% test SR within the framework
figure(3);clf;

% from TestSR
% g = [2.5670    1.1465    2.2885    1.2476];

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
%%
clf;
simulateForceLengthEstim([], params0, true);

% writeParamsToMFile('ModelParams_SR.m', params0);

%% try slack
figure(1);clf;
ModelParams_SR;
params0.justPlotStateTransitionsFlag = false;
params0.FudgeVmax = true;
params0.FudgeA = 0;
params0.FudgeB = 121.922;
params0.FudgeC = -209.18;

sd = @(kx, alphaL, alphaR, dr,eL, eR) min(1e4, kx*(exp((alphaL*(s-dr)).^eL).*(s<dr) + exp((alphaR*(s-dr)).^eR).*(s>=dr)));


% alpha2_L, params.alpha2_R, params.dr2, params.e2L

params0.e2l = 1.8;
params0.k2_L = 190.15;
params0.alpha2_L = -192.115;

params0.e2R = 2;
params0.k2_R = 8000;
params0.alpha2_R = 25;


params0.alpha0 = 0*36.1309;
params0.kd = 4;

params0.ka = 20;

params0.xrate = 6;

params0.kstiff1 = 1.4e4*1.3;
params0.kstiff2 = 1.4e4*1.3;
% params0.dr = 

% params 
% params0.sigma2 = 
% params: ksr, kmsr, sigma2, kstiff1 = kstiff2, 4 fwd rate consts, 1 rwd rate const, alpha2, dr, 1to2?
params0.mods = {'ksr', 'kmsr', 'sigma2', 'kstiff2', 'ka', 'kd', 'k1', 'k2', 'kah', 'alpha2_R', 'alpha2_L', 'e2R', 'e2L'};

% params0.k2 = params0.k2/exp((params0.alpha2_L*(0-params0.dr2)).^params0.e2L)

% have the detachment in the ballpark, otherwise unstable
params0.RunSlackSegments = 'Fourth-rampuponly';
% simulateThat([], params0, true);
params0.RunForceLengthEstim = true;

tic
RunBakersExp;
toc


%%
params0.vmax1 = 20;
params0.FudgeA = 0;
params0.FudgeB = 121.922;
params0.FudgeC = -209.18;

SL = 1.85:0.05:2.05;
velHS = -(params0.FudgeA*SL.^2 +params0.FudgeB*SL + params0.FudgeC);
clf;plot(SL, -velHS, 'x-', LineWidth=2);
    title('Fudged velocity');xlabel('SL');ylabel('Velocity (um/s)');
    
