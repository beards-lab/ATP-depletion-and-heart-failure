clear;
figure(2);
clf; 

% initialize parameters
params0 = getParams();
ModelParamsOptim_fudgeSlackVelocities

% new modifiers - to be overwritten
params0.UseOverlapFactor = true;
% params0.ksr = params0.ksr0;

g = [2.5832    1.5169    1.7629    1.4585]; % better with actual estimated data

params0.kmsr = 5*g(1);
params0.sigma2 = 15*g(2);
params0.ksr = 100*g(3);
params0.sigma1 = 1e6;


% params0.alpha0 = 0;
% params0.alpha1 = 0;
% params0.alpha2_L = 0;
% params0.alpha2_R = 0;
% params0.alpha3 = 0;

params0.ka = 100;
% params0.kd = 0;
% 
params0.k1 = 25;
% params0.k_1 = 0;
params0.alpha1 = 19;
% 
params0.k2 = 10;
% params0.k_2 = 0;
% 
params0.kah = 100;
% params0.kadh = 0;
% 
params0.dr = 0.01;
% 
% params0.kSE = 5e3; 
% 
% params0.kstiff1 = 1e4*g(4);
% params0.kstiff2 = 1e4*g(4);

params0.UseOverlapFactor = true;
params0.ShowStatePlots = true;
params0.drawForceOnset = true;
params0.UseTitinInterpolation = false;
params0.UsePassive = true;
 
params0.justPlotStateTransitionsFlag = false;
params0.dS = 0.005;

LoadData; 
tic
RunBakersExp;
toc

%%

% all possible
params0.mods = {'alpha0', 'alpha1', 'alpha2', 'alpha2_L', 'alpha2_R', 'dr0', 'dr1', 'dr2', 'dr_1', 'e2L', 'e2R', 'k1', 'k2', 'kSE', 'k_1', 'ka', 'kah', 'kd', 'kmsr', 'ksr', 'kstiff1', 'kstiff2', 'alpha_1', 'sigma1', 'sigma2'};
params0.g = ones(size(params0.mods));
saSet = 1:length(params0.mods);

%% RUN driverSA on that, selecting 11 most influential, ie. 
% {'sigma2'}    {'alpha2_R'}    {'e2R'}    {'kSE'}    {'alpha2_L'}    {'alpha0'}    {'e2L'}    {'alpha1'}    {'kstiff2'}    {'ksr'}  {'dr2'}
% running the best optim candidate
g = [1.0265    1.0938    1.1562    1.2040    0.3318    1.0419    1.0189 0.8591    0.9833    1.1513    1.1437];
% the best money can buy?
g = [1.0345    1.1061    1.1532    1.2557    0.3307    1.0592    1.0439    0.8015    0.9718    1.1849    1.1407];
params0.g = g;
RunBakersExp;

