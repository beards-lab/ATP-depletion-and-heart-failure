% driver code
g = ones(30, 1);
g = x;


params = getParams();
params.g = x;
% optimized
params.mods = {"kstiff1", "kstiff2", "kstiff3", "k1", "k2", "k_2", "k3", "s3", "alpha3", 'ksr0', 'sigma0', 'kmsr'};

% optimized for force-velocity only
params.g = [1.3581    1.0578    0.0299    0.4113    1.1269    0.5243    3.3260    1.0322    0.0274    1.7252    1.6593    0.0212];

% g = x2;
% g = p_OptimGA;
% g = ones(30, 1);
% params.g = x;

params.SL0 = 2.2;
% params.Slim = 0.18;
params.Slim = 0.3;
params.dS = 10e-3;
% params.N = 30;
params.MgATP = 8;
params.SimTitle = '';
figure(12);clf;

ghostSave = '';
% ghostSave = 'ones';
% ghostSave = 'ga2';
% ghostSave = 'prev';
% ghostSave = 'mybase';
% ghostSave = 'beardsOrig_passive';
% ghostSave = 'ShiftingStrain40_Slim0_3';
% ghostSave = 'ShiftingStrain80_Slim0_3';
% ghostSave = 'operatorSplittingPU020';
% ghostSave = 'operatorSplittingPU0';
% ghostSave = 'operatorSplittingPU080';
% ghostSave = 'ShiftingStrain160_Slim0_3';
% ghostSave = 'beardsOrig_all60';
% ghostSave = 'beardsOrig_OV_Pas_SS';
% ghostSave = 'beardsOrig';
% ghostSave = 'bO_Fp31';
% ghostSave = 'beardsOrig50'; % strain shifting dS 50 nm
% ghostSave = 'beardsOrig5'; % strain shifting dS 5 nm
% ghostSave = 'beardsOrig5_BW';%bells and whistles :)


ghostLoad = 'ones';
% ghostLoad = 'InterpolateBins';
% ghostLoad = 'beardsOrig';
% ghostLoad = 'ShiftingStrain40';
% ghostLoad = 'operatorSplitting'
% ghostLoad = 'operatorSplittingPU020';
% ghostLoad = 'operatorSplittingPU0';% N = 40
% ghostLoad = 'operatorSplittingPU080';
% ghostLoad = 'ShiftingStrain40_Slim0_3';
% ghostLoad = 'ShiftingStrain80_Slim0_3';
% ghostLoad = 'ShiftingStrain160_Slim0_3';
% ghostLoad = 'ShiftingStrainTest20';
% ghostLoad = 'ShiftingStrainTest40';
% ghostLoad = 'ShiftingStrainTest80';
% ghostLoad = 'ShiftingStrainTest160';
% ghostLoad = 'ShiftingStrainTest240';
% ghostLoad = 'beardsOrig_all30';
% ghostLoad = 'beardsOrig_all60';
% ghostLoad = 'operatorSplittingPU0';% N = 40
% ghostLoad = 'beardsOrig5';% N = 40
% ghostLoad = 'beardsOrig5_BW';% N = 40
% ghostLoad = 'bO_Fp31';

% testing setup
% params.UseOverlap = true;
% params.UsePassive = true;

params.UseTORNegShift = false;
params.UseMutualPairingAttachment = false;
params.UseSlack = false;
params.UseOverlap = true;
params.UsePassive = true;
params.UseSerialStiffness = false;
params.UseKstiff3 = true;

% params.alpha1 = 0;
% params.alpha2 = 0;

% need a ksttiff1 and kstiff2 parameter retune
params.F_act_UseP31 = true;
params.UseP31Shift = true;
params.kstiff1 = 1000;
params.kstiff2 = 14000;
params.kstiff3 = 1400;
% params.kstiff2 = 7000;
% params.kstiff1 = 700000;
params.k1 = 4000;
params.k2 = 10000;
params.k3 = 100;
params.dr = 0.01;
params.mu = 1e-6;

% set zero transition slopes
params.alpha1 = 0;
params.alpha2 = 0;
params.alpha3 = 300;
% zero reverse-flows
% params.k_1 = 10000;
% params.k_2 = 0;

params.PlotEachSeparately = true;
params.PlotFullscreen = false;
params.ghostLoad = ghostLoad;
params.ghostSave = ghostSave;
% params.UseSpeedHalving = true;

plotTransitions = true;
params.PlotFullscreen = true;
% init
clf;

% save as default
params0 = getParams(params, params.g, false, true);

% params0.mods = {"kstiff1", "kstiff2", "k1", "k2", "k3"};
% params0.mods = {"kstiff1", "kstiff2", "k1", "k2", "k_2", "k3", "s3", "alpha3"};
% params0.g = load('p_OptimGa.mat').p_OptimGA;
% params0.g(1) = params0.g(1)*2;
% params0.g(2) = 40;
% params0.dr = 0.01;
params0.UseSpaceDiscretization = false;
params0.UseSpaceInterpolation = true;
% g = x;

LoadData;
% vel = vel(1:6);
t_ss = [0 1];
t_sl0 = [0 0.1];
tic

% params0 = load('gaOutput2_params.mat').params;
% params0 = load('bestParams.mat').params;
params0.PlotEachSeparately = true;
% params0.UseKstiff3 = false;
params0.dS = 0.01;
params0.ksr0 = params0.ksr0;
params0.sigma0 = params0.sigma0;

RunBakersExp;
toc
E
xlim('auto')
ylim('auto')

return
%% SA
E = [];
params0.mods = {"kstiff1", "kstiff2", "kstiff3", "k1", "k2", "k_2", "k3", "s3", "alpha3", 'ksr0', 'sigma0', 'kmsr'};
paramsfn = params0.mods;
% paramsfn = fieldnames(params)
params_d = params0;
RunBakersExp; e0 = E;
for i = 1:length(paramsfn)
    disp(['Computing ' char(paramsfn{i}) '..'])
    params0 = params_d;
    params0.(paramsfn{i}) = params0.(paramsfn{i})*0.95;
    RunBakersExp;
    e_d = E; % E -- delta
    params0 = params_d;
    params0.(paramsfn{i}) = params0.(paramsfn{i})*1.05;
    RunBakersExp;
    epd = E;% e + delta 

    err_1(i) = min(e_d, epd);
end
params0 = params_d;
figure;bar(err_1);hold on;plot([1 length(paramsfn)], [e0 e0])    

%% OPTIM
% return
params0.mods = {"kstiff1", "kstiff2", "kstiff3", "k1", "k2", "k_2", "k3", "s3", "alpha3", 'ksr0', 'sigma0', 'kmsr'};
g = ones(size(params0.mods))
% g = x;

params0.PlotEachSeparately = false;
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.1, 'PlotFcns', @optimplotfval, 'MaxIter', 500);
% g = [1, 1, 1, 1, 1, 1, 1, 1];
% g = [1.2539    0.4422];
optimfun = @(g)evaluateBakersExp(g, params0);
x = fminsearch(optimfun, g, options)
% g = x;
%% GA

%% Attempt on GA
% parpool
ga_Opts = optimoptions('ga', ...
    'PopulationSize',64, ...            % 250
    'Display','iter', ...
    'MaxStallGenerations',4, ...  % 10
    'UseParallel',true);


params0.mods = {"kstiff1", "kstiff2", "k1", "k2", "k_2", "k3", "s3", "alpha3", "sigma0", "ksr0"};
Ng = length(params0.mods);
params0.UseKstiff3 = false;
params0.SaveBest = true;
params0.ghostSave = false;


params0.PlotEachSeparately = false;

% evaluateBakersExp(g, params0)

optimfun = @(g)evaluateBakersExp(g, params0);

% default bounds struct - based on the mods
upp = 100; lwr = 0.01;
p_lb = params0; for i = 1:length(params0.mods), p_lb.(params0.mods{i}) = params0.(params0.mods{i})*lwr; end
p_ub = params0; for i = 1:length(params0.mods), p_ub.(params0.mods{i}) = params0.(params0.mods{i})*upp; end

% specific bounds - this is universal and non-blocking. 
% Might not be used, if not within 'mods' - this is the whole reason to make all this mess
p_lb.dr = 0.005; p_ub.dr = 0.015;
p_lb.s3 = 0.005; p_ub.s3 = 0.015;

% convert back to vector of ratios
for i = 1:length(params0.mods)
    lb(i) = p_lb.(params0.mods{i})/params0.(params0.mods{i});
    ub(i) = p_ub.(params0.mods{i})/params0.(params0.mods{i});
end

[p_OptimGA,Res_OptimGA,~,~,FinPopGA,FinScoreGA] = ...
    ga(optimfun,Ng, ...
    [],[],[],[],lb,ub,[],ga_Opts);

save env;

optimfun(p_OptimGA)

x2 = fminsearch(optimfun, p_OptimGA, options)
save x2
optimfun(x2)

%% Calculation of max speed
dr = params.dr; % um
v = 10; % ML/s = um/s in half-sarcomere
% to produce F we need to have sum < dr, so we need to unbind the rest
% thus, zero force at sum = dr and all binded are unbinded
% tor = ; 

