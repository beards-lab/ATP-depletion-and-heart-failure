% driver code
g = ones(30, 1);
% g = x;


params = getParams();
% params.g = x;
% optimized
% params.mods = {"kstiff1", "kstiff2", "kstiff3", "k1", "k2", "k_2", "k3", "s3", "alpha3", 'ksr0', 'sigma0', 'kmsr'};

% optimized for force-velocity only
% params.g = [1.3581    1.0578    0.0299    0.4113    1.1269    0.5243    3.3260    1.0322    0.0274    1.7252    1.6593    0.0212];

% g = x2;
% g = p_OptimGA;
% g = ones(30, 1);
% params.g = x;

params.SL0 = 2.2;
% params.Slim = 0.18;
params.Slim = 0.3;
params.dS = 1e-3;
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


ghostLoad = '';
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

params.UseTORNegShift = true;
params.UseMutualPairingAttachment = false; % just breaks the sim
params.UseSlack = true;
params.UsePassive = true;
params.UseSerialStiffness = true;
params.kSE = 10e3;
params.UseKstiff3 = true;

% params.alpha1 = 0;
% params.alpha2 = 0;

% need a ksttiff1 and kstiff2 parameter retune
params.F_act_UseP31 = true;
params.UseP31Shift = true;
params.kstiff1 = 1000*0.7;
params.kstiff2 = 10;
params.kstiff3 = 100000;
% params.kstiff2 = 7000;
% params.kstiff1 = 700000;
params.ka = 10;
params.kd = 500;%.1*5000.2;
params.k1 = 5000;
% params.k_1 = 170*1e6; % NA
params.k2 = 1000;
% params.k_2 = 140*1e6; % NA
params.k3 = 50;
params.dr = 0.01;
params.mu = 1e-4;
params.TK = 0;
params.TK0 = 0.01;

% set zero transition slopes
params.alpha1 = 0*0.500;
params.alpha2 = 0;
params.alpha3 = 148;
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

params.WindowsOverflowStepCount = 4;
params.dS = 0.004;
params.Slim_l = 1.6;
params.Slim_r = 2.25;

% init
clf;

% % optimized for force-velocity only
% params.g = [1.3581*1    1.05788*0.9    0.0299    0.4113    1.1269    0.5243    3.3260    1.0322    0.0274    1.7252    1.6593    0.0212];

% save as default, applying the modifiers
params0 = getParams(params, params.g, false, true);
params0.ss
%
% % add some more modifiers, optimized for ATP
% params0.mods = {"K_T1", "K_T3"}; g = [1.4054    0.7373];
% params0 = getParams(params0, g, false, true);


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
params0.PlotFullscreen = false;
% params0.UseKstiff3 = false;
% params0.dS = 0.01;
params0.ksr0 = params0.ksr0;
params0.sigma0 = params0.sigma0;
params0.EvalAtp = [1];
% params0.UseAtpOnUNR = true;
% params0.kstiff3 = params0.kstiff2;
params0.UseTitinModel = false;
params0.UsePassive = false;

params0.RunForceVelocity = false;
params0.RunKtr = false;
params0.RunSlack = true;
params0.RunStairs = false;
params0.UseOverlap = true;

params0.UseSerialStiffness  = true;
params0.UseSlack = true;
params0.UseKstiff3 = false;
params0.F_act_UseP31 = false;


RunBakersExp;
% plot(out.t, out.XB_TORs, 'g--')
% xlim([-0.05, 0.15])
toc
E
% xlim('auto')
% ylim('auto')
6.3/6.04
%%
plot(out.t, out.p1_0, out.t, out.p2_0, out.t, out.p3_0, out.t, out.p2_1, out.t, out.p3_1, out.t, 1-out.NR)
legend('1_0' , '2_0', '3_0', '2_1', '3_1', 'SR')

return
%%
StatesInTime;
%% SA
E = [];
params0.mods = {"kstiff1", "kstiff2", "kstiff3", "k1", "k2", "k_2", "k3", "s3", "alpha3", 'ksr0', 'sigma0', 'kmsr'};
paramsfn = params0.mods;
% paramsfn = fieldnames(params)
params_d = params0;
RunBakersExp; e0 = sum(E);
for i = 6:length(paramsfn)
    disp(['Computing ' char(paramsfn{i}) '..'])
    params0 = params_d;
    params0.(paramsfn{i}) = params0.(paramsfn{i})*0.95;
    RunBakersExp;
    e_d(i) = sum(E); % E -- delta
    params0 = params_d;
    params0.(paramsfn{i}) = params0.(paramsfn{i})*1.05;
    RunBakersExp;
    epd(i) = sum(E);% e + delta 

    % err_1(i) = min(e_d, epd);
end
params0 = params_d;
%%
err_1 = min(e_d, epd)
figure;bar([e_d; epd]');hold on;plot([1 length(paramsfn)], [e0 e0])    
xticklabels(params0.mods)

%% OPTIM
% return
params0.mods = {"kstiff2", "kstiff3", "ka", "kd", "k1", "k2", "k3", "alpha3", 'sigma0', 'kmsr'}; % tune everything
% params0.g = ones(1, length(params0.mods))
params0.g = g;
% params0.mods = {"K_T1", "K_T3"};
% g = ones(size(params0.mods))
% optimized
% g = [1.4054    0.7373];

% modtbl = readtable("modifierstbl.csv");
% params0.mods = modtbl.Properties.VariableNames;
% params.g = modtbl(1, :).Variables;
%%

params0.PlotEachSeparately = false;
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.1, 'PlotFcns', @optimplotfval, 'MaxIter', 1500);
% g = [1, 1, 1, 1, 1, 1, 1, 1];
% g = [1.2539    0.4422];
optimfun = @(g)evaluateBakersExp(g, params0);
x = fminsearch(optimfun, params0.g, options)
% g = exp(x);
g = x;
%% Fminsearch optim results - not able to hit the fast down and slow up
clf
params0.PlotEachSeparately = true;
params0.mods = {'kah', 'kadh', 'ka', 'kd', 'k1', 'k_1', 'k2', 'kstiff2', 'alpha2', 'alpha3', 'kSE'};
params0.g = g;
RunBakersExp;
plot(datatable(i_0:i_e, 1), e)
%%
params0.RunForceVelocity = false;
params0.RunSlack = true;
%% saving the params
% writeParamsToFile('params.csv', params)
figure(222);clf;
% params0 = UpdateValuesFromFile('params.csv', params);
params0.g = [];
params0.mods = [];
params0.kah = 100;
params0.kadh = 10;
params0.ka = 500;
params0.kd = 100;
params0.k1 = 1000;
params0.k_1 = 10;
params0.k2 = 200;
params0.k_2 = 1;
params0.k3 = 10000;
params0.ksr0 = 9.0853;
params0.sigma0 = 33.125;
params0.kmsr = 250;
params0.kstiff1 = 1e1;
params0.kstiff2 = 1e4;
% params0.kstiff3 = 1e2;
params0.alpha0 = 1;
params0.dr0 = 0.001;
params0.alpha1 = 1;
params0.dr1 = 0.001;
params0.alpha_1 = 1;
params0.dr_1 = 0.001;
params0.alpha2 = 1;
params0.dr2 = 0.01;
params0.alpha3 = 1;
% params0.dr3 = 0.01;
params0.kSE = 10000;
params0.k_pas = 200;
params0.dr = 0.01;
params0.UsePassive = true;
params0.UseOverlap = false;


params0.PlotEachSeparately = true;
RunBakersExp;
hold on;plot(out.t, out.XB_TORs)
nexttile;
plot(out.t, out.p1_0, '-', out.t, out.p2_0, '-', out.t, out.p3_0, '-',out.t, out.PuATP, '-',out.t, out.PuR, '-', LineWidth=1.5, LineStyle='-')
legend('P1','P2','P3','PuATP','PuR')

xlim([2.6 3.05])
title(sprintf('XB TOR %g', out.XB_TORs(end)))
% allProbs = sum([out.p1_0; out.p2_0; out.p3_0; out.PuATP; out.PuR]);
% plot(out.t, allProbs)
% StatesInTime;
%%
s = -0.05:params.dS:0.05;
plot(s, -1 + min(params.alpha2, (exp(abs(params.alpha2*(s-0.5*params.dr).^params.alpha3)))));
%%
v = diff(out.LXB);
plot(out.t(1:end-1), v)
%%
s = -0.1:0.01:0.1;clf;
% plot(s, (exp(abs(1000*(s+params.dr).^2.2))))
plot(s, min(params.alpha2, (exp(abs(params.alpha2*(s+params.dr).^params.alpha3))))); % P2 to P3)
ylim([1, 10])
%%
params0.PlotEachSeparately = true;
params0.Slim_l = 1.65;
params0.Slim_r = 2.205;
params0.dS = 0.005;
tic
optimfun = @(g)evaluateBakersExp(exp(g), params0);
optimfun(g)
toc
%% GA optim result - pretty bad
% params0.g = [745.1430  808.6788  500.7677  513.4525  734.6267  244.5222   31.5362    1.4373  347.2581  859.2466];
params0.PlotEachSeparately = true;
tic
% evaluateBakersExp(g, params0)
RunBakersExp;
toc
%% new try
params0.mods = {'kah', 'kadh', 'ka', 'kd', 'k1', 'k_1', 'k2', 'kstiff2', 'alpha2', 'alpha3', 'kSE'};
params0.g = ones(1, length(params0.mods));
%% Attempt on GA
% parpool
params0.mods = {'kah', 'ka', 'kd', 'k1', 'k_1', 'k2', 'kstiff1', 'kstiff2', 'alpha0', 'dr0', 'alpha1', 'dr1', 'alpha_1', 'dr_1', 'alpha2', 'dr2', 'alpha3'};
ga_Opts = optimoptions('ga', ...
    'PlotFcn','gaplotbestf', ...
    'PopulationSize',16*4, ...            % 250
    'MaxStallGenerations',10, ...  % 10
    'UseParallel',true, 'Display', 'iter' ...
    );


% params0.mods = {"kstiff1", "kstiff2", "k1", "k2", "k_2", "k3", "s3", "alpha3", "sigma0", "ksr0"};
Ng = length(params0.mods);
% params0.UseKstiff3 = false;
% params0.SaveBest = true;
% params0.ghostSave = false;


params0.PlotEachSeparately = false;

% evaluateBakersExp(g, params0)

optimfun = @(g)evaluateBakersExp(g, params0);

% default bounds struct - based on the mods
upp = 1e4; lwr = 0;
p_lb = params0; for i = 1:length(params0.mods), p_lb.(params0.mods{i}) = params0.(params0.mods{i})*lwr; end
p_ub = params0; for i = 1:length(params0.mods), p_ub.(params0.mods{i}) = params0.(params0.mods{i})*upp; end

% specific bounds - this is universal and non-blocking. 
% Might not be used, if not within 'mods' - this is the whole reason to make all this mess
p_lb.alpha0 = -1e3;
p_lb.dr0 = -10;
p_lb.alpha1 = -1e3;
p_lb.dr1 = -10;
p_lb.alpha_1 = -1e3;
p_lb.dr_1 = -10;
p_lb.alpha2 = -1e3;
p_lb.dr2 = -10;
p_lb.alpha3 = -1e3;

% convert back to vector of ratios
for i = 1:length(params0.mods)
    lb(i) = p_lb.(params0.mods{i})/params0.(params0.mods{i});
    ub(i) = p_ub.(params0.mods{i})/params0.(params0.mods{i});
end

[p_OptimGA,Res_OptimGA,~,~,FinPopGA,FinScoreGA] = ...
    ga(optimfun,Ng, ...
    [],[],[],[],lb,ub,[],ga_Opts);

save env2;

optimfun(p_OptimGA)

x2 = fminsearch(optimfun, p_OptimGA, options)
save x3
optimfun(x2)

%%
% params0.mods = {'kah', 'ka', 'kd', 'k1', 'k_1', 'k2', 'kstiff1', 'kstiff2', 'alpha0', 'dr0', 'alpha1', 'dr1', 'alpha_1', 'dr_1', 'alpha2', 'dr2', 'alpha3'};
% params0.g = x2;
% params0 = getParams(params0, params0.g, true, true);
params0.PlotEachSeparately = true;
figure(123); clf;
RunBakersExp;
nexttile;
plot(out.t, out.p1_0, '-', out.t, out.p2_0, '-', out.t, out.p3_0, '-',out.t, out.PuATP, '-',out.t, out.PuR, '-', LineWidth=1.5, LineStyle='-')
legend('P1','P2','P3','PuATP','PuR')

%% ga params outcome
% writeParamsToFile('gaOutparams.csv', params0, {'kah', 'ka', 'kd', 'k1', 'k_1', 'k2', 'kstiff1', 'kstiff2', 'alpha0', 'dr0', 'alpha1', 'dr1', 'alpha_1', 'dr_1', 'alpha2', 'dr2', 'alpha3'})
params0.kah = 215.70;
params0.ka = 138083.63;
params0.kd = -36144.22;
params0.k1 = 576357.92;
params0.k_1 = 83864.91;
params0.k2 = 5641.24;
params0.kstiff1 = -1003339;
params0.kstiff2 = 107586;
params0.alpha0 = NaN
params0.dr0 = -2.21;
params0.alpha1 = 0;
params0.dr1 = 0.04;
params0.alpha_1 = 0;
params0.dr_1 = -0.17;
params0.alpha2 = -0;
params0.dr2 = 13.43;
params0.alpha3 = -24.44;
figure(423); clf;
RunBakersExp;
hold on;plot(out.t, out.XB_TORs)
nexttile;
plot(out.t, out.p1_0, '-', out.t, out.p2_0, '-', out.t, out.p3_0, '-',out.t, out.PuATP, '-',out.t, out.PuR, '-', LineWidth=1.5, LineStyle='-')
legend('P1','P2','P3','PuATP','PuR')

xlim([2.6 3.05])
title(sprintf('XB TOR %g', out.XB_TORs(end)))
% allProbs = sum([out.p1_0; out.p2_0; out.p3_0; out.PuATP; out.PuR]);
% plot(out.t, allProbs)
StatesInTime;



%% Calculation of max speed
dr = params.dr; % um
v = 10; % ML/s = um/s in half-sarcomere
% to produce F we need to have sum < dr, so we need to unbind the rest
% thus, zero force at sum = dr and all binded are unbinded
% tor = ; 

%% read all, write selected only
function [params, row_names] = UpdateValuesFromFile(filename, params)
    ti = readtable(filename, 'ReadRowNames',true, 'ReadVariableNames',true, 'CommentStyle', '#',NumHeaderLines=0);
    
    % update from file
    row_names = ti.Properties.RowNames;
    for i_row = 1:length(ti.Properties.RowNames)
        params.(row_names{i_row}) = ti.Value(i_row);
    end
    
    % reset the modifiers - already updated
    params.mods = [];
    params.g = [];
end
%% write to file
function writeParamsToFile(filename, params, modNames)
    modNames_exclude = [];
    if nargin < 3
        modNames = fieldnames(params);
    end
    
    row_vals = {};
    row_names = {};
    i_rowc = 0; % cleared index - without the excluded
    for i_row = 1:length(modNames)
        if any(strcmp(modNames_exclude, modNames{i_row})) || length(params.(modNames{i_row})) > 1
            fprintf("Skipping %s\n", modNames{i_row});
            continue;
        end
        i_rowc = i_rowc + 1;
        row_names{i_rowc} = modNames{i_row}; 
        row_vals{i_rowc} = params.(modNames{i_row});
    end
    
    to = cell2table(row_vals', 'VariableNames', {'Value'}, 'RowNames',string(row_names));
    writetable(to, filename, 'WriteRowNames', true)
end
% to.Properties.VariableNames = {'ParamName', 'value'}

