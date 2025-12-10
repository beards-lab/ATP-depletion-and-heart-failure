%% Handtune alternative model, based on mant-ATP knowledge

clf;
LoadData;
params0 = getParams();
ModelParamsFVOptimSurro
params0.UseA1AttachmentKernel = false;
params0.RunForceVelocity = true;
params0.RunForceVelocityTime = false;
params0.BreakOnODEUnstable = false;
params0.RunSlack = false;

params0.FV_velocities = -[0.5 2 3 5];
% params0.FV_velocities = [-0.5];
params0.UseSuperRelaxed = true;
params0.UseSuperRelaxedADP = true;

% new params
params0.kamh = 0.1*params0.kah;
params0.ksr2srd = params0.kah;
params0.ksrd2sr = params0.kamh;

params0.ksr0 = 10;
params0.kmsr = 1;
params0.kmsrd = 1;
params0.ksrd = 1;
params0.sigma1 = 20;
params0.sigma2 = 40;
params0.sigma_srd1 = 20;
params0.sigma_srd2 = 40;

params0.PieceWiseStrainDepParams = [50 1.61509116782026 2.08654906258156 50 50 50 50];
params0.PieceWiseStrainDepX = [0.1 0.01 0 -0.005 -0.0075 -0.01 -1];


params0.justPlotStateTransitionsFlag = false;
params0.A1AttachmentWidth = 3*params0.dS;
params0.UseA1AttachmentKernel = false;
tic
RunBakersExp
toc

xlim([-inf inf ])

%% 
figure(41);clf;
out = outs(end);
StatesInTime

%% plot the PWSD
figure(42);clf; 

params = getParams(params0, x, false, true);
params.PieceWiseStrainDepParams(4) = params.PieceWiseStrainDepParams(3);
params = getParams(params, x, false, true);
xxx = linspace(-params.dr*3, params.dr*3, 100);
plot(xxx, ppval(params.PieceWiseStrainDep,xxx), '-', params.PieceWiseStrainDepX, params.PieceWiseStrainDepParams, 'o', LineWidth=1.5, MarkerSize=12)

%% parametrization 
disp('params:')
params.dr2 - params.dr
params.sigma1
params.sigma2
params.sigma_srd1
params.sigma_srd2
%%
i = 2;
(outs(i).SRD(end) + outs(i).SR(end))*100
out = outs(i)

%% prepare optim
paramNames = {'kstiff2','k_pas','ka','kah','dr2','k2','k1','PieceWiseStrainDepParams__3','PieceWiseStrainDepParams__2','PieceWiseStrainDepParams__4','mu','kstiff1','kd'};
paramNames = {'PieceWiseStrainDepParams__3','PieceWiseStrainDepParams__4','PieceWiseStrainDepParams__5','PieceWiseStrainDepX__3','PieceWiseStrainDepX__4','PieceWiseStrainDepX__5'};
paramNames = {'dr2'                        ,
'gamma'                      ,
'estiff'                     ,
'k2'                         ,
'k_pas'                      ,
'ekSE'                       ,
'PieceWiseStrainDepParams__2',
'kstiff2'                    ,
'PieceWiseStrainDepX__4'     ,
'PieceWiseStrainDepParams__4',
'PieceWiseStrainDepParams__3',
'ka'                         };
% params0.PieceWiseStrainDepParams__3 = params0.PieceWiseStrainDepParams(3);
% params0.PieceWiseStrainDepParams__4 = params0.PieceWiseStrainDepParams(4);
% params0.PieceWiseStrainDepParams__5 = params0.PieceWiseStrainDepParams(5);
% params0.PieceWiseStrainDepX__3 = params0.PieceWiseStrainDepX(3);
% params0.PieceWiseStrainDepX__4 = params0.PieceWiseStrainDepX(4);
% params0.PieceWiseStrainDepX__5 = params0.PieceWiseStrainDepX(5);

params0.mods = paramNames;
params0.g = ones(1, length(paramNames));

% g0 = x;
g0 = params0.g;
params0.PlotEachSeparately = false;
params0.justPlotStateTransitionsFlag = false;
params0.BreakOnODEUnstable = true;

params0.MaxRunTime = 10;

optimfun = @(g)evaluateBakersExp(g, params0, false);
%% RUN fminsearch optim
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.1, 'PlotFcns', @optimplotfval, 'MaxIter', 1500);

% optimfun = @(pw)evalSDP(params0, pw, pwsel, vel);
% p = parpool('threads', 4);
x = fminsearch(optimfun, g0, options)
% writeParamsToMFile('ModelOptParams/ModelParamsFVOptimBaseline.m', params0, '', 'Optim of force velocity curve. Good but oscillating a bit.');
%% RUN surrogate optim

options = optimoptions('surrogateopt','Display','iter', 'MaxTime', 60*60, 'UseParallel',false, 'PlotFcn', 'surrogateoptplot', 'InitialPoints', g0, MaxFunctionEvaluations=600);
% p = parpool('threads', 4)
% p = parpool('processes', 4)
% delete(p)
[x,fval,exitflag,output,trials] = surrogateopt(optimfun, g0*0.1,g0*10, options);
save env
% writeParamsToMFile('../ModelOptParams/ModelParamsFVOptimSRXTDSurro.m', params0, '', 'Optim of force velocity curve. Good but oscillating.');

%% Test after optim
clf;
params0.FV_velocities = -[0 0.5 1 2 3 4 5 6];
params0.mods
params0.g = x;
params0.PlotEachSeparately = true;
RunBakersExp;

out = outs(1);
StatesInTime

%% test
LoadData
params0 = getParams();
ModelParamsFVOptimSRXTDSurro
figure(553);
clf;
% params0.g = x;
% params0 = params;
% params0.ghostSave = '';
% params0.ghostLoad = 'FVSurro_base';
params0.FV_velocities = -[0, 0.5, 2, 5];
% params0.FV_velocities = -[0, 0.5, 1, 2, 3, 4, 5];
% params0.g = g0;
params0.PlotEachSeparately = true;
params0.justPlotStateTransitionsFlag = false;
params0.RunForceVelocity = false;
params0.RunForceVelocity = true;
params0.RunForceVelocityTime = false;
% params0.RunForceVelocityTime = true;

params0.RunSlack = false;
% params0.RunSlack = false;
params0.ghostSave = '';
params0.ghostLoad = 'oujaDoubleBaseNiceFit';
params0.UseSuperRelaxed = true;
params0.UseSuperRelaxedADP = true;
params0.dS = 0.0011*2*2;
params0.MaxStrainArraySize = 40;
params0.A1AttachmentWidth = 0.003;
params0.MaxSpaceExtensionCount = inf;
% params = getParams(params0);
% params.PieceWiseStrainDep
% getParams
params0.PieceWiseStrainDepX =      [0.1000    0.0100       0   -0.0005   -0.01   -0.100];
params0.PieceWiseStrainDepParams = [50.0000    1.6*2    3    0.1      50      50.0000];
params0.PieceWiseStrainDepX =      [0.1000    0.0100  0   -0.01   -0.100];
params0.PieceWiseStrainDepParams = [50.0000    3      3    50      50.0000];
params0.kd = 10;
params0.dr2 = 2*params0.dr;

pwsdX2 = [0.1000    0.0100       0    -0.01   -0.100] - 0.02;
pwsdP2 = [50    1.6    3.537    50      50.0000];

% params0.PieceWiseStrainDepX = pwsdX2 + 0.01;
% params0.PieceWiseStrainDepParams = pwsdP2;
params0.PieceWiseStrainDepX = pwsdX2 + 0.02;
params0.PieceWiseStrainDepParams = pwsdP2;

% params0.PieceWiseStrainDep2 = pchip(pwsdX2, pwsdP2);
params0.justPlotStateTransitionsFlag = false;
params0.UseStrainDep4R1D = false;
% params0.kd_sigma = 0.002;
params0.EvalFeatures = false;

params0.ksr0 = 40;
params0.sigma2 = 100;
params0.kmsr = 1;
params0.sigma1 = 12;


params0.ksrd = 40;
params0.sigma_srd2 = 100;
params0.kmsrd = 1;
params0.sigma_srd1 = 12;

params0.OptimizeFVInit = false;

params0.MaxStrainArraySize = 40;
params0.dS = 0.00021;

params0.RunForceVelocity = true;
params0.RunSlack = false;
params0.FV_velocities = -[0];

tic
RunBakersExp;
% StatesInTime
toc
%%
out = outs(end);
figure(808);clf;StatesInTime
%%

nexttile();
plot(outs(1).t, outs(1).Force, outs(2).t, outs(2).Force, outs(3).t, outs(3).Force)
% writeParamsToMFile('ModelOptParams/ModelParamsFVOptimValley.m', params0, '', ['Optim of force velocity curve. Good but oscillating, hard to explain valley in R12. ' newline ... 
%     'pwsdX2 = [0.1000    0.0100       0    -0.01   -0.100] - 0.02; pwsdP2 = [50    1.6    3.537    50      50.0000];' newline ...
%     'params0.PieceWiseStrainDep2 = pchip(pwsdX2, pwsdP2);']);

%%
figure(182); clf;
params0.RunForceVelocityTime = false;
params0.RunForceVelocity = false;
params0.RunSlack = true;
params0.MaxStrainArraySize = 120; 
params0.MaxRunTime = 30;
params0.ShowArrayShiftWarnings = true;
params0.UseSuperRelaxed = true;
params0.UseSuperRelaxedADP = true;
RunBakersExp;
%%
out = outs(1);
StatesInTime


% plot(out.t, out.p1_0 + out.p2_0 + out.PuATP + out.PuR -1)
%% Fit a FV metric
FVfit_fun = @(a, b, c, x) a +b*exp(-x*c);
FV_t = 0:.1:6;

clf; hold on;
legstr = [];

win = 2:8;
win = [2 3 4 5 6];
scaleTo = 2;

FV_x = -params0.FV_velocities';FV_y = F_active'/F_active(scaleTo);
[fvfit fvgood] = fit(FV_x(win), FV_y(win), FVfit_fun, 'StartPoint', [1, 100, 0.5], 'Lower', [0 0 0]);
plot(FV_x, FV_y, 'x', FV_t, fvfit(FV_t), '-', LineWidth=1.5, MarkerSize=12);
% plot(FV_x(win), FV_y(win), 'x', LineWidth=3, MarkerSize=18); 
legstr = [legstr, "Simulated points", sprintf("Fit: a=%.2f,b=%0.1f,c=%0.2f", fvfit.a, fvfit.b, fvfit.c)];


FV_x = Data_ATP(:, 1);FV_y = Data_ATP(:, 2)/Data_ATP(scaleTo, 2);
[fvfit fvgood] = fit(FV_x(win), FV_y(win), FVfit_fun, 'StartPoint', [1, 100, 0.5], 'Lower', [0 0 0]);
plot(FV_x, FV_y, 'x', FV_t, fvfit(FV_t), '-', LineWidth=1.5, MarkerSize=12);
% plot(FV_x(win), FV_y(win), 'x', LineWidth=3, MarkerSize=18); 
legstr = [legstr, "Data points", sprintf("Fit: a=%.2f,b=%0.1f,c=%0.2f", fvfit.a, fvfit.b, fvfit.c)];

% FV_x = Data_ATP(:, 1);FV_y = Data_ATP(:, 3)/Data_ATP(scaleTo, 3);
% [fvfit fvgood] = fit(FV_x(win), FV_y(win), FVfit_fun, 'StartPoint', [1, 100, 0.5], 'Lower', [0 0 0]);
% plot(FV_x, FV_y, 'x', FV_t, fvfit(FV_t), '-.', LineWidth=1.5, MarkerSize=12);
% % plot(FV_x(win), FV_y(win), 'x', LineWidth=3, MarkerSize=18); 
% legstr = [legstr, "Data points 4 mM", sprintf("Fit: a=%.2f,b=%0.1f,c=%0.2f", fvfit.a, fvfit.b, fvfit.c)];

FV_x = Data_ATP(:, 1);FV_y = Data_ATP(:, 4)/Data_ATP(scaleTo, 4);
[fvfit fvgood] = fit(FV_x(win), FV_y(win), FVfit_fun, 'StartPoint', [1, 100, 0.5], 'Lower', [0 0 0]);
plot(FV_x, FV_y, 'x', FV_t, fvfit(FV_t), '--', LineWidth=1.5, MarkerSize=12);
% plot(FV_x(win), FV_y(win), 'x', LineWidth=3, MarkerSize=18); 
legstr = [legstr, "Data points 2 mM", sprintf("Fit: a=%.2f,b=%0.1f,c=%0.2f", fvfit.a, fvfit.b, fvfit.c)];

FV_x = Data_ATP_AB2021(:, 1);FV_y = Data_ATP_AB2021(:, 2)/Data_ATP_AB2021(scaleTo, 2);
[fvfit fvgood] = fit(FV_x(win), FV_y(win), FVfit_fun, 'StartPoint', [1, 100, 0.5], 'Lower', [0 0 0]);
plot(FV_x, FV_y, '^', FV_t, fvfit(FV_t), '-', LineWidth=1.5);
% plot(FV_x(win), FV_y(win), 'x', LineWidth=3, MarkerSize=18);
legstr = [legstr, "Data 2021 points", sprintf("Fit: a=%.2f,b=%0.1f,c=%0.2f", fvfit.a, fvfit.b, fvfit.c)];

FV_x = Data_ATP_AB2021(:, 1);FV_y = Data_ATP_AB2021(:, 3)/Data_ATP_AB2021(scaleTo, 3);
[fvfit fvgood] = fit(FV_x(win), FV_y(win), FVfit_fun, 'StartPoint', [1, 100, 0.5], 'Lower', [0 0 0])
plot(FV_x, FV_y, '>', FV_t, fvfit(FV_t), '--', LineWidth=1.5);
% plot(FV_x(win), FV_y(win), 'x', LineWidth=3, MarkerSize=18);
legstr = [legstr, "Data 2021 ATP 2 mM points", sprintf("Fit: a=%.2f,b=%0.1f,c=%0.2f", fvfit.a, fvfit.b, fvfit.c)];

legend(legstr);

metric = (fvfit(0)-FV_y(1))/FV_y(1);
%% construct feature metric - test this!

% construct parameter set to vary
paramNames = {'kstiff2','k_pas','ka','kah','dr2','k2','k1'};
% paramNames = {'baseline_dummy', 'kstiff2','k_pas'};
paramNames = {'baseline_dummy', 'PieceWiseStrainDepX__3','PieceWiseStrainDepX__4','PieceWiseStrainDepX__5','kstiff2','k_pas','ka','kah','kamh','dr2','k2','k1','PieceWiseStrainDepParams__3','PieceWiseStrainDepParams__2','PieceWiseStrainDepParams__4','mu','kstiff1','kd','sigma1','sigma2','sigma_srd1','sigma_srd2','dr','dS','estiff','k_pas','gamma','kSE','ekSE'};
% construct feature set to evaluate
% fn = {'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'SLslack|t0|0', 'Am|SLslack|', 'peak1_y|SLslack|0', 'peak1_y|v_restretch|0','peak1_dSL', 'peak2|v_restretch', 'vall_t|v_restretch|0.1', 'vall_y', 'vall2_dy|v_restretch|0', 'ovrsht_dy|_|0', 'steady'};
% reduced set


% simulate force-isovelocity
% take only non-zero velocities with > 10 kPa
% params0.FV_velocities = -[0.5, 1, 3, 4];
% params0.RunForceVelocity = false;
% params0.RunSlack = false;
% params0.RunForceLengthEstim = false;
params0.EvalFeatures = false;
params0.PlotEachSeparately = false;

% params0.RunForceVelocity = true;
% params0.RunSlack = true;
params0.MaxRunTime = 30;
params0.baseline_dummy = 0;

% params0.EvalFeatures = true;
% params0.RunForceLengthEstim = true;
fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y|SLslack', 'peak1_dSL', 'peak2', 'ovrsht_dy|_|0', 'steady', 'XTOR'};
fn = {'FV_f|FV_v'};

featureMatrixPlus = zeros(length(paramNames), length(fn));
% params_base = params0;
for i_param = 1:length(paramNames)
    try
        fprintf('Running %s..', paramNames{i_param});
        params0.mods = paramNames(i_param);
        params0.g = 1.01;
        
        figure(333); clf;
        
        RunBakersExp;
        % % check the init
        % [costs, weights, cost] = evalFeatureCost(features_data, features_model, fn, 1);
        costs = E;
        featureMatrixPlus(i_param, :) = costs;
        % savedFeats{i_param} = features_model;
        fprintf('Done \n');
    catch e
        fprintf('failed. [%s] \n', e.message);
    end
end
%%


%% RUN THE OPTIM SCRIPT FOR LSQNON USING CUSTOM JACOBIAN

% main_lsqnonlin_script.m

% 1. Setup
paramNames = {'kstiff2','k_pas','ka','kah','dr2','k2','k1'};
params0.g = ones(1, length(paramNames));
params0.mods = paramNames;
params0.FV_velocities = -[0.5, 1, 2, 4];
params0.RunSlack = false;
params0.RunForceVelocity = true;
% matchStructFields(params0, 'Run*', true)
params0.MaxStrainArraySize = 40;
params0.MaxRunTime = 30;
params0.EvalFeatures = true;
params0.BreakOnODEUnstable = true;
params0.PlotEachSeparately = false;

tic
RunBakersExp;
toc
%%
[Residuals, weights, cost] = evalFeatureCost(features_data, features_model, params0.fn, 1);

%%
P0 = ones(1, length(paramNames)); 
LB = P0*0.1; UB = P0*10;

params0.fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y|SLslack', 'peak1_dSL', 'peak2', 'ovrsht_dy|_|0', 'steady', 'XTOR'};

params0.RunSlack = false;
params0.MaxStrainArraySize = 40;
params0.fn = {'FV_f|FV_v'};


% 2. Options
options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true);

% CRITICAL: Tell lsqnonlin your function provides the Jacobian
% options.SpecifyObjectiveAndJacobian = true; 
options.Jacobian = 'on';


% Set MaxFunevals carefully. You only need ~1 + M model calls per outer update.
% If you want 10 Jacobian updates, you need 10 * (1 + 15) = 160 calls.
options.MaxFunctionEvaluations = 200; 
options.Display = 'iter';

% 3. Run the Solver
fprintf('Starting Decoupled lsqnonlin Optimization...\n');

% lsqnonlin handles the trust-region/Levenberg-Marquardt logic using
% the Residuals and the returned Jacobian (S).
[P_optimized, Resnorm, Residuals] = lsqnonlin(@ResidualAndJacobian, P0, LB, UB, options, params0);

fprintf('\nOptimized Parameters: %s\n', mat2str(P_optimized, 4));
fprintf('Minimum Total Cost (Sum of Squares): %.4f\n', Resnorm);
%%


figure(4);
featureMatrixPlusN = featureMatrixPlus - featureMatrixPlus(1, :);

[vals, paramPos] = sort(abs(featureMatrixPlusN), 1, "desc")
featureMatrixPlusNSorted = featureMatrixPlusN(paramPos, :)./featureMatrixPlusN(paramPos(1), :)*100;

imagesc(featureMatrixPlusNSorted);
colorbar;
axis ij tight;

set(gca, 'XTick', 1:length(fn), ...
         'XTickLabel', fn, ...
         'YTick', 1:length(paramNames), ...
         'YTickLabel', paramNames(paramPos), 'TickLabelInterpreter', 'none');

xtickangle(45);
title('Parameter Sensitivities');

% Print numeric values in each cell
for i = 1:size(featureMatrixPlusNSorted,1)
    for j = 1:size(featureMatrixPlusNSorted,2)
        text(j, i, sprintf('%.3g', featureMatrixPlusNSorted(i,j)), ...
            'HorizontalAlignment', 'center', ...
            'Color', 'w');
    end
end

paramNames(paramPos)'

%%


clf;
plotFeatures(features_data, savedFeats{1}, [], fn);

%%
params0.A1AttachmentWidth = 2*params0.dS;
params0 = getParams(params0);
params0.UseA1AttachmentKernel = true;
params0.MaxStrainArraySize = 40;
params0.UseSuperRelaxed = false;
params0.UseSuperRelaxedADP = false;
params0.justPlotStateTransitionsFlag = true;
tic
RunBakersExp; 
toc
StatesInTime;

%%


clf;
pp = params.PieceWiseStrainDep;

pwsdX2 = params.PieceWiseStrainDepX;
pwsdP2 = params.PieceWiseStrainDepParams;

    % params.PieceWiseStrainDepFun = @(x)pchip(pwsdX, pwsdP, x);   
pwsdX2(3) = -0.001;
params.PieceWiseStrainDep2 = pchip(pwsdX2, pwsdP2);
s = -2*params.dr:params.dr/100:params.dr*2;
plot(params.PieceWiseStrainDepX, params.PieceWiseStrainDepParams, 'o', s, ppval(params.PieceWiseStrainDep, s), ...
    pwsdX2, pwsdP2,'x',s, ppval(params.PieceWiseStrainDep2, s));
xlim([min(s) max(s)])
%% atp 8 new data
FV_winm = 8;
FV_winp = 0;
FV_max = length(datatable)

FV_y = arrayfun(@(k) ...
    mean(datatable(k-FV_winm:min(FV_max, k+FV_winp), 3)), ...
    smashed2_0_data);

clf; plot(datatable(:, 1), datatable(:, 3), '-', datatable(:, 1), datatable(:, 2), '|-');
hold on;
% Draw thick horizontal lines showing the averaging window
for idx = smashed2_0_data(:)'   % ensure row vector
    span = idx-FV_winm : min(FV_max, idx+FV_winp);
    FV_ms = mean(datatable(span,3));
    plot([datatable(idx-FV_winm,1), datatable(idx,1)], [FV_ms, FV_ms], 'x-', 'LineWidth', 4);
end
%%
[FV_x, i_FVs] = sort(speeds) ;
FV_ys = FV_y(i_FVs);

fvfit = fit(FV_x(2:end), FV_ys(2:end), FVfit_fun, 'StartPoint', [1, 100, 0.5])

clf;plot(FV_x, FV_ys, 'x', FV_t, fvfit(FV_t), '--', LineWidth=1.5);
metric = (fvfit(0)-FV_ys(1))/FV_ys(1)

%% 
figure(222);clf;
LoadData
ModelParamsFVOptimSRXTDSurro
params0.MaxSpaceExtensionCount = inf;
params0.RunForceVelocity = true;
params0.RunForceVelocityTime = false;
params0.RunSlack = false;
params0.BreakOnODEUnstable = false;
params0.justPlotStateTransitionsFlag = false;
RunBakersExp

legend('Data', 'Model')

StatesInTime

%% test: accessing struct IS costly
N = 1e7;
tic
for i = 1:N
    % a1 = params.dS;
    % a2 = params.dS;
    % a = max(0, a1);
    a = a1;
end
toc

%% Compare force velocity and slack datasets

fv_data = load('data/bakers_isovelocity.mat', 'datatable', 'velocitytable').datatable;
fv_win2_2 = fv_data(:, 1) > 3.45 & fv_data(:, 1) < 3.6;
% fv_win2_2 = fv_data(:, 1) > 0 & fv_data(:, 1) < 0.3;
fv_win2_0 = fv_data(:, 1) > 3.9 & fv_data(:, 1) < 4;
slack_data = load('data/bakers_slack8mM_all.mat').datatable;
slack_win2_2 = slack_data(:, 1) > 2.72 & slack_data(:, 1) < 2.75;
% slack_win2_2 = slack_data(:, 1) > 0 & slack_data(:, 1) < 1.15;
slack_win2_0 = slack_data(:, 1) > 3.04 & slack_data(:, 1) < 3.6;

figure(181);clf; 
plot(fv_data(:, 1),fv_data(:, 3), slack_data(:, 1), slack_data(:, 3), '-', LineWidth=1.5);hold on;
plot(fv_data(fv_win2_2 | fv_win2_0, 1),fv_data(fv_win2_2 | fv_win2_0, 3), slack_data(slack_win2_2 | slack_win2_0, 1), slack_data(slack_win2_2 | slack_win2_0, 3), '-', LineWidth=1.5);

R1 = mean(fv_data(fv_win2_0, 3))/mean(fv_data(fv_win2_2, 3));
text(max(fv_data(:, 1)), mean(fv_data(fv_win2_0, 3)) + 5, sprintf('%0.0f%%', R1*100), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontWeight','bold')
R2 = mean(slack_data(slack_win2_0, 3))/mean(slack_data(slack_win2_2, 3));
text(max(slack_data(:, 1)), mean(slack_data(slack_win2_0, 3)) + 5, sprintf('%0.0f%%', R2*100), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontWeight','bold')
xlabel('Time (s)');ylabel('Force (kPa)');
legend('Isovelocity', 'Slack');