%% Handtune alternative model, based on mant-ATP knowledge

clf;
LoadData;
params0 = getParams();
ModelParamsFVOptimSurro
params0.RunForceVelocity = true;
params0.RunForceVelocityTime = false;
params0.BreakOnODEUnstable = false;
params0.RunSlack = false;

params0.FV_velocities = [0 -0.5 -2 -3];
% params0.FV_velocities = [-0.5];
params0.UseSuperRelaxed = true;
params0.UseSuperRelaxedADP = true;

% new params
params0.kamh = 0.1*params0.kah;
params0.ksr2srd = params0.kah;
params0.ksrd2sr = params0.kamh;

params0.ksr0 = 1;
params0.kmsr = 0;
params0.kmsrd = 1;
params0.ksrd = 0;
params0.sigma_srd1 = 10

params0.justPlotStateTransitionsFlag = false;
RunBakersExp

xlim([-inf inf ])

%% 
figure(41);clf;
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
params0.PieceWiseStrainDepParams__3 = params0.PieceWiseStrainDepParams(3);
params0.PieceWiseStrainDepParams__4 = params0.PieceWiseStrainDepParams(4);
params0.PieceWiseStrainDepParams__5 = params0.PieceWiseStrainDepParams(5);
params0.PieceWiseStrainDepX__3 = params0.PieceWiseStrainDepX(3);
params0.PieceWiseStrainDepX__4 = params0.PieceWiseStrainDepX(4);
params0.PieceWiseStrainDepX__5 = params0.PieceWiseStrainDepX(5);

params0.mods = paramNames;
params0.g = ones(1, length(paramNames));

g0 = x;
params0.PlotEachSeparately = false;
params0.justPlotStateTransitionsFlag = false;
params0.BreakOnODEUnstable = true;

params0.MaxRunTime = 10;

optimfun = @(g)evaluateBakersExp(g, params0, false);
%% RUN optim

options = optimoptions('surrogateopt','Display','iter', 'MaxTime', 60*60, 'UseParallel',false, 'PlotFcn', 'surrogateoptplot', 'InitialPoints', g0, MaxFunctionEvaluations=600);
% p = parpool('threads', 4)
% p = parpool('processes', 4)
% delete(p)
[x,fval,exitflag,output,trials] = surrogateopt(optimfun, g0*0.1,g0*10, options);
save env
% writeParamsToMFile('../ModelOptParams/ModelParamsFVOptimSRXTDSurro.m', params0, '', 'Optim of force velocity curve. Good but oscillating.');

%% test

figure(553);
clf;
% params0.g = x;
% params0 = params;
% params0.ghostSave = '';
% params0.ghostLoad = 'FVSurro_base';
params0.FV_velocities = -[0, 0.5, 1, 2, 3, 4, 5];
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
params0.ghostLoad = 'ouja';
params0.UseSuperRelaxed = false;
params0.UseSuperRelaxedADP = false;
% params0.dS = 0.0011;
% params0.MaxStrainArraySize = 60;
% params0.A1AttachmentWidth = 0.003;
params0.MaxSpaceExtensionCount = inf;
params = getParams(params0);
params.PieceWiseStrainDep

RunBakersExp;
%%
out = outs(end);
StatesInTime

% plot(out.t, out.p1_0 + out.p2_0 + out.PuATP + out.PuR -1)
%% Fit a FV metric
FVfit_fun = @(a, b, c, x) a +b*exp(-x*c);
FV_x = -params0.FV_velocities';FV_y = F_active';
FV_t = 0:.1:10;
[fvfit fvgood] = fit(FV_x(2:end), FV_y(2:end), FVfit_fun, 'StartPoint', [1, 100, 0.5], 'Lower', [0 0 0])

clf;plot(FV_x, FV_y, 'x', FV_t, fvfit(FV_t), '--', LineWidth=1.5);
metric = (fvfit(0)-FV_y(1))/FV_y(1);
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