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

%% Extend optim

params0.RunSlack = true;
params0.RunForceVelocity = true;
params0.EvalFeatures = true;
params0.fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y|SLslack', 'peak1_dSL', 'peak2', 'ovrsht_dy|_|0', 'steady', 'XTOR'};


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
paramNames = {'estiff'                    , ...
'dr2'                        , ...
'ekSE'                       , ...
'PieceWiseStrainDepParams__2', ...
'kstiff2'                    , ...
'PieceWiseStrainDepX__2'     , ...
'kSE'                        , ...
'kah'                        , ...
'kstiff1'                    , ...
'PieceWiseStrainDepX__4'     , ...
'kmsrd'                      , ...
'sigma_srd1'                 , ...
'kd'                         , ...
'ka'                         , ...
'sigma1'                     , ...
'PieceWiseStrainDepParams__3', ...
'ksr0'                       , ...
'kamh'                       , ...
'PieceWiseStrainDepX__5'     , ...
'sigma2'                     , ...
'PieceWiseStrainDepParams__4'};

paramNames = {'k_pas'                      ,...
'gamma'                      ,...
'sigma_srd1'                 ,...
'PieceWiseStrainDepX__4'     ,...
'ksrd'                       ,...
'ksr0'                       ,...
'k2'                         ,...
'kmsrd'                      ,...
'sigma1'                     ,      ...
'PieceWiseStrainDepParams__5'};

paramNames = {
 'PieceWiseStrainDepX__2', ... 
'PieceWiseStrainDepParams__2', ... 
'PieceWiseStrainDepParams__3', ... 
'PieceWiseStrainDep2X__2', ... 
'PieceWiseStrainDep2X__4', ... 
'PieceWiseStrainDep2Params__2', ... 
'PieceWiseStrainDep2Params__3', ... 
'PieceWiseStrainDepR1DX__2', ... 
'PieceWiseStrainDepR1DX__3', ... 
'PieceWiseStrainDepR1DParams__2', ... 
'PieceWiseStrainDepR1DParams__3', ... 
'PieceWiseStrainDepR21X__2', ... 
'PieceWiseStrainDepR21X__3', ... 
'PieceWiseStrainDepR21Params__2', ... 
'PieceWiseStrainDepR21Params__3'};

paramNames = {'kd'                            ,...
'PieceWiseStrainDepX__2'        ,...
'sigma_srd1'                    ,...
'PieceWiseStrainDepR21Params__2',...
'ksr2srd'                       ,...
'k1'                            ,...
'sigma_srd2'                    ,...
'kstiff1'                       ,...
'k2'                            ,...
'gamma'                         ,...
'PieceWiseStrainDep2Params__3'  ,...
'PieceWiseStrainDepR1DX__3'     ,...
'ksrd2sr'                       ,...
'PieceWiseStrainDepR1DParams__3',...
'sigma1'                        };

paramNames = {'estiff'                       , ...}
'kstiff2'                       , ...
'ekSE'                          , ...
'kah'                           , ...
'PieceWiseStrainDep2X__2'       , ...
'PieceWiseStrainDep2Params__3'  , ...
'k2'                            , ...
'L_thin'                        , ...
'PieceWiseStrainDep2X__3'       , ...
'PieceWiseStrainDep2Params__2'  , ...
'ka'                            , ...
'PieceWiseStrainDep2Params__4'  , ...
'kSE'                           , ...
'L_thick'                       , ...
'gamma'                         , ...
'kmsrd'                         , ...
'k_pas'                         , ...
'k_1'                           , ...
'ksr2srd'                       , ...
'kd'                            , ...
'k1'                            , ...
'PieceWiseStrainDepR21Params__2', ...
'PieceWiseStrainDepParams__3'   , ...
'L_hbare'                       , ...
'kmsr'                          , ...
'PieceWiseStrainDepParams__2'   , ...
'ksrd2sr'                       , ...
'sigma_srd1'                    , ...
'sigma1'                        , ...
'PieceWiseStrainDepX__2'        };

paramNames = {'k2', 'xrate', 'estiff'                       , ...}
'kstiff2'                       , ...
'ekSE'                          , ...
'kSE'                           , ...
'PieceWiseStrainDep2X__2'       , ...
'PieceWiseStrainDep2X__3'       , ...
'PieceWiseStrainDepR21X__2', ...
'PieceWiseStrainDepR21X__3', ...
'PieceWiseStrainDepX__2'   , ...
'PieceWiseStrainDepX__3'   , ...
'PieceWiseStrainDepX__2'        };

paramNames = {
'kstiff2'                     ,
'k2'                          ,
'dr2'                         ,
'ekSE'                        ,
'kah'                         ,
'estiff'                      ,
'L_thin'                      ,
'sigma_srd1'                  ,
'PieceWiseStrainDep2X__6'     ,
'kmsr'                        ,
'PieceWiseStrainDep2Params__5',
'k_pas'                       ,
'PieceWiseStrainDepParams__4' ,
'sigma_srd2'                  ,
'kSE'                         ,
'kstiff2_n',
'kstiff1_n'
};

paramNames = {
's_threshold_R',
'slope',
'dr2'                        , 
'xrate'                      , 
'k_pas'                      , 
'kstiff2'                    , 
'L_thick'                    , 
'ka'                         , 
'k2'                         , 
'sigma2'                     , 
'L_thin'                     , 
'kah'                        , 
'ksr0'                       , 
'kamh'                       , 
'PieceWiseStrainDepR21X__1'  , 
'k1'                         , 
'kSE'                        , 
'PieceWiseStrainDepX__4'     , 
'PieceWiseStrainDep2X__5'    , 
'kstiff2_n'                  , 
'PieceWiseStrainDepParams__4', 
'k2d',
'sigma1'                     
}

% params0.kstiff1_n = params0.kstiff1;
% params0.kstiff2_n = params0.kstiff2;
% params0.PieceWiseStrainDepParams__3 = params0.PieceWiseStrainDepParams(3);
% params0.PieceWiseStrainDepParams__4 = params0.PieceWiseStrainDepParams(4);
% params0.PieceWiseStrainDepParams__5 = params0.PieceWiseStrainDepParams(5);
% params0.PieceWiseStrainDepX__3 = params0.PieceWiseStrainDepX(3);
% params0.PieceWiseStrainDepX__4 = params0.PieceWiseStrainDepX(4);
% params0.PieceWiseStrainDepX__5 = params0.PieceWiseStrainDepX(5);
%%
% params0.xrate = 1;
params0.mu_neg = params0.mu;
params0.mods = paramNames;
params0.g = ones(1, length(paramNames));

% g0 = x;
params0.RunForceVelocity = true;
g0 = params0.g;
params0.PlotEachSeparately = false;
params0.justPlotStateTransitionsFlag = false;
params0.BreakOnODEUnstable = true;

params0.MaxRunTime = 40;
% params0.fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y|SLslack', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1'};

optimfun = @(g)sum(ResidualAndJacobian(g, params0, true));

%% RUN fminsearch optim
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.1, 'PlotFcns', @optimplotfval, 'MaxIter', 15000);

% optimfun = @(pw)evalSDP(params0, pw, pwsel, vel);
% p = parpool('threads', 5);
% history = [];
x = fminsearch(optimfun, g0, options)
save env_fmin9
params0.g = x;
% writeParamsToMFile('ModelOptParams/ModelParamsFVOptimBaseline.m', params0, '', 'Optim of force velocity curve. Good but oscillating a bit.');
% writeParamsToMFile('../ModelOptParams/ModelParamsFeats_FVSlackUpdate.m', params0, '', 'Optim of force velocity curve AND the slack. Good starting vehicle.');
% writeParamsToMFile('../ModelOptParams/ModelParamsFeats_FVSlackUpdate_.m', params0, '', 'Optim of force velocity curve AND the slack. Good starting vehicle.');
% writeParamsToMFile('../ModelOptParams/ModelParamsFeats_FVSlackUpdateValley.m', params0, '', 'Optim of force velocity curve AND the slack. Take valley feat aboard, did not help much.');
% writeParamsToMFile('ModelOptParams/ModelParamsFeats_FVSlackUpdateValley2.m', params0, '', 'Second round optim of force velocity curve AND the slack. Take valley feat aboard, did not help much.');
% writeParamsToMFile('../ModelOptParams/ModelParamsFeats_ReasonableStartingPoint.m', params0, '', 'Based ona better estim for kSE');
% writeParamsToMFile('../ModelOptParams/ModelParamsFeats_ReasonableStartingPointB.m', params0, '', 'Based ona better estim for kSE and piecewise nonlinear kstiff2');
%% RUN surrogate optim

options = optimoptions('surrogateopt','Display','iter', 'MaxTime', 60*60*8, 'UseParallel',false, 'PlotFcn', 'surrogateoptplot', 'InitialPoints', g0, MaxFunctionEvaluations=600);
% p = parpool('threads', 4)
% p = parpool('processes', 4)
% delete(p

[x,fval,exitflag,output,trials] = surrogateopt(optimfun, g0*0.1,g0*10, options);
save env_surro
% writeParamsToMFile('../ModelOptParams/ModelParamsFVOptimSRXTDSurro.m', params0, '', 'Optim of force velocity curve. Good but oscillating.');
% writeParamsToMFile('../ModelOptParams/ModelParamsCorrBackwardProfiles.m', params0, '', 'Fixed the inverted backward rate strain sensitivity');
%% Test after optim
figure(802);clf;
% params0.g = x;
% 
% params0bak = params0;
% ModelOptParamsFeaturesOvernight
% params0 = getParams(params0, x, false, true);
params0.FV_velocities = -[0.5 2 3 4];
% params0.FV_velocities = -[0 0.5 1 2 3 4 5 6];
% params0.mods = [];
% params0.g = x;
params0.PlotEachSeparately = true;
% params0.PieceWiseStrainDep2Params = 
% params0.PieceWiseStrainDep2X
% params0.k2
% params0.PieceWiseStrainDep2X = [0.1000  0.03  0.0108         0   -0.0055   -0.0080 -0.0100];
% params0.PieceWiseStrainDep2Params   = [50.0000   2.5 1.1066    1.0383  50.6847   50.0000  50.0000];
params0.RunForceVelocity = true;
params0.RunSlack = true;
params0.RunSlackSegments = 'AllPar';
% params0.RunSlackSegments = 'All';
% params0.RunForceVelocityTime = false;
params0.MaxRunTime = 60;
% parpool('Threads', 5);

% params0.kstiff1 = 0.5*100*733;
% params0.kstiff2 = 0.5*1e5;
% params0.kSE = 833*3.7;
% params0.PieceWiseStrainDepR1DParams = [0.0000    1.1474    1.5953   1500.0000 1500.0000];
% params0.PieceWiseStrainDepR1DX = [-0.0500   -0.0050    0.0063    0.0200 0.05];
% params0.PieceWiseStrainDepR21Params = [0.0000    10    10   0.0000];
% params0.PieceWiseStrainDepR21X = [-0.0500   -0.0049    0.0062    0.0500] - 0.015;
% params0.k2d = 1000;
% params0.drmr = 0.01*0.8;
% params0.k2 = 132*0.6;
% params0.mu = 0.0019;
% params0.mu = 0.0019*100*10;
% params0.mu_neg = params0.mu*0.10;
% params0.mu_neg = params0.mu;
% params0.UseNegativeKstiff = true;
% params0.kstiff2 = 8e4;
% params0.kstiff1 = params0.kstiff2*0.2;

params0.justPlotStateTransitionsFlag = false;
params0.ghostLoad = '';
% params0.UseA2MechanicalRecocking = true;
% params0.xrate = 1;
% params0.kstiff1_n = params0.kstiff1*0.1;
% params0.kstiff2_n = params0.kstiff2*0.1;
% params0.kstiff2_n = params0.kstiff2*10;
% params0.kstiff2_n = 9.1413e+03;
% params0.mu = 
% params0.PieceWiseStrainDep2X__1 = 0.1;

params0.UseA2MechanicalRecocking = false;
params0.k2d = 10;
params0.drmr = params0.dr;
params0.ShowStatePlots = false;

params0.RunForceVelocity = true;
params0.RunSlack = true;
params0.RunMinStretch = false;
% params0.kSE = params0.kSE;

params0.PlotFeatureFitting = false;
params0.mu2 = 1e-9;
params0.mu = 1.9;
params0.ghostSave = '';
% params0.ghostLoad = 'mu1';
params0.ghostLoad = 'RA';
params0.UseA2AttachmentShift = true;
params0.d_actin = 5.5e-3; % 5.5 nm
params0.s_threshold_L = 1000*1.0*5.5e-3;    
params0.s_threshold_R = 1*1.5*5.5e-3;    
params0.slope = 20*4*2e2;
params0.xrate = 0.8;
params0.kstiff1 = 16000;
params0.RunSlackSegments = 'AllPar';


% testing Maxwell dashpot;
params0.kSE_M = 1000;
params0.eta_M = 10;
params0.UseMaxwellDashpot = true;
params0.RunSlackSegments = 'First';
params0.RunForceVelocity = false;


tic
RunBakersExp;
toc

%% figure(808);clf;
% params0.fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'restretchSlopeEnd'};
params0.fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};
% params0.fn = {'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};

% params0.fn = {'restretchSlopeStart'};
plotFeatures(features_data, features_model, [], params0.fn);

sum(evalFeatureCost(features_data, features_model, params0.fn, 1))

%% 
extractSlackAttributes(out.t, out.Force, out.SL, velocitytable, struct(), out, true)

%% tradeoff visualization

h = evalin('base', 'optim_history');
% h = optim_history(15:end, :);
corr_matrix = corr(h);
figure(201);clf;
imagesc(corr_matrix);
colormap(redbluecmap)
colorbar;
caxis([-1 1]);
title('Feature Correlation (Blue = Tradeoff, Red = Redundancy)');
nFeat = numel(fn);

xticks(1:nFeat);
yticks(1:nFeat);
xticklabels(fn);
yticklabels(fn);

xtickangle(45);
set(gca,'TickLabelInterpreter','none');
axis square;



% Let's say Feature 1 (Peak Force) and Feature 4 (Relaxation Time) are the issue
featA = h(:, 1); 
featB = h(:, 4);

figure(202);clf;
scatter(featA, featB, 10, 'filled', 'MarkerFaceAlpha', 0.3);
hold on;
plot(0, 0, 'rp', 'MarkerSize', 15); % The "Origin" is the perfect experimental fit
xlabel('Residual Feature A');
ylabel('Residual Feature B');
grid on;
title('Tradeoff Frontier');

%% 1. Load data from workspace
clf;
% 2. Identify the 4 "costliest" features (highest Mean Squared Error)
mse_per_feature = mean(h.^2, 1);
[~, sorted_idx] = sort(mse_per_feature, 'descend');
top_4_idx = sorted_idx(1:4);

h_top = h(:, top_4_idx);
% feature_names = cellstr("Feature " + top_4_idx);
feature_names = params0.fn(top_4_idx);

% 3. Create the Scatter Plot Matrix
figure('Color', 'w', 'Name', 'Pareto Tradeoff Matrix');
[S, AX, BigAx, H, HAx] = plotmatrix(h_top);

% 4. Formatting for clarity
title(BigAx, 'Pareto Tradeoff Analysis (Top 4 Costly Features)');

% Add labels to the diagonal and adjust axes
for i = 1:4
    ylabel(AX(i,1), feature_names{i}, 'FontWeight', 'bold', Interpreter='none');
    xlabel(AX(4,i), feature_names{i}, 'FontWeight', 'bold', Interpreter='none', Rotation=15);
    
    % Style the histograms on the diagonal
    H(i).FaceColor = [0.2 0.5 0.8]; 
    H(i).EdgeColor = 'w';
    
    for j = 1:4
        % Add a red star at (0,0) - the target experimental fit
        % grid(AX(i,j), 'on');
        
        % Color the dots by iteration (darker = later in the search)
        % This helps see if the optimizer is "drifting" or "stuck"
        if i ~= j
            hold(AX(i,j), 'on');
            % plot(AX(i,j), 0, 0, 'r.', 'MarkerSize', 12, 'LineWidth', 2);
            % xlim()
            % S(i,j).MarkerEdgeColor = 'none';
            % S(i,j).MarkerFaceColor = 'k';
            % S(i,j).MarkerFaceAlpha = 0.1; % Transparency helps see density
        end
    end
end

maxAbs = max(abs(h_top(:)));

for i = 1:4
    for j = 1:4
        if i ~= j
            xlim(AX(i,j), [0 maxAbs]);
            ylim(AX(i,j), [0 maxAbs]);
            hold(AX(i,j),'on');
            plot(AX(i,j),0,0,'r.','MarkerSize',12,'LineWidth',2);
        end
    end
end


fprintf('Top 4 costly features identified: %s\n', strjoin(feature_names, ', '));
%%
figure(4004);clf;hold on;
plot(datatable(:, 1), datatable(:, 3)/2,datatable(:, 1), datatable(:, 3)./datatable(:, 2), '|-')

%%

% original
params0.PieceWiseStrainDepX =      [0.1000     0.0107         0   -0.0052   -0.0080   -0.0100];
params0.PieceWiseStrainDepParams = [50.0000    1.1554    2.1802   50.6847   50.0000   50.0000];
% modified
params0.PieceWiseStrainDepX =      [0.1000     0.0107         0   -0.0052   -1.0000];
params0.PieceWiseStrainDepParams = [0.5          1.1554        2.1802   50    50.0000];

% original
params0.PieceWiseStrainDep2X =      [0.1000     0.0107         0   -0.0052   -0.0080   -0.0100];
params0.PieceWiseStrainDep2Params = [50.0000    1.1554    2.1802   50.6847   50.0000   50.0000];
% modified
params0.PieceWiseStrainDep2X =      [0.1000     0.0107         0   -0.0052   -0.0080   -0.0100];
params0.PieceWiseStrainDep2Params = [50.0000    1.1554    2.1802   50.6847   50.0000   50.0000];
% params0.PieceWiseStrainDep2Params = [];


params0.PieceWiseStrainDepR1DX = [-0.05 -0.005 0.005 0.05];
params0.PieceWiseStrainDepR1DParams = [50 1 1 50];
params0.kd = 10;

params0.PieceWiseStrainDepR21X = [-0.05 -0.005 0.005 0.05];
params0.PieceWiseStrainDepR21Params = [50 1 1 50];
% params0.PieceWiseStrainDepR21Params = [];




tic
RunBakersExp;
toc

[costs, weights, cost] = evalFeatureCost(features_data, features_model, params0.fn(1:end), 1);
sum(costs)
%% features_model_g = features_model;
% plotFeatures(features_data, features_model, features_model_g, params0.fn(1:end))
plotFeatures(features_data, features_model, [], params0.fn(1:end))
% out = outs(1);
% StatesInTime
% plotFeatures

%% 
figure;
StatesInTime

%%
% savedParams = params0;
params0.justPlotStateTransitionsFlag = false;
figure(8483);clf;
% params0.RunForceVelocityTime = false;
% params0.FV_velocities = [-0.5000   -2.0000   -3.0000   -4.0000];
% params0.FV_velocities = [0 -0.5000   -1 -2.0000   -3.0000   -4.0000 -5 -6];
% params0.PieceWiseStrainDepParams__2 = 1.1554;
params0.RunForceVelocity = true;
params0.RunSlack = true;
params0.PlotEachSeparately = true;
params0.mods = [];
params0.g = [];
% params0.k2 = 131.3860*1;
params0.UseOverlapFactor = false;
params0.kSE = 1e3*3.5;
% params0
% params0.fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y'};
params0.UseA2MechanicalRecocking = false;
params0.k2d = 2000;
params0.drmr = params0.dr;
params0.mu_neg = params0.mu;


% params0.g(2) = 0;
params0.RunSlackSegments = 'AllPar';
tic
RunBakersExp
%

% features_model = extractSlackAttributes(out.t, out.Force, out.SL, velocitytable, features_model, out, false);

toc
[costs, weights, cost] = evalFeatureCost(features_data, features_model, params0.fn(1:end), 1);
sum(costs)

plotFeatures(features_data, features_model, [], params0.fn(1:end));

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
%%
params0 = getParams(params0, x, false, true);
%% construct feature metric - test this!
% params0_bak = params0;
% params0 = getParams(params0, params0.g, false, true);
% construct parameter set to vary
paramNames = {'kstiff2','k_pas','ka','kah','dr2','k2','k1'};
% paramNames = {'baseline_dummy', 'kstiff2','k_pas'};
paramNames = {'baseline_dummy', 'PieceWiseStrainDepX__2', 'PieceWiseStrainDepX__3','PieceWiseStrainDepX__4','PieceWiseStrainDepX__5','kstiff2','k_pas','ka','kah','kamh','dr2','k2','k1','PieceWiseStrainDepParams__2','PieceWiseStrainDepParams__3','PieceWiseStrainDepParams__4','PieceWiseStrainDepParams__5','mu','kstiff1','kd','sigma1','sigma2','sigma_srd1','sigma_srd2','estiff','gamma','kSE','ekSE'};
% construct feature set to evaluate
% fn = {'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'SLslack|t0|0', 'Am|SLslack|', 'peak1_y|SLslack|0', 'peak1_y|v_restretch|0','peak1_dSL', 'peak2|v_restretch', 'vall_t|v_restretch|0.1', 'vall_y', 'vall2_dy|v_restretch|0', 'ovrsht_dy|_|0', 'steady'};
% Full set
paramNames =  {'baseline_dummy', ... 
'dr2', ... 
'ekSE', ... 
'estiff', ... 
'gamma', ... 
'k_pas', ... 
'k1', ... 
'k2', ... 
'ka', ... 
'kah', ... 
'kamh', ... 
'kd', ... 
'ksr0', ... 
'ksrd', ... 
'ksrd2sr', ... 
'ksr2srd', ... 
'kmsrd', ... 
'kmsr', ... 
'kSE', ... 
'kstiff1', ... 
'kstiff2', ... 
'MaxSlackNegativeForce', ... 
'mu', ... 
'PieceWiseStrainDepParams__2', ... 
'PieceWiseStrainDepParams__3', ... 
'PieceWiseStrainDepParams__4', ... 
'PieceWiseStrainDepParams__5', ... 
'PieceWiseStrainDepParams__6', ... 
'PieceWiseStrainDepX__2', ... 
'PieceWiseStrainDepX__3', ... 
'PieceWiseStrainDepX__4', ... 
'PieceWiseStrainDepX__5', ... 
'PieceWiseStrainDepX__6', ... 
'sigma_srd1', ... 
'sigma_srd2', ... 
'sigma1', ... 
'sigma2'};

paramNames = {'baseline_dummy', ... 
 'PieceWiseStrainDepX__2', ... 
'PieceWiseStrainDepParams__2', ... 
'PieceWiseStrainDepParams__3', ... 
'PieceWiseStrainDep2X__2', ... 
'PieceWiseStrainDep2X__4', ... 
'PieceWiseStrainDep2Params__2', ... 
'PieceWiseStrainDep2Params__3', ... 
'PieceWiseStrainDepR1DX__2', ... 
'PieceWiseStrainDepR1DX__3', ... 
'PieceWiseStrainDepR1DParams__2', ... 
'PieceWiseStrainDepR1DParams__3', ... 
'PieceWiseStrainDepR21X__2', ... 
'PieceWiseStrainDepR21X__3', ... 
'PieceWiseStrainDepR21Params__2', ... 
'PieceWiseStrainDepR21Params__3'};

paramNames = {'baseline_dummy', ... 
    'MaxSlackNegativeForce', ... 
'PieceWiseStrainDep2Params__2', ... 
'PieceWiseStrainDep2Params__3', ... 
'PieceWiseStrainDep2Params__4', ... 
'PieceWiseStrainDep2X__2', ... 
'PieceWiseStrainDep2X__3', ... 
'PieceWiseStrainDep2X__4', ... 
'PieceWiseStrainDepParams__2', ... 
'PieceWiseStrainDepParams__3', ... 
'PieceWiseStrainDepR1DParams__2', ... 
'PieceWiseStrainDepR1DParams__3', ... 
'PieceWiseStrainDepR1DX__2', ... 
'PieceWiseStrainDepR1DX__3', ... 
'PieceWiseStrainDepR21Params__2', ... 
'PieceWiseStrainDepR21Params__3', ...
'PieceWiseStrainDepR21X__2', ... 
'PieceWiseStrainDepR21X__3', ... 
'PieceWiseStrainDepX__2', ... 
'PieceWiseStrainDepX__3', ... 
'dr2', ... 
'ekSE', ... 
'estiff', ... 
'gamma', ... 
'k1', ... 
'k2', ... 
'kSE', ... 
'k_pas', ... 
'ka', ... 
'kah', ... 
'kamh', ... 
'kd', ... 
'kmsr', ... 
'kmsrd', ... 
'ksr0', ... 
'ksr2srd', ... 
'ksrd', ... 
'ksrd2sr', ... 
'kstiff1', ... 
'kstiff2', ... 
'mu', ... 
'sigma1', ... 
'sigma2', ... 
'sigma_srd1', ... 
'sigma_srd2', ...
'k_1',...
'L_thick', ...
'L_hbare',...
'L_thin' };

tp = tunableParams(params0);
paramNames = fieldnames(tp);
paramNames  = ['baseline_dummy'; paramNames];
% paramNames = paramNames(selectedParamsForResim)
% selectedParamsForResim = [1, 3, 4, 6, 7, 45];
% simulate force-isovelocity
% take only non-zero velocities with > 10 kPa
% params0.FV_velocities = -[0.5, 1, 3, 4];
% params0.RunForceVelocity = false;
% params0.RunSlack = false;
% params0.RunForceLengthEstim = false;
% params0.EvalFeatures = false;
params0.PlotEachSeparately = true;

% params0.RunForceVelocity = true;
% params0.RunSlack = true;
params0.MaxRunTime = 60;
params0.baseline_dummy = 0;
% params0 = getParams(params0, params0.g, false, true);
params0.FV_velocities = -[0.5, 2, 3, 4];
params0.PlotEachSeparately = false;

% params0.EvalFeatures = true;
% params0.RunForceLengthEstim = true;
% fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y|SLslack', 'peak1_dSL', 'peak2', 'ovrsht_dy|_|0', 'steady', 'XTOR'};
% fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y|SLslack', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1'};
% fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y'};

fn = params0.fn;

% fn = {'FV_f|FV_v'};
savedFeats = {};

RunDeltaPlus = true;
RunDeltaMinus = true;

if RunDeltaMinus    
    featureMatrixMinus = zeros(length(paramNames), length(fn));
end
if RunDeltaPlus
    featureMatrixPlus = zeros(length(paramNames), length(fn));
end

delta = 0.01;
% params_base = params0;
tic
% params0.
for i_param = 1:length(paramNames)
% for i_param = selectedParamsForResim
    try
        fprintf('Running %s..', paramNames{i_param});
        params0.mods = paramNames(i_param);
        
        if RunDeltaPlus
            params0.g = 1 + delta;
            fprintf('+...');
            figure(333); clf;        
            RunBakersExp;
            % % check the init
            [costs, weights, cost] = evalFeatureCost(features_data, features_model, fn, 1);
            % costs = E;
            featureMatrixPlus(i_param, :) = costs;
            savedFeats{i_param, 1} = features_model;
        end
        if RunDeltaMinus
            if all(featureMatrixPlus(i_param, :) == featureMatrixPlus(1, :)) 
                % featureMatrixMinus(i_param, :) = costs;
                fprintf('x...');
            else
                params0.g = 1 - delta;
                fprintf('-...');
                figure(334); clf;        
                RunBakersExp;
                % % check the init
                [costs, weights, cost] = evalFeatureCost(features_data, features_model, fn, 1);
                % costs = E;
                featureMatrixMinus(i_param, :) = costs;
                savedFeats{i_param, 2} = features_model;
            end
        end
                
        fprintf('Done \n');
    catch e
        fprintf('failed. [%s] \n', e.message);
    end
end
% save('SA_pwsd', "savedFeats", "featureMatrixMinus", "featureMatrixPlus");
save env_SA
toc
%%


%% RUN THE OPTIM SCRIPT FOR LSQNON USING CUSTOM JACOBIAN

% main_lsqnonlin_script.m
ModelParamsFVOptimBaseline
% 1. Setup
% paramNames = {'kstiff2','k_pas','ka','kah','dr2','k2','k1'};
% params0.g = ones(1, length(paramNames));
% params0.mods = paramNames;
params0.FV_velocities = -[0.5, 2, 3, 4];
params0.RunSlack = true;
params0.RunForceVelocity = true;
% matchStructFields(params0, 'Run*', true)
params0.MaxStrainArraySize = 40;
params0.MaxRunTime = 30;
params0.EvalFeatures = true;
params0.BreakOnODEUnstable = true;
params0.PlotEachSeparately = true;

tic
RunBakersExp;
toc
%%
fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y|SLslack', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1'};
[Residuals, weights, cost] = evalFeatureCost(features_data, features_model, fn, 1);
plotFeatures(features_data, features_model, [], fn)

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
%% Visualize the double-sided 


% figure(4);
% same values, so did not run - indexes
sameNotRun = (featureMatrixPlus  > 0 & featureMatrixMinus == 0);
featureMatrixMinus(sameNotRun) = featureMatrixPlus(sameNotRun);

featureMatrixPlusDiff = (featureMatrixPlus - featureMatrixPlus(1, :));
featureMatrixMinusDiff = (featureMatrixMinus - featureMatrixMinus(1, :));
if ~RunDeltaMinus
    featureMatrixMinusDiff   = featureMatrixPlusDiff;
end
if ~RunDeltaPlus
    featureMatrixPlusDiff = featureMatrixMinusDiff;
end
featureMatrixPlusDiff(featureMatrixPlus == 0) = nan;
featureMatrixMinusDiff(featureMatrixMinus == 0) = nan;


featureMatrixPlusN = featureMatrixPlusDiff./featureMatrixPlus(1, :)*100;
featureMatrixMinusN = featureMatrixMinusDiff./featureMatrixMinus(1, :)*100;
%%

featureMatrixN = min(featureMatrixMinusDiff, featureMatrixPlusDiff);

paramCost = sum(featureMatrixN, 2);
[vals, paramPos] = sort(abs(paramCost), 1, "desc");

% cutoff too bads
tooBadCost = 1;
tooBad = featureMatrixN > tooBadCost;
featureMatrixN(tooBad) = tooBadCost;

% cutoff too good to be true
tgtbtCost = -1;
tgtbt = featureMatrixN < tgtbtCost;
featureMatrixN(tgtbt) = tgtbtCost;

% paramPos = 1:length(paramNames);
featureMatrixNSorted = featureMatrixN(paramPos, :);%./featureMatrixPlusN(paramPos(1), :)*100;

validRows = ~isnan(sum(featureMatrixNSorted, 2)) & sum((featureMatrixNSorted), 2) ~= 0;
paramPos = paramPos(validRows);
featureMatrixNSorted = featureMatrixNSorted(validRows, :);

imagesc(featureMatrixNSorted);
colorbar;
axis ij tight;

set(gca, 'XTick', 1:length(fn), ...
         'XTickLabel', fn, ...
         'YTick', 1:length(paramNames), ...
         'YTickLabel', paramNames(paramPos), 'TickLabelInterpreter', 'none');

xtickangle(45);
title('Parameter Sensitivities');

% Print numeric values in each cell
for i = 1:size(featureMatrixNSorted,1)
    for j = 1:size(featureMatrixNSorted,2)
        text(j, i, sprintf('%.3g', featureMatrixNSorted(i,j)), ...
            'HorizontalAlignment', 'center', ...
            'Color', 'k');
    end
end

paramNames(paramPos(1:20))

%% gemini udpated version

posInfluential = 0.1; % how much do we consider worsening the fit an influential parameter

% 1. Calculate the total cost for each direction (Row-based)
sumMinus = sum(featureMatrixMinusDiff, 2);
sumPlus  = sum(featureMatrixPlusDiff, 2);

% 2. Determine the "Winner" for each row
% This finds which direction has the lower total sum for the whole row
[paramCost, bestDirIdx] = min([sumMinus, sumPlus], [], 2);

% 3. Build the final matrix based on the winning direction for that row
featureMatrixN = zeros(size(featureMatrixMinusN));
for r = 1:size(featureMatrixMinusN, 1)
    if bestDirIdx(r) == 1
        featureMatrixN(r, :) = featureMatrixMinusN(r, :);
    else
        featureMatrixN(r, :) = featureMatrixPlusN(r, :);
    end
end

% 3. Filter out NaNs or all-zero rows before sorting
validMask = ~any(isnan(featureMatrixN), 2) & any(featureMatrixN ~= 0, 2);

% 4. Sort based on the COST (Ascending usually makes sense for "lowest cost")
% If you want highest sensitivity at the top, use 'descend' on abs(paramCost)
[~, sortedIdx] = sort(paramCost.*(paramCost > 0)*posInfluential + abs(paramCost).*(paramCost < 0), 'descend');

% Intersect sorting with validity
paramPos = sortedIdx(validMask(sortedIdx));

% 5. Apply sorting to the matrix and create labels
featureMatrixNSorted = featureMatrixN(paramPos, :);
newYLabels = cell(length(paramPos), 1);

for k = 1:length(paramPos)
    rowIdx = paramPos(k);
    
    if bestDirIdx(rowIdx) == 1
        dirStr = '[-]';
    else
        dirStr = '[+]';
    end
    
    % Construct label: Name [+/-] (TotalCost)
    newYLabels{k} = sprintf('%s %s (%.3f)', ...
        paramNames{rowIdx}, dirStr, paramCost(rowIdx));
end

% --- Visualization ---
featureMatrixNSortedCO = featureMatrixNSorted;
cutOffMax = 10;
cutOffMin = -10;
featureMatrixNSortedCO(featureMatrixNSorted > cutOffMax) = cutOffMax;
featureMatrixNSortedCO(featureMatrixNSorted < cutOffMin) = cutOffMin;

imagesc(featureMatrixNSortedCO);
colorbar;
colormap(turbo);
axis ij tight;

set(gca, 'XTick', 1:length(fn), ...
         'XTickLabel', fn, ...
         'YTick', 1:length(paramPos), ...
         'YTickLabel', newYLabels, ...
         'TickLabelInterpreter', 'none');

xtickangle(45);
title('Sorted Parameter Sensitivities');

% Print numeric values in each cell
for i = 1:size(featureMatrixNSorted,1)
    for j = 1:size(featureMatrixNSorted,2)
        text(j, i, sprintf('%.3g', featureMatrixNSorted(i,j)), ...
            'HorizontalAlignment', 'center', ...
            'Color', 'k');
    end
end
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
