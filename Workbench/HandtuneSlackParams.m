%% Compare simulation versions

%% Extract data features 
datastruct = load('../data/bakers_slack8mM_all_fix.mat');
velocitytable =    datastruct.velocitytable;
datatable = datastruct.datatable;

figure(221); clf;
nexttile;hold on; 
feats_data = extractSlackAttributes(datatable(:, 1), datatable(:, 3), datatable(:, 2), velocitytable);
title('data');

%% baseline
params0 = getParams();
ModelParamsInitNiceSlack;
params0.modelFcn = 'dPUdTCaSimpleAlternative2State';
params0.modelFcn = 'dPUdT_CombinedTransitions';
ksr0 = params0.kmsr;
params0.kmsr = params0.ksr0;
params0.ksr0 = ksr0;
% params0.dS = 0.002;
% params0.dS = 0.004;

% params0.drp3 = 0.01;
% params0.dr = 0.007;
params0.drp3 = params0.dr;

params0.k3m = 0;
params0.k3 = 1500; 

params0.RunStairs = false;    
params0.RunForceVelocity = false;
    params0.RunSlack = true;
    params0.EvalFitSlackOnset = true;
    params0.RunForceLengthEstim = true;
    params0.RunSlackSegments = 'All';
    params0.drawForceOnset = true;


    if params0.alpha0 ~= 0 && params0.alpha0_R == 0 && params0.alpha0_R == 0
        params0.alpha0_R = params0.alpha0 ;
        params0.alpha0_R = params0.alpha0 ;
        fprintf('Updating alpha... \n');
    end

    params0.NumberOfStates = 2;



% params0.ghostSave = 'testdS';
params0.ghostLoad = 'testdS';
params0.dS = 0.0027;
params0.Slim_r = 2.3;
params0.dr = 0.01;
% params0.A1AttachmentWidth = 3*params0.dS;
% params0.kstiff1
% params0.kstiff2 = params0.kstiff2*2.2;
% params0.dr = 0.01;
% params0.xrate = 2;
params0.k2_L = params0.k2_L*5;
params0.alpha2_R = 1.2;
params0.kstiff1 = 12000;
params0.kstiff2  = params0.kstiff1;
figure(10);
LoadData;
params0.ShowStatePlots = true;
params0.EvalFitSlackOnset = true;
params0.UseNegativeForceRip = false;
params0.MaxSlackNegativeForce = -6;

params0.RunForceLengthEstim = false;
params0.justPlotStateTransitionsFlag = false;

params0.ksr0 = params0.ksr0 /10;
params0.kmsr = params0.kmsr/10;
params0.kd = 20;
params0.ka = 50;
params0.UseOverlap = false;
params0.UseSuperRelaxed = false;

params0.ksr0 = 4.3;
params0.kmsr = 6.4e-3;
params0.UsePassiveForSR = true;
params0.sigma1 = 750;

params0.A1AttachmentWidth = 1*params0.dS;
params0.kstiff1 = 12000;
params0.kstiff2  = params0.kstiff1;

clf;
tic
params0.SL0 = 2.2;
params0.kSE = 500;
params0.ekSE = 1.5;

params0.kSE = 50000;
params0.ekSE = 1;

params0.MaxSlackNegativeForce = -5;

params0.ResetSRat = [velocitytable(velocitytable(:, 2) < 0, 1), linspace(20, 80, 5)'];
params0.velocitytableonfile = 'bakers_slack8mM_all_extendedSlack.mat';
% params0.ksr0 = 0;
RunBakersExp;
% plot(out.t, out.RD1);
toc


% writeParamsToMFile('ModelParamsInitNiceSlack_prescribedSR_var1.m', params0, [], "Pretty nice fit, but the SR is fudged");
%
% StatesInTime
%% prescribed SR
clf;
clear params0;
ModelParamsInitNiceSlack_prescribedSR
params0.velocitytableonfile = 'bakers_slack8mM_all_fix.mat';
params0.velocitytableonfile = 'bakers_slack8mM_all_extendedSlack.mat';

sr_dist = (linspace(0, 1, 5)').^0.5;
sr_mi = 0.2; sr_ma = 0.7;
sr_dist = sr_dist*(sr_ma - sr_mi) + sr_mi;
params0.ResetSRat = [velocitytable(velocitytable(:, 2) < 0, 1), sr_dist];
params0.kmsr = .01;
params0.ksr0 = 5;
% params0.kmsr = 0;
% params0.ksr0 = 0;
params0.kstiff1 = 13000;
params0.kstiff2 = params0.kstiff1; 
params0.UseForceOnsetShift = false;
params0.UseSuperRelaxed = true;
params0.UseOverlap = true;
params0.UseNegativeForceRip = true;
% params0.k1 = 1500;
params0.justPlotStateTransitionsFlag = false;
RunBakersExp;

%%
matchStructFields(params0, 'Use*', [0], true)   % print values
%%
clear params0
ModelParamsInitNiceSlack_prescribedSR_var1;

params0.UseOverlap = false;
params0.UseOverlapFactor = true;

RunBakersExp;
%% Extract simulated features
% datastruct = load('data/bakers_slack8mM_all_extendedSlack.mat');
% datastruct = load('data/bakers_slack8mM_all_fix.mat');
datastruct = load(['data/' params0.velocitytableonfile]);
velocitytable =    datastruct.velocitytable;

figure(221); clf;hold on;
feats_sim = extractSlackAttributes(out.t, out.Force, out.SL, velocitytable, out);
title('Sim');   

% plot the features
fn = {'ktr|SLdiff', 'A|SLdiff', 't0|SLslack', 'SLslack|t0|0', 'Am|SLslack|', 'peak1_y|SLslack|0', 'peak1_y|v_restretch|0','peak1_dSL', 'peak2|v_restretch', 'vall_t|v_restretch|0.1', 'vall_y', 'vall2_dy|v_restretch|0', 'ovrsht_dy|_|0', 'steady'};
% figure
plotFeatures(feats_data, feats_sim, feats_ghost, fn);
%%
F = 0:1:100;
params0.kSE = 500;
LSE = F.^(1/params0.ekSE)./params0.kSE;
plot(F, LSE, F, F/1e4)

%% call it a ghost
feats_ghost = feats_sim;

%% Manipulate with 
figure(220);clf; hold on;
params0.UseForceOnsetShift = false;
params0.RunStairs = false;
params0.dr = 0.02;

RunBakersExp;
figure(221); nexttile;hold on;
feats_sim = extractSlackAttributes(out.t, out.Force, out.SL, velocitytable);
title('Sim');   

%% plot the features
fn = {'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'SLslack|t0|0', 'Am|SLslack|', 'peak1_y|SLslack|0', 'peak1_y|v_restretch|0','peak1_dSL', 'peak2|v_restretch', 'vall_t|v_restretch|0.1', 'vall_y', 'vall2_dy|v_restretch|0', 'ovrsht_dy|_|0', 'steady'};
plotFeatures(feats_data, feats_sim, feats_ghost, fn);
%%
figure(222)
% plotRates(out)
runStatesInTime(out, params);

%% test  fit

sld = -[feats_data.SLdiff];
t = [feats_data.t0];

[t',sld']
% USING MYCURVEFIT, predicts -0.125 at t = 0
t_ = linspace(0, t(end)*1.5, 100);
y = -0.4667455 + 0.3412689*exp(-54.36518*t_);
plot(t, sld, 'x-', t_, y, '--');

%% test load 
datastruct = load('../data/bakers_slack2mM.mat');
velocitytable =    datastruct.velocitytable;
datatable = datastruct.datatable;

figure(221); clf;
nexttile;hold on; 
feats_data2 = extractSlackAttributes(datatable(:, 1), datatable(:, 3), datatable(:, 2), velocitytable);
title('data');

sld2 = -[feats_data2.SLdiff];
t2 = [feats_data2.t0];

[t',sld']

%%
plotFeatures(feats_data, feats_data2, feats_ghost, fn);

%%
clf;
% USING MYCURVEFIT, predicts -0.169 at t = 0 for 2 ATP
t_ = linspace(0, t(end)*1.5, 100);
y8 = -0.4667455 + 0.3412689*exp(-54.36518*t_);
y2 = -0.4384684 + 0.2694423*exp(-66.72785*t_);

k8 = [feats_data.steady]./y8(1);
k2 = [feats_data2.steady]./y2(1);

plot(t, sld, 'x-', t_, y8, '--', t2, sld2, 'x-', t_, y2, '--', LineWidth=2);
ylim([-inf, 0]);

