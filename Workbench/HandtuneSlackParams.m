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

% params0.ResetSRat = [velocitytable(velocitytable(:, 2) < 0, 1), linspace(20, 80, 5)'];
params0.velocitytableonfile = 'bakers_slack8mM_all_extendedSlack.mat';

% params0.velocitytableonfile = 'bakers_isovelocity.mat';
% params0.ksr0 = 0;
RunBakersExp;
% plot(out.t, out.RD1);
% writeParamsToMFile('ModelParamsInitNiceSlack_prescribedSrxD.m', params0, [], "Pretty nice fit, but the SrxD is fudged");
%
% StatesInTime
%% prescribed SR
clf;
clear params0;
ModelParamsInitNiceSlack_prescribedSR
params0.velocitytableonfile = 'bakers_slack8mM_all_fix.mat';
% params0.velocitytableonfile = 'bakers_slack8mM_all_extendedSlack.mat';
params0.velocitytableonfile = 'bakers_isovelocity.mat';

params0.RunSlackSegments = 'All';

sr_dist = (linspace(0, 1, 5)').^.6;
sr_mi = 0.2; sr_ma = 0.6;
sr_dist = sr_dist*(sr_ma - sr_mi) + sr_mi;
params0.ResetSRat = [velocitytable(velocitytable(:, 2) < 0, 1), sr_dist];

params0.UseSafeSRReset = true;

params0.ResetSRat = [];
params0.kmsr = .01;
params0.ksr0 = 5;
% params0.kmsr = 0;
% params0.ksr0 = 0;
params0.kstiff1 = 13000;
params0.kstiff2 = params0.kstiff1; 
params0.UseForceOnsetShift = false;
params0.UseSuperRelaxed = true;
params0.UseOverlap = true;
params0.UseOverlapFactor = true;
params0.UseNegativeForceRip = true;
params0.mu = 1e-1;
% params0.k1 = 1500;
params0.justPlotStateTransitionsFlag = false;


params0.UseSuperRelaxed = false
params0.UseSuperRelaxedADP = true;

params0.kmsrd = 40;
params0.sigma_srd1 = 1e6;
params0.ksrd = .07;
params0.sigma_srd2 = 1e6;


RunBakersExp;
%% prestretch
clf;
params0 = getParams();
ModelParamsInitNiceSlack
% params0.SL0 = 2.0;
% params0.RunForceVelocity = false;
% params0.RunKtr = false;
params0.kSE = 600;
params0.rsl0 = findLSEPreStretch(2.0, 57, params0);

RunBakersExp;
StatesInTime
%% Test Srx_ADP

clf;
params0 = getParams();
ModelParamsInitNiceSlack_dr01

params0.UseSuperRelaxed = false;
params0.UseSuperRelaxedADP = true;

params0.ksrd =  0.0038016;
params0.sigma_srd1 = 1e6;

params0.kmsrd =  0.0073827;
params0.sigma_srd2 = 1e6;

% params0.SL0 = 2.0;
% params0.RunForceVelocity = false;
% params0.RunKtr = false;
params0.kSE = 5000;
params0.ShowStatePlots  = true;
% params0.rsl0 = findLSEPreStretch(2.0, 57, params0);
params0.velocitytableonfile = 'bakers_slack8mM_all_fix.mat';

RunBakersExp;


%%
matchStructFields(params0, 'Use*', [0], true)   % print values
%%
clf;
clear params0
% ModelParamsInitNiceSlack_prescribedSR_var1;
ModelParamsInitNiceSlack_prescribedSR_var2;

params0.RunSlackSegments = 'First';
% params0.RunSlackSegments = 'All';
params0.UseSafeSRReset = true;
params0.kmsr = 0.02;
params0.k1 = 85;
params0.kSE = 600;

% params0.UseNegativeForceRip = false;
% params0.UseOverlap = false;
% params0.UseOverlapFactor = true;
params0.justPlotStateTransitionsFlag = false;
% params0.UsePassiveForSR = false;
% params0.MaxSlackNegativeForce = -2;

RunBakersExp;
%% 

%% read bakers titin data to read the slack

base = 'data/PassiveCaSrc2/20241220/';
base = 'data/PassiveCaSrc2/20241219/';
base = 'data/PassiveCaSrc2/20241217/';
base = 'data/PassiveCaSrc2/20241212/';
base = 'data/PassiveCaSrc2/20241121/';
fn = [base  '01_05_0.1s_RelaxS2F.txt'];
fn = [base  '03_05_0.1s_Relax_PNB_Mava.txt'];
fn = [base  '03_10_0.1s_pCa_4.4_PNB_Mava.txt'];


fn = [base  '01_05_0.1s_RelaxS2F.txt'];
fn = [base  '03_02_100s_Relax_PNB_Mava.txt'];
fn = [base  '01_09_100s_RelaxF2S.txt'];
fn = [base  '03_10_0.1s_pCa_4.4_PNB_Mava.txt'];

fn = [base  '01_01_Log_Relax.txt'];
fn = [base  '02_01_Log_Fmax.txt'];
fn = [base  '03_01_Log_PNB_Mava.txt'];

tb = readtable(fn, "NumHeaderLines",3);
tb.Properties.VariableNames = {'t', 'ML', 'F', 'SL'};

figure;clf;hold on;
clf;hold on;
nexttile;
plot(tb.t, tb.ML*2, tb.t, tb.SL)
nexttile;
plot(tb.t, tb.F)
% plot(tb.t, tb.SL)

% figure
% plot(tb.SL, tb.F)

%% for time ~ 1814 and 20241220/03_01_Log_PNB_Mava.txt, the force to 2*ML -
% SL:
FL = [1.4, +0.0560 - 1.9;6.424,+ 0.0679-2.35;27.663,+0.15-2.35];
clf; plot(FL(:, 2), FL(:, 1), '*-'); hold on;

fun = @(a, SL0, c, x) a*(x-SL0).^c;

% fitoptions('init')
f = fit(FL(:, 2), FL(:, 1), fun, 'StartPoint', [0.1, 1.5, 1]);
sl = 1.9:0.05:2.5;
plot(sl,f(sl))
%%

% Parameters
S = 30000;            % max samples to show ("snake" length)
R = 30000;             % new samples per second
T = 10;             % total duration in seconds

dt = 1/R;           % time step
N = T * R;          % total number of samples

% Example signal (circle + noise)
t = linspace(0, T, N);
xData = tb.SL;
yData = tb.F;

% Initialize plot
figure;clf;
hPlot = plot(nan, nan, '-o');
% xlim([-1.5 1.5]); ylim([-1.5 1.5]);
axis equal;
grid on;
title('Snake XY Trace');

% Animate
for k = 1:100:N
    idxStart = max(1, k-S+1);   % sliding window start
    set(hPlot, 'XData', xData(idxStart:k), 'YData', yData(idxStart:k));
    xlim([1.5, 2.5])
    drawnow;
    pause(dt);  % control update rate
end



%% further testing srx d
ModelParamsInitNiceSlack_prescribedSrxD;
clf;
params0.justPlotStateTransitionsFlag = false;
params0.RunSlackSegments = 'All';
params0.UsePassiveForSR = false;

RunBakersExp;
%% Extract simulated features
% datastruct = load('data/bakers_slack8mM_all_extendedSlack.mat');
% datastruct = load('data/bakers_slack8mM_all_fix.mat');
datastruct = load(['data/' params0.velocitytableonfile]);
velocitytable =    datastruct.velocitytable;

figure(221); clf;hold on;
feats_sim = extractSlackAttributes(out.t, out.Force, out.SL, velocitytable, out);
% feats_sim.params = params0;
title('Sim');   

% plot the features
fn = {'ktr|SLdiff', 'A|SLdiff', 't0|SLslack', 'SLslack|t0|0', 'Am|SLslack|', 'peak1_y|SLslack|0', 'peak1_y|v_restretch|0','peak1_dSL', 'peak2|v_restretch', 'vall_t|v_restretch|0.1', 'vall_y', 'vall2_dy|v_restretch|0', 'ovrsht_dy|_|0', 'steady'};
figure(80085)
plotFeatures(feats_data, feats_sim, feats_ghost, fn);
%%

matchStructFields(params0, 'Use*', [], true);
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


[t',sld']

%%
plotFeatures(feats_data, feats_data2, feats_ghost, fn);

%%
clf;

sld = -[feats_data.SLdiff];
t = [feats_data.t0];

sld2 = -[feats_data2.SLdiff];
t2 = [feats_data2.t0];


% USING MYCURVEFIT, predicts -0.169 at t = 0 for 2 ATP
t_ = linspace(0, t(end)*1.5, 100);
y8 = -0.4667455 + 0.3412689*exp(-54.36518*t_);
y2 = -0.4384684 + 0.2694423*exp(-66.72785*t_);

k8 = [feats_data.steady]./y8(1);
k2 = [feats_data2.steady]./y2(1);

plot(t, sld, 'x-', t_, y8, '--', t2, sld2, 'x-', t_, y2, '--', LineWidth=2);
ylim([-inf, 0]);
xlim([0 0.1])

hold on;

% clf
% t_seg = find(velocitytable(:, 2) < 0, 1);
isegs = out.t>velocitytable(19, 1);
t = out.t(isegs);
lxb = out.LXB(isegs);

plot(t - t(1) - 0.015, lxb - 2.2)

%% try to get the best params

%% baseline
figure(13);clf;
params0 = getParams();
ModelParamsInitOptim_slackAll
% params0.modelFcn = 'dPUdTCaSimpleAlternative2State';
% params0.modelFcn = 'dPUdT_CombinedTransitions_old';
params0.velocitytableonfile = 'bakers_isovelocity.mat';

 % adjsut ksr0 and kmsr bullshiut
    if params0.kmsr > params0.ksr0
        ksr0 = params0.kmsr;
        params0.kmsr = params0.ksr0;
        params0.ksr0 = ksr0;        
    end

    if params0.alpha0 ~= 0 && params0.alpha0_R == 0 && params0.alpha0_R == 0
        params0.alpha0_R = params0.alpha0 ;
        params0.alpha0_R = params0.alpha0 ;
        fprintf('Updating alpha... \n');
    end
params0.Slim_l = 1.85 + 0.01;
% params0.Slim_l = 1.80 + 0.0027*20 + 0.0;
params0.Slim_l = 1.80;
params0.Slim_r = 2.31;
% params0.kSE = 1000;
params0.dS = 0.0027/2;
params0.dS = 0.00;
% params0.Slim_r = 2.3;
% params0.dr = 0.01;
params0.RunSlackSegments = 'All';
params0.A1AttachmentWidth = 0.012;
params0.dS = 0.003;
params0.kstiff1 = 1.2e4;
params0.kstiff2 = params0.kstiff1;
params0.kSE  = 1000;
params0.ka = 30;
params0.k1 = 200;
% params0.dr
params0.ShowStatePlots = true;
tic
RunBakersExp;
toc
params.dS
params.ss
%%
figure(11);
StatesInTime
%%
% ksr0 = params0.kmsr;
% params0.kmsr = params0.ksr0;
% params0.ksr0 = ksr0;
% % params0.dS = 0.002;
% % params0.dS = 0.004;
% 
% % params0.drp3 = 0.01;
% % params0.dr = 0.007;
% params0.drp3 = params0.dr;
% 
% params0.k3m = 0;
% params0.k3 = 1500; 

params0.RunStairs = false;    
params0.RunForceVelocity = false;

params0.RunSlack = true;
params0.EvalFitSlackOnset = false;
params0.RunForceLengthEstim = false;
params0.RunSlackSegments = 'All';
params0.drawForceOnset = true;


if params0.alpha0 ~= 0 && params0.alpha0_R == 0 && params0.alpha0_R == 0
    params0.alpha0_R = params0.alpha0 ;
    params0.alpha0_R = params0.alpha0 ;
    fprintf('Updating alpha... \n');
end

% params0.NumberOfStates = 2;
params0.UseTitinInterpolation = 0;
params0.RunSlackSegments = 'Last';

% params0.ghostSave = 'testdS';
% params0.ghostLoad = 'testdS';
params0.dS = 0.0027;
params0.Slim_r = 2.3;
params0.dr = 0.01;
% params0.A1AttachmentWidth = 3*params0.dS;
% params0.kstiff1
% params0.kstiff2 = params0.kstiff2*2.2;
% params0.dr = 0.01;
% params0.xrate = 2;
params0.k2_L = params0.k2_L*5;
params0.k2_R = params0.k2_R;
params0.alpha2_R = 1;
params0.RunSlackSegments = 'First';
params0.kstiff1 = 15000;
params0.kstiff2  = params0.kstiff1;

figure(10);
LoadData;
params0.ShowStatePlots = true;
params0.EvalFitSlackOnset = true;
params0.UseNegativeForceRip = false;
params0.MaxSlackNegativeForce = -6;

params0.RunForceLengthEstim = false;
params0.justPlotStateTransitionsFlag = false;


% params0.ksr0 = params0.ksr0 /10;
% params0.kmsr = params0.kmsr/10;
% params0.kd = 20;
% params0.ka = 50;
params0.UseOverlap = true;
params0.UseSuperRelaxed = true;

% params0.ksr0 = 4.3;
% params0.kmsr = 6.4e-3;
% params0.UsePassiveForSR = true;
% params0.sigma1 = 750;
% 
% params0.A1AttachmentWidth = 1*params0.dS;
% params0.kstiff1 = 12000;
% params0.kstiff2  = params0.kstiff1;
% 
% clf;
% tic
% params0.SL0 = 2.2;
% params0.kSE = 500;
% params0.ekSE = 1.5;
% 
% params0.kSE = 50000;
% params0.ekSE = 1;

params0.MaxSlackNegativeForce = -5;

% params0.ResetSRat = [velocitytable(velocitytable(:, 2) < 0, 1), linspace(20, 80, 5)'];
% params0.velocitytableonfile = 'bakers_slack8mM_all_extendedSlack.mat';
params0.velocitytableonfile = 'bakers_slack8mM_all_fix.mat';

% params0.velocitytableonfile = 'bakers_isovelocity.mat';
% params0.ksr0 = 0;
RunBakersExp;
% plot(out.t, out.RD1);
% writeParamsToMFile('ModelParamsInitNiceSlack_updated.m', params0, [], "Pretty nice fit, but the SrxD is fudged");
%
% StatesInTime
%% Soo...
% matchStructFields(params0, 'Use*', [], true);
clf;
% Simplest model
% add complexities

params0 = getParams();
ModelParamsInitOptim_slackAll

if params0.alpha0 ~= 0 && params0.alpha0_R == 0 && params0.alpha0_R == 0
    params0.alpha0_R = params0.alpha0 ;
    params0.alpha0_R = params0.alpha0 ;
    fprintf('Updating alpha... \n');
end

params0.UseSerialStiffness = false;


params0.RunSlackSegments = 'Last';
params0.A1AttachmentWidth = 0.002;
params0.dS = 0.001;
params0.Slim_r = 2.2;
params0.Slim_l = 2.1;
params0.WindowsOverflowStepCount = -1;
params0.BreakOnODEUnstable = false;

params0.UseSuperRelaxed = false;
params0.UseTitinInterpolation = false;
params0.ShowStatePlots = true;


params0.k2 = 50;

params0.k_1 = 0;
params0.UseStrictDetachmentAt = 0.02;
params0.justPlotStateTransitionsFlag = false;

params0.Velocity = -1;
t_sl0 = [0 0.1];
params0 = getParams(params0, [], true);

tic
[F out] = evaluateModel(modelFcn, t_sl0, params0);
% RunBakersExp
toc

%% 

ModelParams_FVfit
params0.Velocity = -1;
params = params0;
params0.Slim_r = 2.2;
params0.Slim_l = 2.1;
params0.WindowsOverflowStepCount = -1;

params0 = getParams(params0, [], true);
[F out] = evaluateModel(modelFcn, t_sl0, params0);

StatesInTime
