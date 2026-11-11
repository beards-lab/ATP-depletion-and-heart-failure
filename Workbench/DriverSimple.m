% simplest 
% clear;
figure(2);
clf; 
% initialize parameters
params0 = getParams();

% this file contains generated parameters for a once good fit of the last slack
ModelParamsInitOptim_slack4

% Use the overlap function and the overlap function factor (Dan's NEW updated overlap)
params0.UseOverlap = true;
params0.UseOverlapFactor = true;

% correct force affecting the SR transition should be half, as they are
% distributed triangle-shaped
params0.useHalfActiveForSR = false;

% Use steady-state passive, that has been identified from Bakers
% experiments. No viscoelasticity
params0.UseTitinIdentifiedPassive = true;

% titin viscoelasticity minimally affect the restretch force overshoot
params0.UseTitinInterpolation = false;

params0.alpha0_L = params0.alpha0;
params0.alpha0_R = params0.alpha0;
%% Parameter fine tuning
% optional, already in params0

% New S_D state with the transitions
params0.UseSuperRelaxedADP = true;
params0.ksrd = 100;
params0.kmsrd = 0.25;
params0.ksr2srd = 5;
params0.sigma2 = 1e2;
params0.sigma_srd1 = params0.sigma1;
params0.sigma_srd2 = 5;

%   DAB fiddling with parameter values
% left and right side of the A1 detachment strain sensitivity
params0.alpha0_L = 15*params0.alpha0;
params0.alpha0_R = 0.85*params0.alpha0;

% Already set in ModelParamsInit
params0.kah = 150;
% params0.kadh = 0;
% params0.ka = 25.2761;
% params0.kd = 1.33821;
params0.k1 = 100; % DAB slowed down
params0.k_1 = 0*17.103; % not used, depends on Pi
% params0.k2 = 23.9605;
% params0.k_2 = 2.7901; 
% 
% params0.ksr0 = 71.3756;
params0.kmsr = 1.5*2.95303; % DAB increase
% params0.sigma1 = 22.1073;
% params0.sigma2 = 999999;
% 
% params0.dr = 0.0201168;
params0.kstiff1 = 10000;
params0.kstiff2 = 10000;
% params0.alpha0 = 36.1309;
% params0.alpha1 = 98.2684;
% params0.alpha_1 = 0;
% params0.alpha2 = 31.5061;
% params0.dr0 = 0;
% params0.dr1 = 0;
% params0.dr_1 = 0;
% params0.dr2 = -0.005;
% params0.dr3 = -0.01;
% params0.mu = 0.001;
% params0.kSE = 10879.7;
% params0.gamma = 3.00138;
% 
% params0.alpha2_L = 25.0348;
% params0.alpha2_R = 0.756688;
% params0.k2_R = 41291.3;
% params0.dr2_R = 0.00074017;

% WHAT TO RUN
figure(3);clf;
% Run the force-velocity profile
params0.RunForceVelocity = false;

% Which slack segment to run - try 'Last', 'All', 'First', 'Fourth' -
% defined in RunBakersExp around lines 300
params0.RunSlack = true;
params0.RunSlackSegments = 'All';

% run the force-SL length profile at steady state
params0.RunForceLengthEstim = false;

% Show all the ramps overlapped. Only for 'All' slacks
params0.EvalFitSlackOnset = false;

params0.ShowStatePlots  = false;
% params0.modelFcn = 'dPUdTCaSimpleAlternative2State';

LoadData;
% Using the new Combined DxDt
params0.modelFcn = 'dPUdT_CombinedTransitions';
% use the older one
params0.UseUniformTransitionFunc = false;

% only plot the strain-rate profile
params0.justPlotStateTransitionsFlag = false;
params0.EvalFitSlackOnset  = true;
params0.drawForceOnset = true;

params0.FudgeVmax = false;
% params0.FudgeA = 0;
% params0.FudgeB = 121.922;
% params0.FudgeC = -209.18;
params0.UseForceOnsetShift = true;
params0.BreakOnODEUnstable = false;

params0.UseA2Reattaching = false;
params0.useHalfActiveForSR = true;

% params0.k_2 = 50;

% params0.ghostSave = '';
params0.ghostLoad = 'NoReattach';


% tic
RunBakersExp;
% toctime = [toctime toc];


% mean(toctime)
% toctime = [];
%% fancy plot
fig = figure(80085);clf;hold on;
t0 = 2.7594;

plot(datatable(:, 1) - t0, datatable(:, 3), 'Color',[1 1 1]*0.6, LineWidth=6)
plot(out.t - t0, out.Force, 'k-', LineWidth=1);
ylabel('Stress (kPa)');
ylim([-1 80]);

yyaxis("right");
plot(datatable(:, 1) - t0, datatable(:, 2), 'k--', LineWidth=1.5);
xlim([-0.05, inf]);
ylim([1.85 2.25]);
legend('data', 'model', 'ML*', 'Location','best','fontsize',14);
% fontsize(14);
yticks([1.9 2 2.1 2.2])
ax = gca;
% This sets background color to black
ax.YColor = 'k';
xlabel('Time (s)');
ylabel('*Muscle length (L/L_0)');
% saveas(fig, 'Figures\proposal\')
%% Show states in time

% StatesInTime;

% further on its just some retarded experiments
return;



gh = load('Ghost_NoReattach_slack.mat')
params0 = gh.params0

%%
% clear;clf;
params0 = getParams();
ModelParamsInitNiceSlack;
params0.modelFcn = 'dPUdTCaSimpleAlternative2State';
params0.modelFcn = 'dPUdT_CombinedTransitions';
ksr0 = params0.kmsr;
params0.kmsr = params0.ksr0;
params0.ksr0 = ksr0;

% params0.drp3 = 0.01;
% params0.dr = 0.007;
params0.drp3 = params0.dr;

params0.k3m = 0;
params0.k3 = 1500; 

params0.RunStairs = true;    
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

    params0.NumberOfStates = 3;



LoadData;
params0.ShowStatePlots = true;
params0.EvalFitSlackOnset = true;
params0.RunForceLengthEstim = true;
RunBakersExp;

%%
params0 = getParams();
ModelParams_V1_Fmin2;
params0.UseForceOnsetShift = 1;

if params0.alpha0 ~= 0 && params0.alpha0_R == 0 && params0.alpha0_R == 0
    params0.alpha0_R = params0.alpha0 ;
    params0.alpha0_R = params0.alpha0 ;
    fprintf('Updating alpha... \n');
end

LoadData;
params0.ShowStatePlots = true;
params0.EvalFitSlackOnset = true;
params0.RunForceLengthEstim = true;
RunBakersExp;
%%
figure(243);clf;
nexttile;hold on;
feats_sim = extractSlackAttributes(out.t, out.Force, out.SL, velocitytable);


%% 
clf;
s = out.p1_0 + out.p2_0 + out.p3_0 + out.PuR + out.PuATP + out.SR + out.SRD;
plot(s)

%%
clf;
plot(out.t, out.RTD);
plot(out.t, out.SL, out.t, out.LXB)
%% test stretch
clf;
params0 = getParams();
ModelParamsInitNiceSlack;
params0 = getParams(params0, [], true);

Ns = 2;

params0.Vums = 0;
ss = params0.ss; % space size (length of the s for each of p1-p3)
SL = PU(Ns*ss + 3);
LSE = PU(Ns*ss + 4); % length of the serial stiffness
s = params0.s' + (-(SL - LSE) + params0.LXBpivot)/2;
% s = params.s' - (-(SL - LSE) + params.LXBpivot)/2;
s = flipud(-s);


for pi = 1:ss

    PU = params0.PU0';
    p1 = PU(1:ss);
    p2 = PU(ss+1:2*ss);
    

    p1(pi) = 1/params0.dS;
    PU(1:ss) = p1;
    [~, out] = dPUdT_CombinedTransitions(0, PU, params0);
    F1(pi) = out(2);

    p2(pi) = 1/params0.dS;
    PU(ss+1:2*ss) = p2;
    [~, out] = dPUdT_CombinedTransitions(0, PU, params0);
    F2(pi) = out(2);

    if params0.NumberOfStates == 3
        p3 = PU(2*ss +1:3*ss);
        p3(pi) = 1/params0.dS;
        PU(2*ss+1:3*ss) = p3;
        [~, out] = dPUdT_CombinedTransitions(0, PU, params0);
        F3(pi) = out(2);
    end
end

plot(s, F1, s, F2)

%% test extracting features
% Plot signal
figure(221); clf;
nexttile;hold on;
feats_data = extractSlackAttributes(datatable(:, 1), datatable(:, 3), datatable(:, 2), velocitytable);
title('data');
%%
% nexttile;hold on;
figure(221); clf;hold on;
feats_sim = extractSlackAttributes(out.t, out.Force, out.SL, velocitytable);
title('Sim');   

%% plot the features
fn = {'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'Am|SLslack|0', 'peak1_y|SLslack|0', 'peak1_y|v_restretch|0','peak1_dSL', 'peak2|v_restretch', 'vall_t|v_restretch|0.1', 'vall_y', 'vall2_dy|_|0', 'vall2_t|_|0', 'ovrsht_dy|_|0', 'steady', 'XTOR'};
figure
plotFeatures(feats_data, feats_sim, feats_ghost, fn);

% Eval features
%%
figure(8008135);clf;
cost = evalFeatureCost(feats_data, feats_sim, fn);
pie(cost(cost > 0), fn(cost > 0))

%% 
runStatesInTime(out, params0);
getParams