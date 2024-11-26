%% simplest 
clear;
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
params0.ksrd = 10;
params0.kmsrd = 0.25;
params0.ksr2srd = 5;
params0.sigma_srd1 = params0.sigma1;

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

%% WHAT TO RUN
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

params0.ShowStatePlots  = true;
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

RunBakersExp;


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


