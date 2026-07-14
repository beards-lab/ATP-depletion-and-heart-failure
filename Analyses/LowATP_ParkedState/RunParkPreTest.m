% RunParkPreTest.m
% =========================================================================
% Direction pre-test for the forced-detachment parked state (UseForcedDetachmentPark).
% Step 1: 8 mM perturbation (parking ON vs OFF). Step 2: does parking add a per-segment ktr SLOPE
% in the RATIO (2mM/8mM) toward data (0.64->0.46) without harming peak1/steady?
% RESULT (2026-07-14): REFUTED — ktr ratio goes UP (wrong way), no slope, cost worsens (see conclusions.md).
% Parking is an AMPLITUDE lever (removes heads), not a RATE lever; ktr needs a rate/timescale lever.
% Run: cd(root); addpath(genpath('.')); Analyses/LowATP_ParkedState/RunParkPreTest
% =========================================================================
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
params0=getParams(); run(fullfile(root,'params','params_2state_lowATP.m')); params0=getParams(params0,[],true,false);
params0.RunForceVelocity=false; params0.RunKtr=false; params0.RunSlack=true; params0.RunStairs=false; params0.RunForceVelocityTime=false;
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false; params0.RunSlackSegments='All'; params0.MaxRunTime=120;
params0.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat'; params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0; params0.ShowStatePlots=0;
params0.K_T_park=4;

f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8d=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
ff={'steady','A','Am','peak2','peak1_y','vall_y'}; kf={'ktr','t0'}; target=struct('scalar',struct(),'vector',struct());
for i=1:numel(ff); target.scalar.(ff{i})=mean(double(f2d.(ff{i}))./double(f8d.(ff{i})),'omitnan'); end
for i=1:numel(kf); target.vector.(kf{i})=double(f2d.(kf{i}))./double(f8d.(kf{i})); end
target.weights=struct('steady',50,'A',50,'Am',10,'peak2',5,'peak1_y',10,'vall_y',10,'ktr',2,'t0',1,'restretchSlopeStart',0.1,'ktr2_overshoot',0.1);
fprintf('DATA ktr ratio (target slope): %s\n', num2str(target.vector.ktr,'%.2f '));

p8=params0; p8.UseForcedDetachmentPark=0; p8.MgATP=8;p8.MgADP=0;p8.Pi=0; [~,~,fm8,~]=runSlackExperiment(p8);
p2=params0; p2.UseForcedDetachmentPark=0; p2.MgATP=2;p2.MgADP=1.6;p2.Pi=0.8; [~,~,fm2,~]=runSlackExperiment(p2);
[c0,d0]=atpRatioCost(fm2,fm8,target);
fprintf('NO parking : cost %.3f | ktr ratio %s | peak1R %.2f steadyR %.2f\n', c0, num2str(d0.ktr.ratio,'%.2f '), d0.peak1_y.ratio, d0.steady.ratio);
for kp=[40 120]
  q8=params0; q8.UseForcedDetachmentPark=1; q8.kpark=kp; q8.kunpark=500; q8.MgATP=8;q8.MgADP=0;q8.Pi=0;
  q2=params0; q2.UseForcedDetachmentPark=1; q2.kpark=kp; q2.kunpark=500; q2.MgATP=2;q2.MgADP=1.6;q2.Pi=0.8;
  [~,~,g8,~]=runSlackExperiment(q8); [~,~,g2,~]=runSlackExperiment(q2);
  [c,d]=atpRatioCost(g2,g8,target);
  fprintf('kpark=%3d : cost %.3f | ktr ratio %s | peak1R %.2f steadyR %.2f\n', kp, c, num2str(d.ktr.ratio,'%.2f '), d.peak1_y.ratio, d.steady.ratio);
end
