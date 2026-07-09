% DriverLowATP_2State.m
% =====================================================================
% Runs the best 2-STATE low-ATP parametrization (params/params_2state_lowATP.m) and
% reports its relative-ratio fit to the 03/27/2026 data (2 mM vs 8 mM).
%
% ONE file serves both conditions:
%   * as-is (MgATP=2, MgADP=1.6, Pi=0.8)  -> low-ATP  (ratio cost ~0.47)
%   * with MgATP=8, MgADP=0, Pi=0         -> high-ATP (= optfull_opt baseline, cost 2.82)
% The ATP mechanism (UseAdpTrap + UsePiForce) is IDENTICALLY INERT at 8 mM, so the
% high-ATP fit is unaffected (verified max|d|=0).
%
% Run:  cd(root); addpath(genpath('.')); DriverLowATP_2State
% =====================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root)); LoadData;
%%
params0 = getParams(); run(fullfile(root,'params','params_2state_lowATP.m'));
params0 = getParams(params0, [], true, false);
params0.RunForceVelocity=false; params0.RunKtr=false; params0.RunSlack=true;
params0.RunStairs=false; params0.RunForceVelocityTime=false;
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false; params0.RunSlackSegments='AllPar';
params0.MaxRunTime=120; params0.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
params0.PlotEachSeparately=1; params0.ShowStatePlots=0; params0.PlotFeatureFitting=0;

tic
[~,~,fm2,~] = runSlackExperiment(params0);                     % low-ATP (baked 2 mM)
toc
ph=params0; ph.MgATP=8; ph.MgADP=0; ph.Pi=0; ph.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
tic
[~,~,fm8,~] = runSlackExperiment(ph);                          % high-ATP (8 mM)
toc

f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8d=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
ff={'steady','A','Am','peak2','peak1_y','vall_y'}; kf={'ktr','t0'};
target=struct('scalar',struct(),'vector',struct());
for i=1:numel(ff); target.scalar.(ff{i})=mean(double(f2d.(ff{i}))./double(f8d.(ff{i})),'omitnan'); end
for i=1:numel(kf); target.vector.(kf{i})=double(f2d.(kf{i}))./double(f8d.(kf{i})); end
target.weights=struct('steady',50,'A',50,'Am',10,'peak2',5,'peak1_y',10,'vall_y',10, ...
                      'ktr',2,'t0',1,'restretchSlopeStart',0.1,'ktr2_overshoot',0.1);
[c,det]=atpRatioCost(fm2,fm8,target);

fprintf('\n==== 2-state low-ATP relative-ratio cost = %.4f ====\n', c);
fprintf('%-10s  %10s  %10s  %8s\n','feature','model 2/8','data 2/8','w*cost');
sf=fieldnames(target.scalar);
for i=1:numel(sf); f=sf{i}; d=det.(f); fprintf('%-10s  %10.3f  %10.3f  %8.3f\n', f, d.ratio, d.target, d.cost); end
vf=fieldnames(target.vector);
for i=1:numel(vf); f=vf{i}; d=det.(f);
  fprintf('%-10s  model %s\n', f, num2str(d.ratio,'%5.2f '));
  fprintf('%-10s  data  %s  (w*cost %.3f)\n', '', num2str(d.target(:)','%5.2f '), d.cost);
end

%% 4-line feature comparison: low/high ATP x data/model, per slack segment.
% Colour = ATP level (red=2mM, blue=8mM); style = source (solid+filled=data,
% dashed+open=model). So the four lines are easily distinguishable in every tile.
feats_cell = {f2d, fm2, f8d, fm8};
labels     = {'2mM data', '2mM model', '8mM data', '8mM model'};
cRed = [0.85 0.16 0.16];  cBlu = [0.12 0.35 0.80];
colors     = [cRed; cRed; cBlu; cBlu];
markers    = {'o', 's', 'o', 's'};        % data = o, model = s
lineStyles = {'-', '--', '-', '--'};      % data = solid, model = dashed
fillMarkers= [true false true false];     % data = filled, model = open
fnPlot     = {'steady','A','Am','ktr','t0','peak1_y','vall_y','peak2'};
figure('Name','Low-ATP vs High-ATP: data vs model','Color','w');
plotMultipleFeatures(feats_cell, labels, colors, markers, fnPlot, lineStyles, fillMarkers);
sgtitle(sprintf('Absolute features per slack segment   (2-state low-ATP ratio cost %.3f)', c));
