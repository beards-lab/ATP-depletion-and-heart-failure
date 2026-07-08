% RunCoupledCleanBaseline.m
% =========================================================================
% PART C: re-validate the coupled ATP mechanism on the CLEAN opt2state_v2
% baseline (high-ATP feature cost 3.22). Overlays the honest coupled ATP
% mechanism (UseAtpDetach + UseAdpTrap + UsePiForce), compensates k3 for the
% 8 mM rigor-detach gate so the 8 mM baseline is preserved, runs 8 mM (inert)
% and 2 mM (coupled), and scores the low-ATP relative-ratio objective
% (atpRatioCost) against the 03/27/2026 data.
%
% Needs data/protocol_03_27_2026_{2,8}mM_slack.mat. Outputs results/lowatp_partC.mat.
% Run: cd(root); addpath(genpath('.')); Analyses/LowATP_Mechanisms/RunCoupledCleanBaseline
% =========================================================================
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

% ---- clean baseline (opt2state_v2) ----
params0 = getParams();
run(fullfile(root,'params','opt2state_v2_opt.m'));
params0 = getParams(params0, [], true, false);
params0.RunForceVelocity=false; params0.RunKtr=false; params0.RunSlack=true;
params0.RunStairs=false; params0.RunForceVelocityTime=false;
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false;
params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0; params0.ShowStatePlots=0;
params0.RunSlackSegments='All'; params0.MaxRunTime=90;
params0.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
k3_base = params0.k3;

% ---- overlay honest coupled ATP mechanism; compensate k3 for the 8 mM detach gate ----
params0.UseAtpDetach=1; params0.UseAdpTrap=1; params0.UsePiForce=1;
params0.UsePiReversal=0; params0.UsePiReverseStroke=0;
params0.k3 = k3_base * (8 + params0.K_T_detach)/8;   % preserve 8 mM effective k3

p8 = params0; p8.MgATP=8; p8.MgADP=0;   p8.Pi=0;      % 8 mM (inert/compensated)
[~,~,fm8,~] = runSlackExperiment(p8);
p2 = params0; p2.MgATP=2; p2.MgADP=1.6; p2.Pi=0.8;    % 2 mM (coupled active)
[~,~,fm2,~] = runSlackExperiment(p2);

% ---- data target (relative ratio 2/8, CompareATPData logic) ----
f2d = load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8d = load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
forceFeats = {'steady','A','Am','peak2','peak1_y','vall_y'};  kinFeats = {'ktr','t0'};
target = struct('scalar',struct(),'vector',struct());
for i=1:numel(forceFeats); f=forceFeats{i}; target.scalar.(f)=mean(double(f2d.(f))./double(f8d.(f)),'omitnan'); end
for i=1:numel(kinFeats);   f=kinFeats{i};   target.vector.(f)=double(f2d.(f))./double(f8d.(f)); end
target.weights = struct('steady',50,'A',50,'Am',10,'peak2',5,'peak1_y',10,'vall_y',10, ...
                        'ktr',2,'t0',1,'restretchSlopeStart',0.1,'ktr2_overshoot',0.1);

[cost, detail] = atpRatioCost(fm2, fm8, target);
fprintf('\n==== LOW-ATP relative-ratio cost (coupled on clean baseline) = %.4f ====\n', cost);
save(fullfile(resdir,'lowatp_partC.mat'),'fm2','fm8','target','detail','cost');
