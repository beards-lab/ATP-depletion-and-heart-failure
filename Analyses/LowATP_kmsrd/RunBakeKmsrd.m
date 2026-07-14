% RunBakeKmsrd.m
% =========================================================================
% Bake + verify the ONE-FILE ATP-dependent SRX-ADP return: params/params_2state_kmsrd.m
% (= params_2state_lowATP + UseAtpKmsrd). kmsrd_eff = kmsrd*gATP,
% gATP = [MgATP/(MgATP+K_srx)]/[MgATP_ref/(MgATP_ref+K_srx)] (=1 at 8 mM => high-ATP byte-identical).
% Verifies: (a) 8 mM max|d steady| = 0 vs optfull, (b) 2 mM cost ~0.135 + ktr slope, from concentrations.
% NB: load snapshots with loadParams (sandboxed) — never run() two snapshots in the base workspace
% (each assigns to `params0` and the second clobbers the first).
% Run: cd(root); addpath(genpath('.')); Analyses/LowATP_kmsrd/RunBakeKmsrd
% =========================================================================
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;

% (re)bake params_2state_kmsrd.m = params_2state_lowATP + the flag
src=fullfile(root,'params','params_2state_lowATP.m'); dst=fullfile(root,'params','params_2state_kmsrd.m');
if ~exist(dst,'file')
  copyfile(src,dst); fid=fopen(dst,'a');
  fprintf(fid,'\r\n%%%% ---- ATP-dependent SRX-ADP return (one-file weld-breaker) ----\r\n');
  fprintf(fid,'params0.UseAtpKmsrd = 1;\r\nparams0.K_srx = 9.3;\r\nparams0.MgATP_ref = 8;\r\n'); fclose(fid);
end

setrun=@(p) setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield( ...
  p,'RunForceVelocity',false),'RunKtr',false),'RunSlack',true),'RunStairs',false),'RunForceVelocityTime',false), ...
  'EvalFeatures',true),'BreakOnODEUnstable',false),'RunSlackSegments','All'),'MaxRunTime',120),'velocitytableonfile','protocol_03_27_2026_8mM_slack.mat');

params0 = setrun(getParams(loadParams('params_2state_kmsrd'), [], true, false));
opt     = setrun(getParams(loadParams('optfull_opt'),        [], true, false));

f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8d=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
ff={'steady','A','Am','peak2','peak1_y','vall_y'}; kf={'ktr','t0'}; target=struct('scalar',struct(),'vector',struct());
for i=1:numel(ff); target.scalar.(ff{i})=mean(double(f2d.(ff{i}))./double(f8d.(ff{i})),'omitnan'); end
for i=1:numel(kf); target.vector.(kf{i})=double(f2d.(kf{i}))./double(f8d.(kf{i})); end
target.weights=struct('steady',50,'A',50,'Am',10,'peak2',5,'peak1_y',10,'vall_y',10,'ktr',2,'t0',1,'restretchSlopeStart',0.1,'ktr2_overshoot',0.1);

opt.MgATP=8;opt.MgADP=0;opt.Pi=0; [~,~,fmO,~]=runSlackExperiment(opt);
ph=params0; ph.MgATP=8;ph.MgADP=0;ph.Pi=0; [~,~,fm8,~]=runSlackExperiment(ph);
fprintf('8mM invariance (params_2state_kmsrd @8mM vs optfull): max|d steady| = %.4g\n', max(abs(double(fm8.steady)-double(fmO.steady))));
[~,~,fm2,~]=runSlackExperiment(params0);   % 2 mM from the file
[c,d]=atpRatioCost(fm2,fm8,target);
fprintf('2mM one-file kmsrd: cost %.4f | ktr %s (data %s) | steadyR %.3f peak1R %.3f\n', ...
  c, num2str(d.ktr.ratio(:)','%.2f '), num2str(target.vector.ktr(:)','%.2f '), d.steady.ratio, d.peak1_y.ratio);
