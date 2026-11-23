% RunSrxDecouple.m
% =========================================================================
% Tests the SRX + trap decoupling: an ATP-dependent SRX-ADP return kmsrd (kept at the 8 mM baseline,
% reduced at 2 mM) crossed with ADP-trap strength K_D. Shows kmsrd moves ktr while holding force ->
% breaks the 2-state force<->ktr weld. Best: K_D=0.7, kmsrd(2mM)=8 -> ratio cost 0.135 (from 0.47).
% Run: cd(root); addpath(genpath('.')); Analyses/LowATP_2State/RunSrxDecouple
% =========================================================================
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;

params0=getParams(); run(fullfile(root,'params','optfull_opt.m')); params0=getParams(params0,[],true,false);
params0.RunForceVelocity=false; params0.RunKtr=false; params0.RunSlack=true; params0.RunStairs=false; params0.RunForceVelocityTime=false;
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false; params0.RunSlackSegments='All'; params0.MaxRunTime=150;
params0.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat'; params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0; params0.ShowStatePlots=0;
params0.UseAdpTrap=1;

f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8d=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
ff={'steady','A','Am','peak2','peak1_y','vall_y'}; kf={'ktr','t0'}; target=struct('scalar',struct(),'vector',struct());
for i=1:numel(ff); target.scalar.(ff{i})=mean(double(f2d.(ff{i}))./double(f8d.(ff{i})),'omitnan'); end
for i=1:numel(kf); target.vector.(kf{i})=double(f2d.(kf{i}))./double(f8d.(kf{i})); end
target.weights=struct('steady',50,'A',50,'Am',10,'peak2',5,'peak1_y',10,'vall_y',10,'ktr',2,'t0',1,'restretchSlopeStart',0.1,'ktr2_overshoot',0.1);

p8=params0; p8.K_D=0.7; p8.MgATP=8; p8.MgADP=0; p8.Pi=0; [~,~,fm8,~]=runSlackExperiment(p8);   % 8 mM ref (kmsrd baseline)
fprintf('targets: force 1.18, ktr 0.54, t0 1.26 | kmsrd_base=%.2f\n', params0.kmsrd);
fprintf('%6s %8s %6s | %6s %6s %6s %6s\n','K_D','kmsrd2','cost','force','ktr','t0','peak1');
for KD=[0.35 0.5 0.7]; for km=[2 8 20.87 60]
  p=params0; p.K_D=KD; p.kmsrd=km; p.MgATP=2; p.MgADP=1.6; p.Pi=0.8; c=10;sR=NaN;ktrR=NaN;t0R=NaN;p1R=NaN;
  try; [~,~,fm2,~]=runSlackExperiment(p); [c,d]=atpRatioCost(fm2,fm8,target);
    sR=d.steady.ratio; ktrR=mean(d.ktr.ratio,'omitnan'); t0R=mean(d.t0.ratio,'omitnan'); p1R=d.peak1_y.ratio; catch; end
  if ~isfinite(c); c=10; end
  fprintf('%6.2f %8.2f %6.3f | %6.3f %6.3f %6.2f %6.3f\n', KD, km, c, sR, ktrR, t0R, p1R);
end; end
