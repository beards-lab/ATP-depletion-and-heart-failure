% RunKmsrdSlope.m
% =========================================================================
% Does the ATP-dependent SRX-ADP return kmsrd produce the DATA's per-segment ktr SLOPE
% (0.64 0.59 0.53 0.47 0.46), not just the right mean? And can its FORCE-sensitivity sigma_srd1
% (RSRD2PD = kmsrd*exp(F_SR/sigma_srd1)*P_SRD) tilt the slope, since restretch force differs per segment?
% Base = params_2state_lowATP (ADP-trap). kmsrd(8mM)=baseline, kmsrd(2mM)=lowered (the decoupling).
% Run: cd(root); addpath(genpath('.')); Analyses/LowATP_kmsrd/RunKmsrdSlope
% =========================================================================
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end
params0=getParams(); run(fullfile(root,'params','params_2state_lowATP.m')); params0=getParams(params0,[],true,false);
params0.RunForceVelocity=false; params0.RunKtr=false; params0.RunSlack=true; params0.RunStairs=false; params0.RunForceVelocityTime=false;
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false; params0.RunSlackSegments='All'; params0.MaxRunTime=120;
params0.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat'; params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0; params0.ShowStatePlots=0;
kmB=params0.kmsrd; ss1B=params0.sigma_srd1;
fprintf('baseline kmsrd=%.2f sigma_srd1=%.2f\n', kmB, ss1B);

f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8d=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
ff={'steady','A','Am','peak2','peak1_y','vall_y'}; kf={'ktr','t0'}; target=struct('scalar',struct(),'vector',struct());
for i=1:numel(ff); target.scalar.(ff{i})=mean(double(f2d.(ff{i}))./double(f8d.(ff{i})),'omitnan'); end
for i=1:numel(kf); target.vector.(kf{i})=double(f2d.(kf{i}))./double(f8d.(kf{i})); end
target.weights=struct('steady',50,'A',50,'Am',10,'peak2',5,'peak1_y',10,'vall_y',10,'ktr',2,'t0',1,'restretchSlopeStart',0.1,'ktr2_overshoot',0.1);
fprintf('DATA ktr ratio (target slope): %s   (want DOWN across segs)\n', num2str(target.vector.ktr,'%.2f '));

cfgs = {  % {label, kmsrd@2mM, sigma_srd1}
  {'ref (no decouple)', kmB, ss1B}, ...
  {'kmsrd2=12',         12,  ss1B}, ...
  {'kmsrd2=8',           8,  ss1B}, ...
  {'kmsrd2=4',           4,  ss1B}, ...
  {'kmsrd2=8 ss1=10',    8,  10}, ...
  {'kmsrd2=8 ss1=80',    8,  80}, ...
  {'kmsrd2=8 ss1=200',   8,  200} };
res=struct([]);
for i=1:numel(cfgs)
  lab=cfgs{i}{1}; km2=cfgs{i}{2}; ss1=cfgs{i}{3};
  p8=params0; p8.MgATP=8;p8.MgADP=0;p8.Pi=0; p8.kmsrd=kmB; p8.sigma_srd1=ss1;
  p2=params0; p2.MgATP=2;p2.MgADP=1.6;p2.Pi=0.8; p2.kmsrd=km2; p2.sigma_srd1=ss1;
  c=NaN; kr=nan(1,5);
  try; [~,~,fm8,~]=runSlackExperiment(p8); [~,~,fm2,~]=runSlackExperiment(p2);
    [c,d]=atpRatioCost(fm2,fm8,target); kr=d.ktr.ratio(:)'; catch e; fprintf('  %s FAILED %s\n',lab,e.message); end
  fprintf('%-20s cost %6.3f | ktr %s\n', lab, c, num2str(kr,'%.2f '));
  res(i).lab=lab; res(i).cost=c; res(i).ktr=kr;
end
save(fullfile(resdir,'kmsrd_slope.mat'),'res','target');
disp('KMSRD SLOPE DONE');
