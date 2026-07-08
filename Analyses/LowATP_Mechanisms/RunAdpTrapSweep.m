% RunAdpTrapSweep.m
% =========================================================================
% PART A: the 2-state low-ATP levers are K_D (ADP-trap strength, via g2) and K_Pi (Pi force trim).
% (UseAtpDetach is INERT in 2-state — no P3.) Sweeps K_D x K_Pi on the clean opt2state_v2 baseline,
% scores the low-ATP relative-ratio objective, and plots the force<->onset weld.
% Reuses fm8/target from RunCoupledCleanBaseline (results/lowatp_partC.mat) if present.
% Run: cd(root); addpath(genpath('.')); Analyses/LowATP_Mechanisms/RunAdpTrapSweep
% =========================================================================
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

params0=getParams(); run(fullfile(root,'params','opt2state_v2_opt.m')); params0=getParams(params0,[],true,false);
params0.RunForceVelocity=false; params0.RunKtr=false; params0.RunSlack=true; params0.RunStairs=false; params0.RunForceVelocityTime=false;
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false; params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0; params0.ShowStatePlots=0;
params0.RunSlackSegments='All'; params0.MaxRunTime=150; params0.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
params0.UseAtpDetach=0; params0.UseAdpTrap=1; params0.UsePiForce=1; params0.UsePiReversal=0;

% fm8 (8 mM baseline) + data target
pc = fullfile(resdir,'lowatp_partC.mat');
if exist(pc,'file'); L=load(pc); fm8=L.fm8; target=L.target;
else
  p8=params0; p8.MgATP=8; p8.MgADP=0; p8.Pi=0; [~,~,fm8,~]=runSlackExperiment(p8);
  f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
  f8d=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
  ff={'steady','A','Am','peak2','peak1_y','vall_y'}; kf={'ktr','t0'}; target=struct('scalar',struct(),'vector',struct());
  for i=1:numel(ff); target.scalar.(ff{i})=mean(double(f2d.(ff{i}))./double(f8d.(ff{i})),'omitnan'); end
  for i=1:numel(kf); target.vector.(kf{i})=double(f2d.(kf{i}))./double(f8d.(kf{i})); end
  target.weights=struct('steady',50,'A',50,'Am',10,'peak2',5,'peak1_y',10,'vall_y',10,'ktr',2,'t0',1,'restretchSlopeStart',0.1,'ktr2_overshoot',0.1);
end

KD_grid=[0.1 0.194 0.35 0.6 1.0]; KPi_grid=[3 4 6]; res=[]; best=struct('cost',inf);
for KD=KD_grid; for KPi=KPi_grid
  p=params0; p.K_D=KD; p.K_Pi=KPi; p.MgATP=2; p.MgADP=1.6; p.Pi=0.8; c=10; sR=NaN;t0R=NaN;p1R=NaN;
  try; [~,~,fm2,~]=runSlackExperiment(p); [c,d]=atpRatioCost(fm2,fm8,target);
    sR=d.steady.ratio; t0R=mean(d.t0.ratio,'omitnan'); p1R=d.peak1_y.ratio; catch; end
  if ~isfinite(c); c=10; end
  fprintf('KD=%5.3f KPi=%3.1f -> cost %.4f | force %.3f  t0 %.2f  peak1 %.3f\n',KD,KPi,c,sR,t0R,p1R);
  res=[res; KD KPi c sR t0R p1R];
  if c<best.cost; best=struct('cost',c,'KD',KD,'KPi',KPi); end
end; end
save(fullfile(resdir,'partA_grid2.mat'),'res','best','KD_grid','KPi_grid');
fprintf('\nbest cost %.4f at KD=%.3f KPi=%.1f\n',best.cost,best.KD,best.KPi);

% figure: force<->onset weld (K_Pi=4 slice)
sel=res(:,2)==4; [KD,ix]=sort(res(sel,1)); t=res(sel,:); t=t(ix,:);
fig=figure('Position',[100 100 820 540]);
yyaxis left; plot(KD,t(:,4),'-o','LineWidth',2.2); hold on; plot(KD,t(:,6),'-^','LineWidth',1.4);
yline(1.18,'--','force target 1.18'); ylabel('force / peak1 ratio (2/8)'); ylim([0.9 2.0])
yyaxis right; plot(KD,t(:,5),'-s','LineWidth',2.2); yline(1.26,':','onset target ~1.26');
ylabel('t0 onset ratio (2/8)'); ylim([1 4.4]); set(gca,'XScale','log'); grid on
xlabel('K_D (mM)  [ADP-trap strength; smaller = stronger]');
title({'2-state low-ATP: force gain and onset delay welded to one knob','force wants K_D\approx0.19, onset wants K_D\approx1.0 \Rightarrow needs a 3rd state'});
legend({'steady (force)','peak1','force target','t0 (onset)','onset target'},'Location','northwest');
exportgraphics(fig, fullfile(resdir,'partA_force_t0_weld.png'),'Resolution',130);
