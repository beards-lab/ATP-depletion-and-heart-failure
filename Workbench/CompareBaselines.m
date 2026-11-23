% CompareBaselines.m — how much does the 08-10 refinement of rskR2_w025_opt actually
% change behaviour vs the 08-04 copy I had been fitting low-ATP against?
% Also re-checks that the low-ATP gates remain 8 mM-inert on the NEW baseline.
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root)); refreshPool(6); LoadData;

cfg=@(s) s;  % readability
snaps={'rskR2_w025_opt_0804','rskR2_w025_opt'};
lbl={'08-04 (what I fitted on)','08-10 (your final best)'};
m=@(x)mean(double(x),'omitnan');
F=cell(1,2);
for i=1:2
  b=getParams(loadParams(snaps{i}),[],true,false);
  b.RunForceVelocity=false; b.RunForceVelocityTime=false; b.RunKtr=false; b.RunStairs=false;
  b.RunSlackPassive=false; b.RunSlack=true; b.EvalFeatures=true; b.BreakOnODEUnstable=false;
  b.RunSlackSegments='AllPar'; b.MaxRunTime=300; b.PlotEachSeparately=0; b.PlotFeatureFitting=0;
  b.UseAtpDetachR2=0; b.UseAtpKmsrd=0; b.UseAtpPowerGates=0; b.UseStrainDepAtpGate=0;
  b.UseR2Ceiling=0; b.UseD0State=0;
  b.MgATP=8; b.MgADP=0; b.Pi=0; b.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
  [~,~,fm,fd]=runSlackExperiment(b);
  fn=b.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
  fnS=fn(~(startsWith(fn,'FV')|startsWith(fn,'PS_')));
  c=sum(evalFeatureCost(fd,fm,fnS,2));
  F{i}=fm;
  fprintf('%-26s | slack cost %7.4f | steady %6.2f  peak1 %6.2f  ktr %6.2f  rsK %6.1f\n', ...
      lbl{i}, c, m(fm.steady), m(fm.peak1_y), m(fm.ktr), m(fm.rsK));
  if i==2; fdRef=fd; end
end
fprintf('%-26s | %18s | steady %6.2f  peak1 %6.2f  ktr %6.2f  rsK %6.1f\n','DATA 8 mM','', ...
    m(fdRef.steady), m(fdRef.peak1_y), m(fdRef.ktr), m(fdRef.rsK));

d=@(a,b_) max(abs(double(a(:))-double(b_(:))));
fprintf('\nmax |delta| between the two baselines: steady %.4g  peak1 %.4g  ktr %.4g  rsK %.4g\n', ...
    d(F{1}.steady,F{2}.steady), d(F{1}.peak1_y,F{2}.peak1_y), ...
    d(F{1}.ktr,F{2}.ktr), d(F{1}.rsK,F{2}.rsK));

% ---- the gates must still be 8 mM-inert on the NEW baseline ----
b=getParams(loadParams('rskR2_w025_opt'),[],true,false);
b.RunForceVelocity=false; b.RunForceVelocityTime=false; b.RunKtr=false; b.RunStairs=false;
b.RunSlackPassive=false; b.RunSlack=true; b.EvalFeatures=true; b.BreakOnODEUnstable=false;
b.RunSlackSegments='AllPar'; b.MaxRunTime=300; b.PlotEachSeparately=0; b.PlotFeatureFitting=0;
b.MgATP_ref=8; b.MgATP=8; b.Pi=0; b.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
b.UseR2Ceiling=0; b.UseD0State=0; b.UseStrainDepAtpGate=0;
off=b; off.UseAtpPowerGates=0; off.UseAtpDetachR2=0; off.UseAtpKmsrd=0;
on =b; on.UseAtpPowerGates=1; on.UseAtpDetachR2=1; on.K_T2=1.8; on.UseAtpKmsrd=1; on.K_srx=9;
P={'k2','kmsrd','kah','kamh','ka','kd','k1','k_1','ksr0','ksr2srd','ksrd2sr','tau_reg','A0','eta_M', ...
   'kstiff2','kSE','kSE_M','mu_neg','dr','dr2','P_bound_max'};
rng(3); v=randn(1,numel(P))*0.7;
for i=1:numel(P); on.(['pAtp_' P{i}])=v(i); end
[~,~,g0,~]=runSlackExperiment(off); [~,~,g1,~]=runSlackExperiment(on);
fprintf('\n8 mM inertness of ALL gates on the NEW baseline: max|dsteady| %.3g  max|drsK| %.3g  (MUST be 0)\n', ...
    d(g0.steady,g1.steady), d(g0.rsK,g1.rsK));
disp('COMPARE DONE');
