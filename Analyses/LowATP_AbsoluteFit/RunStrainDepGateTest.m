% RunStrainDepGateTest.m
% =========================================================================
% Does STRAIN-DEPENDENT ATP gating break the rsK level<->ratio weld?
%
% Hypothesis (physiology): positive strain inhibits ADP release, so for a strained
% post-stroke head the ATP-binding step sits in series behind a strain-slowed step and
% becomes more rate-controlling => apparent Km(ATP) for detachment RISES with strain.
%   Km_eff(s) = K_T2 * (1 + kAtpStrain * max(0, s - sAtpStrain0))
% Still MgATP_ref-normalised => 8 mM must be BIT-IDENTICAL for any setting.
%
% The separation the data demands:  rsK ratio 0.36  vs  ktr ratio 0.55
%   separation = rsKratio / ktrratio = 0.65   (data).  A strain-UNIFORM gate cannot
%   make this < 1 by much, because it gates both populations equally.
% Run: matlab -batch "run('Analyses/LowATP_AbsoluteFit/RunStrainDepGateTest.m')"
% =========================================================================
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); refreshPool(4); LoadData;
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

b=getParams(loadParams('rskR2_w025_opt'),[],true,false);
b.RunForceVelocity=false; b.RunForceVelocityTime=false; b.RunKtr=false; b.RunStairs=false;
b.RunSlackPassive=false; b.RunSlack=true; b.EvalFeatures=true; b.BreakOnODEUnstable=false;
b.RunSlackSegments='AllPar'; b.MaxRunTime=300; b.PlotEachSeparately=0; b.PlotFeatureFitting=0;
b.UseAtpDetachR2=1; b.UseAtpKmsrd=1; b.MgATP_ref=8; b.K_T2=1.80; b.K_srx=9.0;
b.UseR2Ceiling=0; b.UseD0State=0; b.UseStrainDepAtpGate=0; b.kAtpStrain=0; b.sAtpStrain0=0.005;

C={};
C{end+1}={'uniform gate (baseline)', struct()};
C{end+1}={'strain gate k=50',        struct('UseStrainDepAtpGate',1,'kAtpStrain',50)};
C{end+1}={'strain gate k=150',       struct('UseStrainDepAtpGate',1,'kAtpStrain',150)};
C{end+1}={'strain gate k=400',       struct('UseStrainDepAtpGate',1,'kAtpStrain',400)};
% the real test: pair the clean LEVEL fixer with the strain gate -> does the weld break?
C{end+1}={'eta0.3 (level fix) only',        struct('eta_M',b.eta_M*0.3)};
C{end+1}={'eta0.3 + strain gate k=150',     struct('eta_M',b.eta_M*0.3,'UseStrainDepAtpGate',1,'kAtpStrain',150)};
C{end+1}={'eta0.3 + strain gate k=400',     struct('eta_M',b.eta_M*0.3,'UseStrainDepAtpGate',1,'kAtpStrain',400)};

f8d=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
m=@(x)mean(double(x),'omitnan');
dR=m(f2d.rsK)/m(f8d.rsK); dK=m(f2d.ktr)/m(f8d.ktr);
fprintf('DATA: rsKratio %.2f  ktrratio %.2f  SEPARATION %.2f | rsK8 %.1f rsK2 %.1f steady2 %.1f\n\n', ...
    dR, dK, dR/dK, m(f8d.rsK), m(f2d.rsK), m(f2d.steady));

fprintf('%-28s | %6s %6s | %6s %6s %6s | %6s %6s %6s | %5s\n', ...
    'config','rsK8','inv8','rsK2','ktr2','std2','rsKrat','ktrrat','SEPAR','pk1_2');
ref8=[];
Rres=struct('name',{},'rsK8',{},'rsK2',{},'rat',{},'sep',{},'st2',{});
for i=1:numel(C)
  nm=C{i}{1}; ov=C{i}{2}; q=b; f=fieldnames(ov);
  for k=1:numel(f); q.(f{k})=ov.(f{k}); end
  q=getParams(q,[],true,false);
  try
    a=q; a.MgATP=8; a.MgADP=0; a.Pi=0; a.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
    [~,~,x8,~]=runSlackExperiment(a);
    c=q; c.MgATP=2; c.MgADP=0; c.Pi=0.8; c.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
    [~,~,x2,~]=runSlackExperiment(c);
    if isempty(ref8) || i==1; ref8=double(x8.steady(:)); end
    inv8=max(abs(double(x8.steady(:))-ref8));   % 8 mM invariance vs the gate-off reference
    rr=m(x2.rsK)/m(x8.rsK); kr=m(x2.ktr)/m(x8.ktr);
    fprintf('%-28s | %6.1f %6.3g | %6.1f %6.1f %6.1f | %6.2f %6.2f %6.2f | %5.1f\n', ...
        nm, m(x8.rsK), inv8, m(x2.rsK), m(x2.ktr), m(x2.steady), rr, kr, rr/kr, m(x2.peak1_y));
    Rres(end+1)=struct('name',nm,'rsK8',m(x8.rsK),'rsK2',m(x2.rsK),'rat',rr,'sep',rr/kr,'st2',m(x2.steady)); %#ok<SAGROW>
  catch ME
    fprintf('%-28s | FAILED %s\n', nm, ME.message);
  end
end
fprintf('%-28s | %6.1f %6s | %6.1f %6.1f %6.1f | %6.2f %6.2f %6.2f | %5.1f\n', ...
    'DATA TARGET', m(f8d.rsK),'-', m(f2d.rsK), m(f2d.ktr), m(f2d.steady), dR, dK, dR/dK, m(f2d.peak1_y));

fprintf('\nNOTE inv8 = max|d steady(8mM)| vs the gate-off run. MUST be 0: the gate is\n');
fprintf('MgATP_ref-normalised, so it cannot touch 8 mM (incl. the ktr overstretch peak).\n');
save(fullfile(resdir,'straindep_gate_test.mat'),'Rres','dR','dK');
disp('STRAIN-DEP GATE TEST DONE');
