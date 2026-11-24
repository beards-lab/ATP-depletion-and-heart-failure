% RunMechanismRatioScan.m
% =========================================================================
% Can any of the OTHER mechanisms push the post-restretch recovery the right way --
% now that refitting 8 mM is allowed?
%
% ../ApoPoolDetachment scored the ceiling / D0 pool at 8 mM ONLY. The question never
% asked: do they move the ATP RATIO? The data's rsK ratio is x0.36 (2 mM 15.9 /
% 8 mM 43.7); the gates alone reach only x0.53. A mechanism is useful here only if it
% moves BOTH the 8 mM level DOWN (114 -> 44) and the ratio DOWN (0.53 -> 0.36),
% without destroying the restretch peak (the guardrail that falsified the ceiling).
%
% Each config is run at 8 mM and 2 mM. Guardrails: peak1_y (data 8 mM ~111, the
% restretch peak that pins R2(s)) and steady (data 77.3 / 91.3).
% Run: matlab -batch "run('Analyses/LowATP_AbsoluteFit/RunMechanismRatioScan.m')"
% =========================================================================
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

base = getParams(loadParams('rskR2_w025_opt'), [], true, false);
base.RunForceVelocity=false; base.RunForceVelocityTime=false; base.RunKtr=false;
base.RunStairs=false; base.RunSlackPassive=false; base.RunSlack=true;
base.EvalFeatures=true; base.BreakOnODEUnstable=false; base.RunSlackSegments='AllPar';
base.MaxRunTime=300; base.PlotEachSeparately=0; base.PlotFeatureFitting=0;
base.UseAtpDetachR2=1; base.UseAtpKmsrd=1; base.MgATP_ref=8; base.K_T2=1.80; base.K_srx=9.0;
base.UseR2Ceiling=0; base.UseD0State=0;

fn = base.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
fnAbs = fn(~(startsWith(fn,'FV')|startsWith(fn,'PS_')));

% ---- configs: name, then a struct of overrides ----
C = {};
C{end+1} = {'baseline (gates only)',      struct()};
C{end+1} = {'k2 x0.6',                    struct('k2',base.k2*0.6)};
C{end+1} = {'eta_M x0.3',                 struct('eta_M',base.eta_M*0.3)};
C{end+1} = {'k2 x0.6 + eta_M x0.3',       struct('k2',base.k2*0.6,'eta_M',base.eta_M*0.3)};
C{end+1} = {'ceiling k2max=1500',         struct('UseR2Ceiling',1,'k2max',1500,'K_T2c',0.15)};
C{end+1} = {'ceiling k2max=800',          struct('UseR2Ceiling',1,'k2max',800,'K_T2c',0.15)};
C{end+1} = {'ceil400+D0+rip3000',         struct('UseR2Ceiling',1,'k2max',400,'K_T2c',0.15, ...
                                                 'UseD0State',1,'k_D0T',40,'K_D0T',0.15,'k2rip',3000,'dr3',-0.025)};
C{end+1} = {'ceil300+D0+rip10000',        struct('UseR2Ceiling',1,'k2max',300,'K_T2c',0.15, ...
                                                 'UseD0State',1,'k_D0T',40,'K_D0T',0.15,'k2rip',10000,'dr3',-0.025)};
% The ATP question: make the ceiling itself [ATP]-sensitive by raising its Km.
C{end+1} = {'ceiling k2max=1500 K_T2c=6', struct('UseR2Ceiling',1,'k2max',1500,'K_T2c',6)};
C{end+1} = {'ceil800 K_T2c=6 +D0+rip3000',struct('UseR2Ceiling',1,'k2max',800,'K_T2c',6, ...
                                                 'UseD0State',1,'k_D0T',40,'K_D0T',0.15,'k2rip',3000,'dr3',-0.025)};

f8d=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
tgt8=[mean(double(f8d.rsK)) mean(double(f8d.peak1_y)) mean(double(f8d.steady))];
tgt2=[mean(double(f2d.rsK)) mean(double(f2d.peak1_y)) mean(double(f2d.steady))];
fprintf('DATA  8mM: rsK %.1f peak1 %.1f steady %.1f | 2mM: rsK %.1f peak1 %.1f steady %.1f | rsK ratio %.2f\n\n', ...
    tgt8, tgt2, tgt2(1)/tgt8(1));

fprintf('%-26s | %6s %6s %5s | %6s %6s %5s | %7s %7s\n', ...
    'config','rsK8','pk1_8','std8','rsK2','pk1_2','std2','rsKrat','absCost2');
R = struct('name',{},'rsK8',{},'rsK2',{},'rat',{},'pk8',{},'st8',{},'c2',{});
for i=1:numel(C)
  nm=C{i}{1}; ov=C{i}{2};
  q=base; f=fieldnames(ov); for k=1:numel(f); q.(f{k})=ov.(f{k}); end
  q=getParams(q,[],true,false);   % REQUIRED: UseD0State changes the state-vector length
  try
    a=q; a.MgATP=8; a.MgADP=0; a.Pi=0;   a.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
    [~,~,m8,~]=runSlackExperiment(a);
    b=q; b.MgATP=2; b.MgADP=0; b.Pi=0.8; b.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
    [~,~,m2,d2]=runSlackExperiment(b);
    r8=mean(double(m8.rsK),'omitnan'); r2=mean(double(m2.rsK),'omitnan');
    p8=mean(double(m8.peak1_y),'omitnan'); s8=mean(double(m8.steady),'omitnan');
    p2=mean(double(m2.peak1_y),'omitnan'); s2=mean(double(m2.steady),'omitnan');
    c2=sum(evalFeatureCost(d2,m2,fnAbs));
    fprintf('%-26s | %6.1f %6.1f %5.1f | %6.1f %6.1f %5.1f | %7.2f %7.1f\n', ...
        nm, r8,p8,s8, r2,p2,s2, r2/r8, c2);
    R(end+1)=struct('name',nm,'rsK8',r8,'rsK2',r2,'rat',r2/r8,'pk8',p8,'st8',s8,'c2',c2); %#ok<SAGROW>
  catch ME
    fprintf('%-26s | FAILED: %s\n', nm, ME.message);
  end
end

fprintf('\nTARGETS                    | %6.1f %6.1f %5.1f | %6.1f %6.1f %5.1f | %7.2f\n', ...
    tgt8, tgt2, tgt2(1)/tgt8(1));
fprintf('\nVERDICT (needs rsK8 down toward %.1f, ratio down toward %.2f, peak1_8 held near %.1f):\n', ...
    tgt8(1), tgt2(1)/tgt8(1), tgt8(2));
for i=1:numel(R)
  dRsk = R(i).rsK8 - R(1).rsK8; dRat = R(i).rat - R(1).rat; dPk = R(i).pk8 - R(1).pk8;
  fprintf('  %-26s rsK8 %+7.1f | ratio %+6.3f | peak1_8 %+6.1f  %s\n', ...
      R(i).name, dRsk, dRat, dPk, verdict(dRsk, dRat, dPk));
end
save(fullfile(resdir,'mechanism_ratio_scan.mat'),'R','tgt8','tgt2');
disp('MECHANISM RATIO SCAN DONE');

function s = verdict(dRsk, dRat, dPk)
  good = (dRsk < -5) && (dRat < -0.02);
  if good && abs(dPk) < 6;  s = '<= HELPS (peak held)';
  elseif good;              s = '<= helps rsK but MOVES THE PEAK';
  elseif dRsk < -5;         s = '   rsK down, ratio not';
  else;                     s = '   no / wrong way';
  end
end
