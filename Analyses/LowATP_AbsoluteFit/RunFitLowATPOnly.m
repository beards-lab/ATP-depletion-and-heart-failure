% RunFitLowATPOnly.m
% =========================================================================
% Fit the 2 mM data with the 8 mM fit FROZEN, using only the ATP mechanisms that
% leave the 8 mM fit bit-identical.
%
% OBSERVE: dPUdT has exactly four [MgATP]-dependent paths. Only two are
%   MgATP_ref-normalised (gATP = 1 at 8 mM for ANY Km), i.e. provably inert at 8 mM:
%     K_T2  (UseAtpDetachR2, :482)  R2 detachment  -> force amplitude
%     K_srx (UseAtpKmsrd,    :539)  SRX-ADP return -> ktr
%   The other two are NOT normalised, so they would move the 8 mM baseline and are
%   disqualified here (and are falsified anyway -- see ../ApoPoolDetachment):
%     k2max/K_T2c (UseR2Ceiling, :472)   k_D0T/K_D0T (UseD0State, :504)
%   => the justifiable low-ATP lever set is exactly {K_T2, K_srx}.
%
% ORIENT: the 2 mM rsK error factorises
%     x3.77 (observed) = x2.61 (inherited 8 mM defect) * x1.47 (low-ATP specific)
%   With 8 mM frozen only the x1.47 is reachable, so rsK cannot be fixed here. Scoring
%   it at full weight would let it dominate (9.7 of 24.96) and drag the two gates into
%   compensating a defect they cannot reach, distorting the features they CAN fit.
%
% DECIDE: sweep the 2-D lever plane once and score EVERY objective from the same runs:
%     absFull   - absolute 2 mM cost, full feature list (rsK included)
%     absReach  - absolute 2 mM cost, post-restretch recovery cluster removed
%     ratio     - the legacy relative 2/8 cost (kept per request)
%   The 8 mM leg is computed ONCE: the gates are inert there by construction.
%
% Run: matlab -batch "run('Analyses/LowATP_AbsoluteFit/RunFitLowATPOnly.m')"
% =========================================================================
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

BASE = 'rskR2_w025_opt';                       % the FROZEN 8 mM fit

p = getParams(loadParams(BASE), [], true, false);
p.RunForceVelocity=false; p.RunForceVelocityTime=false; p.RunKtr=false;
p.RunStairs=false; p.RunSlackPassive=false; p.RunSlack=true;
p.EvalFeatures=true; p.BreakOnODEUnstable=false; p.RunSlackSegments='AllPar';
p.MaxRunTime=300; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
p.UseAtpDetachR2=1; p.UseAtpKmsrd=1; p.MgATP_ref=8;
p.UseR2Ceiling=0; p.UseD0State=0;              % disqualified: not 8 mM-inert

% ---- feature lists (slack-only day: no FV, no passive) ----
fn = p.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
fnAbs = fn(~(startsWith(fn,'FV') | startsWith(fn,'PS_')));
% the post-restretch RECOVERY cluster == what the inherited 8 mM defect dominates
RECOV = {'rsK','rsR2','coolDownLS','ovrsht_dy','vall2_dy'};
isRecov = @(s) any(cellfun(@(r) startsWith(s,r), RECOV));
fnReach = fnAbs(~cellfun(isRecov, fnAbs));
fprintf('objectives: absFull %d terms | absReach %d terms (dropped %s)\n', ...
    numel(fnAbs), numel(fnReach), strjoin(RECOV,','));

% ---- 8 mM reference: ONE run (gates inert at 8 mM for any Km) ----
q8=p; q8.MgATP=8; q8.MgADP=0; q8.Pi=0; q8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
[~,~,fm8,fd8]=runSlackExperiment(q8);
c8Abs  = sum(evalFeatureCost(fd8,fm8,fnAbs));
c8Reach= sum(evalFeatureCost(fd8,fm8,fnReach));
fprintf('8 mM FROZEN reference: absFull %.3f  absReach %.3f  (rsK %.1f vs data %.1f)\n', ...
    c8Abs, c8Reach, mean(double(fm8.rsK)), mean(double(fd8.rsK)));
tgt = buildAtpRatioTarget(load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data, fd8, 'sparse');

% ---- sweep the lever plane ----
K_T2v  = [0.8 1.4 1.97 3.0 4.5];
K_srxv = [5 9 13 20 30];
nT=numel(K_T2v); nS=numel(K_srxv);
AF=nan(nT,nS); AR=nan(nT,nS); RC=nan(nT,nS); RSK=nan(nT,nS); STE=nan(nT,nS); KTR=nan(nT,nS);
fprintf('\n%6s %6s | %8s %8s %8s | %7s %7s %7s\n','K_T2','K_srx','absFull','absReach','ratio','rsK','steady','ktr');
for i=1:nT
  for j=1:nS
    q=p; q.MgATP=2; q.MgADP=0; q.Pi=0.8;
    q.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
    q.K_T2=K_T2v(i); q.K_srx=K_srxv(j);
    try
      [~,~,fm2,fd2]=runSlackExperiment(q);
      AF(i,j)=sum(evalFeatureCost(fd2,fm2,fnAbs));
      AR(i,j)=sum(evalFeatureCost(fd2,fm2,fnReach));
      RC(i,j)=atpRatioCost(fm2,fm8,tgt);
      RSK(i,j)=mean(double(fm2.rsK),'omitnan');
      STE(i,j)=mean(double(fm2.steady),'omitnan');
      KTR(i,j)=mean(double(fm2.ktr),'omitnan');
      fprintf('%6.2f %6.1f | %8.3f %8.3f %8.4f | %7.1f %7.1f %7.1f\n', ...
        K_T2v(i),K_srxv(j),AF(i,j),AR(i,j),RC(i,j),RSK(i,j),STE(i,j),KTR(i,j));
    catch ME
      fprintf('%6.2f %6.1f | FAILED %s\n', K_T2v(i),K_srxv(j), ME.message);
    end
  end
end

f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
fprintf('\nDATA 2 mM: rsK %.1f  steady %.1f  ktr %.1f\n', ...
  mean(double(f2d.rsK)), mean(double(f2d.steady)), mean(double(f2d.ktr)));

% ---- optima per objective ----
pick = @(M) deal_min(M, K_T2v, K_srxv);
[bAF_T,bAF_S,bAF] = pick(AF);
[bAR_T,bAR_S,bAR] = pick(AR);
[bRC_T,bRC_S,bRC] = pick(RC);
fprintf('\n=== OPTIMA on the frozen 8 mM fit ===\n');
fprintf('  absFull  best %.3f  at K_T2=%.2f K_srx=%.1f\n', bAF, bAF_T, bAF_S);
fprintf('  absReach best %.3f  at K_T2=%.2f K_srx=%.1f   <-- the reachable fit\n', bAR, bAR_T, bAR_S);
fprintf('  ratio    best %.4f at K_T2=%.2f K_srx=%.1f   (legacy relative metric)\n', bRC, bRC_T, bRC_S);
fprintf('  do the three objectives agree on the levers?  %s\n', ...
  ternstr(bAF_T==bAR_T && bAF_S==bAR_S && bAR_T==bRC_T && bAR_S==bRC_S, 'YES', 'NO -- see table'));

% ---- local refinement of the reachable optimum ----
refT = bAR_T*[0.9 1.0 1.12]; refS = bAR_S*[0.78 1.0 1.22];
fprintf('\n--- refine around (%.2f, %.1f) ---\n', bAR_T, bAR_S);
for T=refT
  for S=refS
    q=p; q.MgATP=2; q.MgADP=0; q.Pi=0.8;
    q.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat'; q.K_T2=T; q.K_srx=S;
    try
      [~,~,fm2,fd2]=runSlackExperiment(q);
      c=sum(evalFeatureCost(fd2,fm2,fnReach));
      fprintf('  K_T2=%.2f K_srx=%.1f -> absReach %.4f (steady %.1f, ktr %.1f)\n', ...
        T,S,c,mean(double(fm2.steady)),mean(double(fm2.ktr)));
      if c < bAR; bAR=c; bAR_T=T; bAR_S=S; end
    catch; end
  end
end
fprintf('REFINED reachable optimum: K_T2=%.2f K_srx=%.1f  absReach=%.4f\n', bAR_T, bAR_S, bAR);

save(fullfile(resdir,'lowatp_only_sweep.mat'), ...
  'K_T2v','K_srxv','AF','AR','RC','RSK','STE','KTR','c8Abs','c8Reach','fnAbs','fnReach', ...
  'bAF','bAF_T','bAF_S','bAR','bAR_T','bAR_S','bRC','bRC_T','bRC_S');

% ---- bake the reachable optimum ----
dst=fullfile(root,'params','lowATP_only_opt.m');
copyfile(fullfile(root,'params',[BASE '.m']), dst);
fid=fopen(dst,'a');
fprintf(fid,'\r\n%%%% ---- Low-ATP fit on the FROZEN %s 8 mM fit (%s) ----\r\n', BASE, datestr(now,'yyyy-mm-dd'));
fprintf(fid,'%%%% Only the two MgATP_ref-normalised gates are used, so 8 mM is bit-identical.\r\n');
fprintf(fid,'params0.UseAtpDetachR2 = 1;\r\nparams0.UseAtpKmsrd = 1;\r\nparams0.MgATP_ref = 8;\r\n');
fprintf(fid,'params0.K_T2 = %.4f;\r\nparams0.K_srx = %.4f;\r\n', bAR_T, bAR_S);
fprintf(fid,'params0.MgATP = 2;\r\nparams0.MgADP = 0;\r\nparams0.Pi = 0.8;\r\n');
fprintf(fid,'%%%% absReach %.3f (8 mM floor %.3f) | absFull %.3f | ratio %.4f\r\n', bAR, c8Reach, bAF, bRC);
fclose(fid);
fprintf('\nSaved %s\n', dst);
disp('LOWATP-ONLY FIT DONE');

% ================= helpers =================
function [bx,by,bv]=deal_min(M,xv,yv)
  [bv,k]=min(M(:)); [i,j]=ind2sub(size(M),k); bx=xv(i); by=yv(j);
end
function s=ternstr(c,a,b); if c; s=a; else; s=b; end; end
