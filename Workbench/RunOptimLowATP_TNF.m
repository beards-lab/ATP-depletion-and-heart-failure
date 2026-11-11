% RunOptimLowATP_TNF.m
% Re-base the low-ATP fit onto the ThursdayNightFever high-ATP baseline (merged from Further-optim).
% Two ATP-driven levers, both MgATP_ref-normalized => 8 mM invariant by construction (verified max|Dsteady|=0):
%   K_T2  : Km(ATP) of R2 (P2 detachment) -> the FORCE-amplitude engine (the "k2" lever, via UseAtpDetachR2)
%   K_srx : Km(ATP) of the SRX-ADP return kmsrd -> the ktr weld-breaker + per-segment slope (UseAtpKmsrd)
% Objective = SPARSE ratio cost (steady/A/Am/peak2/peak1_y/vall_y means + ktr/t0 per-seg) -- the robust
% metric that defines the 0.1188 milestone. Enriched per-segment breakdown reported as diagnostic only,
% because two amplitude/kinetic levers cannot (and should not) chase baseline-shape residuals.
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root)); LoadData;

setrun=@(p) setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield( ...
  p,'RunForceVelocity',false),'RunKtr',false),'RunSlack',true),'RunStairs',false),'RunForceVelocityTime',false), ...
  'EvalFeatures',true),'BreakOnODEUnstable',false),'RunSlackSegments','AllPar'),'MaxRunTime',180), ...
  'velocitytableonfile','protocol_03_27_2026_8mM_slack.mat');

base0 = setrun(getParams(loadParams('ThursdayNightFever'),[],true,false));
% ---- apply the ATP-driven low-ATP mechanism (TNF ships with flags OFF) ----
base0.UseAdpTrap=0; base0.UsePiReversal=0; base0.UsePiForce=0; base0.UsePiReverseStroke=0;
base0.UseAtpDetachR2=1; base0.UseAtpKmsrd=1; base0.MgATP_ref=8;
base0.K_T2=2.0527; base0.K_srx=7.3503;          % seed = old-baseline optimum
base0.MgATP=2; base0.MgADP=0; base0.Pi=0.8;     % 2 mM condition (the run being scored)

% ---- 8 mM reference (levers inert at 8 mM) ----
p8=base0; p8.MgATP=8; p8.MgADP=0; p8.Pi=0;
[~,~,fm8,~]=runSlackExperiment(p8);
fprintf('8mM reference computed (TNF baseline).\n');

f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8d=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
targetFull   = buildAtpRatioTarget(f2d,f8d,'full');
targetSparse = buildAtpRatioTarget(f2d,f8d,'sparse');

% ---- baseline (seed) ----
[~,~,fm2_0,~]=runSlackExperiment(base0);
c0s=atpRatioCost(fm2_0,fm8,targetSparse);
c0f=atpRatioCost(fm2_0,fm8,targetFull);
fprintf('\nSEED  K_T2=%.3f K_srx=%.3f :  sparse=%.4f   enriched=%.4f\n',base0.K_T2,base0.K_srx,c0s,c0f);
printBreak(fm2_0,fm8,targetFull,'seed (enriched breakdown)');

% ---- optimize (K_T2,K_srx) in log-space against the SPARSE cost ----
objfun=@(x) evalObj(x,base0,fm8,targetSparse);
x0=log([base0.K_T2, base0.K_srx]);
opts=optimset('MaxFunEvals',80,'MaxIter',50,'TolFun',5e-5,'TolX',5e-4,'Display','iter');
[xopt,copts]=fminsearch(objfun,x0,opts);
K_T2o=exp(xopt(1)); K_srxo=exp(xopt(2));

p_opt=base0; p_opt.K_T2=K_T2o; p_opt.K_srx=K_srxo;
[~,~,fm2_o,~]=runSlackExperiment(p_opt);
coptf=atpRatioCost(fm2_o,fm8,targetFull);
fprintf('\n=== OPTIMUM ===  K_T2=%.4f  K_srx=%.4f :  sparse=%.4f (was %.4f)   enriched=%.4f (was %.4f)\n', ...
  K_T2o,K_srxo,copts,c0s,coptf,c0f);
printBreak(fm2_o,fm8,targetFull,'optimum (enriched breakdown)');

show={'steady','A','peak1_y','vall_y','peak2','ktr','t0'};
fprintf('\nfeature      modelR(opt) 2/8         dataR 2/8\n');
for i=1:numel(show); f=show{i}; if ~isfield(fm2_o,f); continue; end
  r=double(fm2_o.(f))./double(fm8.(f)); rd=double(f2d.(f))./double(f8d.(f));
  fprintf('  %-10s %-22s | %s\n', f, num2str(r(:)','%.2f '), num2str(rd(:)','%.2f '));
end

% ---- save the re-based one-file low-ATP snapshot: TNF baseline + ATP mechanism + optimized Km ----
dst=fullfile(root,'params','lowATP_TNF_optim.m');
copyfile(fullfile(root,'params','ThursdayNightFever.m'), dst);
fid=fopen(dst,'a');
fprintf(fid,'\r\n%%%% ---- Low-ATP mechanism re-based onto ThursdayNightFever (%s) ----\r\n', datestr(now,'yyyy-mm-dd'));
fprintf(fid,'%%%% High-ATP = this file with MgATP=8 (byte-identical; both gates normalize to 1 at MgATP_ref=8).\r\n');
fprintf(fid,'params0.UseAdpTrap = 0;\r\nparams0.UseAtpDetachR2 = 1;\r\nparams0.UseAtpKmsrd = 1;\r\nparams0.MgATP_ref = 8;\r\n');
fprintf(fid,'params0.K_T2 = %.4f;  %%%% R2 ATP-detach Km -> force amplitude\r\n', K_T2o);
fprintf(fid,'params0.K_srx = %.4f; %%%% SRX kmsrd return Km -> ktr magnitude + slope\r\n', K_srxo);
fprintf(fid,'params0.MgATP = 2;\r\nparams0.MgADP = 0;\r\nparams0.Pi = 0.8;\r\n');
fprintf(fid,'%%%% sparse ratio cost %.4f, enriched %.4f (seed K_T2=2.0527 K_srx=7.3503 -> sparse %.4f)\r\n', copts, coptf, c0s);
fclose(fid);
fprintf('\nSaved %s\n', dst);

% ---- ratio-by-feature figure ----
resdir=fullfile(root,'Analyses','LowATP_kmsrd','results'); vf=fieldnames(targetFull.vector);
fig=figure('Position',[40 40 1500 850],'Color','w'); tl=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
for i=1:numel(vf); f=vf{i}; nexttile(tl); hold on; box on;
  rd=double(targetFull.vector.(f)); r0=double(fm2_0.(f))./double(fm8.(f)); ro=double(fm2_o.(f))./double(fm8.(f));
  plot(1:numel(rd),rd,'ko-','LineWidth',1.4,'DisplayName','data ratio');
  plot(1:numel(r0),r0,'b.--','LineWidth',1.0,'DisplayName','model (seed)');
  plot(1:numel(ro),ro,'rx-','LineWidth',1.2,'DisplayName','model (opt)');
  title(sprintf('%s (w=%g)',f,targetFull.weights.(f)),'Interpreter','none'); xlabel('segment');
  if i==1; legend('Location','best','FontSize',7); end
end
sgtitle(sprintf('TNF baseline: 2mM/8mM ratio per feature (sparse %.4f -> %.4f)',c0s,copts));
exportgraphics(fig,fullfile(resdir,'TNF_OPTIM_ratio_by_feature.png'),'Resolution',140);
save(fullfile(resdir,'TNF_OPTIM_result.mat'),'fm2_o','fm8','targetFull','targetSparse','K_T2o','K_srxo','copts','coptf');
disp('TNF OPTIM DONE');

% ================= helpers =================
function cost=evalObj(x,base0,fm8,target)
  K_T2=exp(x(1)); K_srx=exp(x(2));
  if K_T2<0.3||K_T2>15||K_srx<1||K_srx>40; cost=1e6; return; end
  p=base0; p.K_T2=K_T2; p.K_srx=K_srx;
  try; [~,~,fm2,~]=runSlackExperiment(p); cost=atpRatioCost(fm2,fm8,target);
    if ~isfinite(cost); cost=1e6; end
  catch; cost=1e6; end
end
function printBreak(fm2,fm8,target,tag)
  [~,d]=atpRatioCost(fm2,fm8,target); vf=fieldnames(target.vector);
  fprintf('  --- per-feature enriched cost (%s) ---\n',tag);
  rows=cell(numel(vf),2);
  for i=1:numel(vf); f=vf{i}; rows{i,1}=f; rows{i,2}=d.(f).cost; end
  [~,ord]=sort(cell2mat(rows(:,2)),'descend');
  for k=1:numel(ord); i=ord(k); fprintf('   %-20s cost %7.4f   modelR %s\n', rows{i,1}, rows{i,2}, num2str(d.(rows{i,1}).ratio(:)','%.2f ')); end
end
