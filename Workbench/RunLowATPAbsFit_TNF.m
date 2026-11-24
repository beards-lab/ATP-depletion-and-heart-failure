% RunLowATPAbsFit_TNF.m
% Fit the low-ATP levers (K_T2 = R2 ATP-detach Km; K_srx = SRX kmsrd-return Km) to the ABSOLUTE
% 2 mM slack features -- the SAME feature-cost philosophy used for the high-ATP fit (evalFeatureCost
% with the slack selector list) -- and render a nice visual fit of the 2 mM data.
% The 8 mM baseline (ThursdayNightFever) stays byte-identical: both gates are inert at MgATP_ref=8.
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root)); LoadData;
resdir=fullfile(root,'Analyses','LowATP_kmsrd','results'); if ~exist(resdir,'dir'); mkdir(resdir); end

% ---- baseline + low-ATP mechanism ----
base = getParams(loadParams('ThursdayNightFever'),[],true,false);
base.RunForceVelocity=false; base.RunKtr=false; base.RunSlack=true; base.RunStairs=false;
base.RunForceVelocityTime=false; base.EvalFeatures=true; base.BreakOnODEUnstable=false; base.MaxRunTime=200;
base.UseAdpTrap=0; base.UsePiReversal=0; base.UsePiForce=0; base.UsePiReverseStroke=0;
base.UseAtpDetachR2=1; base.UseAtpKmsrd=1; base.MgATP_ref=8;
base.K_T2=1.90; base.K_srx=7.72;                 % seed = ratio optimum on TNF

% ---- high-ATP feature selectors (absolute; same philosophy as the 8 mM fit) ----
% ktr up-weighted (30) so the SRX lever K_srx is chosen for ktr precision (per user).
fnabs = {'steady|50','A|50','peak1_y|10','peak1_dSL|1','vall_y|10','peak2|5','t0|2', ...
         'ktr|30','restretchSlopeStart|1','vall2_dy|1','ovrsht_dy|0.3'};

% ---- 8 mM reference (drive with 8 mM protocol): absolute cost for context ----
p8=base; p8.RunSlackSegments='All'; p8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
p8.MgATP=8; p8.MgADP=0; p8.Pi=0;
[~,out8,fm8,fd8]=runSlackExperiment(p8);
c8=sum(evalFeatureCost(fd8,fm8,fnabs,2));
fprintf('8 mM absolute feature cost (reference) = %.3f\n', c8);

% ---- 2 mM run driven by the 2 mM protocol ----
p2=base; p2.RunSlackSegments='AllPar'; p2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
p2.MgATP=2; p2.MgADP=0; p2.Pi=0.8;

objfun=@(x) evalAbs(x,p2,fnabs);
c2seed=objfun(log([base.K_T2, base.K_srx]));
fprintf('2 mM absolute feature cost  seed (K_T2=%.3f K_srx=%.3f) = %.3f\n', base.K_T2, base.K_srx, c2seed);

opts=optimset('MaxFunEvals',110,'MaxIter',80,'TolFun',3e-4,'TolX',3e-4,'Display','iter');
[xopt,c2opt]=fminsearch(objfun,log([base.K_T2, base.K_srx]),opts);
K_T2o=exp(xopt(1)); K_srxo=exp(xopt(2));
fprintf('\n=== 2 mM ABSOLUTE optimum ===  K_T2=%.4f  K_srx=%.4f :  cost %.3f (seed %.3f; 8mM ref %.3f)\n', ...
  K_T2o,K_srxo,c2opt,c2seed,c8);

% ---- final render run (serial 'All' -> one concatenated trace) ----
p2o=p2; p2o.RunSlackSegments='All'; p2o.K_T2=K_T2o; p2o.K_srx=K_srxo;
[~,out2,fm2,fd2]=runSlackExperiment(p2o);
[cc,~,cu]=evalFeatureCost(fd2,fm2,fnabs,2);
fprintf('final 2 mM abs cost = %.3f\n', sum(cc));

% per-feature model-vs-data (absolute)
show={'steady','A','peak1_y','vall_y','peak2','ktr','t0','restretchSlopeStart'};
fprintf('\nfeature        model(2mM)                  data(2mM)\n');
for i=1:numel(show); f=show{i}; if ~isfield(fm2,f)||~isfield(fd2,f); continue; end
  fprintf('  %-14s %-26s | %s\n', f, num2str(double(fm2.(f)(:))','%6.3g '), num2str(double(fd2.(f)(:))','%6.3g '));
end

% ================= FIGURES =================
% (1) 2 mM slack-restretch trace: model vs data
tmin=min(out2.datatable(:,1)); tmax=max(out2.datatable(:,1));
fig1=figure('Position',[40 40 1500 520],'Color','w');
plot(out2.datatable(:,1), out2.datatable(:,3),'k-','LineWidth',1.2,'DisplayName','data (2 mM)'); hold on; box on;
plot(out2.t, out2.Force,'r-','LineWidth',1.6,'DisplayName','model (2 mM)');
xlim([tmin-0.5, tmax+0.5]); xlabel('time (s)'); ylabel('Force (rel.)');
title(sprintf('Low-ATP (2 mM) slack-restretch: model vs data on ThursdayNightFever baseline  (abs feature cost %.2f; 8 mM ref %.2f)', sum(cc), c8));
legend('Location','best'); grid on;
exportgraphics(fig1, fullfile(resdir,'TNF_2mM_best_trace.png'),'Resolution',140);

% (2) 2 mM feature-by-feature: model vs data
figF=figure('Position',[40 40 1500 820],'Color','w');
plotFeatures(fd2, fm2, [], {'steady','A','Am','ktr','t0','peak1_y','vall_y','peak2','restretchSlopeStart'});
sgtitle(sprintf('Low-ATP (2 mM): absolute features, data (o) vs model (x)   [abs cost %.2f]', sum(cc)));
exportgraphics(figF, fullfile(resdir,'TNF_2mM_best_features.png'),'Resolution',140);

% ---- save the re-based one-file low-ATP snapshot ----
dst=fullfile(root,'params','lowATP_TNF_best.m');
copyfile(fullfile(root,'params','ThursdayNightFever.m'), dst);
fid=fopen(dst,'a');
fprintf(fid,'\r\n%%%% ---- Low-ATP fit to ABSOLUTE 2 mM features on ThursdayNightFever (%s) ----\r\n', datestr(now,'yyyy-mm-dd'));
fprintf(fid,'params0.UseAdpTrap = 0;\r\nparams0.UseAtpDetachR2 = 1;\r\nparams0.UseAtpKmsrd = 1;\r\nparams0.MgATP_ref = 8;\r\n');
fprintf(fid,'params0.K_T2 = %.4f;\r\nparams0.K_srx = %.4f;\r\n', K_T2o, K_srxo);
fprintf(fid,'params0.MgATP = 2;\r\nparams0.MgADP = 0;\r\nparams0.Pi = 0.8;\r\n');
fprintf(fid,'%%%% 2 mM absolute feature cost %.3f (8 mM ref %.3f)\r\n', sum(cc), c8);
fclose(fid);
fprintf('\nSaved %s\n', dst);
save(fullfile(resdir,'TNF_best_result.mat'),'fm2','fd2','fm8','fd8','K_T2o','K_srxo','c2opt','c8','cc');
disp('ABS FIT DONE');

% ================= helpers =================
function c = evalAbs(x, p2, fnabs)
  K_T2=exp(x(1)); K_srx=exp(x(2));
  if K_T2<0.3||K_T2>15||K_srx<1||K_srx>40; c=1e6; return; end
  p=p2; p.K_T2=K_T2; p.K_srx=K_srx;
  try; [~,~,fm2,fd2]=runSlackExperiment(p); c=sum(evalFeatureCost(fd2,fm2,fnabs,2));
    if ~isfinite(c); c=1e6; end
  catch; c=1e6; end
end
