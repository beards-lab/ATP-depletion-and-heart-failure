% RunOptimLowATP.m
% Re-optimize the low-ATP levers (K_T2, K_srx) against the ENRICHED per-segment ratio cost
% (buildAtpRatioTarget 'full' -> every well-behaved slack feature, per segment, mirroring the
% high-ATP selector list). Both levers are MgATP_ref-normalized => 8 mM invariant by construction,
% so the 8 mM reference is computed ONCE.
root = 'C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root)); LoadData;

setrun = @(p) setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield( ...
  p,'RunForceVelocity',false),'RunKtr',false),'RunSlack',true),'RunStairs',false),'RunForceVelocityTime',false), ...
  'EvalFeatures',true),'BreakOnODEUnstable',false),'RunSlackSegments','AllPar'),'MaxRunTime',120), ...
  'velocitytableonfile','protocol_03_27_2026_8mM_slack.mat');

base0 = setrun(getParams(loadParams('params_2state_final'), [], true, false));
p8 = base0; p8.MgATP=8; p8.MgADP=0; p8.Pi=0;
[~,~,fm8,~] = runSlackExperiment(p8);
fprintf('8mM reference computed.\n');

f2d = load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8d = load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
targetFull   = buildAtpRatioTarget(f2d, f8d, 'full');
targetSparse = buildAtpRatioTarget(f2d, f8d, 'sparse');

reportBreakdown = @(fm2,tag) printBreak(fm2, fm8, targetFull, tag);

%% ---- baseline (current best K_T2=2.0, K_srx=7.0) ----
[~,~,fm2_0,~] = runSlackExperiment(base0);
c0f = atpRatioCost(fm2_0, fm8, targetFull);
c0s = atpRatioCost(fm2_0, fm8, targetSparse);
fprintf('\nBASELINE  K_T2=%.3f K_srx=%.3f :  enriched=%.4f   sparse(legacy)=%.4f\n', base0.K_T2, base0.K_srx, c0f, c0s);
reportBreakdown(fm2_0, 'baseline');

% ---- optimize (K_T2, K_srx) in log-space against the enriched cost ----
objfun = @(x) evalObj(x, base0, fm8, targetFull);
x0 = log([base0.K_T2, base0.K_srx]);
opts = optimset('MaxFunEvals',140,'MaxIter',90,'TolFun',5e-5,'TolX',5e-4,'Display','iter');
[xopt, coptf] = fminsearch(objfun, x0, opts);
K_T2o = exp(xopt(1)); K_srxo = exp(xopt(2));

p_opt = base0; p_opt.K_T2=K_T2o; p_opt.K_srx=K_srxo;
[~,~,fm2_o,~] = runSlackExperiment(p_opt);
copts = atpRatioCost(fm2_o, fm8, targetSparse);
fprintf('\n=== OPTIMUM ===  K_T2=%.3f  K_srx=%.3f :  enriched=%.4f (was %.4f)   sparse=%.4f (was %.4f)\n', ...
  K_T2o, K_srxo, coptf, c0f, copts, c0s);
reportBreakdown(fm2_o, 'optimum');

% ---- save optimized snapshot ----
dst = fullfile(root,'params','params_2state_optim.m');
copyfile(fullfile(root,'params','params_2state_final.m'), dst);
fid = fopen(dst,'a');
fprintf(fid,'\r\n%%%% enriched-ratio fminsearch 2026-07-15: K_T2=%.4f K_srx=%.4f (enriched %.4f, sparse %.4f)\r\n', K_T2o, K_srxo, coptf, copts);
fprintf(fid,'params0.K_T2 = %.4f;\r\nparams0.K_srx = %.4f;\r\n', K_T2o, K_srxo);
fclose(fid);
fprintf('Saved %s\n', dst);

% ---- ratio-space timecourse figure: model vs data per feature ----
resdir = fullfile(root,'Analyses','LowATP_kmsrd','results');
vf = fieldnames(targetFull.vector);
fig = figure('Position',[40 40 1500 850],'Color','w');
tl = tiledlayout(fig, 3, 4, 'TileSpacing','compact','Padding','compact');
for i=1:numel(vf); f=vf{i}; nexttile(tl); hold on; box on;
  rd = double(targetFull.vector.(f));
  r0 = double(fm2_0.(f))./double(fm8.(f));
  ro = double(fm2_o.(f))./double(fm8.(f));
  plot(1:numel(rd), rd, 'ko-','LineWidth',1.4,'DisplayName','data ratio');
  plot(1:numel(r0), r0, 'b.--','LineWidth',1.0,'DisplayName','model (base)');
  plot(1:numel(ro), ro, 'rx-','LineWidth',1.2,'DisplayName','model (opt)');
  title(sprintf('%s  (w=%g)', f, targetFull.weights.(f)),'Interpreter','none'); xlabel('segment');
  if i==1; legend('Location','best','FontSize',7); end
end
sgtitle(sprintf('2mM/8mM ratio per feature: data vs model  (enriched cost %.4f -> %.4f)', c0f, coptf));
exportgraphics(fig, fullfile(resdir,'OPTIM_ratio_by_feature.png'),'Resolution',140);
save(fullfile(resdir,'OPTIM_result.mat'),'fm2_o','fm8','targetFull','K_T2o','K_srxo','coptf','copts');
disp('OPTIM DONE');

% ================= helpers =================
function cost = evalObj(x, base0, fm8, target)
  K_T2 = exp(x(1)); K_srx = exp(x(2));
  if K_T2<0.3 || K_T2>15 || K_srx<1 || K_srx>40; cost = 1e6; return; end
  p = base0; p.K_T2 = K_T2; p.K_srx = K_srx;
  try
    [~,~,fm2,~] = runSlackExperiment(p);
    cost = atpRatioCost(fm2, fm8, target);
    if ~isfinite(cost); cost = 1e6; end
  catch
    cost = 1e6;
  end
end
function printBreak(fm2, fm8, target, tag)
  [~,d] = atpRatioCost(fm2, fm8, target);
  vf = fieldnames(target.vector);
  fprintf('  --- per-feature cost (%s) ---\n', tag);
  rows = cell(numel(vf),2);
  for i=1:numel(vf); f=vf{i}; rows{i,1}=f; rows{i,2}=d.(f).cost; end
  [~,ord] = sort(cell2mat(rows(:,2)),'descend');
  for k=1:numel(ord); i=ord(k); f=rows{i,1};
    fprintf('   %-20s cost %7.4f   modelR %s\n', f, rows{i,2}, num2str(d.(f).ratio(:)','%.2f '));
  end
end
