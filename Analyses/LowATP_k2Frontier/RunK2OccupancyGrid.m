%% RunK2OccupancyGrid.m
% =====================================================================
% SECOND LEVER: pair k2 reduction (the ATP detachment axis) with the global
% occupancy-saturation cap P_bound_max (the isometric-force ceiling), and ask:
% does any (k2, P_bound_max) pair land FORCE at +18% AND ktr at x0.54 together?
%
% Occupancy saturation (dPUdT_CombinedTransitions, UseGlobalOccupancySaturation):
%   attachment flux *= max(0, 1 - P_bound / P_bound_max)   [linear form]
% Current best leaves it ~off (P_bound_max=1.024). Lowering it caps the attached
% pool, biting hardest in the high-attachment low-k2 (2 mM) condition.
%
% Both knobs are STRUCTURAL/kinetic and applied to BOTH conditions; the ATP
% contrast is k2 (8 mM = base k2, 2 mM = k2*scale). Relative-ratio scoring.
%
% Run:  run Analyses/LowATP_k2Frontier/RunK2OccupancyGrid
% Out:  results/k2_occupancy_grid.mat + results/k2_occupancy_landing.png
% =====================================================================
clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
addpath(genpath(root));
resdir = fullfile(here, 'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

params0 = getParams();
run(fullfile(root, 'params', 'params_reseeded_regavail_opt2.m'));
params0 = getParams(params0, [], true, false);
params0.RunForceVelocity=false; params0.RunKtr=false; params0.RunStairs=false;
params0.RunForceVelocityTime=false; params0.RunForceLengthEstim=false;
params0.RunSlack=true; params0.EvalFeatures=true; params0.recalculateDataFeats=false;
params0.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
params0.PlotEachSeparately=0; params0.PlotProbsOnFig=0; params0.PlotFeatureFitting=false;
params0.justPlotStateTransitionsFlag=false; params0.BreakOnODEUnstable=false; params0.MaxRunTime=120;
params0.ghostLoad=''; params0.ghostSave=''; params0.EvalPeaks=false; params0.EvalFitSlackOnset=false;
try; if isempty(gcp('nocreate')); parpool('Threads',5); end; params0.RunSlackSegments='AllPar';
catch; params0.RunSlackSegments='All'; end
baseK2 = params0.k2;

% data targets (2 mM / 8 mM)
f2 = load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8 = load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
core = {'steady','A','peak2','peak1_y','ktr','t0'};
dR = zeros(numel(core),1);
for i=1:numel(core); dR(i)=mean(double(f2.(core{i})),'omitnan')/mean(double(f8.(core{i})),'omitnan'); end

%% 2D grid
caps    = [1.024 0.30 0.22];     % P_bound_max (1.024 = effectively off)
scales  = [1.0 0.50 0.40 0.30];  % k2 multiplier (1.0 = 8 mM reference per cap)
nC=numel(caps); nS=numel(scales); nf=numel(core);
V = nan(nC,nS,nf);
fprintf('2D grid (caps x k2scales = %dx%d runs)...\n',nC,nS); tic
for ic=1:nC
  for is=1:nS
    p=params0; p.P_bound_max=caps(ic); p.k2=baseK2*scales(is);
    p=getParams(p,[],true,false);
    [~,~,fm,~]=runSlackExperiment(p);
    for k=1:nf; V(ic,is,k)=mean(double(fm.(core{k})),'omitnan'); end
  end
  fprintf('  cap=%.3f done (%.0fs)\n',caps(ic),toc);
end

%% ratios vs per-cap 8 mM reference (scale=1) + cost
R = nan(nC,nS,nf); cost = nan(nC,nS);
for ic=1:nC
  for is=1:nS
    r = squeeze(V(ic,is,:))./squeeze(V(ic,1,:));
    R(ic,is,:)=r; cost(ic,is)=sum((log(r)-log(dR)).^2);
  end
end
% best over the 2 mM candidates (scale<1)
cc=cost; cc(:,1)=Inf; [cmin,idx]=min(cc(:)); [ibc,ibs]=ind2sub(size(cc),idx);
fprintf('\nCORE-COST surface (rows=cap, cols=k2scale; col1 is the 8mM ref):\n');
fprintf('         k2x:  %s\n', num2str(scales,'%6.2f '));
for ic=1:nC; fprintf('cap=%.3f      %s\n',caps(ic),num2str(cost(ic,:),'%6.3f ')); end
fprintf('\nBEST 2-knob: P_bound_max=%.3f, k2x%.2f -> coreCost=%.3f  (k2-alone best ~0.088)\n',caps(ibc),scales(ibs),cmin);
rb=squeeze(R(ibc,ibs,:));
fprintf('  model/data: steady %.2f/%.2f  A %.2f/%.2f  peak2 %.2f/%.2f  peak1 %.2f/%.2f  ktr %.2f/%.2f  t0 %.2f/%.2f\n',...
  rb(1),dR(1),rb(2),dR(2),rb(3),dR(3),rb(4),dR(4),rb(5),dR(5),rb(6),dR(6));

%% money figure: (steady ratio, ktr ratio) landing, target = (1.18, 0.54)
is_=strcmp(core,'steady'); ik=strcmp(core,'ktr');
fig=figure('Position',[100 100 720 560]); hold on; box on;
cols=lines(nC);
for ic=1:nC
  sx=squeeze(R(ic,:,is_)); ky=squeeze(R(ic,:,ik));
  plot(sx,ky,'o-','Color',cols(ic,:),'LineWidth',1.8,'MarkerFaceColor',cols(ic,:),...
       'DisplayName',sprintf('P_{bound max}=%.3f',caps(ic)));
  for is=2:nS; text(sx(is)+0.005,ky(is),sprintf('%.2f',scales(is)),'FontSize',7,'Color',cols(ic,:)); end
end
plot(dR(is_),dR(ik),'p','MarkerSize',20,'MarkerFaceColor',[0.95 0.8 0.1],'MarkerEdgeColor','k',...
     'LineWidth',1.2,'DisplayName','DATA target (2mM/8mM)');
xline(dR(is_),':k'); yline(dR(ik),':k');
xlabel('steady-force ratio  (model 2mM / 8mM)'); ylabel('ktr ratio  (model 2mM / 8mM)');
title({'Two-lever landing: k2 (labels) x occupancy cap (colors)','target star = +18% force, x0.54 ktr'});
legend('Location','northeast','FontSize',8);
exportgraphics(fig, fullfile(resdir,'k2_occupancy_landing.png'),'Resolution',130);

save(fullfile(resdir,'k2_occupancy_grid.mat'),'caps','scales','V','R','cost','core','dR','baseK2','ibc','ibs');
fprintf('\nSaved results/k2_occupancy_grid.mat + k2_occupancy_landing.png\n');
