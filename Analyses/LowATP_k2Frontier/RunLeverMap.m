%% RunLeverMap.m
% =====================================================================
% Capstone: map every candidate low-ATP lever in the (steady-force ratio,
% ktr ratio) plane, relative to the 8 mM baseline, with peak2 ratio annotated.
% Shows WHY no single 2nd lever added to k2 reaches the data target
% (steady x1.18, ktr x0.54, peak2 x1.18):
%   * k2 down     -> down-RIGHT (force overshoots)
%   * occupancy   -> LEFT-UP  (trims force, speeds ktr)
%   * Pi block(k_1)-> LEFT-UP  (trims force, speeds ktr, inflates peak2)
%   * ka down     -> LEFT      (trims force & holds ktr best, but peak2 CRASHES)
%   * SRX         -> ~null at max-Ca (see conclusions; not plotted)
%
% Out: results/lever_map.png + results/lever_map.mat
% =====================================================================
clear; clc;
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..'); addpath(genpath(root));
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

params0=getParams(); run(fullfile(root,'params','params_reseeded_regavail_opt2.m'));
params0=getParams(params0,[],true,false);
params0.RunForceVelocity=false;params0.RunKtr=false;params0.RunStairs=false;
params0.RunForceVelocityTime=false;params0.RunForceLengthEstim=false;params0.RunSlack=true;
params0.EvalFeatures=true;params0.recalculateDataFeats=false;
params0.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
params0.PlotEachSeparately=0;params0.PlotProbsOnFig=0;params0.PlotFeatureFitting=false;
params0.justPlotStateTransitionsFlag=false;params0.BreakOnODEUnstable=false;params0.MaxRunTime=120;
params0.ghostLoad='';params0.ghostSave='';params0.EvalPeaks=false;params0.EvalFitSlackOnset=false;
try;if isempty(gcp('nocreate'));parpool('Threads',5);end;params0.RunSlackSegments='AllPar';catch;params0.RunSlackSegments='All';end
bK2=params0.k2; bka=params0.ka; bk_1=params0.k_1; bPb=params0.P_bound_max;

% data targets
f2=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
tg=@(f) mean(double(f2.(f)),'omitnan')/mean(double(f8.(f)),'omitnan');
T=[tg('steady') tg('ktr') tg('peak2')];

% each row: {label, fieldoverrides as struct}
runs = {
 'REF',        struct();
 'k2x0.5',     struct('k2',bK2*0.5);
 'k2x0.4',     struct('k2',bK2*0.4);
 'k2x0.3',     struct('k2',bK2*0.3);
 'k2.4+occ.30',struct('k2',bK2*0.4,'P_bound_max',0.30);
 'k2.4+occ.22',struct('k2',bK2*0.4,'P_bound_max',0.22);
 'k2.4+Pi(k_1x3)',struct('k2',bK2*0.4,'k_1',bk_1*3);
 'k2.4+Pi(k_1x6)',struct('k2',bK2*0.4,'k_1',bk_1*6);
 'k2.4+ka.7',  struct('k2',bK2*0.4,'ka',bka*0.7);
 'k2.4+ka.5',  struct('k2',bK2*0.4,'ka',bka*0.5);
};
n=size(runs,1); M=nan(n,3);  % [steady ktr peak2] absolute
tic;
for i=1:n
  p=params0; ov=runs{i,2}; fn=fieldnames(ov);
  for k=1:numel(fn); p.(fn{k})=ov.(fn{k}); end
  p=getParams(p,[],true,false);
  [~,~,fm,~]=runSlackExperiment(p);
  M(i,:)=[mean(fm.steady,'omitnan') mean(fm.ktr,'omitnan') mean(fm.peak2,'omitnan')];
  fprintf('%-16s steady=%.1f ktr=%.1f peak2=%.1f\n',runs{i,1},M(i,1),M(i,2),M(i,3));
end
fprintf('done %.0fs\n',toc);
R=M./M(1,:);  % ratios vs REF

%% figure: (steadyR, ktrR) with peak2R as label; arms colored by lever family
fam = [0 1 1 1 2 2 3 3 4 4];  % 0 ref,1 k2,2 occ,3 Pi,4 ka
cols=[0 0 0; 0.1 0.35 0.85; 0.9 0.4 0.05; 0.15 0.6 0.2; 0.55 0.2 0.65];
fnm={'ref','k2 (ATP axis)','+occupancy','+Pi block (k_1)','+ka'};
fig=figure('Position',[100 100 780 600]); hold on; box on;
for g=1:4
  idx=find(fam==g);
  if g==1; idx=[1 idx]; end  % anchor k2 arm at ref
  plot(R(idx,1),R(idx,2),'-o','Color',cols(g+1,:),'LineWidth',1.7,'MarkerFaceColor',cols(g+1,:),'DisplayName',fnm{g+1});
end
for i=2:n
  text(R(i,1)+0.006,R(i,2),sprintf('pk2=%.2f',R(i,3)),'FontSize',7.5,'Color',cols(fam(i)+1,:));
end
plot(T(1),T(2),'p','MarkerSize',22,'MarkerFaceColor',[0.95 0.8 0.1],'MarkerEdgeColor','k','LineWidth',1.2,'DisplayName',sprintf('DATA target (pk2=%.2f)',T(3)));
xline(T(1),':k'); yline(T(2),':k');
xlabel('steady-force ratio (2mM/8mM)'); ylabel('ktr ratio (2mM/8mM)');
title({'Low-ATP lever map: only k2 reaches the data force/ktr region,','but every force-trim breaks ktr (occ/Pi: up) or peak2 (ka: pk2<1)'});
legend('Location','northeast','FontSize',8);
exportgraphics(fig,fullfile(resdir,'lever_map.png'),'Resolution',130);
save(fullfile(resdir,'lever_map.mat'),'runs','M','R','T');
fprintf('Saved results/lever_map.png + lever_map.mat\n');
