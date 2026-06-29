%% RunThirdState.m
% =====================================================================
% Phase 2 result: does a 3rd cross-bridge state break the 2-state wall?
%
% 2-state wall (Phase 1): cost 0.70; residuals = restretch transients (vall_y,
% peak2) + ktr. The 3-state adds a low-force, slow-detaching post-stroke state
% p3 (p1->p2->p3->PT). ATP lever:
%   k2 down  -> more p2 (high-force isometric)  -> steady up
%   k3 down  -> p3 accumulates (force only when STRAINED during restretch)
%               -> peak2/vall_y up + cycle slows (ktr down)
% kstiff3 lowered to 7000 so p3 is genuinely low-force (deeper k3 without vall
% overshoot). p3 is transient at 8 mM (fast k3=600) so the baseline is sane.
%
% Relative-ratio scoring vs per-segment target (../LowATP_k2Frontier/results/atp_target.mat).
% Out: results/thirdstate_compare.png + results/thirdstate.mat
% NB: the regavail param file assigns to `params0` internally; loaded here via
%     loadParams() (Auxiliary/loadParams.m), which isolates that assignment in
%     its own scope so it cannot silently clobber/be clobbered by this script's
%     own params0.
% =====================================================================
clear; clc;
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..'); addpath(genpath(root));
refreshPool(5);
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end
target=load(fullfile(root,'Analyses','LowATP_k2Frontier','results','atp_target.mat')).target;

params0 = getParams(loadParams('params_reseeded_regavail_opt2'), [], true, false);
params0.RunForceVelocity=false;params0.RunKtr=false;params0.RunStairs=false;
params0.RunForceVelocityTime=false;params0.RunForceLengthEstim=false;params0.RunSlack=true;
params0.EvalFeatures=true;params0.recalculateDataFeats=false;
params0.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
params0.PlotEachSeparately=0;params0.PlotProbsOnFig=0;params0.PlotFeatureFitting=false;
params0.justPlotStateTransitionsFlag=false;params0.BreakOnODEUnstable=false;params0.MaxRunTime=120;
params0.ghostLoad='';params0.ghostSave='';params0.EvalPeaks=false;params0.EvalFitSlackOnset=false;
try;if isempty(gcp('nocreate'));parpool('Threads',5);end;params0.RunSlackSegments='AllPar';catch;params0.RunSlackSegments='All';end
bK2=params0.k2;
runfeat=@(p) getFeat(getParams(p,[],true,false));

%% --- 2-STATE wall best (k2x0.62, ka x1.0) ---
fm8_2 = runfeat(params0);
p2m=params0; p2m.k2=bK2*0.62;            fm2_2 = runfeat(p2m);
[c2,d2]=atpRatioCost(fm2_2,fm8_2,target);

%% --- 3-STATE best (kstiff3=7000, k3=600 baseline; k2x0.65 k3x0.50 lever) ---
p3b=params0; p3b.NumberOfStates=3; p3b.k3=600; p3b.kstiff3=7000;
fm8_3 = runfeat(p3b);
p3m=p3b; p3m.k2=bK2*0.65; p3m.k3=600*0.50;
fm2_3 = runfeat(p3m);
[c3,d3]=atpRatioCost(fm2_3,fm8_3,target);

fprintf('\n2-STATE wall cost = %.3f ;  3-STATE cost = %.3f\n',c2,c3);

%% --- comparison figure ---
feats={'steady','A','Am','peak2','peak1_y','vall_y','ktr','t0'};
getr=@(d,f) mean(d.(f).ratio);
tgt=zeros(1,numel(feats)); r2=tgt; r3=tgt;
for i=1:numel(feats); f=feats{i};
    tgt(i)=mean(d3.(f).target); r2(i)=getr(d2,f); r3(i)=getr(d3,f);
end
fig=figure('Position',[100 100 920 480]);
b=bar([tgt;r2;r3]','grouped'); hold on; box on;
b(1).FaceColor=[0.95 0.8 0.1]; b(2).FaceColor=[0.85 0.4 0.3]; b(3).FaceColor=[0.15 0.5 0.8];
set(gca,'XTickLabel',feats,'XTickLabelRotation',20);
ylabel('ratio (2 mM / 8 mM)'); yline(1,'k:');
legend({sprintf('DATA target'),sprintf('2-state best (cost %.2f)',c2),sprintf('3-state best (cost %.2f)',c3)},'Location','northwest');
title('3rd state fixes the restretch-transient wall (peak2, vall\_y); residual = ktr/t0 kinetics');
exportgraphics(fig,fullfile(resdir,'thirdstate_compare.png'),'Resolution',130);

cfg3=struct('NumberOfStates',3,'k3_base',600,'kstiff3',7000,'k2_scale',0.65,'k3_scale',0.50);
save(fullfile(resdir,'thirdstate.mat'),'c2','c3','d2','d3','feats','tgt','r2','r3','cfg3');
fprintf('Saved results/thirdstate_compare.png + thirdstate.mat\n');

function fm=getFeat(p); [~,~,fm,~]=runSlackExperiment(p); end
