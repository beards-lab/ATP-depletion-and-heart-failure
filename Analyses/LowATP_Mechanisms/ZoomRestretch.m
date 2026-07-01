%% ZoomRestretch.m
% Zoomed model-vs-data force over single slack segments, to see the
% post-restretch peak (peak2 >> steady) that has no counterpart in the data.
% Out: results/zoom_restretch.png
clear; clc;
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..'); addpath(genpath(root)); refreshPool(5);
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

p0=getParams(loadParams('params_reseeded_regavail_opt2'),[],true,false);
p0.RunForceVelocity=false;p0.RunKtr=false;p0.RunStairs=false;p0.RunForceVelocityTime=false;p0.RunForceLengthEstim=false;
p0.RunSlack=true;p0.EvalFeatures=true;p0.recalculateDataFeats=false;
p0.PlotEachSeparately=0;p0.PlotProbsOnFig=0;p0.PlotFeatureFitting=false;p0.justPlotStateTransitionsFlag=false;
p0.BreakOnODEUnstable=false;p0.MaxRunTime=120;p0.ghostLoad='';p0.ghostSave='';p0.EvalPeaks=false;p0.EvalFitSlackOnset=false;
p0.RunSlackSegments='All';   % serial -> contiguous trace
rb=p0; rb.NumberOfStates=3; rb.UseKstiff3=0; rb.kstiff3=p0.kstiff2; rb.drp3=p0.dr; rb.UseAtpDetach=1; rb.UseAdpTrap=1; rb.UsePiReversal=1; rb.k3=1200;

p8=rb;p8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';p8.MgATP=8;p8.MgADP=0;p8.Pi=0;
[~,o8,~,~]=runSlackExperiment(getParams(p8,[],true,false));
p2=rb;p2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';p2.MgATP=2;p2.MgADP=1.5;p2.Pi=3;
[~,o2,~,~]=runSlackExperiment(getParams(p2,[],true,false));

wins=[74.45 75.05; 76.84 77.55];   % segment 1 (fast restretch) and segment 5 (slow)
fig=figure('Position',[100 100 1050 520]); tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
O={o8,o2}; lab={'8 mM baseline','2 mM coupled'};
for r=1:2
  o=O{r}; dt=o.datatable;
  for c=1:2
    nexttile; hold on; box on;
    plot(dt(:,1), dt(:,3), 'k-','LineWidth',1.2);
    plot(o.t, o.Force, 'b-','LineWidth',1.4);
    xlim(wins(c,:)); ylabel('Force (kPa)');
    if r==1 && c==1; legend({'data','model'},'Location','northwest','FontSize',8); end
    title(sprintf('%s — segment %d', lab{r}, (c==1)*1+(c==2)*5));
  end
end
title(tl,'Restretch zoom: model peak2 overshoots steady (post-restretch peak); data peak2 \approx steady');
exportgraphics(fig, fullfile(resdir,'zoom_restretch.png'),'Resolution',130);
fprintf('Saved results/zoom_restretch.png\n');
