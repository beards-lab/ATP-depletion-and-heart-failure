%% PlotRawOutputs.m
% Raw force-time slack traces: model vs experimental data, for the 8 mM baseline
% and the 2 mM COUPLED metabolic state (rigor + ADP-trap + Pi). Shows the actual
% fit (and where the kinetic-timing residual lives) rather than summary ratios.
% Out: results/raw_trace_8mM.png, results/raw_trace_2mM_coupled.png
clear; clc;
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..'); addpath(genpath(root));
refreshPool(5);
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

base = getParams(loadParams('params_reseeded_regavail_opt2'), [], true, false);
base.RunForceVelocity=false;base.RunKtr=false;base.RunStairs=false;base.RunForceVelocityTime=false;base.RunForceLengthEstim=false;
base.RunSlack=true;base.EvalFeatures=true;base.recalculateDataFeats=false;
base.PlotProbsOnFig=0;base.justPlotStateTransitionsFlag=false;base.BreakOnODEUnstable=false;base.MaxRunTime=120;
base.ghostLoad='';base.ghostSave='';base.EvalPeaks=false;base.EvalFitSlackOnset=false;
base.RunSlackSegments='All';            % serial -> clean contiguous trace for plotting
base.PlotEachSeparately=1; base.ShowStatePlots=0; base.ShowResidualPlots=0; base.drawForceOnset=0; base.PlotFeatureFitting=0;

% rigor + coupled structure (physiological: kstiff3 = kstiff2, drp3 = dr)
rb = base; rb.NumberOfStates=3; rb.UseKstiff3=0; rb.kstiff3=base.kstiff2; rb.drp3=base.dr;
rb.UseAtpDetach=1; rb.k3=1200; rb.UseAdpTrap=1; rb.UsePiReversal=1;

% --- 8 mM baseline (MgADP=0, Pi=0) vs 8 mM data ---
p8=rb; p8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat'; p8.MgATP=8;p8.MgADP=0;p8.Pi=0;
figure(101);clf;tiledlayout('flow','TileSpacing','compact','Padding','compact');
runSlackExperiment(getParams(p8,[],true,false));
sgtitle('8 mM: model (blue) vs data (black) — slack-restretch force');
exportgraphics(figure(101),fullfile(resdir,'raw_trace_8mM.png'),'Resolution',130);

% --- 2 mM COUPLED (MgATP=2, MgADP=1.5, Pi=3) vs 2 mM data ---
p2=rb; p2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat'; p2.MgATP=2;p2.MgADP=1.5;p2.Pi=3;
figure(102);clf;tiledlayout('flow','TileSpacing','compact','Padding','compact');
runSlackExperiment(getParams(p2,[],true,false));
sgtitle('2 mM coupled (vATP+^ADP+^Pi): model (blue) vs data (black)');
exportgraphics(figure(102),fullfile(resdir,'raw_trace_2mM_coupled.png'),'Resolution',130);
fprintf('Saved raw_trace_8mM.png + raw_trace_2mM_coupled.png\n');
