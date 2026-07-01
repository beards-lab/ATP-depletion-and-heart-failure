% DriverLowATP_Coupled.m
% =====================================================================
% Run the coupled LOW-ATP mechanism and report how it reproduces the measured
% 2 mM-vs-8 mM cross-bridge signature (relative-ratio scoring).
%
% NB there is an older Workbench/DriverLowATP.m — a DIFFERENT approach (the
% two-Km UseAtpK2/UseAtpKah FV-based scheme). This driver is the current coupled
% rigor + ADP-trap + Pi mechanism (Analyses/LowATP_Mechanisms).
%
%   PARAM FILE : params/ModelParams_coupled_lowATP.m
%                (the coupled structure on the regavail_opt2 baseline; this file
%                IS the 8 mM baseline. Relative-ratio cost ~0.56.)
%
% Mechanism (all physiologically wired into the active path; low ATP at the
% myofibril = simultaneous  vATP + ^ADP + ^Pi):
%   ADP-trap  (UseAdpTrap)    : ^ADP traps force-bearing P2   -> isometric force ^
%   rigor     (UseAtpDetach)  : vATP slows P3 detachment      -> restretch transients ^
%   Pi        (UsePiReversal) : ^Pi inhibits the power stroke -> tempers the force gain
% The low-ATP condition is set purely by CONCENTRATIONS (below).
%
% Run:  cd(root); addpath(genpath('.')); DriverLowATP_Coupled
% See also: DriverHighATP, Analyses/LowATP_Mechanisms/{REPORT,SYNTHESIS}.md
% =====================================================================
% clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root));
refreshPool(5);

PARAM_FILE = 'ModelParams_coupled_lowATP';        % <-- the coupled low-ATP config
MgATP_low = 2;  MgADP_low = 1.5;  Pi_low = 3;      % the 2 mM condition (fit concentrations)

base = getParams(loadParams(PARAM_FILE), [], true, false);
base.RunForceVelocity=true; base.RunKtr=false; base.RunStairs=false;
base.RunForceVelocityTime=false; base.RunForceLengthEstim=false; base.RunSlack=true;
base.EvalFeatures=true; base.recalculateDataFeats=false;
base.PlotEachSeparately=true; base.PlotProbsOnFig=0; base.PlotFeatureFitting=false; base.justPlotStateTransitionsFlag=false;
base.BreakOnODEUnstable=false; base.MaxRunTime=120; base.ghostLoad=''; base.ghostSave='';
base.EvalPeaks=false; base.EvalFitSlackOnset=false; base.RunSlackSegments='AllPar';

%% --- 8 mM baseline (MgADP=0, Pi=0) and 2 mM low-ATP (the concentrations above) ---
p8 = base; p8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat'; p8.MgATP=8; p8.MgADP=0; p8.Pi=0;
params0 = p8;
% params0.kSE = 3.2e3;
% [~,~,fm8,~] = runSlackExperiment(getParams(p8,[],true,false));
figure(8);clf;
RunBakersExp;
fm8 = features_model;
fd8 = features_data;

%% --- 8 mM baseline (MgADP=0, Pi=0) and 2 mM low-ATP (the concentrations above) ---
p2 = base; p2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat'; p2.MgATP=MgATP_low; p2.MgADP=MgADP_low; p2.Pi=Pi_low;
params0 = p2;
% [~,~,fm8,~] = runSlackExperiment(getParams(p8,[],true,false));
RunBakersExp;
fm2 = features_model;
fd2 = features_data;

% p2 = base; p2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat'; p2.MgATP=MgATP_low; p2.MgADP=MgADP_low; p2.Pi=Pi_low;
% [~,~,fm2,~] = runSlackExperiment(getParams(p2,[],true,false));

%%
% fg8 = fm8;
if ~exist('fg8', 'var')
    fg8 = fm8;
end
sum(plotFeatures(fd8, fm8, fg8, base.fn))
%%

plotMultipleFeatures({fd8, fm8, fd2, fm2}, {'Data8', 'model8', 'Data2', 'Model2'}, lines(4), {'x', 'o', 'x', 'o'}, params0.fn)

%% --- relative-ratio score vs the per-segment data target ---
target = load(fullfile(root,'Analyses','LowATP_k2Frontier','results','atp_target.mat')).target;
[c, d] = atpRatioCost(fm2, fm8, target);
fprintf('Coupled low-ATP = params/%s.m\n', PARAM_FILE);
fprintf('  MgATP %g, MgADP %g, Pi %g  ->  relative-ratio cost = %.3f\n', MgATP_low, MgADP_low, Pi_low, c);
fprintf('\n  %-14s %10s %10s\n','feature','model r','data r');
fn = fieldnames(d);
for k = 1:numel(fn)
    dd = d.(fn{k});
    fprintf('  %-14s %10.2f %10.2f\n', fn{k}, mean(dd.ratio), mean(dd.target));
end
fprintf(['\n  Force amplitudes (steady/A/peak2/vall) reproduced; residuals are kinetic:\n' ...
         '  ktr level+slope and t0 (see SYNTHESIS.md - the force<->ktr iron law).\n']);

% --- figure: per-feature 2mM/8mM ratios, model vs data ---
mr = cellfun(@(f) mean(d.(f).ratio),  fn);
tg = cellfun(@(f) mean(d.(f).target), fn);
figure(3); clf; bar([tg(:) mr(:)]); box on;
set(gca,'XTickLabel',fn,'XTickLabelRotation',25); ylabel('ratio (2 mM / 8 mM)'); yline(1,'k:');
legend('data','model','Location','northwest');
title(sprintf('Coupled low-ATP mechanism: 2mM/8mM ratios (cost %.3f)', c));
