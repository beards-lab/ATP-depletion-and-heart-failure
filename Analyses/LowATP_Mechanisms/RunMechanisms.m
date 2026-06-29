%% RunMechanisms.m
% =====================================================================
% Test physiologically-justifiable nucleotide mechanisms for the low-ATP
% (2 mM vs 8 mM) cross-bridge signature. Relative-ratio scoring vs the
% per-segment target (../LowATP_k2Frontier/results/atp_target.mat).
%
% Low ATP at the myofibril = simultaneous  vATP + ^ADP + ^Pi  (coupled).
% Driven by CONCENTRATIONS (MgATP 8->2; MgADP, Pi fit). Mechanisms (all wired
% into the active pchip path of dPUdT_CombinedTransitions, baseline-preserving
% at MgADP=0,Pi=0; flags in getParams):
%   ADP-trap  (UseAdpTrap)   : ^ADP -> R2 *= g2 -> traps force-bearing P2 (force^, ktr v)
%   Pi-revers (UsePiReversal): ^Pi  -> R12 *= f2, R21 *= (1+Pi/K_Pi) -> force v (Cooke&Pate)
%   Rigor     (UseAtpDetach) : vATP -> R3 *= MgATP/(MgATP+K_T_detach) -> rigor P3 piles up
%                              P3 = STIFF rigor: kstiff3 = kstiff2 (physiological >=), drp3 = dr.
% NB physiology (user): NO ad-hoc `ka` lever, and kstiff3 NOT < kstiff2.
%
% Out: results/mechanisms.mat + results/mechanisms_compare.png
% =====================================================================
clear; clc;
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..'); addpath(genpath(root));
refreshPool(5);
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end
target=load(fullfile(root,'Analyses','LowATP_k2Frontier','results','atp_target.mat')).target;

p0 = getParams(loadParams('params_reseeded_regavail_opt2'), [], true, false);
p0.RunForceVelocity=false;p0.RunKtr=false;p0.RunStairs=false;p0.RunForceVelocityTime=false;p0.RunForceLengthEstim=false;
p0.RunSlack=true;p0.EvalFeatures=true;p0.recalculateDataFeats=false;
p0.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
p0.PlotEachSeparately=0;p0.PlotProbsOnFig=0;p0.PlotFeatureFitting=false;p0.justPlotStateTransitionsFlag=false;
p0.BreakOnODEUnstable=false;p0.MaxRunTime=120;p0.ghostLoad='';p0.ghostSave='';p0.EvalPeaks=false;p0.EvalFitSlackOnset=false;
p0.RunSlackSegments='AllPar';
run8 = @(p) getFeat(setfield(getParams(p,[],true,false),'MgATP',8));   %#ok<*SFLD>
runC = @(p) getFeat(getParams(p,[],true,false));

% rigor 3-state structural base (physiological stiffness)
rb = p0; rb.NumberOfStates=3; rb.UseKstiff3=0; rb.kstiff3=p0.kstiff2; rb.drp3=p0.dr; rb.UseAtpDetach=1; rb.k3=1200;

% mechanism table: {name, base8mM struct, perturb2mM fn(base)->struct}
M = {};
M(end+1,:) = {'ADP-trap (MgADP=2)',          setflags(p0,'UseAdpTrap',1),  @(b) setfields(b,'MgADP',2)};
M(end+1,:) = {'Pi only (Pi=4)',              setflags(p0,'UsePiReversal',1), @(b) setfields(b,'Pi',4)};
M(end+1,:) = {'Rigor only (ATP-detach)',     rb,                            @(b) setfields(b,'MgATP',2)};
M(end+1,:) = {'COUPLED (ADP+Pi+rigor)',      setflags(rb,'UseAdpTrap',1,'UsePiReversal',1), @(b) setfields(b,'MgATP',2,'MgADP',1.5,'Pi',3)};

feats={'steady','A','Am','peak2','peak1_y','vall_y','ktr','t0'};
nM=size(M,1); R=nan(nM,numel(feats)); costs=nan(nM,1); tgt=zeros(1,numel(feats));
fprintf('\n%-26s %6s | key ratios\n','mechanism','cost');
for m=1:nM
  b8 = M{m,2}; b8.MgATP=8; b8.MgADP=0; b8.Pi=0;
  fm8 = runC(b8);
  fm2 = runC(M{m,3}(b8));
  [c,d]=atpRatioCost(fm2,fm8,target); costs(m)=c;
  for i=1:numel(feats); R(m,i)=mean(d.(feats{i}).ratio); tgt(i)=mean(d.(feats{i}).target); end
  fprintf('%-26s %6.3f | st=%.2f pk2=%.2f vall=%.2f ktr=%.2f t0=%.2f\n',M{m,1},c,R(m,strcmp(feats,'steady')),R(m,strcmp(feats,'peak2')),R(m,strcmp(feats,'vall_y')),R(m,strcmp(feats,'ktr')),R(m,strcmp(feats,'t0')));
end
fprintf('\n(refs: 2-state wall 0.70; unphysical wobbly-p3+ka 0.39)\n');

%% figure: per-feature ratios, data vs each mechanism
fig=figure('Position',[100 100 1000 520]); b=bar([tgt;R]','grouped'); hold on; box on;
cols=[0.95 0.8 0.1; 0.6 0.6 0.6; 0.4 0.7 0.9; 0.9 0.6 0.3; 0.15 0.5 0.2];
for k=1:numel(b); b(k).FaceColor=cols(min(k,size(cols,1)),:); end
set(gca,'XTickLabel',feats,'XTickLabelRotation',20); ylabel('ratio (2 mM / 8 mM)'); yline(1,'k:');
legend([{'DATA target'}, M(:,1)'],'Location','northwest','FontSize',8);
title('Physiological low-ATP mechanisms vs data: COUPLED nails force+restretch; residual = t0/peak1/ktr kinetics');
exportgraphics(fig,fullfile(resdir,'mechanisms_compare.png'),'Resolution',130);
save(fullfile(resdir,'mechanisms.mat'),'M','R','tgt','costs','feats');
fprintf('Saved results/mechanisms_compare.png + mechanisms.mat\n');

function s=setflags(s,varargin); for i=1:2:numel(varargin); s.(varargin{i})=varargin{i+1}; end; end
function s=setfields(s,varargin); for i=1:2:numel(varargin); s.(varargin{i})=varargin{i+1}; end; end
function fm=getFeat(p); [~,~,fm,~]=runSlackExperiment(p); end
