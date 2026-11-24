% RunPowerGateOptim2.m — round 2: extend the 8 mM-inert gate set and report fully.
% Seeds from round 1 (results/powergate_state.mat) and adds a second tier of gates.
%
% TIER A (round 1, rate constants -- ATP dependence is mechanistically natural):
%   k2 kmsrd kah kamh ka kd k1 k_1 ksr0 ksr2srd ksrd2sr tau_reg A0 eta_M
% TIER B (added here, SHAPE levers -- weaker physiological claim, reported separately):
%   kstiff2 kSE kSE_M mu_neg dr dr2 P_bound_max
%   Rationale, honestly stated: at low ATP a larger fraction of attached heads sit in
%   long-lived post-stroke/rigor-like states, which are stiffer and differently strained,
%   so an EFFECTIVE change in per-head stiffness and mean stroke offset is arguable --
%   but these are not rate constants and the claim is weaker than for tier A. The run
%   reports how much tier B buys so the trade is explicit rather than hidden.
%
% Every gate remains (MgATP/MgATP_ref)^p => 8 mM BIT-IDENTICAL for any exponent.
% Run: matlab -batch "BUDGET_H=3; run('Analyses/LowATP_AbsoluteFit/RunPowerGateOptim2.m')"
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end
if isempty(gcp('nocreate')); parpool('local',6); end
if ~exist('BUDGET_H','var'); BUDGET_H = 3.0; end
T0=tic; BUDGET_S=BUDGET_H*3600;
STATE2=fullfile(resdir,'powergate2_state.mat');

TIERA = {'k2','kmsrd','kah','kamh','ka','kd','k1','k_1','ksr0','ksr2srd','ksrd2sr','tau_reg','A0','eta_M'};
TIERB = {'kstiff2','kSE','kSE_M','mu_neg','dr','dr2','P_bound_max'};
P = [TIERA TIERB];
LB = -1.5*ones(1,numel(P)); UB = 1.5*ones(1,numel(P));
% shape levers get a tighter box: they are not rate constants
ib = numel(TIERA)+1:numel(P); LB(ib) = -0.8; UB(ib) = 0.8;

b=getParams(loadParams('rskR2_w025_opt'),[],true,false);
b.RunForceVelocity=false; b.RunForceVelocityTime=false; b.RunKtr=false; b.RunStairs=false;
b.RunSlackPassive=false; b.RunSlack=true; b.EvalFeatures=true; b.BreakOnODEUnstable=false;
b.RunSlackSegments='AllPar'; b.MaxRunTime=300;
b.PlotEachSeparately=0; b.PlotFeatureFitting=0;
b.UseAtpDetachR2=0; b.UseAtpKmsrd=0; b.UseStrainDepAtpGate=0;
b.UseR2Ceiling=0; b.UseD0State=0; b.MgATP_ref=8; b.UseAtpPowerGates=1;
b.MgATP=2; b.MgADP=0; b.Pi=0.8;
b.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
fn=b.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
fn=cellfun(@(s) regexprep(s,'^rsK(\|.*)?$','rsK|0.1'), fn, 'UniformOutput',false);
FN=fn(~(startsWith(fn,'FV')|startsWith(fn,'PS_'))); b.fn=FN;

% ---- the 8 mM FLOOR on the identical feature list (what the gates cannot beat) ----
q8=b; q8.MgATP=8; q8.Pi=0; q8.UseAtpPowerGates=0;
q8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
[~,~,fm8,fd8]=runSlackExperiment(q8);
FLOOR=sum(evalFeatureCost(fd8,fm8,FN,2));
fprintf('8 mM FLOOR on the same %d-term list: %.4f  <- the 2 mM fit cannot do much better\n', numel(FN), FLOOR);

% ---- seed from round 1 ----
x0=zeros(1,numel(P));
if isfile(fullfile(resdir,'powergate_state.mat'))
    S=load(fullfile(resdir,'powergate_state.mat'));
    for i=1:numel(S.P); j=find(strcmp(P,S.P{i}),1); if ~isempty(j); x0(j)=S.best.x(i); end; end
    fprintf('seeded from round 1 (cost %.4f)\n', S.best.c);
end
best.x=x0; best.c=objPG2(x0,b,P,LB,UB); nEval=1;
fprintf('round-2 start cost: %.4f\n\n', best.c);
save(STATE2,'best','P','FLOOR','nEval');

STEPS={[-0.4 -0.2 -0.1 0.1 0.2 0.4],[-0.15 -0.07 0.07 0.15],[-0.05 0.05]};
for sweep=1:numel(STEPS)
    st=STEPS{sweep}; fprintf('--- sweep %d ---\n', sweep); anyimp=false;
    for i=1:numel(P)
        if toc(T0)>BUDGET_S; break; end
        bl=best.c; bv=best.x(i);
        for d=st
            xt=best.x; xt(i)=best.x(i)+d;
            if xt(i)<LB(i)||xt(i)>UB(i); continue; end
            c=objPG2(xt,b,P,LB,UB); nEval=nEval+1;
            if c<bl-1e-4; bl=c; bv=xt(i); end
        end
        if bl<best.c-1e-4
            best.x(i)=bv; best.c=bl; anyimp=true;
            fprintf('  [%5.1f min|%4d] %-12s -> %+.3f  cost %.4f\n', toc(T0)/60,nEval,P{i},bv,best.c);
            save(STATE2,'best','P','FLOOR','nEval');
        end
    end
    if toc(T0)>BUDGET_S||~anyimp; break; end
end
if toc(T0)<BUDGET_S
    o=optimset('Display','iter','MaxFunEvals',6000,'MaxIter',6000,'TolFun',1e-4,'TolX',1e-4, ...
               'OutputFcn',@(x,ov,st) toc(T0)>BUDGET_S);
    [xo,co]=fminsearch(@(x) objPG2(x,b,P,LB,UB), best.x, o);
    if co<best.c; best.x=xo; best.c=co; end
end
save(STATE2,'best','P','FLOOR','nEval');

fprintf('\n=== ROUND 2 BEST: %.4f   (8 mM floor %.4f) ===\n', best.c, FLOOR);
for i=1:numel(P)
    if abs(best.x(i))>1e-3
        tier=' A'; if i>numel(TIERA); tier=' B'; end
        fprintf('  [%s] pAtp_%-12s %+7.3f  -> x%.3f at 2 mM\n', tier, P{i}, best.x(i), (2/8)^best.x(i));
    end
end

% ---- report + figures ----
q=applyP2(b,P,best.x);
[~,out,fm,fd]=runSlackExperiment(q);
cc=evalFeatureCost(fd,fm,FN,2); [~,od]=sort(cc,'descend');
fprintf('\nper-feature (2 mM absolute):\n');
for k=1:min(12,numel(od)); if cc(od(k))<1e-3; break; end
    fprintf('   %-24s %8.3f\n', FN{od(k)}, cc(od(k))); end
sh={'steady','A','Am','peak1_y','vall_y','peak2','ktr','rsK','t0'};
g=@(x)num2str(double(x(:))','%7.3g');
fprintf('\n%-10s %-40s | %s\n','feature','model 2 mM','data 2 mM');
for i=1:numel(sh); f=sh{i}; if ~isfield(fm,f); continue; end
    fprintf('%-10s %-40s | %s\n', f, g(fm.(f)), g(fd.(f))); end

q8b=q; q8b.MgATP=8; q8b.Pi=0; q8b.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
r8=q8b; r8.UseAtpPowerGates=0;
[~,~,h1,~]=runSlackExperiment(q8b); [~,~,h0,~]=runSlackExperiment(r8);
fprintf('\n8 mM invariance: max|dsteady| %.3g | max|drsK| %.3g  (MUST be 0)\n', ...
  max(abs(double(h1.steady(:))-double(h0.steady(:)))), max(abs(double(h1.rsK(:))-double(h0.rsK(:)))));

try
  tmax=max(out.datatable(:,1));
  f1=figure('Position',[40 40 1500 520],'Color','w');
  plot(out.datatable(:,1),out.datatable(:,3),'k-','LineWidth',1.2,'DisplayName','data (2 mM)'); hold on; box on; grid on;
  plot(out.t,out.Force,'r-','LineWidth',1.6,'DisplayName','model (2 mM)');
  xlim([tmax-4.6 tmax+0.3]); xlabel('time (s)'); ylabel('Force (rel.)');
  title(sprintf('Low-ATP (2 mM) power-gate fit: cost %.2f (8 mM floor %.2f)', best.c, FLOOR));
  legend('Location','southeast');
  exportgraphics(f1,fullfile(resdir,'powergate_2mM_trace.png'),'Resolution',140);
  f2=figure('Position',[40 40 1400 800],'Color','w');
  plotFeatures(fd,fm,[],{'steady','A','Am','ktr','t0','peak1_y','vall_y','peak2','rsK'});
  sgtitle(sprintf('2 mM power-gate fit, absolute features (cost %.2f)', best.c));
  exportgraphics(f2,fullfile(resdir,'powergate_2mM_features.png'),'Resolution',140);
catch ME; fprintf('figure skipped: %s\n', ME.message); end

dst=fullfile(root,'params','lowATP_powergate2_opt.m');
copyfile(fullfile(root,'params','rskR2_w025_opt.m'), dst);
fid=fopen(dst,'a');
fprintf(fid,'\r\n%%%% ==== Low-ATP fit, 8 mM-inert ATP power gates round 2 (%s) ====\r\n', datestr(now,'yyyy-mm-dd'));
fprintf(fid,'%%%% 2 mM absolute cost %.4f (8 mM floor on same list %.4f). 8 mM bit-identical.\r\n', best.c, FLOOR);
fprintf(fid,'params0.UseAtpPowerGates = 1;\r\nparams0.MgATP_ref = 8;\r\n');
fprintf(fid,'params0.UseAtpDetachR2 = 0;\r\nparams0.UseAtpKmsrd = 0;\r\n');
for i=1:numel(P); fprintf(fid,'params0.pAtp_%s = %.6g;\r\n', P{i}, best.x(i)); end
fprintf(fid,'params0.MgATP = 2;\r\nparams0.MgADP = 0;\r\nparams0.Pi = 0.8;\r\n');
fclose(fid);
save(STATE2,'best','P','FLOOR','nEval','cc','FN');
fprintf('\nSaved %s\nPOWERGATE2 DONE\n', dst);

function q=applyP2(q,P,x)
    q.UseAtpPowerGates=1;
    for i=1:numel(P); q.(['pAtp_' P{i}])=x(i); end
end
function c=objPG2(x,b,P,LB,UB)
    if any(x<LB)||any(x>UB); c=1e6; return; end
    q=applyP2(b,P,x);
    try
        [~,~,fm,fd]=runSlackExperiment(q);
        c=sum(evalFeatureCost(fd,fm,b.fn,2));
        if ~isfinite(c); c=1e6; end
    catch; c=1e6; end
end
