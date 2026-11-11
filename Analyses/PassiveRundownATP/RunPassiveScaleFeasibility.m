% RunPassiveScaleFeasibility.m
% =========================================================================
% CAN the model's passive be scaled up -- and by how much is DEFENSIBLE?
%
% Three constraints are compared on the same axis (the passive scale k_pas):
%   (A) the DIRECT measurement: the PNB+Mava passive protocol (PS_* features).
%       This is what passive was fitted to. It is the hard constraint.
%   (B) the INDIRECT demand: peak2 / vall2_dy / ovrsht_dy in the 2 mM slack run,
%       which want MORE passive (peak2 model 82.3 vs data 92.4).
%   (C) what RUNDOWN can justify.
%
% If (A) and (B) point to different scales, the model cannot satisfy both with a
% single passive parameter -- which is itself the finding, and is exactly the
% "passive during the active runs != passive measured at the end" situation.
%
% ACQUISITION ORDER (raw folder "03 27 2026 M"), which sets (C):
%   01_Relax -> 02_8mM_Active -> 03_2mM_Active -> 04_8mM_Active_repeat
%   -> [05 MISSING] -> 06_8mM_Active_PNB_Mava
% So the passive protocol sits after FOUR activations, not two. The bracket
% (02 vs 04) measured -17.5 % active force, of which the parallel-material
% component is ~16 % at full dose. Passive should therefore be ~16-20 % below the
% first 8 mM run by the time PNB+Mava was recorded -- MORE than the 7.6 % used in
% the previous pass, but still nowhere near the ~x3 that peak2 wants.
%
% NOTE ON SL: the raw export carries only  Time (ms) | L in (Lo) | F in (kPa).
% There is NO measured sarcomere length -- the SL column in the curated .mat is
% derived from motor length. So the passive-vs-SL relation cannot be established
% independently of the model's own series-compliance assumption.
% =========================================================================
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

% Branch-portable baseline: the low-ATP snapshot exists only on lowatp-2state, so
% fall back to the shared high-ATP baseline elsewhere (e.g. on Further-optim).
BASE='lowATP_powergate2_opt';
if ~isfile(fullfile(root,'params',[BASE '.m'])); BASE='rskR2_w025_opt'; end
fprintf('baseline: %s\n', BASE);
b=getParams(loadParams(BASE),[],true,false);
b.EvalFeatures=true; b.BreakOnODEUnstable=false; b.MaxRunTime=400;
b.PlotEachSeparately=0; b.PlotFeatureFitting=0; b.RunSlackSegments='AllPar';
b.passive_velocitytableonfile='protocol_03_27_2026_ActivePNBMava_slack.mat';
KP=b.k_pas;
m=@(x)mean(double(x),'omitnan');
f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
PSFN={'PS_restretchPeak','PS_rampupRMSE|0.005','PS_holdDecayRMSE|0.01','PS_steady22|0.5','PS_steady20|0.5'};

sc=[0.5 0.75 1.0 1.25 1.5 2.0 3.0 4.0];
PS=nan(size(sc)); P2=nan(size(sc)); C2=nan(size(sc));
fn=b.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
fn=cellfun(@(s) regexprep(s,'^rsK(\|.*)?$','rsK|0.1'), fn, 'UniformOutput',false);
FN=fn(~(startsWith(fn,'FV')|startsWith(fn,'PS_')));

fprintf('%7s | %10s | %8s %8s | %s\n','k_pas x','PASSIVE PS','peak2_2','cost2','(data peak2 = 92.43)');
for i=1:numel(sc)
    % (A) direct passive protocol
    qp=b; qp.k_pas=KP*sc(i);
    qp.RunForceVelocity=false; qp.RunForceVelocityTime=false; qp.RunKtr=false;
    qp.RunStairs=false; qp.RunSlack=false; qp.RunSlackPassive=true;
    qp.MgATP=8; qp.MgADP=0; qp.Pi=0;
    try
        [fmp,fdp,~]=runPassiveExperiment(qp);   % NB: model first, then data
        PS(i)=sum(evalFeatureCost(fdp,fmp,PSFN,2));
    catch ME
        fprintf('   (passive %0.2f failed: %s)\n', sc(i), ME.message);
    end
    % (B) 2 mM slack
    q2=b; q2.k_pas=KP*sc(i);
    q2.RunForceVelocity=false; q2.RunForceVelocityTime=false; q2.RunKtr=false;
    q2.RunStairs=false; q2.RunSlackPassive=false; q2.RunSlack=true;
    q2.MgATP=2; q2.MgADP=0; q2.Pi=0.8;
    q2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
    try
        [~,~,g2,d2]=runSlackExperiment(q2);
        P2(i)=m(g2.peak2); C2(i)=sum(evalFeatureCost(d2,g2,FN,2));
    catch; end
    fprintf('%7.2f | %10.4f | %8.2f %8.4f |\n', sc(i), PS(i), P2(i), C2(i));
end

[~,iA]=min(PS); [~,iB]=min(C2);
[~,iP2]=min(abs(P2-m(f2d.peak2)));
fprintf('\n(A) DIRECT passive protocol prefers  k_pas x%.2f   (PS cost %.4f)\n', sc(iA), PS(iA));
fprintf('(B) 2 mM slack total cost prefers     k_pas x%.2f   (cost %.4f)\n', sc(iB), C2(iB));
fprintf('(B*) peak2 alone would need           k_pas x%.2f   (peak2 %.2f vs data %.2f)\n', ...
    sc(iP2), P2(iP2), m(f2d.peak2));
fprintf('(C) RUNDOWN can justify at most       k_pas x%.2f   (~16-20%% over 4 activations)\n', 1.20);
fprintf('\nPS cost at the rundown-justifiable x1.20 : %.4f  (vs %.4f at x1.00) -> %+.4f\n', ...
    interp1(sc,PS,1.20), interp1(sc,PS,1.0), interp1(sc,PS,1.20)-interp1(sc,PS,1.0));
if sc(iP2) > 2.0 && sc(iA) < 1.6
    fprintf('\nCONFLICT: peak2 demands a passive scale the DIRECT measurement rejects.\n');
    fprintf('=> the model cannot satisfy both with one passive parameter; peak2 is NOT\n');
    fprintf('   simply an under-scaled-passive problem.\n');
end
save(fullfile(resdir,'passive_scale_feasibility.mat'),'sc','PS','P2','C2');
try
  fg=figure('Position',[60 60 1200 450],'Color','w');
  subplot(1,2,1); plot(sc,PS,'o-','LineWidth',1.6); grid on; box on; hold on;
  xline(1,'k--'); xline(1.2,'g--','rundown limit');
  xlabel('k_{pas} \times'); ylabel('PNB+Mava passive cost'); title('(A) DIRECT passive measurement');
  subplot(1,2,2); yyaxis left; plot(sc,P2,'s-','LineWidth',1.6); ylabel('peak2 (2 mM)');
  hold on; yline(m(f2d.peak2),'r--','data'); yyaxis right; plot(sc,C2,'^-'); ylabel('2 mM cost');
  grid on; box on; xlabel('k_{pas} \times'); title('(B) what the 2 mM slack wants');
  exportgraphics(fg,fullfile(resdir,'passive_feasibility.png'),'Resolution',160);
catch; end
disp('PASSIVE FEASIBILITY DONE');
