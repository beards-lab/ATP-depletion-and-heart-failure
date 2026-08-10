% RunPassiveRundownTest.m
% =========================================================================
% HYPOTHESIS (user): the PNB+Mava passive run was recorded LAST, so it measures passive
% in the MOST rundown state. If rundown removes force-generating material in parallel,
% those torn myofibrils also stop contributing PASSIVE (titin) force, so
%     passive(8 mM run)  >  passive(2 mM run)  >  passive(PNB+Mava, as fitted)
% The model fits passive to the PNB+Mava level, so it under-assigns passive to BOTH
% active runs, and more so at 8 mM. That is a second, independent handle on the
% high-vs-low ATP difference.
%
% STEP 1 -- consistency with ../RundownCorrection/conclusions.md:
%   * kSE: CONFIRMED. Rundown adds ~35 % series compliance.
%   * SL0 : REVISED. The length-shift evidence is "provisional and resting on a single
%           measurement" (§4); on five points the loss is a SCALE, not a shift.
%   * The revised lesion is "loss of force-generating material in parallel (~16 %) plus
%     added series compliance (~35 %)" -- the "tearing the cell" picture.
%   * PNB+Mava "is charged zero damage" and ran last => it reads the most-rundown state.
%   * TENSION FOUND: §4 asserts in passing that "passive does not run down", but the
%     REVISED lesion (torn myofibrils lost in parallel) logically implies passive DOES
%     fall with rundown, by the same parallel-material fraction. The hypothesis is
%     consistent with the revised mechanism and inconsistent only with that older aside.
%
% STEP 2 -- QUANTITATIVE PREDICTION (before touching the model):
%   passive ~ 5 kPa (the PNB+Mava level, §4); parallel loss ~16 % at FULL bracket dose;
%   03/27 dose lambda = 0.17  =>  material lost between the 8 mM and 2 mM runs
%   ~ 0.17*16 % = 2.7 %  =>  predicted passive difference ~ 5 kPa * 0.027 = 0.14 kPa.
%   Total force ~ 90 kPa, and the ATP force effect is x1.36 (~25 kPa).
%   So the predicted passive effect is ~0.15 % of force and ~180x SMALLER than the
%   signal it would correct. PREDICTION: real in direction, negligible in magnitude.
%
% STEP 3 -- THE TEST (falsifiable):
%   Free the passive scale at the 2 mM condition only and ask what the data wants.
%     optimum ~0.97 (-2.7 %)  => consistent with rundown, and negligible  [H supported]
%     optimum far from 0.97   => the knob is absorbing a DIFFERENT model error, so a
%                                passive difference is not the explanation  [H rejected]
%     no cost improvement     => the data does not constrain passive here  [H untestable]
% Run: matlab -batch "run('Analyses/PassiveRundownATP/RunPassiveRundownTest.m')"
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
b.RunForceVelocity=false; b.RunForceVelocityTime=false; b.RunKtr=false; b.RunStairs=false;
b.RunSlackPassive=false; b.RunSlack=true; b.EvalFeatures=true; b.BreakOnODEUnstable=false;
b.RunSlackSegments='AllPar'; b.MaxRunTime=400; b.PlotEachSeparately=0; b.PlotFeatureFitting=0;
fn=b.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
fn=cellfun(@(s) regexprep(s,'^rsK(\|.*)?$','rsK|0.1'), fn, 'UniformOutput',false);
FN=fn(~(startsWith(fn,'FV')|startsWith(fn,'PS_'))); b.fn=FN;
m=@(x)mean(double(x),'omitnan');

% ---- how big is passive in the model, and is it the same at both ATP? ----
q8=b; q8.MgATP=8; q8.Pi=0; q8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
q2=b; q2.MgATP=2; q2.Pi=0.8; q2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
[~,o8,f8,d8]=runSlackExperiment(q8);
[~,o2,f2,d2]=runSlackExperiment(q2);
% FXBPassive is the parallel (titin) passive force trace.
pk=@(o) max(o.FXBPassive);                       % peak, i.e. during the restretch
ss=@(o) median(o.FXBPassive(round(0.9*numel(o.FXBPassive)):end));   % settled level
pas8=pk(o8); pas2=pk(o2);
fprintf('\n--- model passive (parallel/titin, FXBPassive) ---\n');
fprintf('PEAK passive : 8 mM %.3f | 2 mM %.3f\n', pas8, pas2);
fprintf('settled level: 8 mM %.3f | 2 mM %.3f  (identical by construction -- passive is not ATP-gated)\n', ss(o8), ss(o2));
fprintf('peak passive as %% of steady force: 8 mM %.2f%% | 2 mM %.2f%%\n', 100*pas8/m(f8.steady), 100*pas2/m(f2.steady));
c2base=sum(evalFeatureCost(d2,f2,FN,2));
fprintf('2 mM baseline cost: %.4f\n', c2base);

% ---- predicted rundown effect ----
PAS_KPA=5; LOSS_FULL=0.16; LAMBDA_0327=0.17;
predDrop=LOSS_FULL*LAMBDA_0327;
fprintf('\n--- prediction ---\nparallel loss between the 8 and 2 mM runs = %.1f%% -> passive difference ~ %.3f kPa\n', ...
    100*predDrop, PAS_KPA*predDrop);
fprintf('predicted passive SCALE the 2 mM run should want: %.4f\n', 1-predDrop);

% ---- THE TEST: what passive scale does the 2 mM data actually want? ----
sc=[0.4 0.6 0.8 0.9 0.973 1.0 1.1 1.3 1.6 2.0 3.0];
fprintf('\n%8s | %9s | %8s %8s %8s %8s\n','k_pas x','cost','peakPas','steady','peak1_y','vall_y');
C=nan(size(sc)); PK=nan(size(sc));
for i=1:numel(sc)
    q=q2; q.k_pas=b.k_pas*sc(i);
    try
        [~,oo,ff,dd]=runSlackExperiment(q);
        C(i)=sum(evalFeatureCost(dd,ff,FN,2)); PK(i)=max(oo.FXBPassive);
        fprintf('%8.3f | %9.4f | %8.3f %8.2f %8.2f %8.2f\n', sc(i),C(i),PK(i), ...
            m(ff.steady),m(ff.peak1_y),m(ff.vall_y));
    catch ME
        fprintf('%8.3f | FAILED %s\n', sc(i), ME.message);
    end
end
[cb,ib]=min(C);
fprintf('\nDATA 2 mM: steady %.2f peak1_y %.2f vall_y %.2f\n', m(d2.steady),m(d2.peak1_y),m(d2.vall_y));
fprintf('\n=== BEST passive scale wanted by the 2 mM data: x%.3f (cost %.4f, from %.4f) ===\n', ...
    sc(ib), cb, c2base);
fprintf('    rundown PREDICTS x%.4f\n', 1-predDrop);
if abs(sc(ib)-(1-predDrop)) < 0.10
    verdict='CONSISTENT with the rundown prediction';
elseif (c2base-cb) < 0.05
    verdict='UNCONSTRAINED -- passive barely affects the 2 mM cost';
else
    verdict='INCONSISTENT -- the knob is absorbing a different model error';
end
fprintf('    VERDICT: %s\n', verdict);
fprintf('    cost bought by the whole passive axis: %.4f of %.4f (%.1f%%)\n', ...
    c2base-cb, c2base, 100*(c2base-cb)/c2base);

save(fullfile(resdir,'passive_rundown_test.mat'),'sc','C','PK','c2base','predDrop','pas8','pas2');
try
  fg=figure('Position',[60 60 1100 420],'Color','w');
  subplot(1,2,1); plot(sc,C,'o-','LineWidth',1.5); hold on; grid on; box on;
  xline(1,'k--'); xline(1-predDrop,'r--','rundown prediction');
  xlabel('passive scale (k_{pas} \times)'); ylabel('2 mM absolute cost'); title('What passive scale does 2 mM want?');
  subplot(1,2,2); plot(sc,PK,'s-','LineWidth',1.5); grid on; box on;
  xlabel('passive scale'); ylabel('peak titin passive'); title('Passive magnitude');
  exportgraphics(fg,fullfile(resdir,'passive_scale_scan.png'),'Resolution',160);
catch; end
disp('PASSIVE RUNDOWN TEST DONE');
