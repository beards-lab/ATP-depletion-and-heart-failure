% MechanismDirectionScan.m
% =========================================================================
% Systematic direction scan: scale each candidate cross-bridge / SRX rate at 8 mM (simulating a
% nucleotide nudge) and measure force/ktr/t0 vs baseline, to see which sarcomeric target — when
% [ATP]/[ADP]/[Pi] acts on it — moves TOWARD the low-ATP signature (force UP, ktr DOWN, t0 UP).
%
% Findings (2026-07-09, optfull_opt baseline steady 77.3 / ktr 53.1 / t0 0.0091):
%   k2  x0.5 -> force 1.26 ktr 0.69 t0 1.83  RIGHT on all three (R2 detachment; = ATP-detach AND ADP-trap)
%   kmsrd x0.5 -> force 0.99 ktr 0.74        DECOUPLES (ktr down, force held) -- the SRX weld-breaker
%   ksr0 x2   -> force 1.00 ktr 0.94         weak decoupler
%   ka  x0.5 -> force 0.71 ktr 1.11          WRONG (recruitment-limit lowers force)
%   kah x0.5 -> force 0.95 ktr 0.96          WRONG (force down)
%   k1  x0.5 -> force 0.91                    WRONG (Pi forward-inhib lowers force)
%   k_1 x2   -> force 0.92 ktr 1.07           WRONG on both (Pi reverse stroke)
%   k3  -> inert (no P3 rigor state in 2-state; so ATP-detach == ADP-trap here, degenerate)
% Conclusion: in 2-state at max Ca the only productive levers are R2-slowing (k2, force+coarse ktr)
% and SRX kmsrd (fine kinetic decoupling). Cooperative thin-filament activation is untestable here
% (needs force-pCa, not max-Ca slack).
%
% NB: the snapshot assigns to params0 -> initialize a var NAMED params0 before run(), else it is ignored.
% Run: cd(root); addpath(genpath('.')); Analyses/LowATP_2State/MechanismDirectionScan
% =========================================================================
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
params0=getParams(); run(fullfile(root,'params','optfull_opt.m')); params0=getParams(params0,[],true,false);
params0.RunForceVelocity=false; params0.RunKtr=false; params0.RunSlack=true; params0.RunStairs=false; params0.RunForceVelocityTime=false;
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false; params0.RunSlackSegments='All'; params0.MaxRunTime=150;
params0.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat'; params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0; params0.ShowStatePlots=0;
base=params0;

[~,~,fm0,~]=runSlackExperiment(base);
S0=mean(fm0.steady,'omitnan'); K0=mean(fm0.ktr,'omitnan'); T0=mean(fm0.t0,'omitnan');
fprintf('baseline 8mM: steady %.1f ktr %.1f t0 %.4f | target: force UP, ktr DOWN, t0 UP\n', S0,K0,T0);
fprintf('%-8s %6s | %7s %7s %7s\n','rate','x','forceR','ktrR','t0R');
rates={'k2','ka','kah','kamh','k1','k_1','kmsrd','ksr0','k3'};
for ri=1:numel(rates); r=rates{ri};
  if ~isfield(base,r) || base.(r)==0; fprintf('%-8s   (skip: =0)\n', r); continue; end
  for f=[0.5 2]
    p=base; p.(r)=base.(r)*f; fR=NaN;kR=NaN;tR=NaN;
    try; [~,~,fm,~]=runSlackExperiment(p);
      fR=mean(fm.steady,'omitnan')/S0; kR=mean(fm.ktr,'omitnan')/K0; tR=mean(fm.t0,'omitnan')/T0; catch; end
    fprintf('%-8s x%4.1f | %7.3f %7.3f %7.3f\n', r, f, fR, kR, tR);
  end
end
