% RunCostNoiseAmplitude.m -- how much of the 2 mM objective is numerical noise?
%
% Three accidental samples said wildly different things about the size of the floor:
%   rate factors perturbed by 3e-16  ->  cost moved 0.016   (the Hill conversion)
%   exponents perturbed by 8e-05     ->  cost moved 0.067   (dropping dust gates)
%   exponents perturbed by 4e-07     ->  cost moved 0.602   (%.6g snapshot rounding)
% Non-monotone in the perturbation size => this is not sensitivity, it is a rough
% surface, and a single sample cannot size it. So sample it properly: nudge ONE
% exponent by a ladder of deltas far too small to be physically meaningful and look at
% the SPREAD of the resulting costs. Repeat at a tighter tolerance: if the spread is
% truncation error it must shrink.
%
% This bounds every cost comparison in the low-ATP work, including the fit that
% produced the incumbent.
%
% Run: matlab -batch "run('Analyses/AtpGateHillForm/RunCostNoiseAmplitude.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end
if isempty(gcp('nocreate')); parpool('local',6); end

BASE  = 'lowATP_powergate2_opt';
KNOB  = 'pAtp_kmsrd';                       % the strongest gate; p = 0.636193
DELTA = [0 1e-9 -1e-9 1e-8 -1e-8 1e-7 -1e-7 1e-6 -1e-6 1e-5 -1e-5];
TOLS  = [1e-2 1e-3];

b0 = getParams(loadParams(BASE), [], true, false);
fprintf('\nperturbing %s (= %.9g) by deltas that change the 2 mM rate factor by <2e-5 relative\n', ...
        KNOB, b0.(KNOB));

C = nan(numel(TOLS), numel(DELTA));
for it = 1:numel(TOLS)
    fprintf('\n--- RelTol %.0e ---\n', TOLS(it));
    for id = 1:numel(DELTA)
        q = iCfg(getParams(loadParams(BASE), [], true, false));
        q.(KNOB) = q.(KNOB) + DELTA(id);
        q.odeRelTol = TOLS(it); q.odeAbsTol = TOLS(it)/100;
        try
            [~,~,fm,fd] = runSlackExperiment(q);
            C(it,id) = sum(evalFeatureCost(fd,fm,q.fn,2));
        catch ME
            fprintf('   delta %+9.0e  FAILED (%s)\n', DELTA(id), regexprep(ME.message,'\s+',' '));
        end
        fprintf('   delta %+9.0e   cost %.6f\n', DELTA(id), C(it,id));
        save(fullfile(resdir,'costnoise.mat'),'C','DELTA','TOLS','KNOB','BASE');
    end
    v = C(it,:); v = v(isfinite(v));
    fprintf('  spread: min %.4f  max %.4f  PEAK-TO-PEAK %.4f  std %.4f  (n=%d)\n', ...
            min(v), max(v), max(v)-min(v), std(v), numel(v));
end

fprintf('\n=== summary ===\n');
for it = 1:numel(TOLS)
    v = C(it,:); v = v(isfinite(v));
    fprintf('RelTol %.0e : peak-to-peak %.4f   std %.4f   mean %.4f\n', ...
            TOLS(it), max(v)-min(v), std(v), mean(v));
end
fprintf(['\nAny reported improvement smaller than the peak-to-peak at the production tolerance is\n' ...
         'not distinguishable from re-rolling the solver.\n']);
fprintf('\nCOST NOISE AMPLITUDE DONE\n');

function q = iCfg(q)
q.RunForceVelocity=false; q.RunForceVelocityTime=false; q.RunKtr=false; q.RunStairs=false;
q.RunSlackPassive=false;  q.RunSlack=true; q.EvalFeatures=true; q.BreakOnODEUnstable=false;
q.RunSlackSegments='AllPar'; q.MaxRunTime=600;
q.PlotEachSeparately=0; q.PlotFeatureFitting=0;
q.UseAtpDetachR2=0; q.UseAtpKmsrd=0; q.UseStrainDepAtpGate=0; q.UseR2Ceiling=0; q.UseD0State=0;
q.MgATP=2; q.MgADP=0; q.Pi=0.8;
q.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
fn = q.fn; [~,iu] = unique(fn,'stable'); fn = fn(iu);
fn = cellfun(@(s) regexprep(s,'^rsK(\|.*)?$','rsK|0.1'), fn, 'UniformOutput', false);
q.fn = fn(~(startsWith(fn,'FV') | startsWith(fn,'PS_')));
end
