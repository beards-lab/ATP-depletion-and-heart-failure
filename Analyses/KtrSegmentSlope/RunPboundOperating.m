% RunPboundOperating.m -- what bound fraction does the model actually run at?
%
% The scan says pulling P_bound_max from 1.208 to 0.12 takes the objective from
% 0.80 to 599.5, essentially all of it F-V. That only makes sense if the model's
% own P_bound is nowhere near 0.15, i.e. if comparing 1.208 against the
% literature's ~0.12-0.15 bound-head fraction was never apples to apples.
% `P_bound = p1_0 + p2_0 (+ p3_0)` here is a strain-integrated probability in
% this model's normalisation, not a measured duty ratio.
%
% Measure it: isometric plateau, during the slack hold, and at the F-V
% velocities, and report f_sat = 1 - P_bound/P_bound_max at each.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunPboundOperating.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_pbo_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 5); end

p = getParams(loadParams('refit5_seed'), [], true, false);
p.RunSlackSegments='All'; p.EvalFeatures=1; p.PlotEachSeparately=0;
p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0; p.BreakOnODEUnstable=0; p.MaxRunTime=600;
[~, os] = runSlackExperiment(p);

fn = fieldnames(os);
fprintf('out_slack fields: %s\n', strjoin(fn', ', '));
if isfield(os,'outputs')
    fprintf('outputs is %d x %d\n', size(os.outputs,1), size(os.outputs,2));
end

% p1_0 / p2_0 are columns 5,6 of `outputs` for the 2-state layout (see the
% OUTPUTS comment at the head of dPUdT_CombinedTransitions).
ns = p.NumberOfStates;
if isfield(os,'outputs') && ~isempty(os.outputs)
    O = os.outputs;
    if ns >= 3 && size(O,2) >= 7
        Pb = O(:,5) + O(:,6) + O(:,7);
    else
        Pb = O(:,5) + O(:,6);
    end
    t = os.t(:);
    fprintf('\nP_bound over the whole slack protocol: min %.4f  median %.4f  max %.4f\n', ...
            min(Pb), median(Pb), max(Pb));
    fprintf('P_bound_max = %.4f  ->  f_sat = 1 - P_bound/P_bound_max\n', p.P_bound_max);
    fprintf('   at the isometric plateau (P_bound %.4f): f_sat %.4f\n', ...
            median(Pb), 1 - median(Pb)/p.P_bound_max);
    fprintf('   at the protocol maximum  (P_bound %.4f): f_sat %.4f\n', ...
            max(Pb), 1 - max(Pb)/p.P_bound_max);
    fprintf('   the ceiling would have to be below %.3f to shut attachment off\n', max(Pb));
    save(fullfile(resdir,'pbound_operating.mat'),'Pb','t');
else
    fprintf(2,'no `outputs` in out_slack -- cannot read P_bound from here\n');
end

fprintf('\nDONE\n');
