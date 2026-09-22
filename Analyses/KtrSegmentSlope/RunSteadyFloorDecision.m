% RunSteadyFloorDecision.m -- how much of the force block's error is LEVEL and
% how much is SHAPE, and what a tightened `steady` deadband would actually buy.
%
% RunSteadySpread established that the force-class features carry a deadband
% built from a 21-35 % between-HEART difference in absolute force, while the
% shape each one describes agrees between hearts to 1.6-2.6 %. That is a ~100x
% error in the floor, because fl goes as s^2.
%
% But a wide floor only matters if the model is actually sitting in the slack
% it creates. This script measures that directly at refit4_opt:
%   * the raw (unfloored) cost of every force term
%   * split into a LEVEL error (one uniform scale factor) and a SHAPE residual
%   * what each term would cost under a shape-based floor
% so the decision to tighten is made on numbers, not on the principle.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunSteadyFloorDecision.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~isfolder(resdir); mkdir(resdir); end

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_sfd_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 5); end

p = getParams(loadParams('refit4_opt'), [], true, false);
p.RunSlackSegments = 'AllPar'; p.EvalFeatures = 1;
p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0;
p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0; p.BreakOnODEUnstable=0; p.MaxRunTime=600;
[~, ~, fm, fd] = runSlackExperiment(p);

S = load(fullfile(root,'Analyses','KtrSegmentSlope','results','steady_spread.mat'));
LAMBDA = 0.0106294;

FEATS = {'steady','A','peak1_y','vall_y','peak2','vall2_dy','restretchSlopeStart'};
fprintf('=== refit4_opt vs the pooled target, force block ===\n\n');
out = {};
for i = 1:numel(FEATS)
    nm = FEATS{i};
    if ~isfield(fm,nm) || ~isfield(fd,nm); fprintf('%s: missing\n',nm); continue; end
    D = double(fd.(nm)(:))'; M = double(fm.(nm)(:))';
    n = min(numel(D),numel(M)); D = D(1:n); M = M(1:n);
    nF = mean(D);
    c_raw = sum(((D - M)/nF).^2);

    % LEVEL / SHAPE split: the best single scale factor a, then what is left
    a = (M*D')/(M*M');                       % D ~ a*M  in least squares
    c_level = sum(((D - M)/nF).^2) - sum(((D - a*M)/nF).^2);
    c_shape = sum(((D - a*M)/nF).^2);

    j = find(strcmp(S.rows(:,1), nm), 1);
    s_used = S.rows{j,2}; s_lev = S.rows{j,4}; s_shp = S.rows{j,5};
    fl_now   = (1/3)*n*s_used^2;
    fl_shape = (1/3)*n*max(s_shp,0.01)^2;
    w_now    = LAMBDA/s_used^2;

    fprintf('%s  (nEl %d)\n', nm, n);
    fprintf('   data  %s\n', sprintf('%9.3g', D));
    fprintf('   model %s\n', sprintf('%9.3g', M));
    fprintf('   resid %s   (relative)\n', sprintf('%9.2f%%', 100*(M-D)./D));
    fprintf('   level factor model/data = %.4f  (%+.1f %%)\n', 1/a, 100*(1/a - 1));
    fprintf('   raw cost %.5f  =  level %.5f  +  shape %.5f\n', c_raw, c_level, c_shape);
    fprintf('   floor now %.5f (s %.3f)  ->  floored %.5f, weighted %.5f\n', ...
            fl_now, s_used, max(0,c_raw-fl_now), w_now*max(0,c_raw-fl_now));
    fprintf('   floor shape-based %.5f (s_shape %.3f)  ->  floored %.5f\n', ...
            fl_shape, s_shp, max(0,c_raw-fl_shape));
    if c_raw < fl_now
        fprintf('   >> currently INSIDE the deadband by %.1fx\n', fl_now/max(c_raw,1e-12));
    end
    fprintf('\n');
    out(end+1,:) = {nm, n, c_raw, c_level, c_shape, s_used, s_lev, s_shp, fl_now, fl_shape, w_now}; %#ok<SAGROW>
end

%% ---- what a tightened steady would cost, at a range of floors ----------
fprintf('=== sweep: what `steady` pays at different deadbands ===\n');
j = find(strcmp(out(:,1),'steady'),1);
c_raw = out{j,3}; w = out{j,11}; nEl = out{j,2};
fprintf('  raw cost %.5f, weight %.4f -> unfloored weighted %.5f\n', c_raw, w, w*c_raw);
fprintf('%-30s %10s %10s %12s\n','deadband basis','s','fl','weighted cost');
BAS = { 'current (between-heart)',      out{j,6}
        'SEM of the 3-heart mean',      out{j,7}/sqrt(3)
        'shape only',                   max(out{j,8},0.005)
        'shape x2',                     2*max(out{j,8},0.005)
        'zero (no deadband)',           0 };
for i = 1:size(BAS,1)
    s = BAS{i,2}; fl = (1/3)*nEl*s^2;
    fprintf('%-30s %10.4f %10.5f %12.5f\n', BAS{i,1}, s, fl, w*max(0,c_raw-fl));
end

save(fullfile(resdir,'steady_floor_decision.mat'),'out','FEATS');
fprintf('\nDONE\n');
