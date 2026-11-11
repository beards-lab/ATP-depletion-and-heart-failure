% Baseline.m — measure the post-restretch force redevelopment rate,
% model vs data, at the current best parameter set.
%
% Uses the SHARED estimator Analyses/RestretchVsKtrRecovery/recoveryWindows.m
% so the numbers are directly comparable to that analysis' Part 1/2 tables.
%
% Writes results/baseline.mat with the per-cycle rates and the feature cost
% breakdown, so later steps can diff against it.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root);                       % model + data paths are root-relative
addpath(genpath(root));
if ~exist(fullfile(here, 'results'), 'dir'); mkdir(fullfile(here, 'results')); end

SNAP = 'params/optfull7_opt_mov.m';

params0 = getParams(loadParams(SNAP), [], true, false);
params0.PlotEachSeparately = 0;
params0.PlotFeatureFitting = 0;
params0.RunForceVelocity   = 1;
params0.RunSlack           = 1;
params0.RunSlackPassive    = 1;
params0.RunKtr             = 0;
params0.RunStairs          = 0;
params0.EvalFeatures       = 1;
params0.MaxRunTime         = 600;   % diagnostic run: no optimiser time pressure
ATP_c = params0.MgATP;
FV_velocities = params0.FV_velocities;

fprintf('=== Baseline: %s ===\n', SNAP);
t_run = tic;
RunBakersExp;
fprintf('RunBakersExp took %.1f s\n', toc(t_run));

fprintf('\nE = %s   (sum %.4f)\n', mat2str(round(E, 4)), sum(E));

%% Feature cost breakdown
[ct, w, c] = evalFeatureCost(features_data, features_model, params0.fn, 1);
fprintf('\n%-52s %10s %8s %10s\n', 'feature', 'raw', 'weight', 'weighted');
for i = 1:numel(params0.fn)
    nm = params0.fn{i};
    if numel(nm) > 50; nm = [nm(1:47) '...']; end
    fprintf('%-52s %10.4f %8.3g %10.4f\n', nm, c(i), w(i), ct(i));
end
fprintf('%-52s %10s %8s %10.4f\n', 'TOTAL(features)', '', '', sum(ct));

%% Post-restretch / post-slack recovery rates, model vs data
ds  = load(fullfile(root, 'data', params0.velocitytableonfile));
vt  = ds.velocitytable;
dt  = ds.datatable;

% data F_iso: median of the 50 ms before the first slack release
i1    = find(vt(:,2) < -1, 1);
t_rel = vt(i1, 1);
Fiso_d = median(dt(dt(:,1) > t_rel-0.05 & dt(:,1) < t_rel, 3));

win_d = recoveryWindows(dt(:,1), dt(:,3), vt, 'slack', Fiso_d);

% model — NOTE: RunBakersExp leaves `out` pointing at the LAST protocol it ran
% (here the PASSIVE run), so the slack trace must be taken from out_slack.
nonrep = makeMonotonous(out_slack.t);
Fiso_m = interp1(out_slack.t(nonrep), out_slack.Force(nonrep), t_rel - 0.005);
win_m  = recoveryWindows(out_slack.t(nonrep), out_slack.Force(nonrep), vt, 'slack', Fiso_m);

fprintf('\nF_iso: data %.2f kPa   model %.2f kPa\n', Fiso_d, Fiso_m);

types = {'postSlack', 'postRestretch'};
res = struct();
for it = 1:2
    ty = types{it};
    sd = win_d(strcmp({win_d.type}, ty));
    sm = win_m(strcmp({win_m.type}, ty));
    n  = min(numel(sd), numel(sm));
    fprintf('\n--- %s ---\n', ty);
    fprintf('%-14s %s\n', 'cycle',         sprintf('%8d', 1:n));
    fprintf('%-14s %s\n', 'DATA k63',      sprintf('%8.1f', [sd(1:n).k63]));
    fprintf('%-14s %s\n', 'MODEL k63',     sprintf('%8.1f', [sm(1:n).k63]));
    fprintf('%-14s %s\n', 'ratio M/D',     sprintf('%8.2f', [sm(1:n).k63]./[sd(1:n).k63]));
    fprintf('%-14s %s\n', 'DATA kfixA',    sprintf('%8.1f', [sd(1:n).kfixA]));
    fprintf('%-14s %s\n', 'MODEL kfixA',   sprintf('%8.1f', [sm(1:n).kfixA]));
    fprintf('%-14s %s\n', 'ratio M/D',     sprintf('%8.2f', [sm(1:n).kfixA]./[sd(1:n).kfixA]));
    fprintf('%-14s %s\n', 'DATA A/Fiso',   sprintf('%8.3f', [sd(1:n).A]));
    fprintf('%-14s %s\n', 'MODEL A/Fiso',  sprintf('%8.3f', [sm(1:n).A]));
    fprintf('%-14s %s\n', 'DATA r2(Afix)', sprintf('%8.2f', [sd(1:n).r2A]));
    fprintf('%-14s %s\n', 'MODEL r2(Afix)',sprintf('%8.2f', [sm(1:n).r2A]));
    res.(ty) = struct('data', sd(1:n), 'model', sm(1:n));
end

save(fullfile(here, 'results', 'baseline.mat'), ...
     'res', 'E', 'ct', 'w', 'c', 'features_data', 'features_model', ...
     'Fiso_d', 'Fiso_m', 'SNAP');
fprintf('\nSaved results/baseline.mat\n');
