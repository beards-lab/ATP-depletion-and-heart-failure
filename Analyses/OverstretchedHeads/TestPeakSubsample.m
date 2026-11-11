%% TestPeakSubsample.m -- does the sub-sample vertex fix peak1_dSL / vall_t?
%
% peak1_dSL = SL(t_peak) - SLslack is peak1_t re-expressed through a known
% ramp (7.79 um/s SL), and t_peak is a SAMPLE: the model is read at 1 kHz, so
% the feature moves in 7.8 nm steps (28 % of its mean) and the cost is a step
% function of the parameters. extractSlackAttributes now also returns
% peak1_tq / peak1_dSLq / vall_tq, the vertex of a local quadratic. Three tests:
%   1. NOISE: a smooth model trace put on the data's 4 kHz grid, plus white noise
%      at the data's own sigma (1.3 kPa; RecoveryOscillation: the residual is
%      90 % white), NREP replicates through the DATA path: jitter and bias of
%      each estimator against the noise-free truth.
%   2. SMOOTHNESS: a fine kstiff2 scan on the MODEL path: sample vs vertex.
%   3. DATA: the three 8 mM preparations, old vs new, and the between-prep spread.
if ~exist('NREP', 'var'), NREP = 100; end
if ~exist('SIGMA', 'var'), SIGMA = 1.3; end
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
res = fullfile(here, 'results');
FT = {'peak1_dSL', 'peak1_dSLq', 'peak1_t', 'peak1_tq', 'vall_t', 'vall_tq'};
pick = @(f, q) f.(q)(1:4);                          % cycles 1-4 (cycle 5 re-stretches less)

%% 1. noise
L = load(fullfile(res, 'slack_serial_refit9.mat'));   % refit9_opt serial slack run (ProbeRestretchDwell)
o = L.oS; D = load(fullfile(root, 'data', L.P.velocitytableonfile)); vt = D.velocitytable;
[tu, iu] = unique(o.t(:)); tg = (ceil(73.5/2.5e-4):floor(tu(end)/2.5e-4))' * 2.5e-4;   % slack window only: 28 s -> ~1 s per extraction
F0 = interp1(tu, o.Force(iu), tg, 'pchip'); SL0 = interp1(tu, o.SL(iu), tg, 'linear');
truth = extractSlackAttributes(tg, F0, SL0, vt, struct(), [], false, true);
rng(1); N = struct(); for q = FT, N.(q{1}) = nan(NREP, 4); end
for r = 1:NREP
    f = extractSlackAttributes(tg, F0 + SIGMA*randn(size(F0)), SL0, vt, struct(), [], false, true);
    for q = FT, N.(q{1})(r, :) = pick(f, q{1}); end
end
fprintf('\n1. NOISE (%d reps, sigma %.2f kPa, 4 kHz): relative to the mean of the noise-free value\n', NREP, SIGMA);
fprintf('%-12s %10s %12s %12s\n', 'feature', 'truth', 'jitter SD', 'bias');
for q = FT
    tr = pick(truth, q{1}); m = mean(tr);
    fprintf('%-12s %10.5f %11.1f%% %11.1f%%\n', q{1}, m, 100*mean(std(N.(q{1}))) / m, 100*mean(mean(N.(q{1})) - tr) / m);
end

%% 2. smoothness on the model path
P0 = getParams(loadParams('refit9_opt'), [], true, false);
P0.RunSlackSegments = 'AllPar'; P0.MaxRunTime = 900; P0.PlotEachSeparately = 0; P0.PlotFeatureFitting = 0;
P0.ghostSave = ''; P0.ghostLoad = ''; P0.SaveBest = 0; P0.mods = {}; P0.g = ones(1,0);
if isempty(gcp('nocreate')), refreshPool(5); end
gS = linspace(0.98, 1.02, 9); S = struct(); for q = FT, S.(q{1}) = nan(numel(gS), 4); end
for j = 1:numel(gS)
    P = P0; P.kstiff2 = P0.kstiff2 * gS(j); P = getParams(P, P.g, true, false);
    evalc('[~, ~, fm] = runSlackExperiment(P);');
    for q = FT, S.(q{1})(j, :) = pick(fm, q{1}); end
end
fprintf('\n2. SMOOTHNESS: kstiff2 x%.2f..x%.2f (9 points), cycle-1..4 mean\n', gS(1), gS(end));
fprintf('%-12s %s\n', 'kstiff2 x', sprintf('%8.3f', gS));
for q = FT, fprintf('%-12s %s\n', q{1}, sprintf('%8.4f', 1e3*mean(S.(q{1}), 2))); end
fprintf('(dSL in um*1e3, times in ms)\n');

%% 3. the data
DAYS = {'03_27', '04_03', '04_10'}; V = struct(); for q = FT, V.(q{1}) = nan(3, 5); end
for d = 1:3
    D = load(fullfile(root, 'data', sprintf('protocol_%s_2026_8mM_slack.mat', DAYS{d})));
    f = extractSlackAttributes(D.datatable(:,1), D.datatable(:,3), D.datatable(:,2), D.velocitytable, struct(), [], false, true);
    for q = FT, v = f.(q{1})(:)'; V.(q{1})(d, 1:numel(v)) = v; end
end
fprintf('\n3. DATA, 8 mM, per preparation (cycles 1-5)\n');
for q = FT
    fprintf('%-12s', q{1}); for d = 1:3, fprintf('  %s', mat2str(round(1e3*V.(q{1})(d, :), 2))); end
    mu = mean(V.(q{1}), 1); sd = std(V.(q{1}), 0, 1);
    fprintf('  | pooled %s  rel SD %.3f\n', mat2str(round(1e3*mu, 2)), sqrt(mean((sd./abs(mu)).^2)));
end
save(fullfile(res, 'peak_subsample_test.mat'), 'N', 'truth', 'S', 'gS', 'V', 'FT', 'SIGMA', 'NREP');
