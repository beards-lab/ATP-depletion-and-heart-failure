%% RunCostTiming.m -- what does each ingredient of the noise fix cost per evaluation?
%
% Full objective (evaluateBakersExp), run the way the campaigns run it
% (RunSlackSegments = 'AllPar' on a 5-worker pool), on PARAMSET. Each variant is
% evaluated NEV times back to back; the first absorbs pool/JIT warm-up, the rest
% are the timing. Variants override the numerical settings only.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
if ~exist('NEV', 'var'), NEV = 3; end
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..')));
res = fullfile(here, 'results');

if ~exist('VAR', 'var')
VAR = { 'old numerics (RelTol 1e-2, step output, hard rectifier)', struct('SmoothRectifierEps', 0, 'odeRelTol', 1e-2, 'odeAbsTol', 1e-4, 'odeOutputDt', 0, 'ClampStrainDepExtrapolation', 0)
        'fix: SiLU + RelTol 1e-3 + dense 0.5 ms (campaign)',     struct('SmoothRectifierEps', 0.01, 'odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 5e-4, 'ClampStrainDepExtrapolation', 1)
        'fix with dense 1 ms',                                   struct('SmoothRectifierEps', 0.01, 'odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 1e-3, 'ClampStrainDepExtrapolation', 1)
        'fix without dense output',                              struct('SmoothRectifierEps', 0.01, 'odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 0, 'ClampStrainDepExtrapolation', 1)
        'fix at RelTol 1e-4 + dense 0.5 ms',                     struct('SmoothRectifierEps', 0.01, 'odeRelTol', 1e-4, 'odeAbsTol', 1e-7, 'odeOutputDt', 5e-4, 'ClampStrainDepExtrapolation', 1) };
end
if ~exist('OUTTAG', 'var'), OUTTAG = PARAMSET; end

P0 = getParams(loadParams(PARAMSET), [], true, false);
P0.RunSlackSegments = 'AllPar'; P0.MaxRunTime = 900; P0.PlotEachSeparately = 0; P0.PlotFeatureFitting = 0;
P0.ghostSave = ''; P0.ghostLoad = ''; P0.SaveBest = 0; P0.mods = {}; P0.g = ones(1,0);
if isempty(gcp('nocreate')), parpool(5); end
T = table();
for v = 1:size(VAR, 1)
    P = P0; ov = VAR{v, 2}; f = fieldnames(ov);
    for k = 1:numel(f), P.(f{k}) = ov.(f{k}); end
    P = getParams(P, P.g, true, false);
    w = nan(NEV, 1); c = nan(NEV, 1);
    for e = 1:NEV
        t0 = tic; c(e) = evaluateBakersExp(1, P, true, false); w(e) = toc(t0);
    end
    T = [T; table(string(VAR{v,1}), w(1), median(w(2:end)), c(end), 'VariableNames', {'variant', 'first_s', 'wall_s', 'cost'})]; %#ok<AGROW>
    fprintf('%-58s first %6.1f s | timed %6.1f s | cost %.4f\n', VAR{v,1}, w(1), median(w(2:end)), c(end));
    save(fullfile(res, ['costtiming_' OUTTAG '.mat']), 'T', 'VAR');
end
T.ratio = T.wall_s / T.wall_s(1);
disp(T);
save(fullfile(res, ['costtiming_' OUTTAG '.mat']), 'T', 'VAR');
