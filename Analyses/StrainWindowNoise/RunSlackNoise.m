%% RunSlackNoise.m -- which numerical setting carries the objective's noise?
%
% Every configuration evaluates the same 1 + NREP parameter points: the
% unperturbed set and NREP sets whose lean-pool parameters are each multiplied
% by (1 + 1e-6*randn), with identical draws across configurations. One serial
% slack run per point (parfor over points). Reported per configuration: the
% relative SD of each ATP-carrying slack feature (mean over its elements) and
% of the slack part of the objective, plus the unperturbed values, so a
% setting that removes noise but moves the answer is visible as well.
%
% Set CFGSEL (indices into CFG) and PARAMSET before calling if needed.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit8_opt'; end
if ~exist('OUTTAG', 'var'), OUTTAG = PARAMSET; end
if ~exist('NREP', 'var'), NREP = 8; end
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
addpath(genpath(root));
res  = fullfile(here, 'results');

if ~exist('CFG', 'var')
CFG = { 'baseline (RelTol 1e-2, solver-step output)', struct()
        'dense output 0.5 ms',                       struct('odeOutputDt', 5e-4)
        'RelTol 1e-3',                               struct('odeRelTol', 1e-3, 'odeAbsTol', 1e-6)
        'RelTol 1e-3 + dense',                       struct('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 5e-4)
        'RelTol 1e-4 + dense',                       struct('odeRelTol', 1e-4, 'odeAbsTol', 1e-7, 'odeOutputDt', 5e-4)
        'wide grid 128 + clamp (no shifts)',         struct('MaxStrainArraySize', 128, 'ClampStrainDepExtrapolation', true)
        'wide 128 + clamp + RelTol 1e-3 + dense',    struct('MaxStrainArraySize', 128, 'ClampStrainDepExtrapolation', true, 'odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 5e-4) };
end
if ~exist('CFGSEL', 'var'), CFGSEL = 1:size(CFG, 1); end
FEAT = {'ktr80', 'ktr80_shape', 'rsK80', 'rsT20', 'steady', 'vall_t', 'peak1_dSL'};

P = getParams(loadParams(PARAMSET), [], true, false);
P.RunForceVelocity = false; P.RunSlack = true; P.RunStairs = false; P.RunKtr = false;
P.RunForceVelocityTime = false; P.RunSlackPassive = true; P.EvalFeatures = true;
P.BreakOnODEUnstable = false; P.PlotEachSeparately = 0; P.PlotFeatureFitting = 0;
P.RunSlackSegments = 'All'; P.MaxRunTime = 900; P.OptimizeOn = 'Feats';
P.ghostSave = ''; P.ghostLoad = ''; P.SaveBest = 0; P.mods = {}; P.g = ones(1,0);

% the same perturbation draws for every configuration
rng(11);
pool = P.leanPool; D = randn(NREP, numel(pool));

if isempty(gcp('nocreate')), parpool(5); end
if exist(fullfile(res, ['slacknoise_' OUTTAG '.mat']), 'file')
    load(fullfile(res, ['slacknoise_' OUTTAG '.mat']), 'OUT');
else
    OUT = struct();
end
for c = CFGSEL
    ov = CFG{c, 2};
    Pc = P; f = fieldnames(ov); for k = 1:numel(f), Pc.(f{k}) = ov.(f{k}); end
    Pc = getParams(Pc, Pc.g, true, false);
    fm = cell(NREP+1, 1); cost = nan(NREP+1, 1); wall = nan(NREP+1, 1); ok = true(NREP+1, 1);
    tic;
    parfor r = 1:NREP+1
        pv = Pc;
        if r > 1
            for i = 1:numel(pool)
                nm = pool{i};
                if ~contains(nm, '__') && isfield(pv, nm) && isscalar(pv.(nm)) && pv.(nm) ~= 0
                    pv.(nm) = pv.(nm) * (1 + 1e-6*D(r-1, i));
                end
            end
            pv = getParams(pv, pv.g, true, false);
        end
        t0 = tic;
        try
            [~, ~, fmr, fdr] = runSlackExperiment(pv);
            fnS = slackTerms(pv.fn, fmr);
            cost(r) = sum(evalFeatureCost(fdr, fmr, fnS));
            fm{r} = fmr;
        catch err
            ok(r) = false; fm{r} = struct(); fprintf('cfg %d rep %d failed: %s\n', c, r, err.message);
        end
        wall(r) = toc(t0);
    end
    O = struct('name', CFG{c,1}, 'cost', cost, 'wall', wall, 'ok', ok);
    fprintf('\n=== %s  (%.0f s total, %.0f s/run) ===\n', CFG{c,1}, toc, mean(wall));
    fprintf('  slack cost: base %.4f | mean %.4f | SD %.4f (%.1f %%)\n', cost(1), mean(cost(ok)), std(cost(ok)), 100*std(cost(ok))/mean(cost(ok)));
    for q = FEAT
        X = cell2mat(cellfun(@(s) getFeat(s, q{1}), fm(ok), 'uni', 0));
        rsd = std(X, 0, 1) ./ abs(mean(X, 1));
        O.feat.(q{1}) = X;
        fprintf('  %-12s base %-34s relSD mean %6.2f %%  max %6.2f %%\n', q{1}, mat2str(round(X(1,:), 3)), 100*mean(rsd), 100*max(rsd));
    end
    OUT.(sprintf('c%d', c)) = O;
    save(fullfile(res, ['slacknoise_' OUTTAG '.mat']), 'OUT', 'CFG');
end

function fnS = slackTerms(fn, fm)
% the objective's terms that the slack run can score (grouped boundary terms excluded)
base = @(e) regexprep(regexprep(e, '\|.*$', ''), '\[.*$', '');
keep = cellfun(@(e) ~contains(e, ',') && isfield(fm, base(e)), fn);
fnS = fn(keep);
end

function v = getFeat(s, name)
if isfield(s, name), v = s.(name)(:)'; else, v = nan; end
if strcmp(name, 'rsK80') || strcmp(name, 'rsT20'), v = v(1:min(4, numel(v))); end
end
