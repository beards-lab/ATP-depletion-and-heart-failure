% RunFVBeard2022.m -- drive the incumbent model with every measured isovelocity
% record of the Beard/Baker 2022 (Biophys J) dataset, at that record's [ATP].
%
% Model (no refit, the Paper-1 pair used by Analyses/PooledFitFigures):
%   8 mM  = robust1_seed (refit9_opt + robust reader), the incumbent core
%   2 mM  = fit A: the same core with the k2 power gate, R2 x0.54 at 2 mM
%   4 mM  = the SAME gate evaluated at 4 mM (power law, R2 x0.54^0.5 = x0.735):
%           a held-out prediction from concentration, nothing tuned at 4 mM.
% Input  = each sheet's own motor length record, as a piecewise-constant
%          velocity table (buildVelocityTable, tol 5e-4 Lo) after a 2 s
%          isometric warm-up at the record's starting length (~1.10 Lo = SL 2.2 um).
% No force scale is applied here; RunFVBeard2022Analysis applies the one
% per-preparation scale (8 mM isometric at 1.00 Lo).
%
%   matlab -batch "run('Analyses/FV_Beard2022/RunFVBeard2022.m')"
%
% Output -> results/sims.mat  (one row per sheet: t, L, F of the model on the
%           data's own time grid, plus the velocity table used)

rootP = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(rootP);
addpath(genpath(rootP));
% never let a worktree copy shadow the model (genpath picks up .claude/worktrees)
wt = strsplit(genpath(fullfile(rootP, '.claude')), pathsep); wt = wt(~cellfun(@isempty, wt));
if ~isempty(wt), rmpath(wt{:}); end
outP = fullfile(rootP, 'Analyses', 'FV_Beard2022', 'results');
if ~isfolder(outP); mkdir(outP); end
if ~exist('CORE', 'var'); CORE = 'robust1_seed'; end
if ~exist('ONLY', 'var'); ONLY = []; end            % subset of sheets (debug)

D = loadBeard2022FV();
pP = getParams(loadParams(CORE), [], true, false);
pP.mods = {}; pP.g = ones(1,0);
pP.EvalFeatures = false; pP.BreakOnODEUnstable = false; pP.MaxRunTime = 1800;
pP.PlotEachSeparately = 0; pP.PlotFeatureFitting = 0; pP.PlotFullscreen = 0;
pP.ghostSave = ''; pP.ghostLoad = ''; pP.SaveBest = 0;
pP.UseAtpPowerGates = true; pP.MgATP_ref = 8; pP.UseAtpDetachR2 = 0; pP.UseAtpKmsrd = 0;
pP.pAtp_k2 = log(0.54)/log(1/4);                     % fit A: R2 x0.54 at 2 mM

ks = 1:numel(D); if ~isempty(ONLY), ks = ONLY; end
VT = cell(numel(D), 1); TQ = cell(numel(D), 1);
for k = ks, VT{k} = buildVelocityTable(D(k).t, D(k).L, 5e-4, 2); TQ{k} = D(k).t; end
atpAll = [D.atp];
modelFcn = str2func(pP.modelFcn);
% one checkpoint file per record, so a crashed worker costs one record, not the run
ckP = fullfile(outP, 'sim_records'); if ~isfolder(ckP), mkdir(ckP); end
ckf = @(k) fullfile(ckP, sprintf('sim_%02d.mat', k));
todo = ks(~arrayfun(@(k) isfile(ckf(k)), ks));
fprintf('%d of %d records to simulate\n', numel(todo), numel(ks));
tAll = tic;
for attempt = 1:3
    if isempty(todo), break; end
    VTt = VT(todo); TQt = TQ(todo); atpT = atpAll(todo); ckT = arrayfun(ckf, todo, 'uni', 0);
    try
        parfor ii = 1:numel(todo)
            vt = VTt{ii};
            p = pP; p.MgATP = atpT(ii);
            p.SL0 = 2 * vt(1, 4); p.Slim_l = 1.4; p.Slim_r = 2.3;
            p = getParams(p, p.g, true);
            p.Velocity = vt(:, 2);
            t0 = tic; r = struct('t', TQt{ii}, 'L', [], 'F', [], 'ok', false, 'msg', '', 'tSim', NaN);
            try
                [~, o] = evaluateModel(modelFcn, vt(:, 1), p);
                tq = TQt{ii} / 1000;
                [tu, iu] = unique(o.t(:));
                r.L = interp1(tu, o.SL(iu) / 2, tq, 'linear');
                r.F = interp1(tu, o.Force(iu), tq, 'linear');
                r.ok = all(isfinite(r.F)) && tu(end) >= tq(end) - 1e-6;
            catch err
                r.msg = err.message;
            end
            r.tSim = toc(t0);
            parsave(ckT{ii}, r);
        end
    catch err
        fprintf('attempt %d: pool failure (%s); restarting the pool\n', attempt, err.message);
        delete(gcp('nocreate'));
        parpool('local', max(2, 5 - attempt));
    end
    todo = ks(~arrayfun(@(k) isfile(ckf(k)), ks));
end
S = struct('k', num2cell(ks(:)), 'tSim', [], 't', [], 'L', [], 'F', [], 'ok', false, 'msg', 'not run');
for ii = 1:numel(ks)
    if ~isfile(ckf(ks(ii))), continue; end
    r = load(ckf(ks(ii))); r = r.r;
    S(ii).t = r.t; S(ii).L = r.L; S(ii).F = r.F; S(ii).ok = r.ok; S(ii).msg = r.msg; S(ii).tSim = r.tSim;
end
fprintf('%d sims done, %.0f s wall this call (%.0f s per record on average)\n', sum([S.ok]), toc(tAll), mean([S.tSim], 'omitnan'));
for ii = 1:numel(S)
    if ~S(ii).ok, fprintf('FAILED sheet %d (%s %d mM): %s\n', S(ii).k, D(S(ii).k).prep, D(S(ii).k).atp, S(ii).msg); end
end
meta = struct('CORE', CORE, 'pAtp_k2', pP.pAtp_k2, 'tol', 5e-4, 'tWarm', 2, 'when', char(datetime('now')));
if isempty(ONLY)
    save(fullfile(outP, 'sims.mat'), 'S', 'VT', 'meta');
    fprintf('wrote %s\n', fullfile(outP, 'sims.mat'));
end

function parsave(f, r) %#ok<INUSD>
save(f, 'r');
end
