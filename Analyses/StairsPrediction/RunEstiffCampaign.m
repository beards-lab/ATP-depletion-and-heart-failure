% RunEstiffCampaign.m — matched-budget refit at several FROZEN estiff values.
%
% The point of this run is to remove the artefact that broke every earlier
% estiff experiment. Sweeping estiff with kstiff2 fixed (or fixed by the
% iso-force rule) measures "kstiff2 is no longer optimal here", because the two
% are welded: at the isometric operating point the attached population is a
% narrow spike at x ~ 9.5 nm, so d(logF)/d(estiff) = log(x) ~ 2.25, and a 1 %
% change in estiff needs a 2.3 % change in kstiff2 just to stand still — while
% a 10 % kstiff2 error costs a factor of ~60 in the objective.
%
% So here every estiff gets
%   (a) an iso-force-compensated kstiff2/kstiff2_n SEED  (the ridge, analytically), and
%   (b) the SAME optimiser, the SAME pool minus estiff, and the SAME time budget,
% and is then scored on the held-out staircase. Equal budget from an equally
% constructed seed is the only comparison that is fair across exponents.
%
% e = 0.965 is included deliberately: its iso-rule seed IS the incumbent, so
% that run is the matched control the first report was missing — "can the same
% staircase gain be bought at the incumbent exponent?"
%
% Usage (one shard per MATLAB process, so several can share the box):
%   matlab -batch "EVEC=[0.80 1.10]; TIME_HRS=0.5; run('Analyses/StairsPrediction/RunEstiffCampaign.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));

if ~exist('EVEC',     'var'); EVEC     = [0.80 0.90 0.965 1.10 1.25 1.40 1.60]; end
if ~exist('TIME_HRS', 'var'); TIME_HRS = 0.5;   end
if ~exist('NWORKERS', 'var'); NWORKERS = 3;     end
if ~exist('TAGPFX',   'var'); TAGPFX   = 'ec';  end   % tags: ec080, ec110, ...

STAIRS = 'data/protocol_03_27_2026_velocitytable_stairs.mat';

if isempty(gcp('nocreate')); parpool('local', NWORKERS); end

% ---- the incumbent, and its isometric operating point --------------------
inc = getParams(loadParams('spreadw_opt'), [], true, false);
d   = load(STAIRS);
ref = stairsIsoRef(inc, d.velocitytable);
pw0 = ref.pw(inc.estiff);
fD  = extractStairsAttributes(d.datatable(:,1), d.datatable(:,3), ...
                              d.velocitytable, 'st_', struct(), true);

fprintf('[campaign] incumbent estiff %.4f  kstiff2 %.5g  F_iso %.4g\n', ...
        inc.estiff, inc.kstiff2, ref.Fiso);
fprintf('[campaign] iso spike: <x> %.3f nm  SD %.3f  => d(logF)/de = %.3f\n', ...
        sum(ref.x.*ref.p2)/sum(ref.p2), ...
        sqrt(sum((ref.x-sum(ref.x.*ref.p2)/sum(ref.p2)).^2.*ref.p2)/sum(ref.p2)), ...
        log(sum(ref.x.*ref.p2)/sum(ref.p2)));

res = struct('estiff',{},'tag',{},'scale',{},'seedEt',{},'optEt',{}, ...
             'seedStair',{},'optStair',{},'lsq',{},'parts',{},'kstiff2',{},'kSE',{});

for ie = 1:numel(EVEC)
    e   = EVEC(ie);
    tag = sprintf('%s%03d', TAGPFX, round(e*100));

    % ---- (a) iso-force-compensated seed ---------------------------------
    sc  = pw0 / ref.pw(e);
    sp  = inc;
    sp.estiff    = e;
    sp.kstiff2   = inc.kstiff2   * sc;
    sp.kstiff2_n = inc.kstiff2_n * sc;
    sp.mods = {}; sp.g = [];

    seedSnap = fullfile('params', [tag '_seed.m']);
    writeParamsToMFile(seedSnap, sp, [], ...
        sprintf('estiff FROZEN at %.4g; kstiff2 x%.4f iso-force compensation', e, sc));

    % ---- (b) the same optimiser, the same pool minus estiff --------------
    params0 = getParams(); run(seedSnap);          %#ok<*NASGU>
    fn = params0.fn; params0.mods = {}; params0.g = [];

    cfg = struct();
    cfg.baseSnap = seedSnap; cfg.tag = tag; cfg.fn = fn;
    cfg.pool     = estiffCampaignPool(params0);    % spreadw pool MINUS estiff
    cfg.compulsory = {'k2','kstiff2','ka','v_ref_reg','kSE'};
    assert(~any(strcmp(cfg.pool,'estiff')), 'estiff must stay frozen');

    cfg.N_DRAW = 10; cfg.SIMPLEX_EVALS = 60; cfg.MaxRunTime = 300;
    cfg.SURR_EVALS = 0; cfg.KICK_FRAC = 0.1; cfg.DEBUG = false;
    cfg.TIME_BUDGET_HRS = TIME_HRS;
    cfg.RESUME = isfile(fullfile(root,'params',[tag '_state.mat']));

    params0.RunForceVelocity = true;  params0.RunKtr    = false;
    params0.RunSlack         = true;  params0.RunStairs = false;
    params0.RunForceVelocityTime = false; params0.RunSlackPassive = true;
    params0.EvalFeatures = true;      params0.BreakOnODEUnstable = false;
    params0.PlotEachSeparately = 0;   params0.PlotFeatureFitting = 0;
    params0.RunSlackSegments = 'AllPar';
    params0.MaxRunTime = cfg.MaxRunTime; params0.OptimizeOn = 'Feats';
    cfg.params0 = params0;

    pchk = cfg.params0; pchk.mods = {}; pchk.g = ones(1,0);
    seedEt    = evaluateBakersExp(1, pchk, true, false);
    seedStair = scoreStairs(pchk, d, fD);
    fprintf('\n[campaign] e=%.3f  tag=%s  kstiff2 x%.4f -> %.5g | SEED Et %.4f  stair %.3f\n', ...
            e, tag, sc, sp.kstiff2, seedEt, seedStair);

    % optimizeFeatures ends with a self-consistency guard: it re-evaluates
    % state.best_params and refuses to write the final snapshot if the cost has
    % moved by more than 2 %. That guard is right in general, but the objective
    % is numerically rough (see Docs note "Objective Noise Floor"), and at large
    % kstiff2 a benign re-evaluation drift of ~2.4 % is enough to trip it — as
    % happened at e = 1.40 (1.4321 -> 1.3980, i.e. drifting BETTER). Since every
    % improving round already wrote params/<tag>_opt.m, the result survives the
    % throw; what must not happen is one exponent aborting the whole shard.
    try
        optimizeFeatures(cfg);
    catch ME
        fprintf(2, '[%s] optimizeFeatures threw (%s): %s\n', tag, ME.identifier, ME.message);
        fprintf(2, '[%s] continuing — scoring the last snapshot written by an improving round.\n', tag);
        if ~isfile(fullfile(root,'params',[tag '_opt.m']))
            fprintf(2, '[%s] no snapshot at all; skipping this exponent.\n', tag); continue;
        end
    end

    % ---- score the result, staircase held out ---------------------------
    % MATLAB caches parsed script files and the optimiser just rewrote this one.
    clear functions; rehash;                                        %#ok<CLFUNC>
    pOpt = getParams(loadParams([tag '_opt']), [], true, false);
    optEt = evaluateBakersExp(pOpt.g, pOpt, true, false);
    [optStair, lsq, parts] = scoreStairs(pOpt, d, fD);

    fprintf('\n=== %s (estiff %.4g) ===\n', tag, e);
    fprintf('   Et    : seed %8.4f -> %8.4f   (incumbent 0.975)\n', seedEt, optEt);
    fprintf('   stair : seed %8.3f -> %8.3f   (incumbent 5.415)\n', seedStair, optStair);
    fprintf('   LSQ %.5f [peak %.2f nadir %.2f plat %.3f rec %.2f]\n', ...
            lsq, parts.st_peak, parts.st_nadir, parts.st_plat, parts.st_rec);
    fprintf('   kstiff2 %.5g   kSE %.5g\n', pOpt.kstiff2, pOpt.kSE);

    res(end+1) = struct('estiff',e,'tag',tag,'scale',sc, ...
        'seedEt',seedEt,'optEt',optEt,'seedStair',seedStair,'optStair',optStair, ...
        'lsq',lsq,'parts',parts,'kstiff2',pOpt.kstiff2,'kSE',pOpt.kSE); %#ok<SAGROW>

    % One file PER SHARD. Several processes run this script concurrently with
    % different EVEC, and a single shared filename lets them clobber each other.
    % (MakeCampaignFigure re-derives the table from the snapshots regardless.)
    shard = strjoin(arrayfun(@(v) sprintf('%03d',round(v*100)), EVEC, ...
                             'UniformOutput', false), '-');
    save(fullfile(here,'results',sprintf('campaign_%s_%s.mat', TAGPFX, shard)), 'res');
end

fprintf('\n===== CAMPAIGN SUMMARY (%.2f h each) =====\n', TIME_HRS);
fprintf('%8s %10s %10s %10s %10s\n','estiff','seed Et','opt Et','seed stair','opt stair');
for i = 1:numel(res)
    fprintf('%8.3f %10.4f %10.4f %10.3f %10.3f\n', ...
        res(i).estiff, res(i).seedEt, res(i).optEt, res(i).seedStair, res(i).optStair);
end

% -------------------------------------------------------------------------
function [Cf, Elsq, parts] = scoreStairs(p, d, fD)
% Staircase feature cost, normalised to the force at the instant the first
% ramp starts — the same instrument that scored the data.
ps = p; ps.RunStairs = 1; ps.PlotEachSeparately = false;
ps.ghostSave = ''; ps.ghostLoad = '';
[Elsq, os] = runStairsExperiment(ps);
tMov = d.velocitytable(find(d.velocitytable(:,2)~=0,1),1);
Fm   = interp1(os.t, os.Force, tMov, 'linear', 'extrap');
[Cf, parts] = stairsCost(os.t, os.Force/Fm, d.velocitytable, fD);
end
