% RunExampleHighATP.m
% =========================================================================
% STARTER DRIVER — the HIGH-ATP (8 mM) protocol battery.
%
% This is the "hello world" of the cross-bridge model. It loads a frozen
% parametrization, points it at one day's 8 mM experimental recordings, runs
% four mechanical protocols against them, and reports how well the model fits
% — both as a scalar cost and as a per-feature plot.
%
% Run it as a whole (F5), or step through it cell by cell (Ctrl+Enter) — the
% cells are numbered [0]..[7] and are meant to be read in order.
%
%   THE FOUR PROTOCOLS
%   ------------------
%   SLACK           Repeated shorten / hold-short / re-stretch / hold cycles.
%                   The richest protocol: most of the fitted features come
%                   from here (force levels, recovery rates, restretch peaks).
%   KTR             Rate of force redevelopment after a shorten-restretch
%                   pair. Probes how fast cross-bridges re-attach.
%   STAIRS          Staircase of small stretches ("ramp up"). Probes the
%                   passive + short-range elastic response.
%   FORCE-VELOCITY  Isotonic-equivalent: steady force at each of a set of
%                   constant shortening velocities. The classic Hill curve.
%
% See Example/README.md for the physiology, the file formats, and the
% troubleshooting list. See Example/RunExampleProtocol.m for how to drive the
% model with a protocol of your own design.
%
% See also: Model/RunBakersExp.m, Model/getParams.m, Auxiliary/loadParams.m
% =========================================================================


%% [0] ENVIRONMENT ---------------------------------------------------------
% Everything in this project is addressed relative to the REPOSITORY ROOT —
% the protocol runners do `load(['data/' filename])`, so the current directory
% must be the root, not Example/.
clear; clc; close all;
root = fullfile(fileparts(mfilename('fullpath')), '..');
cd(root);
addpath(genpath(root));
fprintf('Repository root: %s\n', pwd);

% NOTE ON THE PARALLEL POOL
% -------------------------
% MATLAB caches compiled code, and parpool workers keep their OWN copy of it.
% So after you edit anything under Model/ or params/, a running pool will
% silently keep using the OLD code. `refreshPool` clears the cache and
% restarts the pool; it is the single most common foot-gun in this project.
%
% This script runs the slack protocol SERIALLY (see cell [3]), so you do not
% need a pool for a first run and we do not start one — it just costs ~30 s.
% Uncomment once you start editing model code, or once you switch the slack
% protocol to the parallel path:
%
%   refreshPool(5);    % 5 workers; refreshPool(0) = clear cache, no pool


%% [1] PARAMETERS ----------------------------------------------------------
% Everything about the model lives in one `params` struct, and `getParams` is
% its single source of truth:
%
%   params = getParams(params, g, updateInit, updateModifiers)
%     params           partially-filled struct; missing fields get defaults
%     g                vector of MULTIPLICATIVE modifiers, applied to the
%                      fields named in params.mods (this is what the optimizer
%                      tunes). Pass [] when you are not optimizing.
%     updateInit       recompute the strain grid params.s and the initial
%                      state vector params.PU0. ALWAYS true after changing SL0.
%     updateModifiers  apply g. False here — the snapshot already has the
%                      modifiers baked in (it ends with mods={}, g=[]).
%
% getParams also resolves "expression" fields: any value that is a string
% starting with '=' is evaluated as MATLAB (e.g. params.kamh = '=0.1*kah'),
% so derived parameters stay derived instead of drifting out of sync.
%
% ALWAYS load a snapshot through `loadParams`, never `run`/`eval` it. loadParams
% executes the file in its own sandboxed workspace; running it in the base
% workspace lets a stale params0 leak fields into the new one.
%
% The path is given EXPLICITLY ('params/...') on purpose. loadParams resolves a
% bare name via which(), and this repo has snapshot names that exist in more
% than one folder — most notably Workbench/ThursdayNightFever.m shadows
% params/ThursdayNightFever.m and the two are NOT identical. An explicit path
% removes all doubt about which file you actually got.
PARAM_FILE = 'params/ExampleHighATP.m';
params0 = getParams(loadParams(PARAM_FILE), [], true, false);
fprintf('Parametrization: %s   (%d-state model, MgATP = %g mM)\n', ...
        PARAM_FILE, params0.NumberOfStates, params0.MgATP);


%% [2] DATASET -------------------------------------------------------------
% Each protocol is driven by a .mat file under data/ holding two tables:
%
%   velocitytable  [ t(s), v(ML/s), v_um(bookkeeping), L(ML) ]
%                  ROW i's VELOCITY APPLIES OVER [t_i, t_i+1]. The model is
%                  driven by COLUMN 2 ONLY — columns 3 and 4 are bookkeeping
%                  and can silently disagree with the integral of column 2.
%                  Row 1 is a warm-start hold that lets the cross-bridge
%                  distribution settle before the protocol proper begins.
%
%   datatable      [ t(s), SL(um), F ]   the measured trace.
%                  F is ABSOLUTE kPa for the slack/passive files, and
%                  NORMALISED by the pre-event isometric force for the
%                  ktr/stairs/FV2 files. The runners know which is which; you
%                  need to know it only if you build your own protocol file.
%
% ---- ACTIVE: HIGH ATP (8 mM), 03/27/2026 Male ---------------------------
params0.MgATP                       = 8;
params0.velocitytableonfile         = 'protocol_03_27_2026_8mM_slack.mat';
params0.ktr_velocitytableonfile     = 'protocol_03_27_2026_velocitytable_ktr.mat';
params0.stairs_velocitytableonfile  = 'protocol_03_27_2026_velocitytable_stairs.mat';
params0.passive_velocitytableonfile = 'protocol_03_27_2026_ActivePNBMava_slack.mat';
DATASET_LABEL = '8 mM ATP, 03/27/2026 Male';

% ---- ALTERNATIVE DATASETS — uncomment ONE block to switch ---------------
%
% (a) LOW ATP (2 mM), same fibre, same day.
%     The signal here is the RATIO to the 8 mM run, not the absolute fit:
%     force goes up ~1.18x while ktr drops to ~0.54x. Reproducing that with
%     the high-ATP parameters alone is not expected — the low-ATP levers
%     (see Analyses/) exist for exactly this.
% params0.MgATP                      = 2;
% params0.velocitytableonfile        = 'protocol_03_27_2026_2mM_slack.mat';
% params0.ktr_velocitytableonfile    = 'protocol_03_27_2026_velocitytable_ktr_2mM.mat';
% DATASET_LABEL = '2 mM ATP, 03/27/2026 Male';
%
% (b) 04/03/2026 Female, 8 mM — a second fibre, same protocol design.
% params0.MgATP                       = 8;
% params0.velocitytableonfile         = 'protocol_04_03_2026_8mM_slack.mat';
% params0.ktr_velocitytableonfile     = 'protocol_04_03_2026_velocitytable_ktr.mat';
% params0.stairs_velocitytableonfile  = 'protocol_04_03_2026_velocitytable_stairs.mat';
% params0.passive_velocitytableonfile = 'protocol_04_03_2026_ActivePNBMava_slack.mat';
% DATASET_LABEL = '8 mM ATP, 04/03/2026 Female';
%
% (c) 04/10/2026 Male 2-8, 8 mM. This day ran the ATP conditions in REVERSED
%     order (2 mM first), which is what lets you separate the real ATP effect
%     on ktr from the rundown that contaminates the force ratio. It is also
%     the only day with the FV2 experiment (see (d)).
% params0.MgATP                       = 8;
% params0.velocitytableonfile         = 'protocol_04_10_2026_8mM_slack.mat';
% params0.ktr_velocitytableonfile     = 'protocol_04_10_2026_velocitytable_ktr.mat';
% params0.stairs_velocitytableonfile  = 'protocol_04_10_2026_velocitytable_stairs.mat';
% params0.passive_velocitytableonfile = 'protocol_04_10_2026_ActivePNBMava_slack.mat';
% DATASET_LABEL = '8 mM ATP, 04/10/2026 Male 2-8';
%
% (d) FV2 — an isovelocity ramp at 2 ML/s embedded in the 04/10 recordings.
%     It is a genuine sub-Vmax force-velocity point, carved into its own file
%     so the slack files keep their strict row layout. Its datatable is
%     normalised by the pre-event isometric force, which is exactly the ktr
%     runner's convention — so FV2 is driven THROUGH runKtrExperiment:
% params0.ktr_velocitytableonfile = 'protocol_04_10_2026_velocitytable_FV2.mat';
%
% (e) LEGACY Baker set — the getParams defaults, kept for continuity with the
%     older analyses. Coarser than the 2026 driving-signal recordings.
%     Expect a LARGE ktr cost here (~700, not ~1). Emptying
%     ktr_velocitytableonfile switches that runner to its hardcoded fallback
%     protocol, which is scored against the single literature value
%     Ktr_mean(1) = 37.8 /s instead of against a measured trace. The current
%     parametrization sits near 56 /s, and the fallback squares that gap. It
%     is a different objective, not a worse fit.
% params0.velocitytableonfile        = 'bakers_slack8mM_all.mat';
% params0.ktr_velocitytableonfile    = '';    % empty -> hardcoded ktr protocol
% params0.stairs_velocitytableonfile = '';    % empty -> data/bakers_rampup8.mat
% DATASET_LABEL = 'Legacy Baker 8 mM';


%% [3] PROTOCOL SELECTION --------------------------------------------------
% Each protocol is switched on by its own boolean. The snapshot ships with ktr
% and stairs OFF (that optimization run scored only slack + FV), so we turn
% them back on here rather than editing the frozen file.
params0.RunSlack         = true;    % shorten / hold / re-stretch cycles
params0.RunKtr           = true;    % force redevelopment rate
params0.RunStairs        = true;    % staircase of small stretches
params0.RunForceVelocity = true;    % steady force vs shortening velocity

% PASSIVE (PNB + mavacamten): the same slack protocol with active attachment
% switched off (ka = 0), so only titin and the serial elastic element bear
% load. Keep this ON. It contributes no cost of its own — it contributes five
% PS_* FEATURES, and the frozen feature list in cell [6] asks for all five.
% Switching it off while the feature list is unchanged makes those five
% features MISSING, at a flat penalty of 100 each: the cost jumps by ~500 and
% looks like a physics change when it is really a bookkeeping one.
params0.RunSlackPassive  = true;

% Diagnostics / long-running extras — off for a starter run.
params0.RunForceVelocityTime  = false;   % re-simulates FV as a timecourse
params0.RunForceLengthEstim   = false;   % force-length relation estimate

% 'All'    — simulate the whole slack protocol serially. One contiguous trace,
%            no parallel pool needed. Slower, but this is the honest default.
% 'AllPar' — split into chunks across parpool workers. Much faster, but you
%            MUST have a live pool (see refreshPool in cell [0]).
params0.RunSlackSegments = 'All';

% Plotting and housekeeping.
params0.PlotEachSeparately = true;    % one tile per protocol
params0.PlotFullscreen     = false;   % false -> use nexttile in our layout
params0.ShowStatePlots     = false;
params0.ShowResidualPlots  = false;
params0.PlotProbsOnFig     = 0;
params0.EvalFeatures       = true;    % extract features + score them
params0.recalculateDataFeats = false; % use the features stored in the .mat
params0.MaxRunTime         = 320;     % s of wall clock per simulation, then abort
params0.BreakOnODEUnstable = false;

% "Ghosts" are saved reference traces the plots overlay for before/after
% comparison. Blank them so this demo neither reads nor writes stray files.
params0.ghostLoad = '';
params0.ghostSave = '';
params0.SaveBest  = false;            % otherwise the FV runner writes *_params.mat


%% [4] RUN THE BATTERY -----------------------------------------------------
% RunBakersExp is a SCRIPT, not a function: it runs in this workspace, reads
% params0 from here, and writes its results back to here. That is deliberate —
% it makes every intermediate available for poking at after the run.
%
% It leaves behind:
%   E               row vector of cost terms (see cell [5])
%   features_model  struct of features extracted from the SIMULATION
%   features_data   struct of the same features measured from the DATA
%   out, outs       raw output struct(s) from the ODE integrator
figure(1); clf;
set(gcf, 'Position', [80 80 1400 900]);   % the flow layout needs the room
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

fprintf('\nRunning the battery on: %s\n', DATASET_LABEL);
fprintf('(black = data, blue = model)\n');
tRun = tic;
RunBakersExp;
fprintf('Battery finished in %.1f s.\n', toc(tRun));

sgtitle(sprintf('%s — model vs data', DATASET_LABEL));


%% [5] READ THE COST -------------------------------------------------------
% E is POSITIONAL, and its layout depends on which protocols ran. RunBakersExp
% appends in a fixed order — FV, ktr, stairs, slack, force-length, then the
% whole feature-cost vector — so switching a protocol off silently shifts the
% index of everything after it. Never hard-code an index into E; rebuild the
% labels from the same flags, as below.
%
% Two things that surprise everyone the first time:
%   * ktr and stairs contribute a scalar MSE only. They produce NO features,
%     so they will not appear in the feature plot in cell [6].
%   * There IS a feature called 'ktr' in that plot — but it is measured from
%     the slack protocol's re-stretch recovery, not from the ktr protocol.
labels = {};
if params0.RunForceVelocity; labels{end+1} = 'FV (trace MSE)';        end
if params0.RunKtr;           labels{end+1} = 'ktr (trace MSE)';       end
if params0.RunStairs;        labels{end+1} = 'stairs (trace MSE)';    end
if params0.RunSlack
    labels{end+1} = 'slack (trace MSE)';
    if params0.EvalFitSlackOnset                % adds 2 more slack terms
        labels{end+1} = 'slack onset dt';
        labels{end+1} = 'slack onset ktr';
    end
end
if params0.RunForceLengthEstim; labels{end+1} = 'force-length';       end
if params0.EvalFeatures && isfield(params0, 'fn')
    % One term per feature-list entry. NB RunBakersExp may have rewritten
    % params0.fn during the run (AssertParams expansion), so read it back.
    for i = 1:numel(params0.fn)
        labels{end+1} = ['feat: ' params0.fn{i}]; %#ok<SAGROW>
    end
end

%
% READ THIS BEFORE JUDGING THE NUMBERS. sum(E) is NOT the quantity this
% parametrization was fitted to. The snapshot sets params0.OptimizeOn =
% 'Feats', which means the optimizer's scalar objective is the FEATURE cost
% alone (cell [6]) — the raw trace-MSE terms below are carried along for
% diagnostics only. The slack trace MSE in particular is large by
% construction: it is a per-sample squared error over a long trace, scaled by
% 20, and nothing ever tried to minimise it. A big sum(E) next to a small
% feature cost is the expected picture, not a broken fit.
fprintf('\n===== COST BREAKDOWN (diagnostic; see note above) =====\n');
if numel(labels) == numel(E)
    [~, ord] = sort(E, 'descend', 'MissingPlacement', 'last');
    for k = ord(:)'
        lbl = labels{k};
        if numel(lbl) > 52; lbl = [lbl(1:49) '...']; end
        fprintf('  %-52s %10.4g\n', lbl, E(k));
    end
else
    % Defensive: if this ever fires, a protocol changed how many terms it
    % contributes and the label list above needs updating.
    warning('Example:costLabelMismatch', ...
        'E has %d terms but %d labels — printing unlabelled.', numel(E), numel(labels));
    disp(E);
end
fprintf('  %-52s %10.4g\n', 'sum(E)  [diagnostic total]', sum(E, 'omitnan'));


%% [6] PLOT THE FEATURES ---------------------------------------------------
% A "feature" is a scalar (or one scalar per protocol cycle) distilled from a
% force trace — a peak height, a recovery rate, a steady level. Fitting
% features rather than raw traces is what keeps the optimizer honest about
% WHICH part of the trace it is getting wrong.
%
% params0.fn is the feature list: which features are scored, and how heavily.
% Its grammar (full version in Model/evalFeatureCost.m):
%
%   'field'                  score this feature, weight 1
%   'field|w'                ... with weight w
%   'field_y|field_x'        scatter y against x (e.g. 'FV_fnorm|FV_v')
%   'field_y|field_x|w'      ... with weight w
%   'field|lb-ub'            BOUNDARY: penalise only if outside [lb,ub].
%                            Scored against the SIMULATION alone — no data
%                            is read. This is how physiological plausibility
%                            is enforced on quantities nobody measured.
%   'field|lb-ub|w'          ... with weight w
%   'f1|w1,f2|lb-ub'         comma at top level = one GROUPED joint cost
%   'field[sel]'             reduce a per-cycle vector: [1] [end] [mean]
%                            [min] [max] [median] [first] [last]
%
% Penalties: each NaN costs 1; a feature the model never produced costs 100.
% A single 100 in the breakdown almost always means a missing protocol, not
% bad physics.
%
% plotFeatures draws figure 80085: one tile per fn entry with data (blue o)
% vs simulation (red x), plus a sorted "cost culprits" bar chart — read that
% chart first, it tells you where the fit is actually losing.
figure(80085); clf;
set(gcf, 'Position', [60 60 1600 950]);   % 26 tiles need the room to be legible
featCost = plotFeatures(features_data, features_model, [], params0.fn);

% THIS is the number that matters — the scalar objective the optimizer
% minimised (evaluateBakersExp returns exactly this when OptimizeOn='Feats').
% For the shipped parametrization on the 8 mM 03/27 dataset it is ~7.5.
% A value near that means your environment reproduces the reference run.
fprintf('\n==> OPTIMIZER OBJECTIVE (feature cost, costExp = 2): %.4g\n', ...
        sum(featCost, 'omitnan'));

% The same information as text, plus a raw data-vs-model table of EVERY
% feature — including the ones not currently being scored.
reportFeatureCost(features_data, features_model, params0.fn);


%% [7] ONE PROTOCOL AT A TIME ----------------------------------------------
% RunBakersExp is just a coordinator. Each protocol is an ordinary function
% you can call directly — useful when you want to iterate on one protocol
% without paying for the other three. Note the return signatures differ, and
% runPassiveExperiment's are in a DIFFERENT ORDER (features first, no cost).
%
% This cell is off by default (the battery in cell [4] already ran all of it).
% To execute it, set the flag in the command window and re-run just this cell:
%   RUN_A_LA_CARTE = true;
if ~exist('RUN_A_LA_CARTE', 'var'); RUN_A_LA_CARTE = false; end

% (The %#ok below tells the Code Analyzer this block is dead ON PURPOSE, so
%  the editor does not flag the whole cell as unreachable.)
if RUN_A_LA_CARTE   %#ok<*UNRCH>
    p = params0;
    p.PlotEachSeparately = true;
    figure(2); clf;
    tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

    % LoadData is a script too; it supplies the force-velocity reference data
    % (ATP_c, Data_ATP) and the literature ktr targets (Ktr_mean).
    LoadData;

    % Slack — cost vector, output struct, model features, data features.
    [E_slack, out_slack, fm_slack, fd_slack] = runSlackExperiment(p);

    % Ktr — scalar cost + output. Ktr_mean is only used by the fallback
    % (hardcoded) protocol; the file-driven path scores the trace instead.
    [E_ktr, out_ktr] = runKtrExperiment(p, Ktr_mean);

    % Stairs — scalar cost + output.
    [E_stairs, out_stairs] = runStairsExperiment(p);

    % Force-velocity — cost, a struct ARRAY (one entry per velocity in
    % p.FV_velocities), and the feature structs.
    [E_fv, outs_fv, fm_fv, fd_fv] = runFVExperiment(p, ATP_c, Data_ATP);

    % Passive — NOTE: features first, output last, and no cost term at all.
    [fm_pas, fd_pas, out_pas] = runPassiveExperiment(p);

    fprintf('\nà la carte costs:  slack %.4g | ktr %.4g | stairs %.4g | FV %.4g\n', ...
            sum(E_slack), E_ktr, E_stairs, E_fv);
end

fprintf('\nDone. Next: Example/RunExampleProtocol.m builds and runs a protocol of your own.\n');
