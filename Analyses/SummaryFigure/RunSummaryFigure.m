% RUNSUMMARYFIGURE  Three-panel summary figure: force-velocity (both [ATP]),
% slack + passive, and the staircase.
%
%   Baseline parameters: params/optfull6_opt.m (the TestRunParams working point).
%
%   Runs, in this order:
%     1. Force-velocity at 8 mM MgATP          -> left panel (solid)
%     2. Force-velocity at 2 mM MgATP          -> left panel (dashed)
%     3. Slack protocol, 8 mM                  -> right-top panel
%     4. Passive (PNB/Mava, ka=0) slack        -> right-top panel, behind the sim
%     5. Staircase protocol, 8 mM              -> right-bottom panel
%
%   Everything the figure needs is cached to results/summary_cache.mat, so the
%   layout can be re-tuned without re-simulating: set RERUN=false (or just call
%   plotSummaryFigure(C, ...) on the loaded struct).
%
%   Two versions of the left panel are produced:
%     'FV' - normalised force vs velocity      -> results/summary_figure_FV.png
%     'PV' - power vs velocity (F_norm * v)    -> results/summary_figure_PV.png
%
%   NOTE ON THE LOW-ATP RUN.  In this parameterisation [MgATP] enters the model
%   only through the two ATP-gated levers ported from the lowatp-2state branch:
%       UseAtpDetachR2 / K_T2  - R2 (strong-bridge detachment) needs ATP  -> force
%       UseAtpKmsrd    / K_srx - the SRX-ADP return needs ATP             -> ktr
%   Both are normalised to 1 at MgATP_ref = 8 mM, so the 8 mM run below is
%   bit-for-bit the untouched optfull6 baseline. The Km values are the
%   lowATP_TNF_best fit, which was identified on the ThursdayNightFever
%   baseline, NOT re-fitted on optfull6 - treat the 2 mM curve as "the low-ATP
%   mechanism applied to this working point", not as a fresh low-ATP fit.

clear; close all;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath('.'));

resDir = fullfile('Analyses', 'SummaryFigure', 'results');
if ~isfolder(resDir); mkdir(resDir); end
cacheFile = fullfile(resDir, 'summary_cache.mat');

RERUN   = false;         % true -> re-simulate; false -> reuse results/summary_cache.mat
FIGSIZE = [3.8 2.9];     % INCHES. 3.8 in = 45% of an 8.5 in page: drop in at
                         % 1/3-1/2 column with no rescaling. Fonts are fixed in
                         % pt, so enlarge here rather than stretching the export.

% ---- Low-ATP levers (see header note) --------------------------------------
LOWATP.MgATP  = 2;         % mM
LOWATP.K_T2   = 1.9676;    % Km(ATP) of the R2 ATP-detachment step -> force amplitude
LOWATP.K_srx  = 13.0;      % Km(ATP) of the SRX-ADP return         -> ktr

% The optimiser's grid -[0 .5 1 2 3 4 5 6] plus intermediate points, so the model
% draws a curve rather than a polyline. Nothing below 0.5 ML/s: there the model is
% NOT monotonic (at 0.25 ML/s it returns 1.23x isometric). That band was never in
% the cost, and the FV protocol's isometric reference is a fresh SL=2.0 run while
% every shortening point is the end of a 2.2->2.0 ramp, so the two are not the
% same state. Plotting it would show an artefact, not a prediction.
FV_VEL = -[0 0.5 0.75 1 1.5 2 2.5 3 4 5 6];

% The staircase datatable stores force NORMALISED by the pre-stretch isometric
% Fss (DataCuration/CreateStairsProtocolVelocityTable.m). To put the panel back
% in kPa we recover Fss from the same raw merged recording, over the same
% window that script used. Keep these two in step with that config table.
STAIRS_RAW     = 'data/03 27 2026 M/02_Merged_8mM_Active.txt';
STAIRS_FSS_WIN = [72.30 72.45];

%% ===================== SIMULATE (or load the cache) ========================
if RERUN || ~isfile(cacheFile)
    refreshPool(5);        % model files were edited -> workers must be refreshed

    params0 = getParams();
    optfull6_opt;                       %#ok<*NASGU> % params/optfull6_opt.m

    params0.PlotEachSeparately = false; % this script owns the plotting
    params0.PlotFullscreen     = false;
    params0.EvalFeatures       = true;
    params0.BreakOnODEUnstable = false;
    params0.MaxRunTime         = 320;
    params0.RunSlackSegments   = 'AllPar';
    params0.FV_velocities      = FV_VEL;
    params0.velocitytableonfile         = 'protocol_03_27_2026_8mM_slack.mat';
    params0.stairs_velocitytableonfile  = 'protocol_03_27_2026_velocitytable_stairs.mat';
    params0.passive_velocitytableonfile = 'protocol_03_27_2026_ActivePNBMava_slack.mat';

    LoadData;    % -> ATP_c = [8 4 2], Data_ATP = [v, F8, F4, F2]

    C = struct();
    C.meta = struct('baseline', 'params/optfull6_opt.m', 'lowatp', LOWATP, ...
                    'FV_velocities', FV_VEL, 'created', datestr(now)); %#ok<TNOW1,DATST>

    %% --- 1. Force-velocity, 8 mM (baseline, ATP gates inert) ---------------
    fprintf('[1/5] FV  8 mM ... '); t0 = tic;
    pHi = params0;
    [~, ~, fmHi] = runFVExperiment(pHi, 8, Data_ATP);
    C.fv.hi = struct('v', fmHi.FV_v(:), 'f', fmHi.FV_f(:), 'fnorm', fmHi.FV_fnorm(:));
    fprintf('%.0f s\n', toc(t0));

    %% --- 2. Force-velocity, 2 mM (ATP gates active) ------------------------
    fprintf('[2/5] FV  2 mM ... '); t0 = tic;
    pLo = params0;
    pLo.UseAtpDetachR2 = true;  pLo.K_T2  = LOWATP.K_T2;
    pLo.UseAtpKmsrd    = true;  pLo.K_srx = LOWATP.K_srx;
    pLo.MgATP_ref      = 8;
    [~, ~, fmLo] = runFVExperiment(pLo, LOWATP.MgATP, Data_ATP);
    C.fv.lo = struct('v', fmLo.FV_v(:), 'f', fmLo.FV_f(:), 'fnorm', fmLo.FV_fnorm(:));
    fprintf('%.0f s\n', toc(t0));

    % Measured FV (Baker): col 1 = velocity (ML/s), col 2 = 8 mM, col 4 = 2 mM.
    C.fvdata.hi = struct('v', Data_ATP(:,1), 'f', Data_ATP(:,2), 'fnorm', Data_ATP(:,2)/Data_ATP(1,2));
    C.fvdata.lo = struct('v', Data_ATP(:,1), 'f', Data_ATP(:,4), 'fnorm', Data_ATP(:,4)/Data_ATP(1,4));

    %% --- 3. Slack protocol, 8 mM ------------------------------------------
    fprintf('[3/5] slack     ... '); t0 = tic;
    [~, out_slack] = runSlackExperiment(params0);
    Sd = load(['data/' params0.velocitytableonfile]);      % full, untrimmed data
    C.slack = struct('t', out_slack.t(:), 'F', out_slack.Force(:), ...
                     'Fpas', out_slack.FXBPassive(:), ...
                     'dt', Sd.datatable(:,1), 'dML', Sd.datatable(:,2), 'dF', Sd.datatable(:,3), ...
                     'tlim', [Sd.velocitytable(2,1)-0.05, Sd.velocitytable(end,1)]);
    fprintf('%.0f s\n', toc(t0));

    %% --- 4. Passive (PNB/Mava, ka = 0) ------------------------------------
    fprintf('[4/5] passive   ... '); t0 = tic;
    [~, ~, out_pas] = runPassiveExperiment(params0);
    C.passive = struct('t', out_pas.t(:), 'F', out_pas.Force(:), ...
                       'dt', out_pas.datatable(:,1), 'dF', out_pas.datatable(:,3));
    fprintf('%.0f s\n', toc(t0));

    %% --- 5. Staircase, 8 mM -----------------------------------------------
    fprintf('[5/5] stairs    ... '); t0 = tic;
    [~, out_stairs] = runStairsExperiment(params0);
    St = load(['data/' params0.stairs_velocitytableonfile]);
    % Same normalisation runStairsExperiment uses for its cost: model isometric
    % force interpolated at the onset of the first commanded length change.
    i_move = find(St.velocitytable(:,2) ~= 0, 1);
    Fm_ss  = interp1(out_stairs.t, out_stairs.Force, St.velocitytable(i_move,1), 'linear', 'extrap');
    if Fm_ss == 0; Fm_ss = 1; end
    % Un-normalise the DATA back to kPa with the Fss the curation divided by.
    Mraw    = readmatrix(STAIRS_RAW);
    Fss_dat = mean(Mraw(Mraw(:,1) >= STAIRS_FSS_WIN(1) & Mraw(:,1) < STAIRS_FSS_WIN(2), 3));
    fprintf('(data Fss = %.1f kPa, model Fss = %.1f kPa) ', Fss_dat, Fm_ss);
    C.stairs = struct('t', out_stairs.t(:), 'F', out_stairs.Force(:), ...
                      'Fnorm', out_stairs.Force(:)/Fm_ss, 'Fss_model', Fm_ss, ...
                      'SL', out_stairs.SL(:), ...
                      'dt', St.datatable(:,1), 'dML', St.datatable(:,2), ...
                      'dF', St.datatable(:,3)*Fss_dat, 'Fss_data', Fss_dat, ...
                      'tlim', [St.datatable(1,1), St.datatable(end,1)]);
    fprintf('%.0f s\n', toc(t0));

    save(cacheFile, 'C');
    fprintf('cached -> %s\n', cacheFile);
else
    load(cacheFile, 'C');
    fprintf('loaded cache %s (created %s)\n', cacheFile, C.meta.created);
end

%% ============================== PLOT =======================================
% PNG at 600 dpi (small figure -> needs the density) plus a vector PDF for print.
for v = {'FV', 'PV'}
    fh = plotSummaryFigure(C, v{1}, [], FIGSIZE);
    exportgraphics(fh, fullfile(resDir, ['summary_figure_' v{1} '.png']), 'Resolution', 600);
    exportgraphics(fh, fullfile(resDir, ['summary_figure_' v{1} '.pdf']), 'ContentType', 'vector');
    savefig(fh, fullfile(resDir, ['summary_figure_' v{1} '.fig']), 'compact');
    
end

fprintf('figures (%.1f x %.1f in) -> %s\n', FIGSIZE(1), FIGSIZE(2), resDir);
