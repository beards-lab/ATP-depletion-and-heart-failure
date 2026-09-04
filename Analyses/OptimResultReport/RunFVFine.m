% RunFVFine.m — re-run the force-velocity protocol on a FINE velocity grid, for
% a smooth model curve against the eight velocities the data actually has.
%
% WHY THIS IS A SEPARATE SCRIPT. The fit evaluates FV at five velocities only
% (`params0.FV_velocities = [0 -0.5 -1 -2 -4]`), because every extra velocity is
% another ODE solve inside every objective evaluation. That is the right economy
% for fitting and the wrong resolution for a figure: four power points cannot
% show the shape of a curve whose peak sits between two of them.
%
% Nothing in runFVExperiment is tied to the data grid — it loops over whatever
% FV_velocities holds, and picks up data only where `ismember` finds a match
% (0, 0.5, 1, 2, 3, 4, 5, 6). So a finer grid needs no model change: the extra
% velocities are simply model-only. What a fine grid DOES break is the feature
% cost, where model and data vectors must be element-comparable — which is why
% this runs the protocol for PLOTTING and never calls evalFeatureCost.
%
% Output -> results/fv_fine.mat and results/slide_fv.png (slide-styled).
% MakeSlideFigures.m picks up fv_fine.mat automatically if it is present.

% The sweep is CACHED in results/fv_fine.mat and reused automatically whenever
% the saved grid matches the requested one, so re-plotting costs nothing.
% FORCE_RUN = true re-simulates anyway (after a refit, or a params change).
clearvars -except TAG_NEW VMAX DV FORCE_RUN VPLOT_MAX VDROP; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results');

if ~exist('TAG_NEW','var');   TAG_NEW   = 'spreadw_opt'; end
if ~exist('VMAX','var');      VMAX      = 7;    end   % simulated range
if ~exist('DV','var');        DV        = 0.25; end   % simulated resolution
% --- plotting knobs; they do NOT change what is simulated or stored ---
if ~exist('VPLOT_MAX','var'); VPLOT_MAX = 6;    end   % x-limit of the figures
% Velocities excluded from the DRAWN curve only. 0.25 ML/s is a protocol
% artifact, not a model prediction: runFVExperiment takes the isometric
% reference at SL 2.0 but starts every shortening run at SL 2.2, so at the
% slowest speed the 0.4 s ramp spends most of its time at longer, stronger
% lengths and the point lands 12 % ABOVE isometric. The value stays in the
% .mat and the .csv — this only keeps it out of a figure where it would read
% as a claim about force-velocity rather than about the protocol.
if ~exist('VDROP','var');     VDROP     = 0.25; end

MATF = fullfile(resDir,'fv_fine.mat');
haveCache = false;
if isfile(MATF) && ~(exist('FORCE_RUN','var') && FORCE_RUN)
    C = load(MATF);
    haveCache = isfield(C,'DV') && isfield(C,'VMAX') && isfield(C,'TAG_NEW') && ...
                C.DV == DV && C.VMAX == VMAX && strcmp(C.TAG_NEW, TAG_NEW);
    if haveCache
        fprintf('cached sweep matches (%s, 0 : %.2f : %g) — reusing, no model evaluated\n', ...
                C.TAG_NEW, C.DV, C.VMAX);
        vFine=C.vFine; F=C.F; Fn=C.Fn; Pw=C.Pw; vD=C.vD; Dn=C.Dn; Dp=C.Dp;
        dMu=C.dMu; dSd=C.dSd; pMu=C.pMu; pSd=C.pSd;
    else
        fprintf('cached sweep does not match the requested grid — re-simulating\n');
    end
end
if ~haveCache

vFine = 0:DV:VMAX;                       % positive = shortening speed
params0 = getParams(loadParams(TAG_NEW), [], true, false);
params0.mods = {}; params0.g = ones(1,0);
params0.FV_velocities = -vFine;          % the model's sign convention: negative
params0.EvalFeatures = true;
params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0;
params0.MaxRunTime = 900;
params0.OptimizeFVInit = true;           % reuse the isometric init across points

if isempty(gcp('nocreate')); parpool('local', 4); end
LoadData;                                                        %#ok<*NASGU>
fprintf('running FV at %d velocities (0 : %.2f : %g ML/s) ...\n', numel(vFine), DV, VMAX);
t0 = tic;
[~, outs, ~, ~] = runFVExperiment(params0, ATP_c, Data_ATP);
F = arrayfun(@(o) o.Force(end), outs);
fprintf('done in %.0f s\n', toc(t0));

if ~isfinite(F(1)) || F(1) <= 0
    error('RunFVFine:badIsometric','isometric force is %.4g — cannot normalise', F(1));
end
Fn = F / F(1);
Pw = Fn .* vFine;
nBad = sum(~isfinite(F));
if nBad; warning('RunFVFine:nonFinite','%d of %d velocities returned non-finite force', nBad, numel(F)); end
fprintf('F_iso = %.2f kPa | F/F0 at %g ML/s = %.4f | peak power %.4f at %.2f ML/s\n', ...
        F(1), VMAX, Fn(end), max(Pw), vFine(find(Pw==max(Pw),1)));
% A force-velocity relation must fall monotonically during shortening. Flag it
% if it does not -- the 5-point fit grid steps straight over the region where
% this model does not.
[mx, im] = max(Fn);
if mx > 1.001
    warning('RunFVFine:nonMonotonic', ...
        ['model F/F0 peaks at %.4f at v = %.2f ML/s, ABOVE isometric. ' ...
         'The fit grid [0 0.5 1 2 4] does not sample this region.'], mx, vFine(im));
end

%% ---- data: the two source curves, at the eight velocities they have ------
vD = [0, 0.5, 1, 2, 3, 4, 5, 6];
S1 = [56.40,   51.8120, 37.4459, 17.8025, 11.4430, 6.2643, 3.2759, 2.2120];  % Baker2022
S2 = [67.3942, 60.9885, 42.7059, 18.8539, 12.5956, 6.4426, 3.5136, 1.8556];  % Isovel.2021
Dn = [S1/S1(1); S2/S2(1)];               % each normalised by its OWN isometric
Dp = Dn .* vD;
dMu = mean(Dn,1); dSd = std(Dn,0,1);
pMu = mean(Dp,1); pSd = std(Dp,0,1);

save(MATF, 'vFine','F','Fn','Pw','vD','Dn','Dp','dMu','dSd','pMu','pSd','TAG_NEW','DV','VMAX');
end

%% ---- store the sweep in a plain-text form too ---------------------------
% The .mat is the working cache; this is the copy that survives MATLAB, can be
% opened in anything, and keeps EVERY simulated velocity including the ones the
% figures leave out.
CSVF = fullfile(resDir,'fv_fine.csv');
Tv = table(vFine(:), F(:), Fn(:), Pw(:), ...
           'VariableNames', {'v_ML_per_s','F_kPa','F_over_F0','norm_power'});
writetable(Tv, CSVF);
fprintf('stored %d simulated velocities -> %s and fv_fine.mat\n', numel(vFine), CSVF);

%% ---- what actually gets drawn -------------------------------------------
mM = vFine <= VPLOT_MAX + 1e-9;                       % model, within the x-limit
for vv = VDROP(:)'; mM = mM & abs(vFine - vv) > 1e-9; end   % ... minus exclusions
mD = vD <= VPLOT_MAX + 1e-9;                          % data, within the x-limit
if ~isempty(VDROP)
    fprintf('plotting %d of %d simulated velocities (<= %g ML/s, excluding %s)\n', ...
        nnz(mM), numel(vFine), VPLOT_MAX, strjoin(compose('%g', VDROP(:)'), ', '));
end

%% ---- figure (same slide style as MakeSlideFigures.m) --------------------
FIG_W = 4.6; FIG_H = FIG_W*2/3; DPI = 300;
FS = 14; FSL = 15; FST = 15; FSLG = 13;
cDat = [0.10 0.12 0.16]; cMod = [0.78 0.10 0.14];

fig = figure('Units','inches','Position',[1 1 FIG_W FIG_H],'Color','w');
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
fill([vD(mD) fliplr(vD(mD))], [pMu(mD)+pSd(mD) fliplr(pMu(mD)-pSd(mD))], cDat, ...
     'FaceAlpha',0.16, 'EdgeColor','none', 'DisplayName',['data ' char(177) ' SD']);
plot(vD(mD), pMu(mD), 'o', 'Color',cDat, 'MarkerFaceColor',cDat, 'MarkerSize',7, ...
     'LineStyle','none', 'DisplayName','data');
plot(vFine(mM), Pw(mM), '-', 'Color',cMod, 'LineWidth',2.8, 'DisplayName','model');
xlim([-0.2 VPLOT_MAX+0.2]); ylim([-0.03 max([Pw(mM) pMu(mD)+pSd(mD)])*1.12]);
xlabel('shortening velocity (ML/s)'); ylabel('normalised power');
title('Force–velocity (power)');
legend('Location','northeast');
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(resDir,'slide_fv.png'), 'Resolution',DPI);

% companion: normalised FORCE on the same grid
fig = figure('Units','inches','Position',[1 1 FIG_W FIG_H],'Color','w');
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
fill([vD(mD) fliplr(vD(mD))], [dMu(mD)+dSd(mD) fliplr(dMu(mD)-dSd(mD))], cDat, ...
     'FaceAlpha',0.16, 'EdgeColor','none', 'DisplayName',['data ' char(177) ' SD']);
plot(vD(mD), dMu(mD), 'o', 'Color',cDat, 'MarkerFaceColor',cDat, 'MarkerSize',7, ...
     'LineStyle','none', 'DisplayName','data');
plot(vFine(mM), Fn(mM), '-', 'Color',cMod, 'LineWidth',2.8, 'DisplayName','model');
yline(1, ':', 'Color',[0.45 0.45 0.45], 'LineWidth',1.2, 'HandleVisibility','off');
xlim([-0.2 VPLOT_MAX+0.2]); ylim([-0.03 max(1.08, max(Fn(mM))*1.06)]);
xlabel('shortening velocity (ML/s)'); ylabel('F / F_0');
title('Force–velocity (force)');
legend('Location','northeast');
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(resDir,'slide_fv_force.png'), 'Resolution',DPI);

fprintf('wrote results/fv_fine.mat, slide_fv.png, slide_fv_force.png\n');

%% ======================= helpers ========================================
function iStyleAx(ax, FS, FSL, FST, FSLG)
    set(ax, 'FontSize',FS, 'LineWidth',1.1, 'GridAlpha',0.16, 'Box','on', ...
            'TickDir','out', 'TickLength',[0.015 0.015]);
    ax.XLabel.FontSize = FSL;  ax.XLabel.FontWeight = 'bold';
    ax.YLabel.FontSize = FSL;  ax.YLabel.FontWeight = 'bold';
    ax.Title.FontSize  = FST;  ax.Title.FontWeight  = 'bold';
    lg = get(ax,'Legend');
    if ~isempty(lg); lg.FontSize = FSLG; lg.Box = 'on'; end
    ax.LooseInset = max(ax.TightInset, [0.02 0.02 0.01 0.01]);
end
