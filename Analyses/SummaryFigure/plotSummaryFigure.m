function fig = plotSummaryFigure(C, leftPanel, figNum, figSize, style)
% PLOTSUMMARYFIGURE  Three-panel summary figure from a pre-simulated cache.
%
%   FIG = PLOTSUMMARYFIGURE(C) plots the normalised force-velocity version.
%   FIG = PLOTSUMMARYFIGURE(C, 'PV') plots the power-velocity version instead.
%   FIG = PLOTSUMMARYFIGURE(C, MODE, FIGNUM) draws into a specific figure.
%   FIG = PLOTSUMMARYFIGURE(C, MODE, FIGNUM, [W H]) sets the size in INCHES.
%
%   Sized for print: the default 3.4 x 2.7 in is 40% of an 8.5 in page width,
%   so the figure can be dropped in at 1/3-1/2 column with no rescaling. All
%   type is 7-8 pt AT THAT SIZE, i.e. legible as printed. Do not enlarge the
%   figure window and re-export - re-call with a bigger figSize instead, so the
%   fonts scale with it.
%
%   Colour code, identical in every panel:
%     black       measured data
%     blue        simulation (dashed = the passive component of the sim)
%     grey        PNB/Mava passive data
%     red         SL, and the right-hand axis it lives on
%   Line style in the left panel: solid = 8 mM MgATP, dashed = 2 mM.
%
%   Layout:
%     left  (full height) - force-velocity OR power-velocity
%     right top           - slack protocol; the passive PNB/Mava protocol sits
%                           behind the dashed blue passive-sim line
%     right bottom        - staircase protocol
%
%   C is the cache struct built by RunSummaryFigure.m.
%
%   See also: RunSummaryFigure

    if nargin < 2 || isempty(leftPanel); leftPanel = 'FV'; end
    if nargin < 3 || isempty(figNum);    figNum    = 900 + strcmpi(leftPanel, 'PV'); end
    if nargin < 4 || isempty(figSize);   figSize   = [3.8 2.9]; end   % inches
    % LEFT-PANEL marker/line style (right panels are unaffected):
    %   'line' - data = open markers + connecting line, sim = line   (default)
    %   'v1'   - data = markers ONLY, no line (8 mM filled, 2 mM open); sim = line
    %   'v2'   - NO lines at all; data and sim each a distinct marker shape
    %            (data = circle, sim = triangle), 8 mM filled / 2 mM open
    if nargin < 5 || isempty(style);     style     = 'line'; end
    % Type sizes (pt) chosen to stay legible at the default figSize.
    FS.axis = 8;    % tick labels
    FS.lab  = 9;    % axis labels
    FS.tit  = 9;    % panel titles
    FS.leg  = 7;    % legends
    % LW.curve is deliberately BELOW LW.data: on the left panel the sim tracks
    % the data closely, so a fat sim line hides it completely.
    LW.data = 0.8;  LW.sim = 2;  LW.curve = 1.1;  LW.sl = 1.1;
    MS      = 3.6;

    CL.data = [0.00 0.00 0.00];   % measured
    CL.sim  = [0.00 0.30 0.85 0.6];   % simulated (blue)
    CL.pas  = [0.45 0.65 0.95];   % lighter blue: passive/PNB sim
    CL.pnb  = [0.55 0.55 0.55];   % grey: PNB/Mava passive data
    CL.sl   = [0.85 0.10 0.10];   % red: SL and its axis

    fig = figure(figNum); clf;
    set(fig, 'Color', 'w', 'Units', 'inches', ...
        'Position', [1 1 figSize], 'PaperUnits', 'inches', ...
        'PaperSize', figSize, 'PaperPosition', [0 0 figSize]);
    % Parent 1x2 layout so the LEFT-RIGHT gap ('compact') is set independently of
    % the vertical gap inside the right column. The right column is a nested 5-row
    % layout, split 3:2 (slack needs the height for five cycles across two force
    % bands; the staircase does not).
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
    rt = tiledlayout(tl, 5, 1, 'TileSpacing', 'tight', 'Padding', 'none');
    rt.Layout.Tile = 2;

    axDefaults = {'FontSize', FS.axis, 'LineWidth', 0.5, 'TickLength', [0.02 0.02], ...
                  'TickDir', 'out', 'Layer', 'top'};

    %% ---------------- LEFT: force-velocity or power-velocity --------------
    ax1 = nexttile(tl, 1); hold(ax1, 'on'); box(ax1, 'on');

    isPV = strcmpi(leftPanel, 'PV');
    if isPV
        % Power = stress x strain rate. Force here is a STRESS (kPa) and velocity
        % is a strain rate (ML/s = 1/s), so P = F*v is a power PER UNIT VOLUME.
        % We report it as P/P0, each SOURCE normalised to its OWN 8 mM peak power
        % (data by the 8 mM data peak, model by the 8 mM model peak) - so both
        % 8 mM curves top out at 1 and the 2 mM deficit reads off directly,
        % without the model's ~15% isometric-stress overshoot contaminating the
        % scale. Absolute units are kW/m^3 (= kPa*1/s exactly, ~0.94 W/kg).
        Pw   = @(s) s.f(:) .* s.v(:);
        P0d  = max(Pw(C.fvdata.hi));
        P0m  = max(Pw(C.fv.hi));
        yDat = @(s) Pw(s) / P0d;
        yMod = @(s) Pw(s) / P0m;
        ylab = 'Power ($P/P_0$, 8 mM peak)';
    else
        yDat = @(s) s.fnorm(:);
        yMod = yDat;
        ylab = 'Force (norm.)';
    end

    % Data first so the simulated lines sit on top. Marker fill encodes [ATP]
    % (filled = 8 mM, open = 2 mM); marker shape in 'v2' encodes source
    % (circle = data, triangle = sim). CL.sim carries an alpha 4th element for
    % the lines; marker face/edge need a plain RGB triplet, so strip it.
    simRGB = CL.sim(1:3);
    switch lower(style)
        case 'v1'   % data markers only (no line); sim lines as in 'line'
            hDhi = plot(ax1, C.fvdata.hi.v, yDat(C.fvdata.hi), 'o', 'Color', CL.data, ...
                        'LineWidth', LW.data + 0.3, 'MarkerSize', MS, 'MarkerFaceColor', CL.data);
            hDlo = plot(ax1, C.fvdata.lo.v, yDat(C.fvdata.lo), 'o', 'Color', CL.data, ...
                        'LineWidth', LW.data + 0.3, 'MarkerSize', MS, 'MarkerFaceColor', 'w');
            hMhi = plot(ax1, C.fv.hi.v, yMod(C.fv.hi), '-',  'Color', CL.sim, 'LineWidth', LW.curve);
            hMlo = plot(ax1, C.fv.lo.v, yMod(C.fv.lo), '--', 'Color', CL.sim, 'LineWidth', LW.curve);

        case 'v2'   % no lines anywhere; distinct marker shape per source
            % Sim is evaluated on a denser velocity grid than the data; as pure
            % markers that grid is too crowded, so subsample the sim down to just
            % the velocities the data has (one sim marker per data marker).
            vHi = C.fv.hi.v; yHi = yMod(C.fv.hi); mHi = ismembertol(vHi, C.fvdata.hi.v, 1e-6);
            vLo = C.fv.lo.v; yLo = yMod(C.fv.lo); mLo = ismembertol(vLo, C.fvdata.lo.v, 1e-6);
            hDhi = plot(ax1, C.fvdata.hi.v, yDat(C.fvdata.hi), 'o', 'LineStyle', 'none', ...
                        'Color', CL.data, 'LineWidth', LW.data + 0.3, 'MarkerSize', MS + 0.4, 'MarkerFaceColor', CL.data);
            hDlo = plot(ax1, C.fvdata.lo.v, yDat(C.fvdata.lo), 'o', 'LineStyle', 'none', ...
                        'Color', CL.data, 'LineWidth', LW.data + 0.3, 'MarkerSize', MS + 0.4, 'MarkerFaceColor', 'w');
            hMhi = plot(ax1, vHi(mHi), yHi(mHi), '^', 'LineStyle', 'none', ...
                        'Color', simRGB, 'LineWidth', LW.data + 0.3, 'MarkerSize', MS + 0.4, 'MarkerFaceColor', simRGB);
            hMlo = plot(ax1, vLo(mLo), yLo(mLo), '^', 'LineStyle', 'none', ...
                        'Color', simRGB, 'LineWidth', LW.data + 0.3, 'MarkerSize', MS + 0.4, 'MarkerFaceColor', 'w');

        otherwise   % 'line' (default, original look)
            hDhi = plot(ax1, C.fvdata.hi.v, yDat(C.fvdata.hi), 'o-', 'Color', CL.data, ...
                        'LineWidth', LW.data + 0.3, 'MarkerSize', MS, 'MarkerFaceColor', 'w');
            hDlo = plot(ax1, C.fvdata.lo.v, yDat(C.fvdata.lo), 's--', 'Color', CL.data, ...
                        'LineWidth', LW.data + 0.3, 'MarkerSize', MS, 'MarkerFaceColor', 'w');
            hMhi = plot(ax1, C.fv.hi.v, yMod(C.fv.hi), '-',  'Color', CL.sim, 'LineWidth', LW.curve);
            hMlo = plot(ax1, C.fv.lo.v, yMod(C.fv.lo), '--', 'Color', CL.sim, 'LineWidth', LW.curve);
    end

    if isPV
        yline(ax1, 0, '-', 'Color', [0.85 0.85 0.85], 'HandleVisibility', 'off');
    end
    xlabel(ax1, 'Velocity (ML/s)', 'Interpreter', 'latex', 'FontSize', FS.lab);
    ylabel(ax1, ylab, 'Interpreter', 'latex', 'FontSize', FS.lab);
    % Default (tex) interpreter with a literal en dash, so this title keeps the
    % same sans face as 'Slack'/'Staircase' instead of switching to LaTeX serif.
    title(ax1, ternary(isPV, ['Power' char(8211) 'velocity'], ...
                             ['Force' char(8211) 'velocity']), 'FontSize', FS.tit);
    set(ax1, axDefaults{:});
    xlim(ax1, [-0.15 6.3]);
    xticks(ax1, 0:2:6);
    % Legend sits top-right where both curves have already decayed; leave a little
    % headroom above the 8 mM peak (=1 for PV) for it.
    if isPV
        ylim(ax1, [0 1.22]);
    else
        ylim(ax1, [0 1.05]);
    end
    lg1 = legend(ax1, [hDhi hMhi hDlo hMlo], ...
        {'8 mM data', '8 mM sim', '2 mM data', '2 mM sim'}, ...
        'Location', 'northeast', 'FontSize', FS.leg, 'Box', 'off');
    lg1.ItemTokenSize = [10 4];

    %% ---------------- RIGHT TOP: slack + passive --------------------------
    ax2 = nexttile(rt, 1, [3 1]); hold(ax2, 'on'); box(ax2, 'on');

    yyaxis(ax2, 'right');
    hSL = plot(ax2, C.slack.dt, C.slack.dML, '-', 'Color', CL.sl, 'LineWidth', LW.sl);
    ylabel(ax2, 'SL ($\mu$m)', 'Interpreter', 'latex', 'FontSize', FS.lab);
    % Squeeze SL into the lower half so it stays clear of the force traces and
    % of the legend band opened up at the top of the left axis.
    set(ax2, 'YColor', CL.sl, 'YTick', [1.9 2.2]);
    ylim(ax2, [1.83 2.62]);

    yyaxis(ax2, 'left');
    % Passive data goes down first so its simulation draws on top of it.
    % NOTE: what is shown as "passive sim" is the separate ka=0 PNB/Mava RUN, not
    % the passive component of the active slack sim (C.slack.Fpas) - that dashed
    % trace was dropped; the two differ slightly and that is discussed in text.
    hPd = plot(ax2, C.passive.dt, C.passive.dF, '-', 'Color', CL.pnb,  'LineWidth', LW.data);
    hPm = plot(ax2, C.passive.t,  C.passive.F,  '-', 'Color', CL.pas,  'LineWidth', LW.sim);
    hSd = plot(ax2, C.slack.dt,   C.slack.dF,   '-', 'Color', CL.data, 'LineWidth', LW.data);
    hSm = plot(ax2, C.slack.t,    C.slack.F,    '-', 'Color', CL.sim,  'LineWidth', LW.sim);

    ylabel(ax2, 'Force (kPa)', 'Interpreter', 'latex', 'FontSize', FS.lab);
    xlabel(ax2, '$t$ (s)', 'Interpreter', 'latex', 'FontSize', FS.lab);
    set(ax2, axDefaults{:}, 'YColor', 'k');
    xlim(ax2, C.slack.tlim);
    ylim(ax2, [-4 138]);          % headroom above the ~99 kPa peaks for the 'north' legend
    yticks(ax2, 0:40:80);
    title(ax2, 'Slack', 'FontSize', FS.tit);
    % Two columns, not four: at this width a single row runs off the axis and
    % into the SL scale. 'north' (auto-anchored) rather than a manual Position:
    % manual figure-normalised placement is computed against the on-screen raster
    % layout and does NOT survive the reflow that vector export (PDF/.fig) does,
    % which mangled the legend there. Location presets are re-anchored by MATLAB
    % in every export path. The ylim headroom below keeps it clear of the peaks.
    lg2 = legend(ax2, [hSd hSm hPd hPm], ...
        {'data', 'sim', 'passive data', 'passive sim'}, ...
        'Location', 'north', 'FontSize', FS.leg - 0.5, 'NumColumns', 2, 'Box', 'off');
    lg2.ItemTokenSize = [6 4];

    %% ---------------- RIGHT BOTTOM: staircase -----------------------------
    ax3 = nexttile(rt, 4, [2 1]); hold(ax3, 'on'); box(ax3, 'on');

    yyaxis(ax3, 'right');
    plot(ax3, C.stairs.dt, C.stairs.dML, '-',  'Color', CL.sl, 'LineWidth', LW.sl);
    % plot(ax3, C.stairs.t,  C.stairs.SL,  '--', 'Color', CL.sl, 'LineWidth', LW.sl);
    ylabel(ax3, 'SL ($\mu$m)', 'Interpreter', 'latex', 'FontSize', FS.lab);
    set(ax3, 'YColor', CL.sl, 'YTick', [2.0 2.2], 'YTickLabel', {'2.0', '2.2'});
    ylim(ax3, [1.98 2.32]);

    yyaxis(ax3, 'left');
    % Normalised on both sides. The cache carries the recording in kPa (undone
    % from the curation's Fss), so divide it back out here to match the model's
    % own normalisation; swap both to C.stairs.dF / C.stairs.F for absolute kPa.
    hFd = plot(ax3, C.stairs.dt, C.stairs.dF/C.stairs.Fss_data, '-', 'Color', CL.data, 'LineWidth', LW.data);
    hFm = plot(ax3, C.stairs.t,  C.stairs.Fnorm,  '-', 'Color', CL.sim,  'LineWidth', LW.sim);
    ylabel(ax3, 'Force (norm.)', 'Interpreter', 'latex', 'FontSize', FS.lab);
    xlabel(ax3, '$t$ (s)', 'Interpreter', 'latex', 'FontSize', FS.lab);
    set(ax3, axDefaults{:}, 'YColor', 'k');
    xlim(ax3, C.stairs.tlim);
    % The staircase rises left-to-right, so its top-LEFT corner is empty: put the
    % legend there rather than reserving a full-width band above the trace. That
    % is what lets this panel run at 2/5 of the column height. Limits come from
    % what is actually inside the time window, not from the whole warm-started run.
    wm = C.stairs.t  >= C.stairs.tlim(1) & C.stairs.t  <= C.stairs.tlim(2);
    lo = min([C.stairs.Fnorm(wm)]);
    hi = max([C.stairs.Fnorm(wm)]);
    ylim(ax3, [lo - 0.06*(hi-lo), hi + 0.1*(hi-lo)]);
    title(ax3, 'Staircase', 'FontSize', FS.tit);
    lg3 = legend(ax3, [hFd hFm], {'data', 'sim'}, ...
        'Location', 'northwest', 'FontSize', FS.leg, 'NumColumns', 2, 'Box', 'off');
    lg3.ItemTokenSize = [6 4];

    % Push each legend flush against the TOP of its panel, into the empty band
    % above the traces (the 'Location' presets leave it sitting down on the
    % peaks). Done after a drawnow so the nested layout has resolved the axis
    % Positions. This is SAFE for the exports because both PNG and PDF are now
    % raster (Resolution, not ContentType 'vector'): they capture the on-screen
    % layout as-is, so manual placement holds. It only broke under vector reflow.
    drawnow;
    set([lg1 lg2 lg3], 'Units', 'normalized');
    lg1.Position(1) = ax1.Position(1) + ax1.Position(3) - lg1.Position(3) - 0.006;
    lg1.Position(2) = ax1.Position(2) + ax1.Position(4) - lg1.Position(4) - 0.006;
    lg2.Position(1) = ax2.Position(1) + 0.5*(ax2.Position(3) - lg2.Position(3));
    lg2.Position(2) = ax2.Position(2) + ax2.Position(4) - lg2.Position(4) - 0.004;
    lg3.Position(1) = ax3.Position(1) + 0.012;
    lg3.Position(2) = ax3.Position(2) + ax3.Position(4) - lg3.Position(4) - 0.004;
end

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end
