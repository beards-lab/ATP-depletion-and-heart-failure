% MakePooledFitFigures.m -- the EMC-poster figure set, redrawn for the POOLED fit.
%
% The poster (Docs/presentations/2026-EMC.pptx, Sep 2026) showed the spreadw fit
% against ONE preparation (03/27), with the other preparations only as an SD
% band. The model is now fitted to the 3-preparation pooled target (mean
% features, UsePooledTarget = 1), so here the DATA LINE IS THE POOLED MEAN and the
% band is the between-preparation SD (n = 3) around it.
%
% Models (no re-fit, no new optimisation):
%   8 mM  robust1_seed  = refit9_opt + robust feature reader (current incumbent)
%   2 mM  fit A         = same core, k2 power gate, R2 x0.54 at 2 mM
%                         (Analyses/LowATP_Refit9K2Sweep; the physiological set)
%
% Inputs (all saved, nothing is simulated here):
%   LowATP_Refit9K2Sweep/results/bestfit_k2only.mat  slack + ktr-protocol traces
%   LowATP_Refit9K2Sweep/results/levers_ref.mat      8 mM model features
%   LowATP_Refit9K2Sweep/results/k2sweep.mat         2 mM fit-A model features
%   PooledFitFigures/results/sims.mat                fine F-V + passive (RunPooledFigureSims.m)
%   data/protocol_<day>_<cond>_slack.mat, ..._ktr*.mat  the 3 preparations
%
% Bands. Slack traces: per-preparation traces put on the pooled kPa scale (the
% same per-prep anchor + rundown lift the target uses), mean +- SD. Features:
% per-preparation values from buildPooledTarget (.Y), mean +- SD, on the
% target's units. Force features are normalised per preparation before pooling,
% so their band is the SHAPE spread; the between-heart LEVEL spread (CV 15 %)
% is not in it. F-V: the poster's construction (Baker2022 vs AB2021, two
% preparations) -- the pooled 3-prep data has only one in-prep F-V point (FV2,
% 2 ML/s, n = 2), drawn as a diamond.
%
% Run:  matlab -batch "run('Analyses/PooledFitFigures/MakePooledFitFigures.m')"
% Out:  figures/*.png, 300 dpi, poster cell = 4.6 x 3.07 in (wide panels 2 cells)

rootF = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(rootF); addpath(genpath(rootF));
here = fullfile(rootF, 'Analyses', 'PooledFitFigures'); outD = fullfile(here, 'figures');
rsK = fullfile(rootF, 'Analyses', 'LowATP_Refit9K2Sweep', 'results');

BF = load(fullfile(rsK, 'bestfit_k2only.mat'));
Lr = load(fullfile(rsK, 'levers_ref.mat'));  fm8 = Lr.fm8L;
K2 = load(fullfile(rsK, 'k2sweep.mat'), 'RK2'); fmA = K2.RK2(abs([K2.RK2.g] - 0.54) < 1e-9).fm;
SIM = load(fullfile(outD, 'sims.mat'));
PT8 = pooledTarget8mM(); PT2 = pooledTarget2mM(); FV2 = fv2Target();
evalc('PTb = buildPooledTarget(struct());');
LoadData;

%% ---- style (poster) --------------------------------------------------------
FIG_W = 4.6; FIG_H = FIG_W*2/3; FIG_W2 = 2*FIG_W + 0.14; DPI = 300;
FS = 14; FSL = 15; FST = 15; FSLG = 11;
C8 = [0.00 0.35 0.70];  C2 = [0.80 0.20 0.10];       % model 8 / 2 mM
K8 = [0.04 0.05 0.07];  K2c = [0.44 0.45 0.48];      % data  8 / 2 mM (pooled mean)
B8 = [0.62 0.74 0.88];  B2 = [0.92 0.72 0.66];       % bands
cLen = [0.42 0.53 0.62];
newFig = @(w) figure('Units', 'inches', 'Position', [1 1 w FIG_H], 'Color', 'w', 'Visible', 'off');

%% ---- pooled slack traces on the full protocol timeline ---------------------
DAYS = {'03_27_2026', '04_03_2026', '04_10_2026'}; CONDS = {'8mM', '2mM'};
Fref = PTb.scale.mean; lift = PTb.scale.lift; alpha = PTb.scale.perPrep(:);
TG = (-0.10 : 2.5e-4 : 3.0)';                          % s from the first release
TR = cell(1, 2); REL = cell(1, 2);
for c = 1:2
    V = nan(3, numel(TG));
    for p = 1:3
        S = load(fullfile('data', sprintf('protocol_%s_%s_slack.mat', DAYS{p}, CONDS{c})));
        rl = S.velocitytable(S.velocitytable(:,2) < -20, 1);
        if p == 1; REL{c} = rl; end
        V(p,:) = interp1(S.datatable(:,1) - rl(1), S.datatable(:,3)*lift(p,c)/alpha(p)*Fref, TG, 'linear', NaN);
        if c == 1 && p == 1
            SV = S.velocitytable; mlT = S.datatable(:,1) - rl(1); mlV = S.datatable(:,2)/2;
        end
    end
    TR{c} = struct('m', mean(V, 1, 'omitnan')', 's', std(V, 0, 1, 'omitnan')', 'n', sum(isfinite(V), 1)');
end
relR = REL{1} - REL{1}(1);
tRest = SV(SV(:,2) > 1, 1) - REL{1}(1); vRest = SV(SV(:,2) > 1, 2);
MOD = {BF.os8, BF.os2};

%% ---- 1. slack protocol, all five cycles (2 cells wide) ---------------------
fig = newFig(FIG_W2); ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
YL = [-5 145];
iLengthAxis(ax, mlT, mlV, cLen, FS, FSL, 0.30, YL);
iTraceBand(ax, TG, TR, {B8, B2}, {K8, K2c}, 2.0);
for c = 1:2
    o = MOD{c}; t = o.t - REL{c}(1); q = t >= TG(1) & t <= TG(end);
    plot(ax, t(q), o.F(q), '-', 'Color', iif(c == 1, C8, C2), 'LineWidth', 1.3);
end
xlim(ax, [TG(1) TG(end)]); ylim(ax, YL);
xlabel(ax, 'time from first release (s)'); ylabel(ax, 'force (kPa)');
title(ax, 'Slack protocol - all five cycles (pooled, n = 3)');
iSlackAnnot(ax, relR, tRest, vRest, mlT, mlV, 1, 11, [0.24 0.30 0.38]);
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(outD, 'slack_all.png'), 'Resolution', DPI); close(fig);

%% ---- 2. slack cycle 4 detail (2 cells wide) -------------------------------
cD = 4; WIN = [-0.100 0.400];
fig = newFig(FIG_W2); ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
mw = mlT - relR(cD) >= WIN(1) & mlT - relR(cD) <= WIN(2);
iLengthAxis(ax, 1000*(mlT(mw) - relR(cD)), mlV(mw), cLen, FS, FSL, 0.30, YL);
h = iTraceBand(ax, 1000*(TG - relR(cD)), TR, {B8, B2}, {K8, K2c}, 1.6);
hm = gobjects(1, 2);
for c = 1:2
    o = MOD{c}; t = o.t - REL{c}(cD); q = t >= WIN(1) & t <= WIN(2);
    hm(c) = plot(ax, 1000*t(q), o.F(q), '-', 'Color', iif(c == 1, C8, C2), 'LineWidth', 2.2);
end
xlim(ax, 1000*WIN); ylim(ax, YL);
xlabel(ax, 'time from slack release (ms)'); ylabel(ax, 'force (kPa)');
title(ax, sprintf('Slack cycle %d (re-stretch %.0f ML/s)', cD, vRest(cD)));
legend(ax, [h(1) hm(1) h(2) hm(2)], {'data 8 mM (mean \pm SD)', 'model 8 mM', 'data 2 mM (mean \pm SD)', 'model 2 mM'}, ...
       'Location', 'northwest', 'NumColumns', 2);
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(outD, 'slack_cycle4.png'), 'Resolution', DPI); close(fig);

%% ---- 3. ktr protocol (new since the poster; pooled, n = 3) -----------------
KT = (-0.03 : 2e-4 : 0.058)'; KTR = cell(1, 2); KD = {'03_27', '04_03', '04_10'}; SUF = {'', '_2mM'};
for c = 1:2
    V = nan(3, numel(KT));
    for p = 1:3
        S = load(fullfile('data', sprintf('protocol_%s_2026_velocitytable_ktr%s.mat', KD{p}, SUF{c})));
        V(p,:) = interp1(S.datatable(:,1) - S.velocitytable(end-1,1), S.datatable(:,3), KT, 'linear', NaN);
    end
    KTR{c} = struct('m', mean(V, 1, 'omitnan')', 's', std(V, 0, 1, 'omitnan')', 'n', sum(isfinite(V), 1)');
end
fig = newFig(FIG_W2); ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
h = iTraceBand(ax, 1000*KT, KTR, {B8, B2}, {K8, K2c}, 1.8);
OK = {BF.ok8, BF.ok2};
for c = 1:2
    o = OK{c}; q = o.t >= KT(1) & o.t <= KT(end);
    hm(c) = plot(ax, 1000*o.t(q), o.F(q), '-', 'Color', iif(c == 1, C8, C2), 'LineWidth', 2.2);
end
xlim(ax, 1000*[KT(1) KT(end)]); ylim(ax, [-0.1 1.45]);
xlabel(ax, 'time from final release (ms)'); ylabel(ax, 'F / F_{iso}');
title(ax, 'k_{tr} protocol (pooled, n = 3)');
legend(ax, [h(1) hm(1) h(2) hm(2)], {'data 8 mM', 'model 8 mM', 'data 2 mM', 'model 2 mM'}, ...
       'Location', 'northwest', 'NumColumns', 2);
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(outD, 'ktr_protocol.png'), 'Resolution', DPI); close(fig);

%% ---- 4. force-velocity + power (1 cell each) ------------------------------
vD = [0 0.5 1 2 3 4 5 6];
Dn8 = [Data_ATP(1:8,2)'; Data_ATP_AB2021(1:8,2)']; Dn8 = Dn8 ./ Dn8(:,1);
Dn2 = [Data_ATP(1:8,4)'; Data_ATP_AB2021(1:8,3)']; Dn2 = Dn2 ./ Dn2(:,1);   % AB2021 has no 4 mM column
v = SIM.vFine; Fn = {SIM.F8/SIM.F8(1), SIM.F2/SIM.F2(1)};
for c = 1:2; Fn{c}(find(Fn{c} <= 0, 1):end) = NaN; end   % each curve ends at its own V_max
for panel = 1:2
    fig = newFig(FIG_W); ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
    if panel == 1; y8 = Dn8; y2 = Dn2; else; y8 = Dn8 .* vD; y2 = Dn2 .* vD; end
    iBand(ax, vD, mean(y8, 1), std(y8, 0, 1), B8); iBand(ax, vD, mean(y2, 1), std(y2, 0, 1), B2);
    d1 = plot(ax, vD, mean(y8, 1), 'o', 'Color', K8, 'MarkerFaceColor', K8, 'MarkerSize', 6);
    d2 = plot(ax, vD, mean(y2, 1), 'o', 'Color', K2c, 'MarkerFaceColor', K2c, 'MarkerSize', 6);
    for c = 1:2
        yy = Fn{c}; if panel == 2; yy = yy .* v; end
        q = v ~= 0.25;                                         % see header of RunFVFineAllRungs.m
        pm(c) = plot(ax, v(q), yy(q), '-', 'Color', iif(c == 1, C8, C2), 'LineWidth', 2.2); %#ok<SAGROW>
    end
    if panel == 1
        dd = errorbar(ax, 2, FV2.m8, FV2.s8, 'd', 'Color', K8, 'MarkerFaceColor', B8, 'MarkerSize', 8, 'LineWidth', 1.4);
        errorbar(ax, 2.1, FV2.m2, FV2.s2, 'd', 'Color', K2c, 'MarkerFaceColor', B2, 'MarkerSize', 8, 'LineWidth', 1.4);
        ylim(ax, [0 1.05]); ylabel(ax, 'F / F_{iso}'); title(ax, 'Force-velocity');
        legend(ax, [d1 pm(1) d2 pm(2) dd], {'data 8 mM', 'model 8 mM', 'data 2 mM', 'model 2 mM', 'pooled in-prep FV2'}, ...
               'Location', 'northeast', 'NumColumns', 1);
        fn = 'fv_force.png';
    else
        ylim(ax, [0 0.75]); ylabel(ax, 'normalised power (s^{-1})'); title(ax, 'Force-velocity (power)');
        fn = 'fv_power.png';
    end
    xlim(ax, [0 6]); xlabel(ax, 'shortening velocity (ML/s)');
    iStyleAx(ax, FS, FSL, FST, FSLG);
    exportgraphics(fig, fullfile(outD, fn), 'Resolution', DPI); close(fig);
end

%% ---- 5. features over the five cycles (1 cell each) ------------------------
% name, title, y label, scale (target units -> plotted units)
FT = { 'steady',   'steady force',               'kPa',      1
       'peak1_y',  're-stretch peak',            'kPa',      1
       't0',       't_0 (k_{tr} onset)',          'ms',       1000
       'ktrfit',   'slack k_{tr}',                's^{-1}',   1
       'A',        'redeveloped force',          'kPa',      1
       'rsK80',    'rsK80 (recovery rate)',      's^{-1}',   1
       'vall_y',   'valley after peak',          'kPa',      1
       'peak2',    'second peak',                'kPa',      1
       'vall2_dy', 'post-peak drop',             'kPa',      1
       'rsTfitE',  'recovery time',              'ms',       1000
       'vall_tq',  'valley time',                'ms',       1000
       'peak1_tq', 'peak time',                  'ms',       1000 };
PTc = {PT8, PT2}; FMc = {fm8, fmA};
for i = 1:size(FT, 1)
    nm = FT{i,1}; sc = FT{i,4};
    fig = newFig(FIG_W); ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
    hh = gobjects(1, 4);
    for c = 1:2
        [mu, sd] = iFeatBand(PTb, PTc{c}, nm, c);
        if isempty(mu); continue; end
        mu = mu*sc; sd = sd*sc; x = 1:numel(mu);
        iBand(ax, x, mu, sd, iif(c == 1, B8, B2));
        hh(2*c-1) = plot(ax, x, mu, '-o', 'Color', iif(c == 1, K8, K2c), 'MarkerFaceColor', iif(c == 1, K8, K2c), ...
                         'LineWidth', 1.2, 'MarkerSize', 6);
        m = iV(FMc{c}, nm)*sc;
        if ~isempty(m); hh(2*c) = plot(ax, 1:numel(m), m, '-', 'Color', iif(c == 1, C8, C2), 'LineWidth', 2.4); end
    end
    set(ax, 'XTick', 1:5, 'XLim', [0.7 5.3]);
    xlabel(ax, 'slack cycle'); ylabel(ax, FT{i,3}); title(ax, FT{i,2});
    if i == 1
        legend(ax, hh, {'data 8 mM', 'model 8 mM', 'data 2 mM', 'model 2 mM'}, 'Location', 'southwest', 'NumColumns', 2);
    end
    iStyleAx(ax, FS, FSL, FST, FSLG);
    exportgraphics(fig, fullfile(outD, ['feat_' nm '.png']), 'Resolution', DPI); close(fig);
end

%% ---- 6. passive, Ca-sensitive (1.5 cells wide; 03/27, the only PNB+Mava set) --
P = SIM.PAS; dt = P.datatable;
rx = readmatrix(fullfile('data', '03 27 2026 M', '01_Merged_Relax.txt'));
fig = figure('Units', 'inches', 'Position', [1 1 1.5*FIG_W FIG_H*0.85], 'Color', 'w', 'Visible', 'off');
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
XLp = [74.5 78];
mlp = dt(:,2); if median(mlp) > 1.5; mlp = mlp/2; end
iLengthAxis(ax, dt(:,1), mlp, cLen, FS, FSL, 0.30, [0 30]);
h1 = plot(ax, rx(:,1), rx(:,3), '-', 'Color', [0.85 0.15 0.10], 'LineWidth', 1.6);
h2 = plot(ax, dt(:,1), dt(:,3), '-', 'Color', K8, 'LineWidth', 1.0);
q = P.t >= XLp(1) & P.t <= XLp(2);
h3 = plot(ax, P.t(q), P.F(q), '-', 'Color', C8, 'LineWidth', 1.8);
xlim(ax, XLp); ylim(ax, [0 30]);
xlabel(ax, 'time (s)'); ylabel(ax, 'force (kPa)'); title(ax, 'Ca-sensitive passive force (03/27)');
legend(ax, [h1 h2 h3], {'data: relaxed', 'data: active + PNB/Mava', 'model'}, 'Location', 'north', 'Orientation', 'horizontal');
iStyleAx(ax, FS, FSL, FST, FSLG);
exportgraphics(fig, fullfile(outD, 'passive.png'), 'Resolution', DPI); close(fig);

%% ---- 7. staircase PREDICTION (3 cells wide; one tile per preparation) --------
% Not poolable: 03/27 is 8 steps x 0.1 s, 04/03 and 04/10 are 4 steps x 0.2 s.
% 8 mM data only; the 2 mM curve is a prediction at an [ATP] never measured.
% Each trace is over its OWN isometric, so only the staircase SHAPE is compared.
fig = figure('Units', 'inches', 'Position', [1 1 3*FIG_W FIG_H], 'Color', 'w', 'Visible', 'off');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
YLS = [0.8 1.75];
for i = 1:numel(SIM.STR)
    s = SIM.STR(i); dts = s.dt; t0 = s.t0;
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
    iLengthAxis(ax, dts(:,1) - t0, dts(:,2)/2, cLen, FS, FSL, 0.14, YLS);
    h1 = plot(ax, dts(:,1) - t0, dts(:,3), '-', 'Color', [K8 0.85], 'LineWidth', 1.6);
    h2 = plot(ax, s.t - t0, s.F8, '-', 'Color', C8, 'LineWidth', 2.0);
    h3 = plot(ax, s.t2 - t0, s.F2, '--', 'Color', C2, 'LineWidth', 2.0);
    xlim(ax, [-0.10 dts(end,1) - t0]); ylim(ax, YLS);
    xlabel(ax, 'time from first step (s)'); if i == 1; ylabel(ax, 'F / F_{iso}'); end
    title(ax, sprintf('%s/%s', s.day(1:2), s.day(4:5)));
    if i == 1
        legend(ax, [h1 h2 h3], {'data 8 mM', 'model 8 mM', 'model 2 mM (prediction)'}, 'Location', 'northwest');
    end
    iStyleAx(ax, FS, FSL, FST, FSLG);
end
title(tl, 'Staircase - PREDICTION (never fitted; no 2 mM data)', 'FontSize', FST, 'FontWeight', 'bold');
exportgraphics(fig, fullfile(outD, 'stairs.png'), 'Resolution', DPI); close(fig);

fprintf('wrote %s\n', outD);

%% ======================= helpers ============================================
function y = iif(c, a, b)
    if c; y = a; else; y = b; end
end
function v = iV(st, nm)
    v = []; if isfield(st, nm); v = double(st.(nm)); v = v(:)'; end
end
function [mu, sd] = iFeatBand(PTb, PT, nm, c)
% Pooled mean as the objective sees it (PT.val, target units) and the
% between-preparation SD, put on the same units by the per-cycle ratio
% val / mean(per-prep values) -- the per-prep values (PTb.Y) are anchor-normalised.
    mu = []; sd = [];
    j = find(strcmp(PT.name, nm), 1); k = find(strcmp(PTb.name, nm), 1);
    if isempty(j) || isempty(k); return; end
    mu = PT.val(j,:); Y = squeeze(PTb.Y(:, c, k, 1:numel(mu)));
    my = mean(Y, 1, 'omitnan'); r = mu ./ my;
    bad = ~isfinite(r) | abs(my) < 1e-9; if any(bad); r(bad) = median(r(~bad)); end
    sd = std(Y, 0, 1, 'omitnan') .* abs(r);
end
function h = iTraceBand(ax, t, TR, BC, LC, lw)
    h = gobjects(1, numel(TR));
    for c = 1:numel(TR)
        m = TR{c}.m; s = TR{c}.s; ok = isfinite(m) & isfinite(s) & TR{c}.n >= 2;
        iBand(ax, t(ok)', m(ok)', s(ok)', BC{c});
        mm = m; mm(TR{c}.n < 2) = NaN;
        h(c) = plot(ax, t, mm, '-', 'Color', LC{c}, 'LineWidth', lw);
    end
end
function iBand(ax, x, mu, sd, col)
    x = x(:)'; mu = mu(:)'; sd = sd(:)';
    g = isfinite(mu) & isfinite(sd); if ~any(g); return; end
    % split into contiguous runs so a NaN gap does not draw a bridge
    d = diff([0 g 0]); st = find(d == 1); en = find(d == -1) - 1;
    for r = 1:numel(st)
        ii = st(r):en(r);
        fill(ax, [x(ii) fliplr(x(ii))], [mu(ii)-sd(ii) fliplr(mu(ii)+sd(ii))], col, ...
             'FaceAlpha', 0.45, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
function iSlackAnnot(ax, tRel, tRest, vRes, tData, mlData, xScale, fs, col)
% ML on each slack plateau + re-stretch speed next to each re-stretch (as on the poster).
    yyaxis(ax, 'right');
    for c = 1:numel(tRel)
        tS = tRel(c); tE = tRest(c); span = tE - tS;
        win = tData >= tS + 0.15*span & tData <= tE - 0.05*span;
        if ~any(win); continue; end
        ml = round(median(mlData(win)), 2);
        text(ax, xScale*(tS + tE)/2, ml + 0.02, sprintf('%.2f', ml), 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'FontSize', fs, 'Color', col, 'FontWeight', 'bold');
        text(ax, xScale*(tE + 0.030), 1.07, sprintf('%.0f ML/s', vRes(c)), 'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'middle', 'FontSize', fs, 'Color', col);
    end
    yyaxis(ax, 'left');
end
function iLengthAxis(ax, t, ml, c, FS, FSL, topFrac, leftYlim)
% Commanded length on the right axis, in the top band of the panel (MakeSlideFigures.m).
    yyaxis(ax, 'right');
    lo = min(ml); hi = max(ml); r = max(hi - lo, eps);
    plot(ax, t, ml, '-', 'Color', c, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    pad = 0.12*r; span = (r + pad)/topFrac; r0 = lo - pad;
    ylim(ax, [r0, r0 + span]); ylabel(ax, 'ML (1)');
    ax.YAxis(2).Color = c; ax.YAxis(2).FontSize = FS - 1;
    ax.YAxis(2).Label.FontSize = FSL - 1; ax.YAxis(2).Label.FontWeight = 'bold';
    ax.YAxis(2).TickValues = unique(round([lo, 1, hi], 2));
    yyaxis(ax, 'left'); ylim(ax, leftYlim);
end
function iStyleAx(ax, FS, FSL, FST, FSLG)
    if numel(ax.YAxis) > 1; yyaxis(ax, 'left'); end
    set(ax, 'FontSize', FS, 'LineWidth', 1.1, 'GridAlpha', 0.16, 'Box', 'on', ...
            'TickDir', 'out', 'TickLength', [0.015 0.015]);
    ax.YAxis(1).Color = [0.15 0.15 0.15];
    ax.XLabel.FontSize = FSL; ax.XLabel.FontWeight = 'bold';
    ax.YAxis(1).Label.FontSize = FSL; ax.YAxis(1).Label.FontWeight = 'bold';
    ax.Title.FontSize = FST; ax.Title.FontWeight = 'bold';
    lg = get(ax, 'Legend'); if ~isempty(lg); lg.FontSize = FSLG; lg.Box = 'on'; end
    ax.LooseInset = max(ax.TightInset, [0.02 0.02 0.01 0.01]);
end
