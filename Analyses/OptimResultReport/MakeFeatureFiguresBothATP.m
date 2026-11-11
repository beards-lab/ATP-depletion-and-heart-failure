% MakeFeatureFiguresBothATP.m -- one slide-sized panel PER FEATURE, both [ATP],
% in the same visual language as the force-velocity panel.
%
% WHAT IS ON A PANEL
%   x            slack cycle 1..5 (protocol order, which is also cost-vector order)
%   band         mean +- SD across the three 2026 preparations, that condition
%   marker+line  03/27, the preparation both models are fitted to
%   line         model
%
% WHY SD AND NOT SE. Both were computed; the figures draw SD. Three reasons, in
% increasing order of force:
%   1. The point drawn is ONE preparation (03/27), not the mean. An SE band is
%      the precision of a mean that is not on the figure -- putting it around a
%      single-preparation point asks the reader to compare incommensurable
%      things.
%   2. The question a reader asks here is "does the model land where real muscle
%      lands", not "is the model distinguishable from the population mean". The
%      first is answered by the spread of preparations; that is SD.
%   3. Decisive: the objective's feature weights are lambda/s^2 built from this
%      same between-preparation SD (Analyses/FeatureWeighting), so a cost of 1.0
%      MEANS one SD. Banding the figure with SE = SD/sqrt(3) would draw a
%      tolerance 42 % tighter than the one the optimiser is actually enforcing,
%      and every fit would look worse on paper than it is in the cost function.
%   Set BAND = 'se' to see the alternative; the label follows.
%
% NORMALISATION. Force-dimension features are divided by that preparation's mean
% 8 mM `steady` -- BOTH conditions by the 8 mM one, so cross-preparation scale
% (77/65/48 kPa: cross-section and damage) is removed while the within-
% preparation ATP effect survives. The models are divided by the 8 mM model's
% own mean steady, the same rule. See RunEnvelopeBothATP.m.
%
% COLOUR. Hue carries the condition (blue 8 mM, red 2 mM), exactly as in
% figs_fvfine/fig_fv_fine_bothatp. Form carries data-vs-model: dark markers on a
% thin line for the measurement (black 8 mM, dark grey 2 mM), a heavy saturated
% line for the model. The BANDS are tinted by condition rather than drawn in the
% data greys -- two greyscale bands at 18 % alpha are the same grey and cannot
% be told apart, which would defeat the point of putting both on one axis.
%
% Inputs:  results/bothatp/envelope_bothatp.mat  (RunEnvelopeBothATP.m)
%          results/fit_report_<tag>.mat          (8 mM model)
%          results/bothatp/<RUNG>/lowatp_runs.mat  (2 mM model, one ladder rung)
% Figures -> results/bothatp/<RUNG>/feat_<name>.png  (300 dpi, 4.6 x 3.07 in)
%            results/bothatp/<RUNG>/feat_contact.png (all on one sheet)

clearvars -except TAG_NEW BAND RUNG RUNGS iR; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results');

if ~exist('TAG_NEW','var'); TAG_NEW = 'spreadw_opt'; end
if ~exist('BAND','var');    BAND    = 'sd';  end     % 'sd' | 'se'
if ~exist('RUNG','var');    RUNG    = 'H4';  end     % ATP-ladder rung for 2 mM

% The DATA envelope is a property of the three preparations, not of any model,
% so it is built once and shared by every rung -- it lives one level up.
envDir = fullfile(resDir,'bothatp');
outDir = fullfile(envDir,RUNG); if ~isfolder(outDir); mkdir(outDir); end

E = load(fullfile(envDir,'envelope_bothatp.mat')); E = E.E;
L = load(fullfile(resDir, sprintf('fit_report_%s.mat', TAG_NEW)));
S = load(fullfile(outDir,'lowatp_runs.mat'));
gateTxt = strjoin(strrep(fieldnames(S.gates),'pAtp_',''), ', ');
if isempty(gateTxt); gateTxt = '(none)'; end
fprintf('%s | 2 mM: %s | gates: %s\n', RUNG, S.snap, gateTxt);

% A band whose n changes along the x axis is not the same band from end to end,
% and nothing on the figure says so. 04/03 has no cycle-5 at 2 mM, so the last
% point of every 2 mM band is a TWO-preparation mean while the rest are three.
% Surfaced here rather than left in the .mat, because a 2-prep SD is a much
% weaker statement than a 3-prep one and the panels draw them identically.
for cc = {'mM8','mM2'}
    Q = E.(cc{1}); fn = fieldnames(Q)';
    bad = fn(cellfun(@(f) numel(unique(Q.(f).n)) > 1, fn));
    for b = bad
        fprintf(2,'  NOTE %s/%s: band n varies along the cycle axis: %s\n', ...
                cc{1}, b{1}, mat2str(Q.(b{1}).n));
    end
end
M8 = L.R(2).fmodel; M2 = S.fmodel;
sc8 = mean(double(M8.steady),'omitnan');              % the models' 8 mM divisor
nC  = E.nC; xc = 1:nC;
iRef = E.iRef;

% which features the objective actually scores -- the name before the first '|'
fnScored = cellfun(@(s) strtrim(extractBefore([s '|'],'|')), L.R(2).fnNames, ...
                   'UniformOutput',false);

XLAB = 'slack cycle';   % the ladder is documented by its own two panels below

%% ---- pretty names, units, and display scaling --------------------------
% {field, title, y-axis unit, multiplier, ylim-group}  -- multiplier 1000 turns
% s into ms; features sharing a ylim-group get identical y limits.
% The y label carries the UNIT only; repeating the title down the side just
% costs width in a 4.6 in panel.
%
% NAMING. `ktr` is the rate constant of the exponential force redevelopment
% after a slack release, and nothing else -- so `A`, which is that fit's
% ASYMPTOTE (ae.A + F0, the force it redevelops TO), cannot sensibly be called a
% "ktr amplitude". It is the redeveloped force. `Am` is the same quantity read
% straight off the trace at the end of the recovery window rather than from the
% fit. And the rate itself is "slack k_tr", to keep it distinct from the k_tr
% PROTOCOL, which is a separate experiment in this project.
%
% GROUPS. `ktr` and `rsK` share limits: both are first-order recovery rates in
% s^-1, one after the slack and one after the re-stretch, and the whole point of
% looking at them together is that low ATP slows them by different factors
% (x0.54 vs x0.34). Different axes would hide that.
LBL = { 'steady',              'steady force',             'rel. 8 mM steady',  1, ''
        'A',                   'redeveloped force',        'rel. 8 mM steady',  1, ''
        'Am',                  'redeveloped force (end)',  'rel. 8 mM steady',  1, ''
        'peak1_y',             're-stretch peak',          'rel. 8 mM steady',  1, ''
        'vall_y',              'post-peak valley',         'rel. 8 mM steady',  1, ''
        'peak2',               'second peak',              'rel. 8 mM steady',  1, ''
        'vall2_dy',            'second valley depth',      'rel. 8 mM steady',  1, ''
        'rsA',                 're-stretch recovery amp',  'rel. 8 mM steady',  1, ''
        'ktr',                 'slack k_{tr}',             's^{-1}',            1, 'rate'
        'rsK',                 'rsK (recovery rate)',      's^{-1}',            1, 'rate'
        'restretchSlopeStart', 're-stretch onset slope',   'kPa s^{-1}',        1, ''
        't0',                  't_0 (k_{tr} onset)',       'ms',             1000, ''
        'peak1_t',             're-stretch peak time',     'ms',             1000, ''
        'peak1_dSL',           're-stretch peak, dSL',     'ML (1)',            1, ''
        'vall_t',              'valley time',              'ms',             1000, ''
        'vall2_t',             'second valley time',       'ms',             1000, ''
        'rsT0',                'rsT_0 (recovery delay)',   'ms',             1000, ''
        'rsR2',                'rsR^2 (recovery fit)',     'R^2 (1)',           1, ''
        'SLslack',             'SL at slack',              '\mum',              1, ''
        'v_restretch',         're-stretch velocity',      'ML s^{-1}',         1, '' };

% The two protocol panels. Their model lines are DASHED: SLslack is read off the
% model's own SL trace and v_restretch off the table that drives it, so the
% curves coincide with the data almost by construction -- a solid line would
% hide the data underneath and claim an agreement that is definitional. The dash
% lets the measured trace show through, while still exposing a model that failed
% to reach the commanded slack.
PROTO = {'SLslack','v_restretch'};

% The one panel that carries the legend for the whole folder. `steady` is the
% right host: it is the first panel, it is the feature the ATP effect is named
% by, and all four series are well separated on it, so a legend parked at the
% top crosses nothing.
LEGEND_ON = 'steady';

%% ---- style -------------------------------------------------------------
FIG_W = 4.6; FIG_H = FIG_W*2/3; DPI = 300;
FS = 14; FSL = 15; FST = 15; FSLG = 10;
C8 = [0.00 0.35 0.70];  C2 = [0.80 0.20 0.10];       % models
K8 = [0.04 0.05 0.07];  K2 = [0.44 0.45 0.48];       % data 03/27
B8 = [0.62 0.74 0.88];  B2 = [0.92 0.72 0.66];       % bands
LW = 2.8; MS = 7;
bandFld = lower(BAND);
bandTxt = upper(BAND);
if ~ismember(bandFld,{'sd','se'})
    error('MakeFeatureFiguresBothATP:badBand','BAND must be ''sd'' or ''se''');
end

if ~isfolder(outDir); mkdir(outDir); end
close all force;

%% ---- pass 1: gather every series, then settle the y limits -------------
% Two passes rather than one, because a ylim-group has to see every member
% before any of them can be drawn.
G = struct('nm',{},'j',{},'d8',{},'d2',{},'b8m',{},'b8s',{},'b2m',{},'b2s',{}, ...
           'm8',{},'m2',{},'yl',{});
for j = 1:size(LBL,1)
    nm = LBL{j,1}; mult = LBL{j,4};
    if ~isfield(E.mM8, nm) || ~isfield(E.mM2, nm)
        fprintf('  skip %s (not in the envelope)\n', nm); continue;
    end
    e8 = E.mM8.(nm); e2 = E.mM2.(nm);
    isF = e8.cls == 'F';
    g = struct('nm',nm,'j',j, ...
        'd8', mult*e8.raw(iRef,:),  'd2', mult*e2.raw(iRef,:), ...
        'b8m',mult*e8.mean, 'b8s',mult*e8.(bandFld), ...
        'b2m',mult*e2.mean, 'b2s',mult*e2.(bandFld), ...
        'm8', iModel(M8, nm, nC, isF, sc8, mult), ...
        'm2', iModel(M2, nm, nC, isF, sc8, mult), 'yl', []);
    % Limits set here rather than left to MATLAB. Only the LEGEND panel needs
    % the tall reserved strip: with four series plus two bands in a 4.6 in box,
    % 'Location','best' has nowhere good to go and lands on the data in about
    % half the features, so the legend gets a fixed parking space instead.
    % Every other panel drops to ordinary padding and spends that height on
    % the data.
    padHi = 0.12; if strcmp(nm, LEGEND_ON); padHi = 0.58; end
    g.yl = iLimits({g.b8m-g.b8s, g.b8m+g.b8s, g.b2m-g.b2s, g.b2m+g.b2s, ...
                    g.d8, g.d2, g.m8, g.m2}, 0.06, padHi);
    G(end+1) = g;                                                           %#ok<SAGROW>
end

% shared limits within a ylim-group
grp = cellfun(@(x) x, {LBL{[G.j],5}}, 'UniformOutput',false);
for u = unique(grp(~cellfun(@isempty,grp)))
    k = find(strcmp(grp, u{1}));
    yl = [min(cellfun(@(v) v(1), {G(k).yl})), max(cellfun(@(v) v(2), {G(k).yl}))];
    for i = k; G(i).yl = yl; end
    fprintf('ylim group ''%s'': %s -> [%.4g %.4g]\n', u{1}, ...
            strjoin({G(k).nm}, ', '), yl(1), yl(2));
end

%% ---- pass 2: draw -----------------------------------------------------
made = {}; rows = {};
for i = 1:numel(G)
    g = G(i); j = g.j; nm = g.nm;
    isProto = any(strcmp(PROTO, nm));
    ls8 = '-'; ls2 = '-';
    if isProto; ls8 = '--'; ls2 = ':'; end

    fig = figure('Units','inches','Position',[1 1 FIG_W FIG_H],'Color','w');
    ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    iBand(xc, g.b8m, g.b8s, B8); iBand(xc, g.b2m, g.b2s, B2);
    h = gobjects(0); lb = {};
    p = plot(xc, g.d8, 'o-','Color',K8,'MarkerFaceColor',K8,'MarkerSize',MS, ...
             'LineWidth',1.2); h(end+1) = p; lb{end+1} = 'data 8 mM';       %#ok<SAGROW>
    if ~all(isnan(g.m8))
        p = plot(xc, g.m8, ls8,'Color',C8,'LineWidth',LW); h(end+1) = p; lb{end+1} = 'model 8 mM';
    end
    p = plot(xc, g.d2, 'o-','Color',K2,'MarkerFaceColor',K2,'MarkerSize',MS, ...
             'LineWidth',1.2); h(end+1) = p; lb{end+1} = 'data 2 mM';       %#ok<SAGROW>
    if ~all(isnan(g.m2))
        p = plot(xc, g.m2, ls2,'Color',C2,'LineWidth',LW); h(end+1) = p; lb{end+1} = 'model 2 mM';
    end
    xlim([0.6 nC+0.4]); xticks(xc);
    xlabel(XLAB);
    ylabel(LBL{j,3});
    title(LBL{j,2});
    ylim(g.yl);
    % ONE legend for the whole folder, on this panel. The colour key is the
    % same on every figure here -- black/grey data, blue/red model, 8/2 mM --
    % so repeating a four-entry legend twenty-odd times would cost a fifth of
    % each panel to restate what the reader learned on the first one.
    if strcmp(nm, LEGEND_ON)
        legend(h, lb, 'Location','north','NumColumns',2);
    end
    iStyleAx(ax, FS, FSL, FST, FSLG);
    f = fullfile(outDir, sprintf('feat_%s.png', nm));
    exportgraphics(fig, f, 'Resolution',DPI);
    made{end+1} = f;                                                        %#ok<SAGROW>

    % z of the model against the band it is being judged by, per condition.
    % NOT for the protocol descriptors: their between-preparation SD is a few
    % parts in 10^4 of a commanded quantity, so a z of 6 there means the model
    % reproduces the commanded slack to within a thousandth of a micron -- the
    % opposite of what a z of 6 means everywhere else on this table.
    z8 = NaN; z2 = NaN;
    if ~isProto
        z8 = iZ(g.m8, g.b8m, g.b8s);
        z2 = iZ(g.m2, g.b2m, g.b2s);
    end
    sc = ternary(any(strcmp(fnScored, nm)), 'scored', 'not scored');
    if isProto; sc = 'protocol'; end
    rows(end+1,:) = {nm, sc, mean(g.b8m,'omitnan'), mean(g.b2m,'omitnan'), z8, z2}; %#ok<SAGROW>
end

%% ---- contact sheet, for reviewing all of them at once ------------------
n = numel(made); nc = 4; nr = ceil(n/nc);
fig = figure('Units','inches','Position',[1 1 nc*FIG_W/1.6 nr*FIG_H/1.6],'Color','w');
tl = tiledlayout(fig, nr, nc, 'Padding','tight','TileSpacing','tight');
for i = 1:n
    nexttile(tl); imshow(imread(made{i})); axis image off
end
title(tl, sprintf('features, both [ATP] - 2 mM = %s (%s) - band = mean %s %s over 3 preparations', RUNG, gateTxt, ...
                  char(177), bandTxt), 'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'feat_contact.png'), 'Resolution',150);

%% ---- the numbers -------------------------------------------------------
fprintf('\n%d feature panels (band = %s) -> %s\n', n, bandTxt, outDir);
fprintf('\n%-22s %-11s %10s %10s %8s %8s\n', ...
        'feature','objective','8 mM mean','2 mM mean','z 8 mM','z 2 mM');
for i = 1:size(rows,1)
    fprintf('%-22s %-11s %10.4g %10.4g %8.2f %8.2f\n', rows{i,:});
end
fprintf(['\nz = mean over cycles of (model - prep mean) / %s. |z| > 1 means the ' ...
         'model sits\noutside the spread of the three preparations.\n'], bandTxt);

%% ======================= helpers ========================================
function v = iModel(F, nm, nC, isF, sc, mult)
% One feature from a model's feature struct, on the cycle grid, normalised by
% the same rule the data used.
    v = nan(1,nC);
    if ~isfield(F, nm); return; end
    x = double(F.(nm)(:))';
    if numel(x) < nC; return; end
    x = x(1:nC);
    if isF; x = x / sc; end
    v = mult*x;
end

function yl = iLimits(sets, padLo, padHi)
% Common limits over every drawn series, with asymmetric padding: a little air
% below, a reserved strip above for the legend.
    v = [];
    for i = 1:numel(sets); v = [v, sets{i}(:)']; end                     %#ok<AGROW>
    v = v(isfinite(v));
    if isempty(v); yl = [0 1]; return; end
    lo = min(v); hi = max(v); r = max(hi-lo, eps);
    yl = [lo - padLo*r, hi + padHi*r];
end

function iBand(x, mu, sd, col)
    g = ~isnan(mu) & ~isnan(sd);
    if ~any(g); return; end
    fill([x(g) fliplr(x(g))], [mu(g)-sd(g) fliplr(mu(g)+sd(g))], col, ...
         'FaceAlpha',0.45,'EdgeColor','none','HandleVisibility','off');
end

function z = iZ(m, mu, sd)
% Mean z of the model against the preparation band. Guarded: for a protocol
% descriptor the between-preparation SD is ~1e-12 and an unguarded ratio
% returns a number like 5e11, which reads as a catastrophic failure when it is
% actually 0/0.
    ok = isfinite(m) & isfinite(mu) & isfinite(sd) & sd > 1e-6*max(abs(mu));
    if ~any(ok); z = NaN; return; end
    z = mean((m(ok) - mu(ok)) ./ sd(ok));
end

function o = ternary(c, a, b)
    if c; o = a; else; o = b; end
end

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
