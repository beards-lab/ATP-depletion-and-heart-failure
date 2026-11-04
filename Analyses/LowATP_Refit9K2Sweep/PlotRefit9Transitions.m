% PlotRefit9Transitions  Slide-grade figure of every cross-bridge transition rate
% of the current best parameter set params/refit9_opt.m.
%
% Rates are computed EXACTLY as Model/dPUdT_CombinedTransitions.m does
% (UsePieceWiseStrainDep branch, ClampStrainDepExtrapolation), from the pchip
% interpolants that getParams Step 6 builds:
%   R1D(s) = kd  * pp_R1D(s)            knots at X_R1D + log(off_R1D)/10
%   R12(s) = k1  * pp_R12(s)            knots at X_R12  (NO log offset applied)
%   R21(s) = k_1 * pp_R21(s)            knots at X_R21 + log(off_R21)/10
%   R2(s)  = k2  * pp_R2(s + dr2 - dr)  knots at X_2 + log(off_2)/10 - (dr2 - dr)
% each held at its edge-knot value beyond the outer knots. The script then
% cross-checks these curves against the RHS itself (justPlotStateTransitionsFlag
% path, which evaluates the rates with p1 = p2 = 1 on s = -20..20 nm).
%
% Output: results/refit9_transitions.png, results/refit9_transitions.txt

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
root = fileparts(fileparts(here));
addpath(genpath(root));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

params = getParams(loadParams('refit9_opt'), [], true, false);
P = params;

%% --- strain-dependent rates (s in um internally, plotted in nm) -------------
cl  = P.ClampStrainDepExtrapolation;
ev  = @(pp, x) ppval(pp, clampx(pp, x, cl));
sG  = P.s(:);                                   % model strain grid (um)
s   = linspace(-0.075, 0.075, 1501)';           % plotting grid, beyond the model grid
sh2 = P.dr2 - P.dr;                             % R2 argument shift

R1D = P.kd  * ev(P.PieceWiseStrainDepR1D, s);
R12 = P.k1  * ev(P.PieceWiseStrainDep,    s);
R21 = P.k_1 * ev(P.PieceWiseStrainDepR21, s);
R2  = P.k2  * ev(P.PieceWiseStrainDep2,   s + sh2);

% SRX routing of the R2 flux above sSRXrip (destination switch, not a rate change)
fSRX = P.SRXFromR2 + P.SRXFromR2HighStrain * min(1, max(0, (s - P.sSRXrip) / max(P.dSRXrip, 1e-9)));
fSRX = min(1, fSRX);

% knot positions exactly as used in getParams (x) and the RHS (R2 shift)
kn.R12 = struct('x', P.PieceWiseStrainDepX(:), 'y', P.k1 * P.PieceWiseStrainDepParams(:));
kn.R1D = struct('x', P.PieceWiseStrainDepR1DX(:) + log(P.PieceWiseStrainDepR1DX_logOffset)/10, ...
                'y', P.kd * P.PieceWiseStrainDepR1DParams(:));
kn.R21 = struct('x', P.PieceWiseStrainDepR21X(:) + log(P.PieceWiseStrainDepR21X_logOffset)/10, ...
                'y', P.k_1 * P.PieceWiseStrainDepR21Params(:));
kn.R2  = struct('x', P.PieceWiseStrainDep2X(:) + log(P.PieceWiseStrainDep2X_logOffset)/10 - sh2, ...
                'y', P.k2 * P.PieceWiseStrainDep2Params(:));

%% --- cross-check against the RHS itself ---------------------------------------
pf = P; pf.justPlotStateTransitionsFlag = true; pf.Vums = 0;
hf = figure('Visible', 'off'); tiledlayout(hf, 'flow');
try
    dPUdT_CombinedTransitions(0, P.PU0(:), pf);
catch me
    if ~strcmp(me.message, 'Quitting after plotting states'), rethrow(me); end
end
ln = findobj(hf, 'Type', 'line');
getLine = @(name) ln(strcmp({ln.DisplayName}, name));
sChk = (-0.02:0.001:0.02)';
chk = {'R1D', P.kd  * ev(P.PieceWiseStrainDepR1D, sChk);
       'R12', P.k1  * ev(P.PieceWiseStrainDep,    sChk);
       'R21', P.k_1 * ev(P.PieceWiseStrainDepR21, sChk);
       'R2',  P.k2  * ev(P.PieceWiseStrainDep2,   sChk + sh2)};
maxRelErr = zeros(1, 4);
for i = 1:4
    L = getLine(chk{i,1});
    maxRelErr(i) = max(abs(L(1).YData(:) - chk{i,2}) ./ abs(chk{i,2}));
end
close(hf);
fprintf('RHS cross-check max rel err  R1D %.2g  R12 %.2g  R21 %.2g  R2 %.2g\n', maxRelErr);
assert(all(maxRelErr < 1e-10), 'Rate formulas disagree with the RHS');

%% --- force-dependent SRX rates -------------------------------------------------
F = linspace(0, 100, 401)';
kPT2SR  = P.ksr0  * exp(-F / P.sigma2);         % PT  -> SR   (SRX-ATP entry)
kSR2PT  = P.kmsr  * exp( F / P.sigma1);         % SR  -> PT   (0 in refit9)
kSRD2PD = P.kmsrd * exp( F / P.sigma_srd1);     % SRD -> PD
kPD2SRD = P.ksrd  * exp(-F / P.sigma_srd2);     % PD  -> SRD
kSR2SRD = P.ksr2srd * ones(size(F));            % SR  -> SRD
kSRD2SR = P.ksrd2sr * ones(size(F));            % SRD -> SR

%% --- figure ------------------------------------------------------------------
cFwd  = [0 39 76]/255;          % navy 00274C: forward
cRev  = [213 94 0]/255;         % vermillion D55E00: reverse / detachment-back
cAtt  = [0 158 115]/255;        % bluish green 009E73: attachment
cHyd  = [230 159 0]/255;        % orange E69F00: hydrolysis
cSky  = [86 180 233]/255;       % sky blue 56B4E9
cPur  = [204 121 167]/255;      % reddish purple CC79A7
cGrey = [0.45 0.45 0.45];
fs = 15; lw = 2.6; ms = 11;
nm = 1000;
xl = [-75 75]; yl = [1 1e4];
gridEdge = [min(sG) max(sG)] * nm;

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 12 6.75]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) attachment, A1 detachment, hydrolysis
ax = nexttile(tl); prepAx(ax, fs, xl, yl, gridEdge);
hA  = yline(ax, P.ka, '-',  'Color', cAtt, 'LineWidth', lw);
hA0 = yline(ax, P.ka * P.A0, '--', 'Color', cAtt, 'LineWidth', lw);
hH  = yline(ax, P.kah,  '-',  'Color', cHyd, 'LineWidth', lw);
hHm = yline(ax, P.kamh, '--', 'Color', cHyd, 'LineWidth', lw);
hD  = plot(ax, s*nm, R1D, '-', 'Color', cRev, 'LineWidth', lw);
plot(ax, kn.R1D.x*nm, kn.R1D.y, 'x', 'Color', cRev, 'MarkerSize', ms, 'LineWidth', 2.2);
title(ax, 'a  Attachment & A1 detachment', 'FontWeight', 'bold');
ylabel(ax, 'Rate (s^{-1})');
lab(ax, -72, P.ka*1.6,        'k_a  PD\rightarrowA1', cAtt, fs);
lab(ax, -72, sqrt(P.kah*P.ka*P.A0), 'k_{ah} hydrolysis', cHyd, fs);
lab(ax, -72, P.ka*P.A0/1.7,   'k_a\cdotA_0 isometric', cAtt, fs);
lab(ax, -72, P.kamh/1.8,     'k_{-ah}  PD\rightarrowPT', cHyd, fs);
lab(ax, 22, 2.2e3,            'R1D  A1\rightarrowPD', cRev, fs);

% (b) power stroke and its reverse
ax = nexttile(tl); prepAx(ax, fs, xl, yl, gridEdge);
h12 = plot(ax, s*nm, R12, '-', 'Color', cFwd, 'LineWidth', lw);
h21 = plot(ax, s*nm, R21, '-', 'Color', cRev, 'LineWidth', lw);
plot(ax, kn.R12.x*nm, kn.R12.y, 'x', 'Color', cFwd, 'MarkerSize', ms, 'LineWidth', 2.2);
plot(ax, kn.R21.x*nm, kn.R21.y, 'x', 'Color', cRev, 'MarkerSize', ms, 'LineWidth', 2.2);
title(ax, 'b  Power stroke A1\leftrightarrowA2', 'FontWeight', 'bold');
lab(ax, 12, 45, 'R12  A1\rightarrowA2', cFwd, fs);
lab(ax, -45, 15, 'R21  A2\rightarrowA1', cRev, fs);

% (c) A2 detachment, with the share routed to SRX
ax = nexttile(tl); prepAx(ax, fs, xl, yl, gridEdge);
h2  = plot(ax, s*nm, R2, '-', 'Color', cFwd, 'LineWidth', lw);
h2s = plot(ax, s*nm, R2 .* fSRX, '--', 'Color', cPur, 'LineWidth', lw);
plot(ax, kn.R2.x*nm, kn.R2.y, 'x', 'Color', cFwd, 'MarkerSize', ms, 'LineWidth', 2.2);
title(ax, 'c  A2 detachment', 'FontWeight', 'bold');
lab(ax, -40, 15, 'R2  A2\rightarrowPT', cFwd, fs);
lab(ax, 20, 150, {'to SRX', sprintf('(%.0f%% of R2)', 100*P.SRXFromR2HighStrain)}, cPur, fs);

% (d) all strain-dependent rates together
ax = nexttile(tl); prepAx(ax, fs, xl, yl, gridEdge);
g1 = plot(ax, s*nm, R1D, '-', 'Color', cRev, 'LineWidth', lw);
g2 = plot(ax, s*nm, R12, '-', 'Color', cFwd, 'LineWidth', lw);
g3 = plot(ax, s*nm, R21, ':', 'Color', cRev, 'LineWidth', lw+0.4);
g4 = plot(ax, s*nm, R2,  '-', 'Color', cSky, 'LineWidth', lw);
g5 = yline(ax, P.ka * P.A0, '--', 'Color', cAtt, 'LineWidth', lw);
title(ax, 'd  All strain-dependent rates', 'FontWeight', 'bold');
xlabel(ax, 'Strain s (nm)'); ylabel(ax, 'Rate (s^{-1})');
lg = legend(ax, [g1 g2 g3 g4 g5], {'R1D', 'R12', 'R21', 'R2', 'k_a\cdotA_0'}, ...
    'Location', 'south', 'NumColumns', 5, 'FontSize', fs-4, 'Box', 'off');
lg.ItemTokenSize = [14 18];

% (e) SRX-ATP exchange vs force
ax = nexttile(tl); prepAxF(ax, fs);
e1 = plot(ax, F, kPT2SR, '-', 'Color', cFwd, 'LineWidth', lw);
if P.kmsr > 0
    e2 = plot(ax, F, kSR2PT, '-', 'Color', cRev, 'LineWidth', lw);
    legE = {'PT\rightarrowSRX', 'SRX\rightarrowPT'}; hE = [e1 e2];
else
    hE = e1; legE = {sprintf('PT\\rightarrowSRX  k_{sr0}e^{-F/%.1f}', P.sigma2)};
    text(ax, 50, 1.6, 'SRX\rightarrowPT direct exit k_{-sr} = 0', 'Color', cRev, ...
        'FontSize', fs-2, 'HorizontalAlignment', 'center');
end
title(ax, 'e  SRX-ATP entry', 'FontWeight', 'bold');
xlabel(ax, 'Force F (kPa)');
legend(ax, hE, legE, 'Location', 'north', 'FontSize', fs-3, 'Box', 'off');

% (f) SRX-ADP exchanges vs force
ax = nexttile(tl); prepAxF(ax, fs);
f1 = plot(ax, F, kPD2SRD, '-',  'Color', cFwd, 'LineWidth', lw);
f2 = plot(ax, F, kSRD2PD, '-',  'Color', cRev, 'LineWidth', lw);
f3 = plot(ax, F, kSR2SRD, '--', 'Color', cSky, 'LineWidth', lw);
f4 = plot(ax, F, kSRD2SR, '--', 'Color', cPur, 'LineWidth', lw);
title(ax, 'f  SRX-ADP exchanges', 'FontWeight', 'bold');
xlabel(ax, 'Force F (kPa)');
lab(ax, 55, 4.5,  'PD\rightarrowSRD', cFwd, fs);
lab(ax, 5, 250,   'SRD\rightarrowPD', cRev, fs);
lab(ax, 55, 240,  'SRX\rightarrowSRD', cSky, fs);
lab(ax, 55, 17,   'SRD\rightarrowSRX', cPur, fs);

outPng = fullfile(resDir, 'refit9_transitions.png');
exportgraphics(fig, outPng, 'Resolution', 200, 'BackgroundColor', 'white');
close(fig);
fprintf('Saved %s\n', outPng);

%% --- caption ------------------------------------------------------------------
at0 = @(pp, x) ppval(pp, clampx(pp, x, cl));
v0 = [P.kd*at0(P.PieceWiseStrainDepR1D, 0), P.k1*at0(P.PieceWiseStrainDep, 0), ...
      P.k_1*at0(P.PieceWiseStrainDepR21, 0), P.k2*at0(P.PieceWiseStrainDep2, sh2)];
cap = {
 sprintf(['refit9_opt transition rates (8 mM MgATP). Strain grid %d bins, %.0f to %.0f nm (grey band edges); ', ...
    'x = pchip knots at their effective positions; rates held at edge-knot values beyond the outer knots (ClampStrainDepExtrapolation).'], ...
    numel(sG), gridEdge(1), gridEdge(2))
 sprintf(['Scale constants (s^-1): ka %.1f, kd %.1f, k1 %.1f, k_1 %.2f, k2 %.1f, kah %.1f, kamh %.2f, ', ...
    'ksr0 %.1f (sigma2 %.1f kPa), kmsr %g, ksr2srd %.1f, ksrd2sr %.2f, kmsrd %.2f (sigma_srd1 %.1f), ksrd %.1f (sigma_srd2 %.1f).'], ...
    P.ka, P.kd, P.k1, P.k_1, P.k2, P.kah, P.kamh, P.ksr0, P.sigma2, P.kmsr, P.ksr2srd, P.ksrd2sr, P.kmsrd, P.sigma_srd1, P.ksrd, P.sigma_srd2)
 sprintf(['Knots: R1D %d, R12 %d (outer knots at %.0f and %.0f nm, off-plot; 3 knots inside the grid), R21 %d, R2 %d (evaluated at s + dr2 - dr = s + %.1f nm). ', ...
    'Value at s = 0 (s^-1): R1D %.0f, R12 %.0f, R21 %.1f, R2 %.1f. Max on the model grid: R12 %.0f, R2 %.0f.'], ...
    numel(kn.R1D.x), numel(kn.R12.x), max(kn.R12.x)*nm, min(kn.R12.x)*nm, numel(kn.R21.x), numel(kn.R2.x), sh2*nm, ...
    v0, max(P.k1*at0(P.PieceWiseStrainDep, sG)), max(P.k2*at0(P.PieceWiseStrainDep2, sG + sh2)))
 sprintf(['Attachment PD->A1 is k_a times gates: registration availability A_reg in [A0 = %.3f, 1] (A0 at isometric steady state), ', ...
    'overlap, lattice spacing, global occupancy (1 - P_bound/%.2f) and load-realignment deficit; new A1 heads land in a %.0f nm triangular kernel at s = 0.'], ...
    P.A0, P.P_bound_max, P.A1AttachmentWidth*nm)
 sprintf(['SRX rates are functions of total force F (F_SR = F_total). A share %.3f of the R2 flux (ramping in over %.0f nm above sSRXrip = %.1f nm) ', ...
    'lands in SRX-ATP instead of PT (panel c, dashed).'], P.SRXFromR2HighStrain, P.dSRXrip*nm, P.sSRXrip*nm)
};
fid = fopen(fullfile(resDir, 'refit9_transitions.txt'), 'w');
fprintf(fid, '%s\n', cap{:});
fclose(fid);
fprintf('%s\n', cap{:});

%% --- helpers --------------------------------------------------------------------
function x = clampx(pp, x, cl)
if cl
    x = min(max(x, pp.breaks(1)), pp.breaks(end));
end
end

function prepAx(ax, fs, xl, yl, gridEdge)
hold(ax, 'on'); box(ax, 'on');
set(ax, 'YScale', 'log', 'FontSize', fs, 'LineWidth', 1, 'Layer', 'top', ...
    'XLim', xl, 'YLim', yl, 'YTick', 10.^(0:4), 'XTick', -75:25:75, ...
    'YMinorTick', 'on', 'TickDir', 'out');
grid(ax, 'on'); ax.GridAlpha = 0.12; ax.MinorGridAlpha = 0.0;
% shade outside the model strain grid
c = [0.90 0.90 0.90];
patch(ax, [xl(1) gridEdge(1) gridEdge(1) xl(1)], [yl(1) yl(1) yl(2) yl(2)], c, 'EdgeColor', 'none', 'HandleVisibility', 'off');
patch(ax, [gridEdge(2) xl(2) xl(2) gridEdge(2)], [yl(1) yl(1) yl(2) yl(2)], c, 'EdgeColor', 'none', 'HandleVisibility', 'off');
xline(ax, 0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

function lab(ax, x, y, str, c, fs)
text(ax, x, y, str, 'Color', c*0.85, 'FontSize', fs-1, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

function prepAxF(ax, fs)
hold(ax, 'on'); box(ax, 'on');
set(ax, 'YScale', 'log', 'FontSize', fs, 'LineWidth', 1, 'Layer', 'top', ...
    'XLim', [0 100], 'YLim', [0.5 2e3], 'YTick', 10.^(0:3), 'XTick', 0:25:100, ...
    'TickDir', 'out');
grid(ax, 'on'); ax.GridAlpha = 0.12;
end
