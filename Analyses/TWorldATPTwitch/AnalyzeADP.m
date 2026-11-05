%% AnalyzeADP.m -- ADP study: how free ADP follows ATP at a given pH, and what it does to the twitch
% Reads results/sa_adp_*.mat (runStandalone('adp')) and results/ad_*.mat (runCoupledADP).
% Writes results/adp_*.csv and figures/figA1_adp_scenarios.png, figA2_traces.png, figA3_split.png.

here = fileparts(mfilename('fullpath'));
res = fullfile(here, 'results');  fig = fullfile(here, 'figures');
ld = @(name) getfield(load(fullfile(res, [name '.mat'])), 'r'); %#ok<GFLD>
pHs = [7.2 6.9 6.6];  As = [8 6 4 2];  paths = {'mild', 'severe'};
fMg = 0.5/(0.68 + 0.5);  K_D = 0.194;  K_T2 = 2.5;

%% 1. the ADP scenarios and what they do to the two ATP/ADP-sensitive steps (analytic)
rows = {};
for q = {'control', 'mild', 'severe'}
    for p = pHs
        for a = As
            [d, f] = adpScenario(a, p, q{1});
            rows(end+1, :) = {q{1}, p, a, f, 1e3*d, fSerca(a, d), gR2(a, fMg*d, 8, fMg*adpScenario(8, 7.2, 'control'), K_T2, K_D), ...
                gR2(a, 0, 8, 0, K_T2, K_D)}; %#ok<SAGROW>
        end
    end
end
T1 = cell2table(rows, 'VariableNames', {'path', 'pH', 'ATP_mM', 'fPCr', 'ADP_uM', 'f_SERCA_vs_ref', 'gR2_ATP_ADP', 'gR2_ATP_only'});
disp(T1);  writetable(T1, fullfile(res, 'adp_scenarios.csv'));

col = containers.Map({7.2, 6.9, 6.6}, {[0.12 0.30 0.55], [0.90 0.55 0.15], [0.75 0.12 0.15]});
lsP = struct('control', ':', 'mild', '--', 'severe', '-');
aa = linspace(2, 8, 61);
fh = figure('Position', [50 50 1200 380], 'Color', 'w');
tl = tiledlayout(fh, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(1, 3);  for j = 1:3, ax(j) = nexttile(tl);  hold(ax(j), 'on'); end
hl = gobjects(0);  lab = {};
for p = pHs
    for q = {'control', 'mild', 'severe'}
        d = arrayfun(@(a) adpScenario(a, p, q{1}), aa);
        h = plot(ax(1), aa, 1e3*d, lsP.(q{1}), 'Color', col(p), 'LineWidth', 1.5);
        plot(ax(2), aa, arrayfun(@(a, x) fSerca(a, x), aa, d), lsP.(q{1}), 'Color', col(p), 'LineWidth', 1.5);
        plot(ax(3), aa, arrayfun(@(a, x) gR2(a, fMg*x, 8, fMg*adpScenario(8, 7.2, 'control'), K_T2, K_D), aa, d), lsP.(q{1}), 'Color', col(p), 'LineWidth', 1.5);
        if strcmp(q{1}, 'severe'), hl(end+1) = h; lab{end+1} = sprintf('pH %.1f', p); end %#ok<SAGROW>
    end
end
plot(ax(1), [4 2], [182 109], 'ko', 'MarkerFaceColor', 'w');  plot(ax(1), [4 2], [284 767], 'ks', 'MarkerFaceColor', 'k');
plot(ax(3), aa, arrayfun(@(a) gR2(a, 0, 8, 0, K_T2, K_D), aa), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 2.5);
set(ax(1), 'YScale', 'log');  ylabel(ax(1), 'free ADP (\muM)');  title(ax(1), 'a  CK-equilibrium ADP', 'FontWeight', 'normal');
ylabel(ax(2), 'f_{SERCA} (TW\_ATP, vs its reference)');  title(ax(2), 'b  SERCA nucleotide factor (Ca side)', 'FontWeight', 'normal');
ylabel(ax(3), 'R2 detachment gate');  title(ax(3), 'c  XB R2 gate, ATP + MgADP competition', 'FontWeight', 'normal');
text(ax(3), 7.8, 0.35, 'thick grey: ATP only, no ADP', 'Color', [0.45 0.45 0.45]);
for j = 1:3, xlabel(ax(j), 'ATP (mM)');  grid(ax(j), 'on');  box(ax(j), 'off');  set(ax(j), 'XDir', 'reverse'); end
legend(ax(1), [hl, plot(ax(1), NaN, NaN, 'k:'), plot(ax(1), NaN, NaN, 'k--'), plot(ax(1), NaN, NaN, 'k-')], ...
    [lab, {'PCr kept (control)', 'mild PCr loss (guinea pig)', 'severe PCr loss (rat)'}], 'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
exportgraphics(fh, fullfile(fig, 'figA1_adp_scenarios.png'), 'Resolution', 150);

%% 2. coupled cases
K = adpCases();
rows = {};
for i = 1:numel(K)
    k = K(i);  r = ld(k.label);
    if ~isempty(r.err), fprintf('%s FAILED: %s\n', k.label, r.err);  continue; end
    S = r.S;  mF = beatMetrics(r.tb, r.Fdev);  mC = beatMetrics(r.tb, S.Ca);
    sa = ld(saLabel(k.A, k.DCa, k.pH));
    rows(end+1, :) = {k.label, k.kind, k.path, k.pH, k.A, 1e3*k.DCa, 1e3*k.DXB, max(r.Fdev), S.Fdia, mF.ttp, mF.d50, mF.d90, ...
        mC.dia, mC.peak, mC.amp, mC.d50, mC.tau, sa.beat(2).APD90, sa.beat(2).CaSR(1), sa.beat(2).f_SERCA_ATP, ...
        trapz(r.tb, S.hyd), max(S.gate), abs(S.drift)/max(r.Fdev), max(abs(S.Ca(:) - r.CaSA(:)))}; %#ok<SAGROW>
end
TC = cell2table(rows, 'VariableNames', {'case', 'kind', 'path', 'pH', 'ATP_mM', 'ADP_Ca_uM', 'ADP_XB_uM', 'Fdev_kPa', 'Fdia_kPa', ...
    'F_TTP_ms', 'F_RT50_ms', 'F_RT90_ms', 'Ca_dia_uM', 'Ca_peak_uM', 'Ca_amp_uM', 'Ca_d50_ms', 'Ca_tau_ms', 'APD90_ms', ...
    'CaSR_dia_mM', 'f_SERCA', 'ATP_per_head_beat', 'gate_max', 'drift_rel', 'maxdCa_vs_standalone_uM'});
writetable(TC, fullfile(res, 'adp_cases.csv'));
fprintf('ADP coupled: %d cases, max drift %.2g %%, max |Ca - standalone| %.2g uM\n', height(TC), 100*max(TC.drift_rel), max(TC.maxdCa_vs_standalone_uM));
disp(TC(:, {'case', 'ADP_Ca_uM', 'Fdev_kPa', 'Fdia_kPa', 'F_TTP_ms', 'F_RT50_ms', 'Ca_dia_uM', 'Ca_amp_uM', 'Ca_d50_ms', 'f_SERCA', 'ATP_per_head_beat'}));

% split relative to the same pH's 8 mM reference
g = @(kind, path, p, a, v) TC.(v)(strcmp(TC.kind, kind) & strcmp(TC.path, path) & abs(TC.pH - p) < 1e-9 & TC.ATP_mM == a);
vars = {'Fdev_kPa', 'F_TTP_ms', 'F_RT50_ms', 'Fdia_kPa', 'Ca_amp_uM', 'Ca_dia_uM', 'Ca_d50_ms', 'ATP_per_head_beat'};
rows = {};
for p = pHs
    for a = As(2:end)
        for q = paths
            for v = vars
                r0 = g("ref", "none", p, 8, v{1});  ra = g("atp", "none", p, a, v{1});
                rc = g("adpCa", q{1}, p, a, v{1});  rf = g("full", q{1}, p, a, v{1});
                rows(end+1, :) = {p, a, q{1}, v{1}, 1e3*adpScenario(a, p, q{1}), 100*(rf/r0 - 1), 100*(ra/r0 - 1), ...
                    100*(rc - ra)/r0, 100*(rf - rc)/r0}; %#ok<SAGROW>
            end
        end
    end
end
TS = cell2table(rows, 'VariableNames', {'pH', 'ATP_mM', 'path', 'metric', 'ADP_uM', 'total_pct', 'ATP_pct', 'ADP_Ca_side_pct', 'ADP_XB_side_pct'});
writetable(TS, fullfile(res, 'adp_split.csv'));
disp(TS(TS.metric == "Fdev_kPa", :));  disp(TS(TS.metric == "F_RT50_ms", :));  disp(TS(TS.metric == "Ca_amp_uM", :));
disp(TS(TS.metric == "Ca_d50_ms", :));  disp(TS(TS.metric == "Fdia_kPa", :));
% the reference row per pH (pH effect at 8 mM with its own CK ADP)
disp(TC(strcmp(TC.kind, "ref"), {'pH', 'ADP_Ca_uM', 'Fdev_kPa', 'Fdia_kPa', 'Ca_dia_uM', 'Ca_amp_uM', 'F_TTP_ms', 'F_RT50_ms', 'ATP_per_head_beat'}));

%% 3. traces: Ca and force at 8/4/2 mM, one column per pH
colA = containers.Map({8, 4, 2}, {[0.12 0.30 0.55], [0.90 0.55 0.15], [0.75 0.12 0.15]});
fh = figure('Position', [50 50 1300 700], 'Color', 'w');
tl = tiledlayout(fh, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
pn = {'unchanged', 'mildly acidic', 'heavily acidic'};
for j = 1:3
    p = pHs(j);  axC = nexttile(tl, j);  hold(axC, 'on');  axF = nexttile(tl, j + 3);  hold(axF, 'on');
    r = ld(sprintf('ad_ref_none_A8_pH%.1f', p));
    plot(axC, 1e3*r.tb, r.S.Ca, '-', 'Color', colA(8), 'LineWidth', 1.6);  plot(axF, 1e3*r.tb, r.Fdev, '-', 'Color', colA(8), 'LineWidth', 1.6);
    for a = [4 2]
        sets = {sprintf('ad_atp_none_A%g_pH%.1f', a, p), ':';  sprintf('ad_full_mild_A%g_pH%.1f', a, p), '--';  sprintf('ad_full_severe_A%g_pH%.1f', a, p), '-'};
        for s = 1:3
            r = ld(sets{s, 1});
            plot(axC, 1e3*r.tb, r.S.Ca, sets{s, 2}, 'Color', colA(a), 'LineWidth', 1.4);
            plot(axF, 1e3*r.tb, r.Fdev, sets{s, 2}, 'Color', colA(a), 'LineWidth', 1.4);
        end
    end
    title(axC, sprintf('pH %.1f (%s): Ca transient', p, pn{j}), 'FontWeight', 'normal');  title(axF, 'XB developed force', 'FontWeight', 'normal');
    ylabel(axC, '[Ca^{2+}]_i (\muM)');  ylabel(axF, 'force (kPa)');  xlabel(axF, 'time after stimulus (ms)');
    for ax = [axC axF], xlim(ax, [0 800]);  grid(ax, 'on');  box(ax, 'off'); end
end
lg = [arrayfun(@(a) plot(nexttile(tl, 1), NaN, NaN, '-', 'Color', colA(a), 'LineWidth', 1.6), [8 4 2]), ...
      plot(nexttile(tl, 1), NaN, NaN, 'k:'), plot(nexttile(tl, 1), NaN, NaN, 'k--'), plot(nexttile(tl, 1), NaN, NaN, 'k-')];
legend(nexttile(tl, 1), lg, {'8 mM', '4 mM', '2 mM', 'ATP only (ADP held)', '+ ADP, mild PCr loss', '+ ADP, severe PCr loss'}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
for a_ = findall(fh, 'Type', 'axes')', a_.Toolbar.Visible = 'off'; end
for j = 4:6, ax = nexttile(tl, j); yl = ylim(ax); ylim(ax, [-0.5 yl(2)]); end
exportgraphics(fh, fullfile(fig, 'figA2_traces.png'), 'Resolution', 150);

%% 4. split bars: ATP effect, Ca-side ADP, XB-side ADP (developed force and RT50)
fh = figure('Position', [50 50 1250 430], 'Color', 'w');
tl = tiledlayout(fh, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for m = {'Fdev_kPa', 'F_RT50_ms'}
    ax = nexttile(tl);  hold(ax, 'on');  Y = [];  xl = {};
    for p = pHs, for a = [4 2], for q = paths
        z = TS(strcmp(TS.metric, m{1}) & abs(TS.pH - p) < 1e-9 & TS.ATP_mM == a & strcmp(TS.path, q{1}), :);
        Y(end+1, :) = [z.ATP_pct z.ADP_Ca_side_pct z.ADP_XB_side_pct]; %#ok<SAGROW>
        xl{end+1} = sprintf('pH %.1f, %g mM, %s', p, a, q{1}); %#ok<SAGROW>
    end, end, end
    bh = bar(ax, Y, 'stacked');  bh(1).FaceColor = [0.6 0.6 0.6];  bh(2).FaceColor = [0.20 0.50 0.75];  bh(3).FaceColor = [0.85 0.45 0.15];
    plot(ax, 1:size(Y, 1), sum(Y, 2), 'kd', 'MarkerFaceColor', 'k');
    set(ax, 'XTick', 1:size(Y, 1), 'XTickLabel', xl, 'XTickLabelRotation', 40, 'FontSize', 8);
    if strcmp(m{1}, 'Fdev_kPa'), ylabel(ax, '\DeltaF_{dev} vs 8 mM, same pH (%)'); else, ylabel(ax, '\DeltaRT50 vs 8 mM, same pH (%)'); end
    grid(ax, 'on');
end
legend(nexttile(tl, 1), {'ATP alone (ADP held)', 'ADP, Ca side (SERCA/NaK)', 'ADP, XB side (R2 competition)', 'total'}, 'Location', 'northeast', 'Box', 'on');
exportgraphics(fh, fullfile(fig, 'figA3_split.png'), 'Resolution', 150);
save(fullfile(res, 'adp_analysis.mat'), 'T1', 'TC', 'TS');

function f = fSerca(A, D)
% TW_ATP / Cortassa 2006 SERCA nucleotide factor, normalised to TW_ATP's reference (7.9875 / 0.0125 mM)
c = @(a, d) 1 ./ ((0.01 ./ a) .* (1 + d/0.14) + 1 + d/5.1);
f = c(A, D) ./ c(7.9875, 0.0125);
end
function g = gR2(A, Dmg, Aref, Dref, KT, KD)
% the XB R2 gate of dPUdT_CombinedTransitions with UseAdpCompR2 (strain-uniform)
g = (A ./ (A + KT*(1 + Dmg/KD))) ./ (Aref ./ (Aref + KT*(1 + Dref/KD)));
end
