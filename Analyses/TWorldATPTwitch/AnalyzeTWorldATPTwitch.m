%% AnalyzeTWorldATPTwitch.m -- tables and figures of the TW_ATP / lowAtp comparison and the coupled twitch
% Reads results/sa_*.mat (runStandalone) and results/cp_*.mat (runCoupled); needs no T-World on the path.
% Writes results/*.csv and figures/fig1_summary.png, fig2_tworld_variants.png, fig3_split.png.

here = fileparts(mfilename('fullpath'));
res = fullfile(here, 'results');  fig = fullfile(here, 'figures');
[ATP, ADP] = atpLadder();  pHs = [7.2 7.0 6.8 6.6];
ld = @(name) getfield(load(fullfile(res, [name '.mat'])), 'r'); %#ok<GFLD>

%% A. standalone: the two ATP-sensitive T-Worlds on the shared ladder (pH 7.2)
sets = {'twatp_A%g_pH7.2',    'TW_ATP, GKATP 0'
        'twatp_A%g_pH7.2_gK', 'TW_ATP, GKATP 0.05 (Michailova)'
        'lowatp_A%g_K1',      'lowAtp, I_KATP on (Mena Tobar)'
        'lowatp_A%g_K0',      'lowAtp, I_KATP off (pumps only)'};
rows = {};
for s = 1:size(sets, 1)
    for i = 1:4
        r = ld(['sa_' sprintf(sets{s, 1}, ATP(i))]);  b = r.beat(2);  b1 = r.beat(1);
        m = beatMetrics(b.t/1e3, b.Ca);  m1 = beatMetrics(b1.t/1e3, b1.Ca);
        fS = getf(b, {'f_SERCA_ATP', 'fATP_SERCA'});  fN = getf(b, {'f_NaK_ATP', 'fATP_NaK'});
        rows(end+1, :) = {sets{s, 2}, ATP(i), 1e3*ADP(i), m.dia, m.peak, m.amp, m.ttp, m.d50, m.d90, m.tau, ...
            b.APD90, b.CaSR(1), b.Na(1), max(b.Ta), max(b.IKATP), fS, fN, abs(m.amp/m1.amp - 1)}; %#ok<SAGROW>
    end
end
TA = cell2table(rows, 'VariableNames', {'model', 'ATP_mM', 'ADP_uM', 'Ca_dia_uM', 'Ca_peak_uM', 'Ca_amp_uM', 'Ca_TTP_ms', ...
    'Ca_d50_ms', 'Ca_d90_ms', 'Ca_tau_ms', 'APD90_ms', 'CaSR_dia_mM', 'Nai_mM', 'Ta_peak_kPa', 'IKATP_peak', ...
    'f_SERCA', 'f_NaK', 'beat_to_beat_amp'});
% relative to the same model's 8 mM row
rel = {'Ca_dia_uM', 'Ca_amp_uM', 'Ca_d50_ms', 'Ca_d90_ms', 'Ca_tau_ms', 'APD90_ms', 'CaSR_dia_mM', 'Nai_mM', 'Ta_peak_kPa'};
for v = rel
    x = TA.(v{1});  y = x;
    for s = 1:size(sets, 1), j = (s-1)*4 + (1:4);  y(j) = 100*(x(j)/x(j(1)) - 1); end
    TA.(['d_' v{1} '_pct']) = y;
end
disp(TA(:, {'model', 'ATP_mM', 'Ca_dia_uM', 'Ca_amp_uM', 'd_Ca_amp_uM_pct', 'd_Ca_d90_ms_pct', 'd_APD90_ms_pct', 'd_CaSR_dia_mM_pct', 'f_SERCA', 'f_NaK', 'IKATP_peak', 'beat_to_beat_amp'}));
writetable(TA, fullfile(res, 'standalone_compare.csv'));

% reference checks
rr = ld('sa_twatp_ref');  rp = ld('sa_twatp_plain');  r8 = ld('sa_twatp_A8_pH7.2');
dRef = maxdiff(rr.beat(2), rp.beat(2));  d8 = maxdiff(r8.beat(2), rr.beat(2));
fprintf('TW_ATP at its own reference vs TW_ATP plain model_TWorld: max|dCa| %.3g uM (peak %.3f), max|dV| %.3g mV\n', dRef.Ca, max(rp.beat(2).Ca), dRef.V);
fprintf('TW_ATP at the ladder 8 mM (ADP 93.7 uM) vs its reference (12.5 uM): max|dCa| %.3g uM, amp %+.2f %%\n', d8.Ca, ...
    100*(beatMetrics(r8.beat(2).t/1e3, r8.beat(2).Ca).amp/beatMetrics(rr.beat(2).t/1e3, rr.beat(2).Ca).amp - 1));
l1 = ld('sa_lowatp_A8_K1');  l0 = ld('sa_lowatp_A8_K0');  dl = maxdiff(l1.beat(2), l0.beat(2));
fprintf('lowAtp at 8 mM, I_KATP on vs off (both at the reference): max|dCa| %.3g uM\n', dl.Ca);
fprintf('lowAtp vs TW_ATP base beat at 8 mM: Ca dia %.4f vs %.4f, peak %.4f vs %.4f uM, APD90 %.1f vs %.1f ms\n', ...
    l1.beat(2).Ca(1), r8.beat(2).Ca(1), max(l1.beat(2).Ca), max(r8.beat(2).Ca), l1.beat(2).APD90, r8.beat(2).APD90);
fs = TA(TA.model == "TW_ATP, GKATP 0" & TA.ATP_mM == 2, :);
fprintf('TW_ATP pH 7.2: amp at 2 mM = %.4f, 8 mM = %.4f\n', fs.Ca_amp_uM, TA.Ca_amp_uM(1));

%% B. standalone TW_ATP over pH (GKATP 0)
rows = {};
for p = pHs
    for i = 1:4
        r = ld(sprintf('sa_twatp_A%g_pH%.1f', ATP(i), p));  b = r.beat(2);  m = beatMetrics(b.t/1e3, b.Ca);
        rows(end+1, :) = {p, ATP(i), m.dia, m.peak, m.amp, m.ttp, m.d50, m.d90, m.tau, b.APD90, b.CaSR(1), b.Na(1), max(b.Ta), ...
            abs(m.amp/beatMetrics(r.beat(1).t/1e3, r.beat(1).Ca).amp - 1)}; %#ok<SAGROW>
    end
end
TB = cell2table(rows, 'VariableNames', {'pH', 'ATP_mM', 'Ca_dia_uM', 'Ca_peak_uM', 'Ca_amp_uM', 'Ca_TTP_ms', 'Ca_d50_ms', ...
    'Ca_d90_ms', 'Ca_tau_ms', 'APD90_ms', 'CaSR_dia_mM', 'Nai_mM', 'LandTa_peak_kPa', 'beat_to_beat_amp'});
disp(TB);  writetable(TB, fullfile(res, 'standalone_twatp_pH.csv'));

%% C. coupled twitch
K = coupledCases();
rows = {};  R = struct();
for i = 1:numel(K)
    k = K(i);  r = ld(k.label);  R.(strrep(strrep(k.label, '.', 'p'), '-', '_')) = r;
    if ~isempty(r.err), fprintf('%s FAILED: %s\n', k.label, r.err);  continue; end
    S = r.S;  mF = beatMetrics(r.tb, r.Fdev);  mC = beatMetrics(r.tb, S.Ca);
    rows(end+1, :) = {k.label, k.core, k.kind, k.aCa, k.aXB, k.pH, max(r.Fdev), S.Fdia, mF.ttp, mF.d50, mF.d90, ...
        mC.dia, mC.peak, mC.amp, mC.ttp, mC.d50, mC.d90, trapz(r.tb, S.hyd), max(S.gate), max(S.bound), ...
        abs(S.drift)/max(r.Fdev), max(abs(S.Ca(:) - r.CaSA(:))), r.wall}; %#ok<SAGROW>
end
TC = cell2table(rows, 'VariableNames', {'case', 'core', 'kind', 'ATP_Ca', 'ATP_XB', 'pH', 'Fdev_kPa', 'Fdia_kPa', 'F_TTP_ms', ...
    'F_RT50_ms', 'F_RT90_ms', 'Ca_dia_uM', 'Ca_peak_uM', 'Ca_amp_uM', 'Ca_TTP_ms', 'Ca_d50_ms', 'Ca_d90_ms', ...
    'ATP_per_head_beat', 'gate_max', 'bound_max', 'drift_rel', 'maxdCa_vs_standalone_uM', 'wall_s'});
disp(TC(:, {'case', 'Fdev_kPa', 'Fdia_kPa', 'F_TTP_ms', 'F_RT50_ms', 'Ca_amp_uM', 'ATP_per_head_beat', 'drift_rel', 'maxdCa_vs_standalone_uM'}));
writetable(TC, fullfile(res, 'coupled_cases.csv'));
fprintf('coupled: %d cases, max beat-to-beat Fdev drift %.2g %%, max |Ca coupled - standalone| %.2g uM\n', ...
    height(TC), 100*max(TC.drift_rel), max(TC.maxdCa_vs_standalone_uM));

% split of each ATP effect at each pH: total = Ca side + XB side + interaction (relative to (8,p|8))
g = @(aCa, aXB, p, v) TC.(v)(TC.core == "t1_isoT" & TC.ATP_Ca == aCa & TC.ATP_XB == aXB & abs(TC.pH - p) < 1e-9);
vars = {'Fdev_kPa', 'F_TTP_ms', 'F_RT50_ms', 'F_RT90_ms', 'Fdia_kPa', 'ATP_per_head_beat'};
rows = {};
for p = pHs
    for a = ATP(2:end)
        for v = vars
            r0 = g(8, 8, p, v{1});  tot = g(a, a, p, v{1})/r0 - 1;  ca = g(a, 8, p, v{1})/r0 - 1;  xb = g(8, a, p, v{1})/r0 - 1;
            rows(end+1, :) = {p, a, v{1}, 100*tot, 100*ca, 100*xb, 100*(tot - ca - xb)}; %#ok<SAGROW>
        end
    end
end
TS = cell2table(rows, 'VariableNames', {'pH', 'ATP_mM', 'metric', 'total_pct', 'Ca_side_pct', 'XB_side_pct', 'interaction_pct'});
disp(TS(TS.metric == "Fdev_kPa" | TS.metric == "F_RT50_ms", :));
writetable(TS, fullfile(res, 'coupled_split.csv'));

% the full grid relative to 8 mM / pH 7.2
F0 = g(8, 8, 7.2, 'Fdev_kPa');  G = zeros(4);  GC = zeros(4);  C0 = g(8, 8, 7.2, 'Ca_amp_uM');
for ip = 1:4, for ia = 1:4
    G(ip, ia) = 100*(g(ATP(ia), ATP(ia), pHs(ip), 'Fdev_kPa')/F0 - 1);
    GC(ip, ia) = 100*(g(ATP(ia), ATP(ia), pHs(ip), 'Ca_amp_uM')/C0 - 1);
end, end
TG = array2table([pHs' G GC], 'VariableNames', [{'pH'}, compose('dFdev_pct_ATP%g', ATP), compose('dCaAmp_pct_ATP%g', ATP)]);
disp(TG);  writetable(TG, fullfile(res, 'coupled_grid.csv'));
rb = TC(TC.core == "s3_free", :);
fprintf('robustness, s3_free core: Fdev %.2f -> %.2f kPa at 2 mM (%+.1f %%) vs t1_isoT %+.1f %%\n', rb.Fdev_kPa(1), rb.Fdev_kPa(2), ...
    100*(rb.Fdev_kPa(2)/rb.Fdev_kPa(1) - 1), 100*(g(2, 2, 7.2, 'Fdev_kPa')/F0 - 1));

%% D. which XB gate carries the XB-side effect (runGateSplit): Ca side at 8 mM, XB at 2 mM
rows = {};
for p = [7.2 6.6]
    c0 = ld(sprintf('cp_t1_isoT_Ca8_XB8_pH%.1f', p));  cb = ld(sprintf('cp_t1_isoT_Ca8_XB2_pH%.1f', p));
    gr = ld(sprintf('gs_R2_pH%.1f', p));  gk = ld(sprintf('gs_Kmsrd_pH%.1f', p));
    lab = {'8 mM (both gates = 1)', '2 mM, both gates', '2 mM, R2 detachment gate only', '2 mM, SRX-return gate only'};
    rr = {c0, cb, gr, gk};
    for j = 1:4
        S = rr{j}.S;
        rows(end+1, :) = {p, lab{j}, max(rr{j}.Fdev), 100*(max(rr{j}.Fdev)/max(c0.Fdev) - 1), S.SRX(1), min(S.SRX), max(S.bound), ...
            beatMetrics(rr{j}.tb, rr{j}.Fdev).ttp, trapz(rr{j}.tb, S.hyd)}; %#ok<SAGROW>
    end
end
TD = cell2table(rows, 'VariableNames', {'pH', 'XB', 'Fdev_kPa', 'dFdev_pct', 'SRX_dia', 'SRX_min', 'bound_max', 'F_TTP_ms', 'ATP_per_head_beat'});
disp(TD);  writetable(TD, fullfile(res, 'gate_split.csv'));

%% Figures
col = containers.Map({8, 6, 4, 2}, {[0.12 0.30 0.55], [0.35 0.60 0.80], [0.90 0.55 0.15], [0.75 0.12 0.15]});
show = [8 4 2];  pShow = [7.2 6.8];  ls = {'-', '--'};

% fig1: Ca transient and XB force, 8/4/2 mM, pH 7.2 (solid) and 6.8 (dashed)
fh = figure('Position', [50 50 720 760], 'Color', 'w');
tl = tiledlayout(fh, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(tl);  hold(ax1, 'on');  ax2 = nexttile(tl);  hold(ax2, 'on');
hl = gobjects(0);  lab = {};
for ip = 1:numel(pShow)
    for a = show
        r = ld(sprintf('cp_t1_isoT_Ca%g_XB%g_pH%.1f', a, a, pShow(ip)));
        h = plot(ax1, 1e3*r.tb, r.S.Ca, ls{ip}, 'Color', col(a), 'LineWidth', 1.6);
        plot(ax2, 1e3*r.tb, r.Fdev, ls{ip}, 'Color', col(a), 'LineWidth', 1.6);
        hl(end+1) = h;  lab{end+1} = sprintf('%g mM ATP, pH %.1f', a, pShow(ip)); %#ok<SAGROW>
    end
end
ylabel(ax1, '[Ca^{2+}]_i (\muM)');  title(ax1, 'a  Ca transient (TW\_ATP, endo, 1 Hz)', 'FontWeight', 'normal');
ylabel(ax2, 'developed force (kPa)');  xlabel(ax2, 'time after stimulus (ms)');
title(ax2, 'b  cross-bridge force (t1\_isoT core, both MgATP gates on)', 'FontWeight', 'normal');
legend(ax1, hl, lab, 'Location', 'northeast', 'Box', 'off', 'NumColumns', 2);
ylim(ax2, [-1 26]);
for ax = [ax1 ax2], xlim(ax, [0 800]);  grid(ax, 'on');  box(ax, 'off');  set(ax, 'FontSize', 11); end
linkaxes([ax1 ax2], 'x');
exportgraphics(fh, fullfile(fig, 'fig1_summary.png'), 'Resolution', 150);
exportgraphics(fh, fullfile(fig, 'fig1_summary.pdf'), 'ContentType', 'vector');

% fig2: the two ATP-sensitive T-Worlds, Ca and AP at 8/4/2 mM, pH 7.2
fh = figure('Position', [50 50 1100 700], 'Color', 'w');
tl = tiledlayout(fh, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
pairs = {'twatp_A%g_pH7.2', 'twatp_A%g_pH7.2_gK', 'TW\_ATP (solid GKATP 0, dashed 0.05)'
         'lowatp_A%g_K0',   'lowatp_A%g_K1',      'TWorld@lowAtp (solid I_{KATP} off, dashed on)'};
for c = 1:2
    axC = nexttile(tl, c);  hold(axC, 'on');  axV = nexttile(tl, c + 2);  hold(axV, 'on');
    for a = show
        for j = 1:2
            b = ld(['sa_' sprintf(pairs{c, j}, a)]).beat(2);
            plot(axC, b.t, b.Ca, ls{j}, 'Color', col(a), 'LineWidth', 1.4);
            plot(axV, b.t, b.V, ls{j}, 'Color', col(a), 'LineWidth', 1.4);
        end
    end
    title(axC, pairs{c, 3}, 'FontWeight', 'normal');  ylabel(axC, '[Ca^{2+}]_i (\muM)');  xlim(axC, [0 800]);
    ylabel(axV, 'V_m (mV)');  xlabel(axV, 'time after stimulus (ms)');  xlim(axV, [0 500]);
    for ax = [axC axV], grid(ax, 'on');  box(ax, 'off'); end
end
lg = arrayfun(@(a) plot(nexttile(tl, 1), NaN, NaN, '-', 'Color', col(a), 'LineWidth', 1.4), show);
legend(lg, compose('%g mM ATP', show), 'Location', 'northeast', 'Box', 'off');
exportgraphics(fh, fullfile(fig, 'fig2_tworld_variants.png'), 'Resolution', 150);

% fig3: grid heat map and the Ca/XB split of the Fdev change
fh = figure('Position', [50 50 1100 420], 'Color', 'w');
tl = tiledlayout(fh, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);  imagesc(ax, G);  colormap(ax, redblue(64));  cl = max(abs(G(:)));  clim(ax, [-cl cl]);  colorbar(ax);
set(ax, 'XTick', 1:4, 'XTickLabel', compose('%g', ATP), 'YTick', 1:4, 'YTickLabel', compose('%.1f', pHs));
xlabel(ax, 'ATP (mM), both sides');  ylabel(ax, 'pH_i (Ca side only)');  title(ax, '\DeltaF_{dev} vs 8 mM / pH 7.2 (%)', 'FontWeight', 'normal');
for ip = 1:4, for ia = 1:4, text(ax, ia, ip, sprintf('%+.0f', G(ip, ia)), 'HorizontalAlignment', 'center', 'FontSize', 10); end, end
ax = nexttile(tl);  hold(ax, 'on');
Y = [];  xl = {};
for p = [7.2 6.6]
    for a = ATP(2:end)
        q = TS(TS.metric == "Fdev_kPa" & TS.ATP_mM == a & abs(TS.pH - p) < 1e-9, :);
        Y(end+1, :) = [q.Ca_side_pct q.XB_side_pct q.interaction_pct]; %#ok<SAGROW>
        xl{end+1} = sprintf('%g mM, pH %.1f', a, p); %#ok<SAGROW>
    end
end
bh = bar(ax, Y, 'stacked');  bh(1).FaceColor = [0.20 0.50 0.75];  bh(2).FaceColor = [0.85 0.45 0.15];  bh(3).FaceColor = [0.6 0.6 0.6];
plot(ax, 1:size(Y, 1), sum(Y, 2), 'kd', 'MarkerFaceColor', 'k');
set(ax, 'XTick', 1:size(Y, 1), 'XTickLabel', xl, 'XTickLabelRotation', 30);  ylabel(ax, '\DeltaF_{dev} vs 8 mM at the same pH (%)');
legend(ax, {'Ca side (TW\_ATP)', 'XB side (MgATP gates)', 'interaction', 'total'}, 'Location', 'best', 'Box', 'off');
title(ax, 'split of the ATP effect on developed force (vs 8 mM, same pH)', 'FontWeight', 'normal');  grid(ax, 'on');
exportgraphics(fh, fullfile(fig, 'fig3_split.png'), 'Resolution', 150);
save(fullfile(res, 'analysis.mat'), 'TA', 'TB', 'TC', 'TS', 'TG', 'TD');

function v = getf(s, names)
v = NaN;
for n = names, if isfield(s, n{1}), v = s.(n{1}); return; end, end
end
function d = maxdiff(a, b)
tq = linspace(0, min(a.t(end), b.t(end)), 4000)';
[ta, ia] = unique(a.t);  [tb, ib] = unique(b.t);
d.Ca = max(abs(interp1(ta, a.Ca(ia), tq) - interp1(tb, b.Ca(ib), tq)));
d.V  = max(abs(interp1(ta, a.V(ia), tq) - interp1(tb, b.V(ib), tq)));
end
function c = redblue(n)
% diverging map: negative -> red, 0 -> white, positive -> blue
x = linspace(-1, 1, n)';
c = [1 - max(x, 0), 1 - abs(x), 1 + min(x, 0)];
end
