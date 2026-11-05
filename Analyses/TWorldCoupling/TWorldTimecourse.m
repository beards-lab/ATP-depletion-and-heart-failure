%% TWorldTimecourse.m -- electrical -> Ca -> force time course of the EHT twitch fits with T-World's Ca.
% For tw1_isoT (adult T-World beat as granted) and tw2_caH (its SERCA/ICaL/NCX freed), reads the
% T-World reference beat stored with each fit (R.S.P.ca.tw.ref: V, Ca, Land Ta and the 93 states),
% recomputes the currents (getCurrentsStructure) and puts membrane voltage, Ca, the fitted XB force
% and the data on one time axis (data time: the model beat starts at R.shift).
% Writes results/tworld_timecourse.png and results/tworld_timecourse.csv (timing ladder).

here = fileparts(mfilename('fullpath'));
res  = fullfile(here, 'results');
fitd = fullfile(here, '..', 'HumanTwitch', 'results');
tiers = {'tw1_isoT', 'tw2_caH'};
lab   = {'T-World as granted (tw1)', 'T-World Ca handling freed (tw2)'};
col   = [0 0.45 0.74; 0.85 0.33 0.1];
R0 = load(fullfile(fitd, 't1_isoT_fit.mat'));  R0 = R0.R;

D = struct();
for i = 1:numel(tiers)
    S = load(fullfile(fitd, [tiers{i} '_fit.mat']));  R = S.R;
    tw = R.S.P.ca.tw;  ref = tw.ref;
    c = getCurrentsStructure({1e3*ref.t}, {ref.X}, tw.param, 0);
    d.t   = 1e3*ref.t + 1e3*R.shift;                 % ms, data time
    d.V   = ref.V;  d.Ca = ref.Ca;  d.Ta = ref.Ta;
    d.ICaL = c.ICaLtot;  d.INaCa = c.INaCa;  d.Jup = c.Jup;  d.Jrel = c.Jrel;
    d.tF  = 1e3*R.tfine;  d.F = R.sigma*R.fm;  d.R = R;
    d.m   = timing(d);
    D.(tiers{i}) = d;
end
tb = 1e3*R0.tb;  Fd = R0.Fd;  dm = R0.mD;

%% timing ladder (ms after the data's stimulus)
names = {'AP upstroke (dV/dt max)', 'APD50 end', 'APD90 end', 'Ca peak', 'Ca 50 % decay', 'Ca 90 % decay', ...
         'force peak', 'force RT50', 'force RT90', 'Land Ta peak'};
T = table(names', 'VariableNames', {'event'});
for i = 1:numel(tiers), m = D.(tiers{i}).m;
    T.(tiers{i}) = [m.up; m.apd50; m.apd90; m.caPk; m.ca50; m.ca90; m.fPk; m.f50; m.f90; m.taPk];
end
T.data = [NaN; NaN; NaN; NaN; NaN; NaN; 1e3*dm.TTP; 1e3*(dm.TTP + dm.RT50); 1e3*(dm.TTP + dm.RT90); NaN];
disp(T);  writetable(T, fullfile(res, 'tworld_timecourse.csv'));
for i = 1:numel(tiers), m = D.(tiers{i}).m;
    fprintf('%s: APD90 %.0f ms, Vrest %.1f mV, Vpeak %.1f mV, Ca %.3f->%.3f uM, AP->Ca %.0f ms, Ca->F %.0f ms\n', ...
        tiers{i}, m.apd90 - m.up, m.Vrest, m.Vpk, m.Cadia, m.CaPkv, m.caPk - m.up, m.fPk - m.caPk);
end

%% figure: A the one-beat cascade (tw2), B action potential and C Ca, baseline T-World vs reparametrized
lab2 = {'T-World baseline (default parameters)', 'T-World reparametrized (SERCA \times0.81, I_{CaL} \times2.76, NCX \times2.01)'};
fh = figure('Color', 'w', 'Position', [40 40 1400 860], 'Visible', 'off');
tl = tiledlayout(fh, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% A: normalised cascade, tw2
ax = nexttile([1 2]);  hold(ax, 'on');  d = D.tw2_caH;
nrm = @(y) (y - min(y))/(max(y) - min(y));
plot(ax, d.t, nrm(d.V), '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.6, 'DisplayName', 'membrane voltage');
plot(ax, d.t, nrm(d.Ca), '-', 'Color', [0.1 0.6 0.3], 'LineWidth', 1.8, 'DisplayName', 'cytosolic Ca');
plot(ax, d.t, nrm(d.Ta), '--', 'Color', [0.55 0.35 0.7], 'LineWidth', 1.2, 'DisplayName', 'T-World''s own Land tension (not used)');
plot(ax, d.tF, d.F/max(Fd), '-', 'Color', col(2, :), 'LineWidth', 2, 'DisplayName', 'XB-model force (fitted)');
plot(ax, tb, Fd/max(Fd), 'k.', 'MarkerSize', 8, 'DisplayName', 'EHT force data');
m = d.m;
xline(ax, m.apd90, ':', sprintf('APD90 %.0f', m.apd90), 'Color', [0.2 0.2 0.2], 'LabelVerticalAlignment', 'middle', 'HandleVisibility', 'off');
xline(ax, m.caPk, ':', sprintf('Ca peak %.0f', m.caPk), 'Color', [0.1 0.6 0.3], 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xline(ax, m.fPk, ':', sprintf('force peak %.0f', m.fPk), 'Color', col(2, :), 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xlim(ax, [0 800]);  ylim(ax, [-0.05 1.1]);  grid(ax, 'on');
legend(ax, 'Location', 'northeast', 'Box', 'off');
xlabel(ax, 'ms after stimulus');  ylabel(ax, 'normalised (min-max)');
title(ax, 'A  one beat, reparametrized T-World: voltage \rightarrow Ca \rightarrow force', 'FontWeight', 'normal');
% B: action potential
ax = nexttile;  hold(ax, 'on');
for i = 1:numel(tiers), d = D.(tiers{i});  plot(ax, d.t, d.V, '-', 'Color', col(i, :), 'LineWidth', 1.6); end
for i = 1:numel(tiers), m = D.(tiers{i}).m;  xline(ax, m.apd90, ':', sprintf('APD90 %.0f', m.apd90), 'Color', col(i, :), 'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off'); end
legend(ax, lab2, 'Box', 'off', 'FontSize', 9, 'Location', 'southwest');
xlim(ax, [0 600]);  grid(ax, 'on');  xlabel(ax, 'ms after stimulus');  ylabel(ax, 'V_m (mV)');
title(ax, 'B  action potential (EP not fitted; APD90 moves with I_{CaL})', 'FontWeight', 'normal');
% C: Ca
ax = nexttile;  hold(ax, 'on');
for i = 1:numel(tiers), d = D.(tiers{i});  plot(ax, d.t, d.Ca, '-', 'Color', col(i, :), 'LineWidth', 1.6); end
plot(ax, 1e3*R0.tfine, R0.S.Ca, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
legend(ax, [lab2, {'de Lange transient (prescribed, t1)'}], 'Box', 'off', 'FontSize', 9);
xlim(ax, [0 800]);  grid(ax, 'on');  xlabel(ax, 'ms after stimulus');  ylabel(ax, '[Ca^{2+}]_i (\muM)');
title(ax, 'C  cytosolic Ca', 'FontWeight', 'normal');
exportgraphics(fh, fullfile(res, 'tworld_timecourse.png'), 'Resolution', 110);
fprintf('written %s\n', fullfile(res, 'tworld_timecourse.png'));

function m = timing(d)
t = d.t;  V = d.V;  Ca = d.Ca;
[~, j] = max(gradient(V, t));  m.up = t(j);
m.Vrest = V(1);  [m.Vpk, jp] = max(V);
amp = m.Vpk - m.Vrest;
k = find(t > t(jp) & V < m.Vpk - 0.5*amp, 1);  m.apd50 = t(k);
k = find(t > t(jp) & V < m.Vpk - 0.9*amp, 1);  m.apd90 = t(k);
m.Cadia = Ca(1);  [m.CaPkv, jc] = max(Ca);  m.caPk = t(jc);
k = find(t > t(jc) & Ca < m.CaPkv - 0.5*(m.CaPkv - m.Cadia), 1);  m.ca50 = t(k);
k = find(t > t(jc) & Ca < m.CaPkv - 0.9*(m.CaPkv - m.Cadia), 1);  m.ca90 = t(k);
[Fp, jf] = max(d.F);  m.fPk = d.tF(jf);
k = find(d.tF > m.fPk & d.F < 0.5*Fp, 1);  m.f50 = d.tF(k);
k = find(d.tF > m.fPk & d.F < 0.1*Fp, 1);  m.f90 = d.tF(k);
[~, jt] = max(d.Ta);  m.taPk = t(jt);
end
