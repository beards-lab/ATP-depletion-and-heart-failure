%% PrepareEHTControls.m -- Plate 6 paced (1 Hz) control EHT twitches -> fitting target
%
% Input:  data/travis/processed/plate6_stim_*.csv  (from extractEHTPlate6.py)
% Output: results/eht_controls.mat   D: per-well beat averages, the control mean
%                                    template and its metrics
%         results/eht_controls.png   wells, control template, well selection
%
% Control wells: the 15 wells that do not carry the isoproterenol signature
% (columns 1-3, plus D4-D6). See EHT_identification_report.md, section 1.
% Beats are cut from each stimulus onset (-50 ms .. +950 ms); the first 3 beats
% of the recording are dropped.

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
res  = fullfile(here, 'results');
src  = fullfile(root, 'data', 'travis', 'processed');

W  = readtable(fullfile(src, 'plate6_stim_waveforms.csv'), 'VariableNamingRule', 'preserve');
S  = readmatrix(fullfile(src, 'plate6_stim_times.csv'));
M  = readtable(fullfile(src, 'plate6_stim_wellmetrics.csv'), 'TextType', 'string');
wells = string(W.Properties.VariableNames(2:end));
t  = W.t_s;  F = W{:, 2:end};
isCtrl = ismember(wells, ["A1" "B1" "C1" "D1" "A2" "B2" "C2" "D2" "A3" "B3" "C3" "D3" "D4" "D5" "D6"]);

dt  = 0.01;  tb = (-0.05:dt:0.95)';                 % grid relative to the stimulus (data are 100 Hz)
ons = S(S + tb(end) <= t(end) & S + tb(1) >= t(1));
ons = ons(4:end);                                   % drop the first 3 beats
nb  = numel(ons);  nw = numel(wells);
beat = zeros(numel(tb), nb, nw);
for w = 1:nw
    for k = 1:nb
        beat(:, k, w) = interp1(t, F(:, w), ons(k) + tb, 'pchip');
    end
end
Fw = squeeze(mean(beat, 2));                        % per-well beat average (uN), numel(tb) x nw
% baseline: Pulse3D exports active force with valleys at ~0; re-zero on the pre-stimulus window
Fw = Fw - mean(Fw(tb < 0, :), 1);

C  = Fw(:, isCtrl);
Cn = C ./ max(C, [], 1);                            % each control well normalised to its peak
D  = struct();
D.tb = tb;  D.wells = wells;  D.isCtrl = isCtrl;  D.Fw = Fw;  D.nBeats = nb;  D.onsets = ons;
D.meanF  = mean(C, 2);   D.sdF  = std(C, 0, 2);     % absolute template (uN)
D.meanFn = mean(Cn, 2);  D.sdFn = std(Cn, 0, 2);    % shape template (-)
D.beatCV = squeeze(std(max(beat(:, :, isCtrl), [], 1), 0, 2) ./ mean(max(beat(:, :, isCtrl), [], 1), 2));
D.metrics = M;
D.m = twitchShape(tb, D.meanF);
D.mWell = arrayfun(@(w) twitchShape(tb, Fw(:, w)), 1:nw);
save(fullfile(res, 'eht_controls.mat'), 'D');

fprintf('%d beats per well, %d control wells\n', nb, nnz(isCtrl));
fprintf('control template: Fpk %.0f uN (well SD %.0f), stim->peak %.0f ms, t10->peak %.0f ms, RT50 %.0f, RT90 %.0f ms, W50 %.0f ms\n', ...
    D.m.Fpk, std(max(C)), 1e3*D.m.TTP, 1e3*D.m.C10toPk, 1e3*D.m.RT50, 1e3*D.m.RT90, 1e3*D.m.W50);
fprintf('Pulse3D (control mean): t10->peak %.0f ms, RT50 %.0f, RT90 %.0f, W50 %.0f ms\n', ...
    1e3*mean(M.C10toPk_s(isCtrl)), 1e3*mean(M.PktoR50_s(isCtrl)), 1e3*mean(M.PktoR90_s(isCtrl)), 1e3*mean(M.W50_s(isCtrl)));

%% figure
Cb = [42 120 214]/255;  Co = [235 104 52]/255;  G = [0.6 0.6 0.58];
fig = figure('Color', 'w', 'Position', [40 40 1250 420], 'Visible', 'off');
tl_ = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_, 'Plate 6, 1 Hz field stimulation: beat-averaged active twitch force', 'FontWeight', 'normal');
ax = nexttile; hold(ax, 'on');
plot(ax, 1e3*tb, Fw(:, ~isCtrl), '-', 'Color', Co);  plot(ax, 1e3*tb, Fw(:, isCtrl), '-', 'Color', Cb);
style(ax, 'A  All 24 wells (blue: control-like, orange: Iso signature)', 'force (\muN)');
ax = nexttile; hold(ax, 'on');
fill(ax, 1e3*[tb; flipud(tb)], [D.meanFn - D.sdFn; flipud(D.meanFn + D.sdFn)], Cb, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(ax, 1e3*tb, Cn, '-', 'Color', [Cb 0.35]);  plot(ax, 1e3*tb, D.meanFn, 'k-', 'LineWidth', 1.8);
style(ax, 'B  Control shape: wells normalised to peak, mean \pm SD', 'F / F_{peak}');
ax = nexttile; hold(ax, 'on');
x = M.spont_freq_Hz;  y = 1e3*M.PktoR50_s;
scatter(ax, x(isCtrl), y(isCtrl), 40, Cb, 'filled');  scatter(ax, x(~isCtrl), y(~isCtrl), 40, Co, 'filled');
text(ax, x, y, "  " + wells', 'FontSize', 7, 'Color', G);
style(ax, 'C  Well selection', 'peak \rightarrow 50% relaxation, paced (ms)');  xlabel(ax, 'spontaneous beating rate (Hz)');
for k = 1:2, xlabel(tl_.Children(end-k+1), 'time from stimulus (ms)'); end
exportgraphics(fig, fullfile(res, 'eht_controls.png'), 'Resolution', 150);

function m = twitchShape(tb, F)
% TTP from stimulus; Pulse3D-like t10->peak, RT50/RT90 from peak, width at 50 %
[Fp, ip] = max(F);  up = F(1:ip);  dn = F(ip:end);  tu = tb(1:ip);  td = tb(ip:end);
cross = @(tt, ff, lvl) interp1(ff(find(ff >= lvl, 1) + [-1 0]), tt(find(ff >= lvl, 1) + [-1 0]), lvl);
crossDn = @(lvl) interp1(dn(find(dn <= lvl, 1) + [-1 0]), td(find(dn <= lvl, 1) + [-1 0]), lvl);
m.Fpk = Fp;  m.TTP = tb(ip);
m.t10 = cross(tu, up, 0.1*Fp);  m.t50 = cross(tu, up, 0.5*Fp);
m.C10toPk = tb(ip) - m.t10;
m.RT50 = crossDn(0.5*Fp) - tb(ip);  m.RT90 = crossDn(0.1*Fp) - tb(ip);
m.W50 = crossDn(0.5*Fp) - m.t50;
end

function style(ax, ttl, yl)
title(ax, ttl, 'FontWeight', 'normal');  ylabel(ax, yl);  box(ax, 'off');  grid(ax, 'on');
set(ax, 'GridAlpha', 0.12, 'XColor', [0.32 0.32 0.31], 'YColor', [0.32 0.32 0.31], 'TickDir', 'out');
end
