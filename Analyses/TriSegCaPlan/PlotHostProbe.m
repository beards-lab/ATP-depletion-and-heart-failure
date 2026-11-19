%% PlotHostProbe  Figures and summary numbers for the #82 open-loop probes.
% Reads results/hostProbe.mat (RunHostProbe.m) and the host fit itself.
here = fileparts(mfilename('fullpath'));  root = fileparts(fileparts(here));  addpath(here);
if ~exist('HOSTROOT', 'var'), HOSTROOT = fullfile(fileparts(root), 'Multiscale-Hx-P21-Model'); end
load(fullfile(here, 'results', 'hostProbe.mat'), 'res');
H = hostSignals('Nx51', HOSTROOT);  hr = H.raw;  tb = hr.t;  T = H.T;
fd = fullfile(here, 'figures');  if ~exist(fd, 'dir'), mkdir(fd); end
ok = @(nm) isfield(res, nm) && ~isfield(res.(nm), 'err');
lastBeat = @(R) R.t >= (round(R.t(end)/T) - 1)*T - 1e-9;
ms = @(R) 1e3*(R.t(lastBeat(R)) - (round(R.t(end)/T) - 1)*T);
set(0, 'DefaultAxesFontSize', 10);

%% summary numbers: waveform metrics of the last beat (active stress)
fprintf('\nWaveform of the active stress, last beat (t from the start of the beat; the Ca transient starts at t = 0):\n');
fprintf('%-22s %8s %8s %8s %8s %8s %10s\n', 'trace', 'peak', 'tPeak', 'RT50', 'RT90', 'min', 'min/peak');
wf('host XB (9 moments)', tb, hr.LV.sigma_XB);
for c = {'E1_ss60', 'E5_kon1', 'E5_kon10', 'E7_x2', 'E7_x5'; 'XB model, k_on x3', 'XB model, k_on x1', 'XB model, k_on x10', 'rates x2', 'rates x5'}
    if ok(c{1}), R = res.(c{1}); k = lastBeat(R); wf(c{2}, ms(R)/1e3, R.Factive(k)); end
end
if ok('E1_ss60')
    R = res.E1_ss60;  k = lastBeat(R);
    fprintf('\nLV, last beat: contractile element shortens at up to %.1f um/s while the host shortens at up to %.1f um/s\n', -min(R.velHS(k)), -min(R.vh(k)));
    tk = R.t(k);  Fk = R.F(k);  dtk = [diff(tk); 0];
    fprintf('   time with total stress < 0.5 kPa: %.1f ms; < 0.05 kPa: %.1f ms; passive at end diastole %.1f kPa (host %.2f kPa)\n', 1e3*sum(dtk(Fk < 0.5)), 1e3*sum(dtk(Fk < 0.05)), R.Fpassive(find(k, 1)), hr.LV.sigma_pas(1));
    fprintf('   A_reg %.3f..%.3f, gate %.3f..%.3f, bound fraction up to %.3f\n', min(R.A_reg(k)), max(R.A_reg(k)), min(R.gate(k)), max(R.gate(k)), max(R.bound(k)));
end

%% Fig 1: host beat and the XB response on it (LV)
if ok('E1_ss60')
    R = res.E1_ss60;  k = lastBeat(R);
    f = figure('Position', [50 50 1150 720], 'Visible', 'off');  tiledlayout(2, 2, 'TileSpacing', 'compact');
    nexttile; hold on;
    plot(1e3*tb, hr.LV.Ls, 'LineWidth', 1.5); plot(1e3*tb, hr.SEP.Ls, 'LineWidth', 1.5); plot(1e3*tb, hr.RV.Ls, 'LineWidth', 1.5);
    ylabel('Ls (\mum)'); xlabel('t (ms)'); legend('LV', 'SEP', 'RV', 'Location', 'best'); title(sprintf('Host sarcomere length (P21 Nx51, %.0f bpm)', H.HR)); grid on;
    nexttile; hold on;
    plot(ms(R), R.F(k), 'k', 'LineWidth', 1.6);  plot(ms(R), R.Factive(k), 'b', 'LineWidth', 1.2);  plot(ms(R), R.Fpassive(k), 'b--', 'LineWidth', 1.2);
    plot(1e3*tb, hr.LV.sigma, 'r', 'LineWidth', 1.4);  plot(1e3*tb, hr.LV.sigma_pas, 'r--', 'LineWidth', 1.2);
    ylabel('stress (kPa)'); xlabel('t (ms)'); title('LV stress on the host''s Ls(t) and Ca(t)'); grid on;
    legend('XB model total', 'XB model active', 'XB model passive', 'host total', 'host passive', 'Location', 'best', 'FontSize', 8);
    nexttile; hold on;
    plot(ms(R), R.vh(k), 'r', 'LineWidth', 1.4); plot(ms(R), R.velHS(k), 'b', 'LineWidth', 1.2);
    yline(-4.6*2.2, 'k:', 'XB model V_{max} 4.6 ML/s'); ylabel('dSL/dt (\mum/s)'); xlabel('t (ms)'); grid on;
    legend('host dLs/dt', 'XB contractile element', 'Location', 'southeast', 'FontSize', 8); title('LV velocity');
    nexttile; hold on;
    yyaxis left;  plot(ms(R), R.gate(k), 'LineWidth', 1.3); plot(ms(R), R.A_reg(k), '--', 'LineWidth', 1.3); ylabel('gate, A_{reg}'); ylim([0 1.05]);
    yyaxis right; plot(1e3*tb, hr.Ca, 'LineWidth', 1.3); ylabel('Ca_i (\muM)');
    xlabel('t (ms)'); legend('Ca gate', 'A_{reg}', 'Ca_i', 'Location', 'best', 'FontSize', 8); title('Activation'); grid on;
    exportgraphics(f, fullfile(fd, 'host_beat_lv.png'), 'Resolution', 140);  close(f);
end

%% Fig 2: grid size, sliding window, approach to the periodic state
f = figure('Position', [50 50 1150 400], 'Visible', 'off');  tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on;
for ss = [60 120 250]
    nm = sprintf('E1_ss%d', ss);  if ~ok(nm), continue; end
    R = res.(nm);  k = lastBeat(R);
    plot(ms(R), R.F(k), 'LineWidth', 1.3, 'DisplayName', sprintf('ss = %d (%.1f s/beat, %d slides in 4 beats)', ss, R.wallPerBeat, R.nShift));
end
if ok('E1_ss60'), R = res.E1_ss60; tsl = R.tShift(R.tShift >= 3*T) - 3*T; for x = 1e3*tsl, xline(x, 'k:', 'HandleVisibility', 'off'); end, end
xlabel('t (ms)'); ylabel('total stress (kPa)'); legend('Location', 'best', 'FontSize', 8); title('LV, last of 4 beats (dotted: window slides at ss = 60)'); grid on;
nexttile; hold on;
if ok('E6'), R = res.E6; plot(R.t/T, R.F, 'k'); plot((1:12) - 0.5, R.peak, 'ro-', 'MarkerFaceColor', 'r'); end
xlabel('beats'); ylabel('total stress (kPa)'); title('Approach to the periodic state (ss = 60)'); grid on;
exportgraphics(f, fullfile(fd, 'grid_size.png'), 'Resolution', 140);  close(f);

%% Fig 3: velocity signal
nm = {'E1_ss60', 'E3_zero', 'E3_filtHS_1ms', 'E3_filtHS_5ms', 'E3_filtHost_1ms', 'E3_filtHost_5ms', 'E3_velHS'};
lab = {'imposed dLs/dt (reference)', 'zero', 'low-pass velHS, 1 ms', 'low-pass velHS, 5 ms', 'low-pass dLs/dt, 1 ms', 'low-pass dLs/dt, 5 ms', 'velHS (algebraic loop)'};
f = figure('Position', [50 50 1150 400], 'Visible', 'off');  tiledlayout(1, 2, 'TileSpacing', 'compact');
ax1 = nexttile; hold on;  ax2 = nexttile; hold on;
for i = 1:numel(nm)
    if ~ok(nm{i}), continue; end
    R = res.(nm{i});  k = lastBeat(R);
    plot(ax1, ms(R), R.F(k), 'LineWidth', 1.2, 'DisplayName', lab{i});  plot(ax2, ms(R), R.A_reg(k), 'LineWidth', 1.2, 'DisplayName', lab{i});
end
xlabel(ax1, 't (ms)'); ylabel(ax1, 'total stress (kPa)'); legend(ax1, 'Location', 'best', 'FontSize', 8); title(ax1, 'LV stress by velocity signal'); grid(ax1, 'on');
xlabel(ax2, 't (ms)'); ylabel(ax2, 'A_{reg}'); title(ax2, 'Registration availability'); grid(ax2, 'on');
exportgraphics(f, fullfile(fd, 'velocity_signal.png'), 'Resolution', 140);  close(f);

%% Fig 4: walls, activation amplitude, rate scaling
f = figure('Position', [50 50 1350 400], 'Visible', 'off');  tiledlayout(1, 3, 'TileSpacing', 'compact');
nexttile; hold on;
nm = {'E1_ss60', 'E4_SEP', 'E4_RV'};  lab = {'LV', 'SEP', 'RV'};  cw = lines(3);
for i = 1:3
    if ~ok(nm{i}), continue; end
    R = res.(nm{i}); k = lastBeat(R); plot(ms(R), R.F(k), 'Color', cw(i, :), 'LineWidth', 1.3, 'DisplayName', ['XB model ' lab{i}]);
    plot(1e3*tb, hr.(lab{i}).sigma, '--', 'Color', cw(i, :), 'LineWidth', 1, 'DisplayName', ['host ' lab{i}]);
end
xlabel('t (ms)'); ylabel('total stress (kPa)'); legend('Location', 'best', 'FontSize', 8); title('All walls (ss = 60, k_{on} \times 3)'); grid on;
nexttile; hold on;
nm = {'E5_kon1', 'E1_ss60', 'E5_kon10'};  lab = {'k_{on} \times 1 (Marzban 2020)', 'k_{on} \times 3', 'k_{on} \times 10'};
for i = 1:3, if ok(nm{i}), R = res.(nm{i}); k = lastBeat(R); plot(ms(R), R.Factive(k), 'LineWidth', 1.3, 'DisplayName', lab{i}); end, end
plot(1e3*tb, hr.LV.sigma_XB, 'k--', 'LineWidth', 1.2, 'DisplayName', 'host XB');
xlabel('t (ms)'); ylabel('active stress (kPa)'); legend('Location', 'best', 'FontSize', 8); title('LV activation amplitude'); grid on;
nexttile; hold on;
nm = {'E1_ss60', 'E7_x2', 'E7_x5'};  lab = {'cycle rates \times 1', '\times 2', '\times 5'};
for i = 1:3, if ok(nm{i}), R = res.(nm{i}); k = lastBeat(R); plot(ms(R), R.Factive(k), 'LineWidth', 1.3, 'DisplayName', lab{i}); end, end
plot(1e3*tb, hr.LV.sigma_XB, 'k--', 'LineWidth', 1.2, 'DisplayName', 'host XB');
xlabel('t (ms)'); ylabel('active stress (kPa)'); legend('Location', 'best', 'FontSize', 8); title('LV: cross-bridge cycle rates scaled'); grid on;
exportgraphics(f, fullfile(fd, 'walls_amplitude_rates.png'), 'Resolution', 140);  close(f);
disp('figures written');

function wf(name, t, y)
% peak, time of peak, time from peak to 50 % and 10 % of (peak - min), minimum
t = t(:);  y = y(:);  [pk, i] = max(y);  mn = min(y);  a = pk - mn;
r50 = find(y(i:end) <= mn + 0.5*a, 1);  r90 = find(y(i:end) <= mn + 0.1*a, 1);
rt = @(r) iff(isempty(r), NaN, 1e3*(t(min(i + max(r, 1) - 1, numel(t))) - t(i)));
fprintf('%-22s %8.1f %8.1f %8.1f %8.1f %8.2f %10.3f\n', name, pk, 1e3*(t(i) - t(1)), rt(r50), rt(r90), mn, mn/pk);
end

function v = iff(c, a, b)
if c, v = a; else, v = b; end
end
