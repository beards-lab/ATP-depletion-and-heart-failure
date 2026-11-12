%% PlotHostProbe  Figures for the #82 open-loop probes (reads results/hostProbe.mat).
here = fileparts(mfilename('fullpath'));  if isempty(here), here = fullfile(pwd, 'Analyses', 'TriSegCaPlan'); end
load(fullfile(here, 'results', 'hostProbe.mat'), 'res');
fd = fullfile(here, 'figures');  if ~exist(fd, 'dir'), mkdir(fd); end
T = res.T;  hr = res.Hraw;  tb = hr.t;
lastBeat = @(R) R.t >= R.t(end) - T + 1e-9;
ms = @(R) 1e3*(R.t(lastBeat(R)) - (R.t(end) - T));

%% Fig 1: host beat and the XB response on it (LV)
f = figure('Position', [50 50 1100 700], 'Visible', 'off');
tiledlayout(2, 2, 'TileSpacing', 'compact');
nexttile; hold on;
plot(1e3*tb, hr.LV.Ls, 'LineWidth', 1.5); plot(1e3*tb, hr.SEP.Ls, 'LineWidth', 1.5); plot(1e3*tb, hr.RV.Ls, 'LineWidth', 1.5);
ylabel('Ls (\mum)'); xlabel('t (ms)'); legend('LV', 'SEP', 'RV', 'Location', 'best'); title('Host sarcomere length (P21 Nx51, 325 bpm)');
nexttile; plot(1e3*tb, hr.Ca, 'k', 'LineWidth', 1.5); ylabel('Ca_i (\muM)'); xlabel('t (ms)'); title('Host Ca transient');
nexttile; hold on;
R = res.E1_ss60;  k = lastBeat(R);
plot(ms(R), R.F(k), 'LineWidth', 1.5);  plot(ms(R), R.Fpassive(k), '--', 'LineWidth', 1.2);
plot(1e3*tb, hr.LV.sigma_pas, ':', 'LineWidth', 1.5);
ylabel('stress (kPa)'); xlabel('t (ms)'); title('LV: XB model on host kinematics');
legend('XB model total', 'XB model passive (Ca-stiffened)', 'host passive', 'Location', 'best');
nexttile; hold on;
plot(ms(R), R.vh(k), 'LineWidth', 1.5); plot(ms(R), R.velHS(k), '--', 'LineWidth', 1.2);
yline(-15, 'k:', 'XB V_{max} \approx 7.5 ML/s'); ylabel('dSL/dt (\mum/s)'); xlabel('t (ms)');
legend('host dLs/dt', 'XB contractile element', 'Location', 'best'); title('LV velocity');
exportgraphics(f, fullfile(fd, 'host_beat_lv.png'), 'Resolution', 150);  close(f);

%% Fig 2: grid size and the sliding window
f = figure('Position', [50 50 1100 400], 'Visible', 'off');
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on;
for ss = [60 120 200]
    R = res.(sprintf('E1_ss%d', ss));  k = lastBeat(R);
    plot(ms(R), R.F(k), 'LineWidth', 1.3, 'DisplayName', sprintf('ss = %d (%.1f s/beat, %d slides)', ss, R.wallPerBeat, R.nShift));
end
R = res.E1_ss60;  tsl = R.tShift(R.tShift >= R.t(end) - T) - (R.t(end) - T);
for x = 1e3*tsl, xline(x, 'k:', 'HandleVisibility', 'off'); end
xlabel('t (ms)'); ylabel('total stress (kPa)'); legend('Location', 'best'); title('LV, last of 4 beats (dotted: ss = 60 slides)');
nexttile; hold on;
for ss = [60 120 200]
    R = res.(sprintf('E1_ss%d', ss));
    plot(R.t/T, R.F, 'DisplayName', sprintf('ss = %d', ss));
end
xlabel('beats'); ylabel('total stress (kPa)'); title('Approach to periodic steady state'); legend('Location', 'best');
exportgraphics(f, fullfile(fd, 'grid_size.png'), 'Resolution', 150);  close(f);

%% Fig 3: velocity signal
f = figure('Position', [50 50 1100 400], 'Visible', 'off');
tiledlayout(1, 2, 'TileSpacing', 'compact');
nm = {'E1_ss60', 'E3_velHS', 'E3_zero'};  lab = {'imposed dLs/dt', 'velHS', 'zero'};
nexttile; hold on;
for i = 1:3, R = res.(nm{i}); k = lastBeat(R); plot(ms(R), R.F(k), 'LineWidth', 1.3, 'DisplayName', lab{i}); end
xlabel('t (ms)'); ylabel('total stress (kPa)'); legend('Location', 'best'); title('LV force by velocity signal');
nexttile; hold on;
for i = 1:3, R = res.(nm{i}); k = lastBeat(R); plot(ms(R), R.A_reg(k), 'LineWidth', 1.3, 'DisplayName', lab{i}); end
xlabel('t (ms)'); ylabel('A_{reg}'); legend('Location', 'best'); title('Registration availability');
exportgraphics(f, fullfile(fd, 'velocity_signal.png'), 'Resolution', 150);  close(f);

%% Fig 4: walls and activation amplitude
f = figure('Position', [50 50 1100 400], 'Visible', 'off');
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on;
nm = {'E1_ss60', 'E4_SEP', 'E4_RV'};  lab = {'LV', 'SEP', 'RV'};
for i = 1:3, R = res.(nm{i}); k = lastBeat(R); plot(ms(R), R.F(k), 'LineWidth', 1.3, 'DisplayName', lab{i}); end
xlabel('t (ms)'); ylabel('total stress (kPa)'); legend('Location', 'best'); title('All walls, ss = 60, k_{on}\times3');
nexttile; hold on;
nm = {'E5_kon1', 'E1_ss60'};  lab = {'k_{on}\times1 (Marzban 2020)', 'k_{on}\times3'};
for i = 1:2, R = res.(nm{i}); k = lastBeat(R); plot(ms(R), R.F(k), 'LineWidth', 1.3, 'DisplayName', lab{i}); end
xlabel('t (ms)'); ylabel('total stress (kPa)'); legend('Location', 'best'); title('LV activation amplitude');
exportgraphics(f, fullfile(fd, 'walls_amplitude.png'), 'Resolution', 150);  close(f);
disp('figures written');
