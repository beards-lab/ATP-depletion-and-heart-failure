%% PlotHostKinetics  Tables and figure for RunHostKinetics (#98).
here = fileparts(mfilename('fullpath'));  ana = fileparts(here);  root = fileparts(fileparts(ana));
addpath(fullfile(root, 'Analyses', 'TriSegCaPlan'));
if ~exist('HOSTROOT', 'var'), HOSTROOT = fullfile(fileparts(root), 'Multiscale-Hx-P21-Model'); end
L = load(fullfile(ana, 'results', 'hostKinetics_Nx51.mat'));  res = L.res;
H = hostSignals('Nx51', HOSTROOT);
logf = fullfile(ana, 'results', 'tables.log');  if isfile(logf), delete(logf); end
diary(logf);

% host reference, LV
th = H.raw.t;  hx = H.raw.LV.sigma_XB;  ej = th >= 0.020 & th <= 0.070;  sp = th >= 0.020 & th <= 0.040;
hPeak = max(hx);  hEj = trapz(th(ej), hx(ej))/0.05;  hDip = min(hx(sp));
v = H.LV.v(linspace(0, H.T, 4001));  Ls = H.LV.Ls(linspace(0, H.T, 4001));
fprintf('Host LV (Nx51): active peak %.1f kPa at %.1f ms, mean 20-70 ms %.2f (%.0f %% of peak), lowest 20-40 ms %.2f; peak shortening %.1f um/s = %.1f lengths/s\n', ...
    hPeak, 1e3*th(find(hx == hPeak, 1)), hEj, 100*hEj/hPeak, hDip, -min(v), -min(v)/max(Ls));

fprintf('\nLV arms. Active stress (kPa): peak, lowest in 20-40 ms (the shortening spike), mean 20-70 ms; retention = mean/peak; slack = time with total stress < 0.05 kPa; neg = time with negative cross-bridge stress\n');
fprintf('%-8s %5s %6s %6s | %6s %7s %6s %9s %9s | %6s %6s | %7s\n', 'arm', 's', 'cycle', 'k_on', 'peak', 'dip', 'mean', 'mean/host', 'retention', 'slack', 'neg', 'vCE um/s');
nm = fieldnames(res);  T = struct('name', {}, 's', {}, 'x', {}, 'k', {}, 'peak', {}, 'dip', {}, 'mean', {}, 'slack', {}, 'neg', {});
for i = 1:numel(nm)
    A = res.(nm{i});  if isfield(A, 'err') || ~strcmp(A.wall, 'LV'), continue; end
    k = A.t >= 0.020 & A.t <= 0.040;
    r = struct('name', nm{i}, 's', A.s, 'x', A.xrate, 'k', A.konMult, 'peak', A.m.FaPeak, 'dip', min(A.Fa(k)), 'mean', A.m.FaEj, 'slack', 1e3*A.m.tSlack, 'neg', 1e3*A.m.tNeg);
    T(end+1) = r; %#ok<SAGROW>
    fprintf('%-8s %5.1f %6.2f %6d | %6.1f %7.2f %6.2f %8.0f%% %8.0f%% | %6.2f %6.2f | %7.1f\n', r.name, r.s, r.x, r.k, r.peak, r.dip, r.mean, 100*r.mean/hEj, 100*r.mean/r.peak, r.slack, r.neg, A.m.vCEmax);
end
fprintf('host retention %.0f %%\n', 100*hEj/hPeak);

fprintf('\nWhere in the ejection window the stress is missing: mean active stress (kPa) over 20-40 ms and 40-70 ms, and the thin-filament gate at 30 / 50 / 70 ms\n');
w1 = th >= 0.020 & th <= 0.040;  w2 = th > 0.040 & th <= 0.070;
fprintf('%-8s %8.2f %8.2f\n', 'host', trapz(th(w1), hx(w1))/0.02, trapz(th(w2), hx(w2))/0.03);
for c0 = {'U1', 'U3', 'U5', 'U6p5', 'U10', 'C5', 'C10', 'K5_10'}
    A = res.(c0{1});  k1 = A.t >= 0.020 & A.t <= 0.040;  k2 = A.t > 0.040 & A.t <= 0.070;
    fprintf('%-8s %8.2f %8.2f   gate %.3f %.3f %.3f (max %.3f at %.0f ms)\n', c0{1}, trapz(A.t(k1), A.Fa(k1))/0.02, trapz(A.t(k2), A.Fa(k2))/0.03, ...
        interp1(A.t, A.gate, 0.03), interp1(A.t, A.gate, 0.05), interp1(A.t, A.gate, 0.07), max(A.gate), 1e3*A.t(find(A.gate == max(A.gate), 1)));
end
fprintf('host Ca (uM) at 0 / 10 / 20 / 30 / 50 / 70 ms: %s\n', mat2str(round(H.Ca([0 .01 .02 .03 .05 .07]), 3)));

fprintf('\nOther walls: active mean 20-70 ms over the host''s, slack (ms)\n');
for w = {'SEP', 'RV'}
    for s = [1 5]
        A = res.(sprintf('W%d_%s', s, w{1}));
        fprintf('%-3s uniform x%d: peak %.1f, mean %.2f of host %.2f (%.0f %%), slack %.2f ms\n', w{1}, s, A.m.FaPeak, A.m.FaEj, A.m.hostXBEj, 100*A.m.FaEj/A.m.hostXBEj, 1e3*A.m.tSlack);
    end
end

fprintf('\nQ10 needed from a rig temperature Trig to 37 C for a uniform factor s: Q10 = s^(10/(37-Trig))\n');
for s = [2 3 5 6.5 10], fprintf('s = %4.1f: Trig 20 C -> %.2f, 22 C -> %.2f, 25 C -> %.2f\n', s, s^(10/17), s^(10/15), s^(10/12)); end
fprintf('factor given by a Q10 from Trig to 37 C:\n');
for q = [2 2.5 3.5 4.6], fprintf('Q10 %.1f: 20 C x%.1f, 22 C x%.1f, 25 C x%.1f\n', q, q^1.7, q^1.5, q^1.2); end
diary off;

%% figure
fig = figure('Visible', 'off', 'Position', [50 50 1500 820], 'Color', 'w');
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
cm = turbo(8);
nexttile; hold on;                                   % uniform scaling, traces
U = {'U1', 'U2', 'U3', 'U5', 'U6p5', 'U10'};
for i = 1:numel(U), A = res.(U{i}); plot(1e3*A.t, A.Fa, 'Color', cm(i+1, :), 'LineWidth', 1.2, 'DisplayName', sprintf('all rates \\times %g', A.s)); end
plot(1e3*th, hx, 'k--', 'LineWidth', 1.6, 'DisplayName', 'host''s own cross-bridges');
yline(0, ':', 'HandleVisibility', 'off'); xlim([0 90]); xlabel('time in the beat (ms)'); ylabel('active stress (kPa)');
title('LV, every model rate scaled'); legend('Location', 'northeast', 'Box', 'off');
nexttile; hold on;                                   % cycle only
C = {'U1', 'C2', 'C3', 'C5', 'C10'};
for i = 1:numel(C), A = res.(C{i}); plot(1e3*A.t, A.Fa, 'Color', cm(i+1, :), 'LineWidth', 1.2, 'DisplayName', sprintf('cycle rates \\times %g', A.xrate)); end
A = res.N10; plot(1e3*A.t, A.Fa, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'DisplayName', 'all but the cycle rates \times 10');
plot(1e3*th, hx, 'k--', 'LineWidth', 1.6, 'DisplayName', 'host');
yline(0, ':', 'HandleVisibility', 'off'); xlim([0 90]); xlabel('time in the beat (ms)'); ylabel('active stress (kPa)');
title('LV, cycle rates only'); legend('Location', 'northeast', 'Box', 'off');
nexttile; hold on;                                   % host kinematics
tt = linspace(0, 0.09, 1000);
yyaxis left; plot(1e3*tt, H.LV.Ls(tt), 'LineWidth', 1.4); ylabel('LV sarcomere length (\mum)');
yyaxis right; plot(1e3*tt, -H.LV.v(tt)./max(Ls), 'LineWidth', 1.4); ylabel('shortening speed (lengths/s)');
xline([20 70], ':'); xlabel('time in the beat (ms)'); title('What the host asks for (fit Nx51)');
nexttile; hold on;                                   % mean ejection stress vs factor
sel = @(c) T(arrayfun(@(r) r.name(1) == c, T));
u = sel('U'); c = sel('C'); n = sel('N');
[~, iu] = sort([u.s]); u = u(iu);
plot([u.s], 100*[u.mean]/hEj, 'o-', 'LineWidth', 1.4, 'DisplayName', 'all rates');
plot([1 c.x], 100*[u(1).mean c.mean]/hEj, 's-', 'LineWidth', 1.4, 'DisplayName', 'cycle rates only');
plot([1 n.s], 100*[u(1).mean n.mean]/hEj, 'd-', 'LineWidth', 1.4, 'DisplayName', 'all but the cycle rates');
yline(50, 'k:', 'DisplayName', 'criterion: half of the host''s');
xlabel('rate factor'); ylabel('mean active stress 20-70 ms (% of host)'); legend('Location', 'northwest', 'Box', 'off'); grid on;
title('Stress carried through ejection');
nexttile; hold on;                                   % retention
plot([u.s], 100*[u.mean]./[u.peak], 'o-', 'LineWidth', 1.4, 'DisplayName', 'all rates');
plot([1 c.x], 100*[u(1).mean c.mean]./[u(1).peak c.peak], 's-', 'LineWidth', 1.4, 'DisplayName', 'cycle rates only');
yline(100*hEj/hPeak, 'k--', 'DisplayName', 'host');
xlabel('rate factor'); ylabel('mean 20-70 ms / peak (%)'); ylim([0 80]); legend('Location', 'east', 'Box', 'off'); grid on;
title('Shape: how much of the peak is kept');
nexttile; hold on;                                   % slack and dip
yyaxis left; plot([u.s], [u.slack], 'o-', 'LineWidth', 1.4); ylabel('time slack (ms)');
yyaxis right; plot([u.s], [u.dip], 's-', 'LineWidth', 1.4); yline(hDip, '--'); ylabel('lowest active stress, 20-40 ms (kPa)');
xlabel('rate factor (all rates)'); grid on; title('Slack and the dip at ejection onset (dashed: host)');
title(tl, 'XB model (robust1\_seed, k_{on} \times 3) on the P21 host''s LV beat, open loop: what rate scaling buys');
exportgraphics(fig, fullfile(ana, 'figures', 'rate_scaling.png'), 'Resolution', 150);
close(fig);
