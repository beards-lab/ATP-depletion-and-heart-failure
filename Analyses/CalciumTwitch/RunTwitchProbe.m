%% RunTwitchProbe.m -- can the skinned-calibrated cross-bridge core make a twitch?
%
% Measured probe behind Analyses/CalciumTwitch/conclusions.md. Does NOT modify
% the model: dPUdT_CaGate wraps dPUdT_CombinedTransitions and gates attachment
% by one thin-filament state (identity to machine precision at gate = 1).
%
%   E0  driving signals     shape metrics of the candidate Ca transients
%   E1  activation clamp    gate(t) prescribed -> the cross-bridge core's own
%                           bandwidth: the fastest twitch ANY Ca model can get
%   E2  steady force-pCa    literature thin filament (Marzban 2020, rat) on this core
%   E3  Ca-driven twitches  Rice 2008 (rat, 24 C, 1 Hz), Marzban 2020 (rat, 6 Hz)
%   E5  k_on scan: is twitch amplitude only a Ca-sensitivity problem?
%   E6  numerical check (tight tolerances, dense output)
%   E4  figures (reads the saved E0, E1, E2, E5 results)
%
% Set RUN (e.g. RUN = 1) before calling to run a subset; default runs everything.
% Params: params/fastgate_best.m (the TriSegCoupling set). Isometric at SL 2.2 um
% with fixed ends, so the series elastic element is in the loop; 8 mM MgATP.

if ~exist('RUN', 'var'), RUN = [0 1 2 3 5 6 4]; end
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..')));
res  = fullfile(here, 'results');

P = getParams(loadParams('fastgate_best'), [], true, false);
P.Vums = 0; P.v = 0;
TF   = struct('k_on', 101.185, 'k_off', 723.852, 'K_coop', 9.6846);  % Marzban 2020 rat
opts = odeset('AbsTol', 1e-4, 'RelTol', 1e-3, 'MaxStep', 1e-3);
dt   = 5e-4;
iSR  = P.NumberOfStates*P.ss + [1 6];          % P_SR, P_SRD (thick-filament OFF pools)

%% E0 -- driving signals ------------------------------------------------------
if any(RUN == 0)
    sigs = {'rice2008', 1; 'fibersim2c', 1; 'tewari_tanh', 1; 'marzban2020', 1/2; 'marzban2020', 1/6};
    tt = (0:1e-4:1)';
    E0 = table();
    for k = 1:size(sigs, 1)
        T  = sigs{k, 2};
        t1 = tt(tt < T);
        Ca = caTransient(sigs{k, 1}, t1, struct('T', T));
        [pk, ip] = max(Ca);  dia = Ca(end);
        j50 = find(Ca(ip:end) <= dia + 0.5*(pk - dia), 1) + ip - 1;
        E0 = [E0; table(string(sigs{k,1}), 1/T, Ca(1), pk, 1e3*t1(ip), 1e3*(t1(j50) - t1(ip)), ...
            'VariableNames', {'signal', 'Hz', 'Ca_onset_uM', 'Ca_peak_uM', 'TTP_ms', 'RT50_ms'})]; %#ok<AGROW>
    end
    disp(E0);
    save(fullfile(res, 'probe_E0.mat'), 'E0', 'sigs');
end

%% E1 -- activation clamp -------------------------------------------------------
if any(RUN == 1)
    Pc = P;  Pc.tf = TF;  Pc.tf.mode = 'clamp';
    tic;
    Pc.tf.gateFcn = @(t) 0;                                    % diastole: nothing attaches
    [t, Y, O] = integrate(Pc, [P.PU0(:); 1], [0 1.5 2], dt, opts);
    rest = Y(end, :)';
    E1.rest = struct('F', O(end,1), 'SRX', sum(Y(end, iSR)), ...
        'drift_last0p5s', max(abs(Y(end,:) - Y(find(t >= 1.5, 1), :))));

    Pc.tf.gateFcn = @(t) 1;                                    % step on, 1 s
    [t, Y, O] = integrate(Pc, rest, [0 1], dt, opts);
    on = Y(end, :)';
    E1.on = riseMetrics(t, O(:,1));
    E1.on.SRX = sum(Y(:, iSR), 2);  E1.on.t = t;  E1.on.F = O(:,1);  E1.on.bound = O(:,5) + O(:,6);

    Pc.tf.gateFcn = @(t) 0;                                    % step off from steady state
    [t, ~, O] = integrate(Pc, on, [0 0.5], dt, opts);
    E1.off = struct('t', t, 'F', O(:,1), 'RT50', crossTime(t, O(:,1), 0.5), ...
        'RT90', crossTime(t, O(:,1), 0.1));

    D = [0.010 0.020 0.050 0.100];                             % twitch-like pulses
    for k = 1:numel(D)
        Pc.tf.gateFcn = @(t) double(t < D(k));
        [t, Y, O] = integrate(Pc, rest, [0 D(k) 0.6], dt, opts);
        m = twitchMetrics(t, O(:,1), 0);
        m.D = D(k);  m.frac = m.Fdev/(E1.on.Fss - E1.on.F0);  m.t = t;  m.F = O(:,1);
        m.SRXmin = min(sum(Y(:, iSR), 2));
        E1.pulse(k) = m;
    end
    E1.wall_s = toc;  E1.sim_s = 2 + 1 + 0.5 + 4*0.6;
    fprintf('E1 rest F %.2f kPa, SRX %.3f | step-on F0 %.2f -> Fss %.2f kPa, t50 %.1f ms, t90 %.1f ms, SRX %.3f -> %.3f\n', ...
        E1.rest.F, E1.rest.SRX, E1.on.F0, E1.on.Fss, 1e3*E1.on.t50, 1e3*E1.on.t90, E1.on.SRX(1), E1.on.SRX(end));
    fprintf('E1 step-off RT50 %.1f ms, RT90 %.1f ms\n', 1e3*E1.off.RT50, 1e3*E1.off.RT90);
    for k = 1:numel(D)
        fprintf('E1 pulse %3.0f ms: Fdev %.2f kPa (%.0f %% of Fss), TTP %.1f ms, RT50 %.1f ms, SRX min %.3f\n', ...
            1e3*D(k), E1.pulse(k).Fdev, 100*E1.pulse(k).frac, 1e3*E1.pulse(k).TTP, 1e3*E1.pulse(k).RT50, E1.pulse(k).SRXmin);
    end
    fprintf('E1 wall %.0f s for %.1f s simulated\n', E1.wall_s, E1.sim_s);
    save(fullfile(res, 'probe_E1.mat'), 'E1', 'rest', 'on');
end

%% E2 -- steady force-pCa with the literature thin filament ---------------------
if any(RUN == 2)
    load(fullfile(res, 'probe_E1.mat'), 'E1', 'rest');
    pCa = [7.0 6.5 6.2 6.0 5.8 5.6 5.4 5.0 4.51];
    Fss = nan(size(pCa));  gss = Fss;  drift = Fss;
    tic;
    for k = 1:numel(pCa)
        Pk = P;  Pk.tf = TF;  Pk.tf.mode = 'ca';
        Pk.tf.CaFcn = @(t) 10^(6 - pCa(k));
        [t, ~, O] = integrate(Pk, rest, [0 1.5], dt, opts);
        Fss(k) = O(end,1);  gss(k) = O(end,13);
        drift(k) = O(end,1) - O(find(t >= 1.4, 1), 1);
    end
    E2 = table(pCa', Fss', gss', drift', 'VariableNames', {'pCa', 'F_kPa', 'gate', 'drift_last100ms'});
    Fn  = (Fss - Fss(1))/(Fss(end) - Fss(1));
    hill = @(x) 1./(1 + 10.^(x(2)*(pCa - x(1))));
    x = fminsearch(@(x) sum((hill(x) - Fn).^2), [5.8 3]);
    E2fit = struct('pCa50', x(1), 'nH', x(2), 'F451_over_clamp', (Fss(end) - E1.on.F0)/(E1.on.Fss - E1.on.F0), ...
                   'wall_s', toc, 'sim_s', 1.5*numel(pCa));
    disp(E2);
    fprintf('E2 Hill: pCa50 %.2f, nH %.2f | developed F(pCa 4.51)/F(gate=1) = %.3f | wall %.0f s\n', ...
        E2fit.pCa50, E2fit.nH, E2fit.F451_over_clamp, E2fit.wall_s);
    save(fullfile(res, 'probe_E2.mat'), 'E2', 'E2fit');
end

%% E3 -- Ca-driven twitches -------------------------------------------------------
if any(RUN == 3)
    load(fullfile(res, 'probe_E1.mat'), 'E1', 'rest');
    cases = {'rice2008', 1, 3; 'marzban2020', 1/6, 6};         % signal, period, beats
    for c = 1:size(cases, 1)
        [sig, T, nb] = cases{c, :};
        Pk = P;  Pk.tf = TF;  Pk.tf.mode = 'ca';
        dia = caTransient(sig, T - 1e-4, struct('T', T));       % diastolic pre-equilibration
        Pk.tf.CaFcn = @(t) dia;
        [~, Y] = integrate(Pk, rest, [0 1], dt, opts);
        Pk.tf.CaFcn = @(t) caTransient(sig, t, struct('T', T));
        knots = unique([(0:nb)*T, (0:nb-1)*T + 0.005*strcmp(sig, 'rice2008')]);
        tic;
        [t, Y, O] = integrate(Pk, Y(end,:)', knots, dt, opts);
        wall = toc;
        m = twitchMetrics(t, O(:,1), (nb-1)*T);
        m.sig = sig;  m.Hz = 1/T;  m.t = t;  m.F = O(:,1);  m.Ca = O(:,15);  m.gate = O(:,13);
        m.SRX = sum(Y(:, iSR), 2);  m.bound = O(:,5) + O(:,6);
        m.frac = m.Fdev/(E1.on.Fss - E1.on.F0);  m.wall_per_beat = wall/nb;
        E3(c) = m; %#ok<SAGROW>
        fprintf('E3 %-12s %4.1f Hz: Fdia %.2f, Fdev %.2f kPa (%.0f %% of Fss), TTP %.1f ms, RT50 %.1f ms, gate peak %.3f, wall %.1f s/beat\n', ...
            sig, 1/T, m.Fd, m.Fdev, 100*m.frac, 1e3*m.TTP, 1e3*m.RT50, max(m.gate), m.wall_per_beat);
    end
    save(fullfile(res, 'probe_E3.mat'), 'E3');
end

%% E5 -- is twitch amplitude just a Ca-sensitivity (k_on/k_off) problem? ---------
if any(RUN == 5)
    load(fullfile(res, 'probe_E1.mat'), 'E1', 'rest');
    mult = [1 3 10];   T = 1;  nb = 3;
    for k = 1:numel(mult)
        Pk = P;  Pk.tf = TF;  Pk.tf.mode = 'ca';  Pk.tf.k_on = TF.k_on*mult(k);
        Pk.tf.CaFcn = @(t) caTransient('rice2008', T - 1e-4, struct('T', T));
        [~, Y] = integrate(Pk, rest, [0 1], dt, opts);
        Pk.tf.CaFcn = @(t) caTransient('rice2008', t, struct('T', T));
        [t, ~, O] = integrate(Pk, Y(end,:)', unique([(0:nb)*T, (0:nb-1)*T + 0.005]), dt, opts);
        b = arrayfun(@(i) twitchMetrics(t, O(:,1), i*T), 0:nb-1);
        m = b(end);  m.k_on_mult = mult(k);  m.frac = m.Fdev/(E1.on.Fss - E1.on.F0);
        m.gate_peak = max(O(:,13));  m.beat_drift_Fdev = b(end).Fdev - b(end-1).Fdev;
        m.t = t;  m.F = O(:,1);  m.gate = O(:,13);
        E5(k) = m; %#ok<SAGROW>
        fprintf('E5 k_on x%-3g: Fdia %.2f, Fdev %.2f kPa (%.0f %% of Fss), TTP %.1f ms, RT50 %.1f ms, gate peak %.3f, beat2->3 dFdev %.2g kPa\n', ...
            mult(k), m.Fd, m.Fdev, 100*m.frac, 1e3*m.TTP, 1e3*m.RT50, m.gate_peak, m.beat_drift_Fdev);
    end
    save(fullfile(res, 'probe_E5.mat'), 'E5');
end

%% E6 -- numerical check: tolerance and output-grid convergence -------------------
if any(RUN == 6)
    load(fullfile(res, 'probe_E1.mat'), 'E1', 'rest');
    tight = odeset('AbsTol', 1e-6, 'RelTol', 1e-5, 'MaxStep', 2e-4);
    Pc = P;  Pc.tf = TF;  Pc.tf.mode = 'clamp';  Pc.tf.gateFcn = @(t) double(t < 0.010);
    [t, ~, O] = integrate(Pc, rest, [0 0.010 0.6], 1e-4, tight);
    m = twitchMetrics(t, O(:,1), 0);
    Pk = P;  Pk.tf = TF;  Pk.tf.mode = 'ca';  Pk.tf.CaFcn = @(t) caTransient('marzban2020', t, struct('T', 1/6));
    [t, ~, O] = integrate(Pk, rest, (0:6)/6, 1e-4, tight);
    m3 = twitchMetrics(t, O(:,1), 5/6);
    load(fullfile(res, 'probe_E3.mat'), 'E3');
    fprintf('E6 10 ms pulse, tight vs probe: Fdev %.3f vs %.3f kPa, TTP %.2f vs %.2f ms, RT50 %.2f vs %.2f ms\n', ...
        m.Fdev, E1.pulse(1).Fdev, 1e3*m.TTP, 1e3*E1.pulse(1).TTP, 1e3*m.RT50, 1e3*E1.pulse(1).RT50);
    fprintf('E6 Marzban 6 Hz beat 6, tight (no diastolic pre-eq) vs probe: Fdev %.3f vs %.3f kPa, TTP %.2f vs %.2f ms, RT50 %.2f vs %.2f ms\n', ...
        m3.Fdev, E3(2).Fdev, 1e3*m3.TTP, 1e3*E3(2).TTP, 1e3*m3.RT50, 1e3*E3(2).RT50);
end

%% E4 -- figures ---------------------------------------------------------------
if any(RUN == 4)
    load(fullfile(res, 'probe_E0.mat'), 'sigs');
    load(fullfile(res, 'probe_E1.mat'), 'E1');
    load(fullfile(res, 'probe_E2.mat'), 'E2', 'E2fit');

    fig = figure('Color', 'w', 'Position', [50 50 1200 800]);
    tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile; hold on;
    tt = (0:1e-4:0.5)';
    for k = 1:size(sigs, 1)
        T = sigs{k, 2};  u = tt(tt < T);
        plot(1e3*u, caTransient(sigs{k,1}, u, struct('T', T)), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s, %.0f Hz', sigs{k,1}, 1/T));
    end
    xlabel('time since stimulus (ms)'); ylabel('[Ca^{2+}]_i (\muM)');
    title('A  Candidate prescribed driving signals'); legend('Interpreter', 'none', 'Box', 'off'); box on;

    nexttile; hold on;
    Fs = E1.on.Fss - E1.on.F0;
    plot(1e3*E1.on.t, (E1.on.F - E1.on.F0)/Fs, 'k', 'LineWidth', 1.5, 'DisplayName', 'gate 0 to 1 (step)');
    plot(1e3*E1.off.t, (E1.off.F - E1.on.F0)/Fs, 'k--', 'LineWidth', 1.5, 'DisplayName', 'gate 1 to 0 from steady');
    for k = 1:numel(E1.pulse)
        plot(1e3*E1.pulse(k).t, (E1.pulse(k).F - E1.on.F0)/Fs, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%.0f ms gate pulse', 1e3*E1.pulse(k).D));
    end
    xlim([0 500]); xlabel('time (ms)'); ylabel('developed force / F_{ss}');
    title('B  Activation clamp: bandwidth of the cross-bridge core'); legend('Box', 'off'); box on;

    nexttile; hold on;
    Fn = (E2.F_kPa - E2.F_kPa(1))/(E2.F_kPa(end) - E2.F_kPa(1));
    pp = linspace(4.4, 7.1, 200);
    plot(E2.pCa, Fn, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'model, steady state');
    plot(pp, 1./(1 + 10.^(E2fit.nH*(pp - E2fit.pCa50))), 'k-', 'DisplayName', ...
        sprintf('Hill: pCa_{50} %.2f, n_H %.2f', E2fit.pCa50, E2fit.nH));
    set(gca, 'XDir', 'reverse'); xlabel('pCa'); ylabel('normalised force');
    title('C  Force-pCa: Marzban 2020 thin filament on this core'); legend('Box', 'off', 'Location', 'northwest'); box on;

    nexttile; hold on;
    load(fullfile(res, 'probe_E5.mat'), 'E5');
    u = (0:1e-4:1)';
    yyaxis right; plot(1e3*u, caTransient('rice2008', u), ':', 'LineWidth', 1.2, 'DisplayName', 'Ca, rice2008 1 Hz');
    ylabel('[Ca^{2+}]_i (\muM)');
    yyaxis left;  cols = [0.00 0.45 0.74; 0.47 0.67 0.19; 0.49 0.18 0.56];
    for k = 1:numel(E5)
        s = E5(k).t >= E5(k).t(end) - 1 - 1e-9;
        plot(1e3*(E5(k).t(s) - E5(k).t(end) + 1), E5(k).F(s), '-', 'Color', cols(k,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('force, k_{on} x%g (%.0f %% of F_{ss})', E5(k).k_on_mult, 100*E5(k).frac));
    end
    ylabel('force (kPa)');  xlim([0 600]); xlabel('time in last beat (ms)');
    title('D  Ca-driven twitch: amplitude is one steep Ca-sensitivity knob'); legend('Box', 'off'); box on;

    exportgraphics(fig, fullfile(res, 'twitch_probe.png'), 'Resolution', 150);
end

%% ---- local functions -----------------------------------------------------------
function [t, Y, O] = integrate(P, Y0, knots, dt, opts)
% Integrate dPUdT_CaGate across KNOTS, restarting at every knot (input
% discontinuities), then evaluate the outputs on the returned grid.
t = [];  Y = [];
for k = 1:numel(knots) - 1
    ts = unique([knots(k):dt:knots(k+1), knots(k+1)]);
    [tk, Yk] = ode15s(@dPUdT_CaGate, ts, Y0, opts, P);
    if k > 1, tk = tk(2:end);  Yk = Yk(2:end, :); end
    t = [t; tk];  Y = [Y; Yk];  Y0 = Yk(end, :)'; %#ok<AGROW>
end
[~, o] = dPUdT_CaGate(t(1), Y(1,:)', P);
O = zeros(numel(t), numel(o));
for i = 1:numel(t)
    [~, O(i,:)] = dPUdT_CaGate(t(i), Y(i,:)', P);
end
end

function m = riseMetrics(t, F)
% Rise from F(1) to the final plateau: times to 10/50/90 % of the change.
F0 = F(1);  Fss = F(end);  x = (F - F0)/(Fss - F0);
m = struct('F0', F0, 'Fss', Fss, 't10', firstCross(t, x, 0.1), ...
           't50', firstCross(t, x, 0.5), 't90', firstCross(t, x, 0.9));
end

function tc = crossTime(t, F, level)
% Time at which a relaxation from F(1) toward the final level crosses LEVEL of its range.
x = (F - F(end))/(F(1) - F(end));
tc = firstCross(t, 1 - x, 1 - level);
end

function tc = firstCross(t, x, level)
j = find(x >= level, 1);
if isempty(j) || j == 1, tc = NaN; return; end
tc = interp1(x([j-1 j]), t([j-1 j]), level);
end

function m = twitchMetrics(t, F, t0)
% Developed force, time to peak and time from peak to 50 % relaxation of the beat at t0.
i0 = find(t >= t0 - 1e-12, 1);
Fd = F(i0);
[Fp, ip] = max(F(i0:end));  ip = ip + i0 - 1;
half = Fd + 0.5*(Fp - Fd);
j = find(F(ip:end) <= half, 1) + ip - 1;
if isempty(j), rt50 = NaN; else, rt50 = interp1(F([j-1 j]), t([j-1 j]), half) - t(ip); end
m = struct('Fd', Fd, 'Fdev', Fp - Fd, 'TTP', t(ip) - t0, 'RT50', rt50);
end
