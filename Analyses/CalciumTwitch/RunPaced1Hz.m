%% RunPaced1Hz.m -- 1 Hz paced isometric twitch, CaModel 'twoState' + analytic Ca
%
% Production Ca layer (calciumActivation.m) on the frozen best parameter set
% params/ExampleHighATP.m (8 mM MgATP, run at SL_TWITCH = 2.0 um, fixed ends: the series
% elastic element is in the loop). Ca: Rice 2008 bi-exponential (rat, 24 C),
% 1 Hz. Thin filament: Marzban 2020 constants (tf_* defaults).
%
%   P1  baseline: 2 s diastolic equilibration, then NB beats; last beat analysed
%   P2  references: saturating Ca (CaModel 'none') and constant Ca = peak Ca
%   P3  illustrative single-parameter perturbations (drug-mechanism directions)
%   P4  figures -> results/paced1Hz_*.png, numbers -> results/paced1Hz.mat
%
% Report: Paced1Hz_report.md. Integrator restarted at every beat onset and at
% the transient's onset kink (+5 ms), as in TestCaActivation.m.

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
res  = fullfile(here, 'results');
opts = odeset('AbsTol', 1e-4, 'RelTol', 1e-3, 'MaxStep', 1e-3);
dt   = 5e-4;
T    = 1;  T_PRE = 2;  NB = 6;
if ~exist('SL_TWITCH', 'var'), SL_TWITCH = 2.0; end        % um; the snapshot's own SL0 is 2.2
sfx  = '';  if SL_TWITCH ~= 2.0, sfx = sprintf('_SL%02.0f', 10*SL_TWITCH); end
base = loadParams(fullfile(root, 'params', 'ExampleHighATP.m'));
base.SL0 = SL_TWITCH;                                        % getParams below rebuilds PU0 at this SL
twitch = {'CaModel', 'twoState', 'CaInput', 'analytic', 'CaTransientKind', 'rice2008', ...
          'CaPeriod', T, 'CaFirstBeat', T_PRE};

%% P1 -- baseline -------------------------------------------------------------------
P  = setParams(base, twitch);
Ns = P.NumberOfStates;  ss = P.ss;  nC = numel(P.PU0) - P.ca.nStates;
ix = struct('SR', Ns*ss+1, 'SL', Ns*ss+3, 'LSE', Ns*ss+4, 'PD', Ns*ss+5, 'SRD', Ns*ss+6, ...
            'Areg', Ns*ss+8, 'Adef', Ns*ss+8+P.UseRegistrationAvailability+P.UseD0State, 'N', P.ca.idx);
tic;
[t, Y, O, R] = integrate(P, P.PU0(:), beatKnots(T_PRE, T, NB), dt, opts);
wall = toc;
B = collect(t, Y, O, R, ix);
beats = arrayfun(@(k) twitchMetrics(t, B.F, T_PRE + k*T, T), 0:NB-1);
last  = beats(end);
fprintf('P1 wall %.0f s for %.0f s simulated (%.1f s/beat)\n', wall, T_PRE + NB*T, wall/(T_PRE/T + NB));
fprintf('P1 per beat Fdev: %s kPa\n', mat2str([beats.Fdev], 5));
fprintf('P1 last beat: Fdia %.3f, Fdev %.3f kPa, TTP %.1f ms, RT50 %.1f ms, RT90 %.1f ms\n', ...
    last.Fd, last.Fdev, 1e3*last.TTP, 1e3*last.RT50, 1e3*last.RT90);

%% P2 -- references ---------------------------------------------------------------------
rest = Y(find(t >= T_PRE - 1e-12, 1), :)';                   % diastolic state at beat-1 onset
Pn = setParams(base, {'CaModel', 'none'});
[tn, ~, On] = integrate(Pn, rest(1:nC), [0 1.5], dt, opts);
Fsat = On(end, 1);
lb = t >= T_PRE + (NB-1)*T - 1e-12;                           % last-beat mask
Ca_pk = max(B.Ca(lb));
Pp = setParams(base, {'CaModel', 'twoState', 'CaInput', 'const', 'pCa', 6 - log10(Ca_pk)});
[~, ~, Op] = integrate(Pp, rest, [0 1.5], dt, opts);
Fpk_ss = Op(end, 1);  gate_pk_ss = Op(end, end);
fprintf('P2 F_rest %.3f | F_sat (none) %.3f | F_ss at peak Ca %.3f uM (pCa %.2f): %.3f kPa, gate %.3f\n', ...
    B.F(find(lb,1)), Fsat, Ca_pk, 6 - log10(Ca_pk), Fpk_ss, gate_pk_ss);

% energetics over the last beat: net hydrolysis flux (PT->PD and SR->SRD), head fraction / s
tl = t(lb);  hyd = B.hyd(lb);
ATP_beat = trapz(tl, hyd);
sys = tl - tl(1) <= last.TTP + last.RT90;
ATP_sys = trapz(tl(sys), hyd(sys));
TTI = trapz(tl, B.F(lb) - last.Fd);                           % developed tension-time integral
fprintf('P2 ATP/head/beat %.4f (systole %.4f, rest-equivalent %.4f) | dev. TTI %.3f kPa s\n', ...
    ATP_beat, ATP_sys, hyd(1)*T, TTI);

%% P3 -- perturbations --------------------------------------------------------------------
V = { 'baseline',            {}
      'tf\_k\_on x3',         {'tf_k_on', 3*P.tf_k_on}
      'tf\_k\_off x2',        {'tf_k_off', 2*P.tf_k_off}
      'k1 x0.5',              {'k1', 0.5*P.k1}
      'ka x0.5',              {'ka', 0.5*P.ka}
      'ksr0 x3',              {'ksr0', 3*P.ksr0}
      'kmsrd x3',             {'kmsrd', 3*P.kmsrd}
      'k2 x0.5',              {'k2', 0.5*P.k2}
      'Ca peak x1.5',         {'CaTransientOpt', struct('Camax', 0.09 + 1.5*(1.45 - 0.09))}
      'Ca decay tau x0.6',    {'CaTransientOpt', struct('tau2', 0.6*0.110)} };
nbV = 4;  preV = 1.5;
PV = struct('name', {}, 'm', {}, 't', {}, 'F', {}, 'gate', {}, 'SRX', {}, 'bound', {});
for k = 1:size(V, 1)
    nv = [twitch, {'CaFirstBeat', preV}, V{k, 2}];
    Pk = setParams(base, nv);
    [tk, Yk, Ok, Rk] = integrate(Pk, Pk.PU0(:), beatKnots(preV, T, nbV), dt, opts);
    Bk = collect(tk, Yk, Ok, Rk, ix);
    bk = arrayfun(@(j) twitchMetrics(tk, Bk.F, preV + j*T, T), 0:nbV-1);
    s  = tk >= preV + (nbV-1)*T - 1e-12;
    m  = bk(end);  m.drift = bk(end).Fdev - bk(end-1).Fdev;
    m.gate_pk = max(Bk.gate(s));  m.SRX_dia = Bk.SRX(find(s,1));  m.SRX_min = min(Bk.SRX(s));
    m.bound_pk = max(Bk.bound(s));  m.ATP = trapz(tk(s), Bk.hyd(s));
    PV(k) = struct('name', V{k,1}, 'm', m, 't', tk(s) - tk(find(s,1)), 'F', Bk.F(s), ...
                   'gate', Bk.gate(s), 'SRX', Bk.SRX(s), 'bound', Bk.bound(s));
    fprintf('P3 %-18s Fdia %.2f Fdev %.2f kPa TTP %.1f RT50 %.1f ms gate pk %.3f SRX %.3f->%.3f ATP %.4f drift %.2g\n', ...
        strrep(V{k,1}, '\', ''), m.Fd, m.Fdev, 1e3*m.TTP, 1e3*m.RT50, m.gate_pk, m.SRX_dia, m.SRX_min, m.ATP, m.drift);
end

save(fullfile(res, ['paced1Hz' sfx '.mat']), 'SL_TWITCH', 't', 'B', 'beats', 'Fsat', 'Ca_pk', 'Fpk_ss', 'gate_pk_ss', ...
     'ATP_beat', 'ATP_sys', 'TTI', 'PV', 'wall');

%% P4 -- figures ------------------------------------------------------------------------------
C = struct('b', [42 120 214]/255, 'o', [235 104 52]/255, 'a', [27 175 122]/255, ...
           'y', [237 161 0]/255, 'g', [0.55 0.55 0.53], 'ink', [0.32 0.32 0.31]);
tb = 1e3*(t(lb) - t(find(lb,1)));                              % ms into last beat
fig = figure('Color', 'w', 'Position', [40 40 1250 1050], 'Visible', 'off');
tl_ = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_, sprintf('1 Hz paced isometric twitch at SL %.1f um, last of %d beats (twoState + Rice 2008, ExampleHighATP)', SL_TWITCH, NB), 'FontWeight', 'normal');

ax = nexttile; plotL(ax, tb, B.Ca(lb), C.b, '[Ca^{2+}]_i');
style(ax, 'A  Ca transient (input)', '[Ca^{2+}]_i (\muM)');

ax = nexttile; plotL(ax, tb, B.F(lb), C.b, 'total (= SE force)');
plotL(ax, tb, B.Fxb(lb), C.o, 'cross-bridge');  plotL(ax, tb, B.Fpas(lb), C.a, 'passive');
style(ax, 'B  Force', 'force (kPa)'); legend(ax, 'Box', 'off', 'Location', 'east');

ax = nexttile; plotL(ax, tb, 1 - B.N(lb), C.b, 'permissive sites 1-N');
plotL(ax, tb, B.gate(lb), C.o, 'attachment gate');
style(ax, 'C  Thin filament', 'fraction'); legend(ax, 'Box', 'off', 'Location', 'northeast');

ax = nexttile; plotL(ax, tb, B.bound(lb), C.b, 'attached p1+p2');
plotL(ax, tb, B.PD(lb), C.o, 'primed PD (M.ADP.Pi)');
plotL(ax, tb, B.PT(lb), C.a, 'PT (M.ATP)');  plotL(ax, tb, B.SRX(lb), C.y, 'SRX (P_{SR}+P_{SRD})');
style(ax, 'D  Myosin head pools', 'fraction of heads'); legend(ax, 'Box', 'off', 'Location', 'east');

ax = nexttile; plotL(ax, tb, B.RD1(lb), C.b, 'attachment PD\rightarrowp1');
plotL(ax, tb, B.R2(lb), C.o, 'detachment p2\rightarrowPT');
plotL(ax, tb, B.srxOut(lb), C.a, 'SRX exit SRD\rightarrowPD');  plotL(ax, tb, B.srxIn(lb), C.y, 'SRX entry PT\rightarrowSR');
style(ax, 'E  Fluxes', 'flux (heads s^{-1})'); legend(ax, 'Box', 'off', 'Location', 'east');

ax = nexttile; plotL(ax, tb, B.Areg(lb), C.b, 'registration A_{reg}');
plotL(ax, tb, B.Adef(lb), C.o, 'realignment deficit A_{def}');
style(ax, 'F  Slow attachment gates', 'fraction'); legend(ax, 'Box', 'off', 'Location', 'east');
for a = findall(fig, 'Type', 'axes')', xlabel(a, 'time since stimulus (ms)'); end
exportgraphics(fig, fullfile(res, ['paced1Hz' sfx '_lastbeat.png']), 'Resolution', 150);

fig = figure('Color', 'w', 'Position', [40 40 1250 800], 'Visible', 'off');
tl_ = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_, 'Whole run: 2 s diastolic equilibration, then 1 Hz pacing', 'FontWeight', 'normal');
ax = nexttile; plotL(ax, t, B.F, C.b, 'force');  style(ax, 'A  Force', 'force (kPa)');
ax = nexttile; plotL(ax, t, B.gate, C.o, 'gate');  style(ax, 'B  Attachment gate', 'fraction');
ax = nexttile; plotL(ax, t, B.SRX, C.y, 'SRX');  plotL(ax, t, B.bound, C.b, 'attached');
style(ax, 'C  SRX and attached heads', 'fraction of heads'); legend(ax, 'Box', 'off', 'Location', 'east');
xlabel(ax, 'time (s)');
exportgraphics(fig, fullfile(res, ['paced1Hz' sfx '_overview.png']), 'Resolution', 150);

fig = figure('Color', 'w', 'Position', [40 40 620 520], 'Visible', 'off');
ax = axes(fig); plotL(ax, B.Ca(lb), B.F(lb), C.b, 'last beat');
[~, ip] = max(B.F(lb));  ca_l = B.Ca(lb);  f_l = B.F(lb);
plot(ax, ca_l(1), f_l(1), 'o', 'Color', C.ink, 'MarkerFaceColor', C.ink, 'MarkerSize', 7);
text(ax, ca_l(1), f_l(1), '  stimulus', 'Color', C.ink, 'VerticalAlignment', 'bottom');
plot(ax, ca_l(ip), f_l(ip), 'o', 'Color', C.ink, 'MarkerFaceColor', C.ink, 'MarkerSize', 7);
text(ax, ca_l(ip), f_l(ip), '  force peak', 'Color', C.ink, 'VerticalAlignment', 'bottom');
style(ax, 'Force-Ca loop (counter-clockwise = force lags Ca)', 'force (kPa)'); xlabel(ax, '[Ca^{2+}]_i (\muM)');
exportgraphics(fig, fullfile(res, ['paced1Hz' sfx '_loop.png']), 'Resolution', 150);

fig = figure('Color', 'w', 'Position', [40 40 1300 620], 'Visible', 'off');
tl_ = tiledlayout(fig, 2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_, 'Illustrative single-parameter perturbations (last of 4 beats); grey = baseline', 'FontWeight', 'normal');
yl = [min(arrayfun(@(p) min(p.F), PV)) max(arrayfun(@(p) max(p.F), PV))];
for k = 1:numel(PV)
    ax = nexttile;
    if k > 1, plotL(ax, 1e3*PV(1).t, PV(1).F, C.g, 'baseline'); end
    plotL(ax, 1e3*PV(k).t, PV(k).F, C.b, PV(k).name);
    style(ax, PV(k).name, '');  ylim(ax, yl + [-1 1]);  xlim(ax, [0 1000]);
    if mod(k-1, 5) == 0, ylabel(ax, 'force (kPa)'); end
    if k > 5, xlabel(ax, 'ms'); end
end
exportgraphics(fig, fullfile(res, ['paced1Hz' sfx '_perturb.png']), 'Resolution', 150);
fprintf('figures written\n');

%% ---- local functions ------------------------------------------------------------------------
function p = setParams(base, nv)
for i = 1:2:numel(nv), base.(nv{i}) = nv{i+1}; end
p = getParams(base, [], true, false);
p.Vums = 0;  p.v = 0;
end

function k = beatKnots(t1, T, nb)
k = unique([0, t1 + (0:nb)*T, t1 + (0:nb-1)*T + 0.005]);
end

function [t, Y, O, R] = integrate(P, Y0, knots, dt, opts)
t = [];  Y = [];
for k = 1:numel(knots) - 1
    ts = unique([knots(k):dt:knots(k+1), knots(k+1)]);
    [tk, Yk] = ode15s(@dPUdT_CombinedTransitions, ts, Y0, opts, P);
    if k > 1, tk = tk(2:end);  Yk = Yk(2:end, :); end
    t = [t; tk];  Y = [Y; Yk];  Y0 = Yk(end, :)'; %#ok<AGROW>
end
[~, o, r] = dPUdT_CombinedTransitions(t(1), Y(1,:)', P);
O = zeros(numel(t), numel(o));  R = zeros(numel(t), numel(r));
for i = 1:numel(t)
    [~, O(i,:), R(i,:)] = dPUdT_CombinedTransitions(t(i), Y(i,:)', P);
end
end

function B = collect(t, Y, O, R, ix)
% outputs (2-state): 1 F, 2 F_active, 3 F_passive, 5 p1_0, 6 p2_0, 9 PT, ..., end-1 Ca, end gate
% rates: 1 RTD 2 RDT 3 RD1 7 R2 10 RPT2SR 11 RSRD2PD 13 RSR2SRD 14 RSRD2SR
B.t = t;  B.F = O(:,1);  B.Fxb = O(:,2);  B.Fpas = O(:,3);
B.bound = O(:,5) + O(:,6);  B.PT = O(:,9);  B.Ca = O(:,end-1);  B.gate = O(:,end);
B.PD = Y(:, ix.PD);  B.SRX = Y(:, ix.SR) + Y(:, ix.SRD);  B.N = Y(:, ix.N);
B.Areg = Y(:, ix.Areg);  B.Adef = Y(:, ix.Adef);  B.LSE = Y(:, ix.LSE);
B.RD1 = R(:,3);  B.R2 = R(:,7);  B.srxIn = R(:,10);  B.srxOut = R(:,11);
B.hyd = (R(:,1) - R(:,2)) + (R(:,13) - R(:,14));
end

function m = twitchMetrics(t, F, t0, T)
s  = t >= t0 - 1e-12 & t <= t0 + T + 1e-12;
ts = t(s);  Fs = F(s);
Fd = Fs(1);  [Fp, ip] = max(Fs);
m = struct('Fd', Fd, 'Fp', Fp, 'Fdev', Fp - Fd, 'TTP', ts(ip) - t0, ...
    'RT50', relax(ts, Fs, ip, Fd, 0.5), 'RT90', relax(ts, Fs, ip, Fd, 0.9));
end

function r = relax(ts, Fs, ip, Fd, frac)
lvl = Fs(ip) - frac*(Fs(ip) - Fd);
j = find(Fs(ip:end) <= lvl, 1) + ip - 1;
if isempty(j), r = NaN; else, r = interp1(Fs([j-1 j]), ts([j-1 j]), lvl) - ts(ip); end
end

function plotL(ax, x, y, c, name)
hold(ax, 'on');
plot(ax, x, y, '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', name);
end

function style(ax, ttl, yl)
title(ax, ttl, 'FontWeight', 'normal', 'Interpreter', 'tex');
ylabel(ax, yl);  box(ax, 'off');  grid(ax, 'on');
set(ax, 'GridAlpha', 0.12, 'XColor', [0.32 0.32 0.31], 'YColor', [0.32 0.32 0.31], 'TickDir', 'out');
end
