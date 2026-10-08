%% NonUniqueEHTTwitch.m -- the EHT twitch fit is not unique: exact symmetries and equivalent refits
%
% Starting point: fit C (results/eht_fitC.mat; FitCEHTTwitch.m).
%   N1  exact symmetries at fit C: (a) Ca levels x lambda with tf_k_on / lambda; (b) tau1 <-> tau2.
%       Force, gate, bound heads and ATP are compared sample by sample.
%   N2  refits of all 6 free parameters (tau1, tau2, t0 = onset offset, tf_k_on, tf_K_coop,
%       tf_k_off; sigma in closed form) under different FIXED assumptions:
%         Cadia 0.05 / 0.20 uM; cross-bridge cycling rates x0.5 / x2 (kah kamh ka kd k1 k_1 k2); K_coop fixed at 10;
%         tf_k_off fixed at 0.5x / 2x fit C (the others refitted)
%       Each: fit quality, parameters, steady pCa50/nH, internal states, ATP, ATP x sigma.
%   N3  figure -> results/eht_nonunique.png; numbers -> results/eht_nonunique.mat
% Fits are cached in results/eht_nonunique_fits.mat (delete to redo, ~40 min).
% Report: EHT_identification_report.md, section 5.

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
res  = fullfile(here, 'results');
L  = load(fullfile(res, 'eht_controls.mat'));  D = L.D;
FC = load(fullfile(res, 'eht_fitC.mat'));
base0 = loadParams(fullfile(root, 'params', 'ExampleHighATP.m'));  base0.SL0 = 2.0;
sel = D.tb >= 0;  tb = D.tb(sel);  Fd = D.meanF(sel);
sem = sqrt(mean(D.sdFn(sel).^2))/sqrt(nnz(D.isCtrl));
tfine = (0:0.0025:1)';
names6 = {'tau1', 'tau2', 't0', 'tf_k_on', 'tf_K_coop', 'tf_k_off'};
lb6 = [0.010 0.080 -0.150 log10(0.1) 1 log10(0.5)];
ub6 = [0.150 0.300  0.060 log10(3000) 20 log10(2000)];
xC = FC.xC;  fixC = struct('Cadia', 0.10, 'Camax', 1.0);
if isempty(gcp('nocreate')), refreshPool(5); end

%% N1 -- exact symmetries ----------------------------------------------------------------------
thC = toTh6(xC, fixC, names6);
[fC, SC] = simEHTTwitch(base0, thC, tfine);
sym = struct('name', {}, 'dF', {}, 'dGate', {}, 'dBound', {}, 'dATP', {});
for lam = [0.5 2]
    t = thC;  t.Cadia = lam*thC.Cadia;  t.Camax = lam*thC.Camax;  t.tf_k_on = thC.tf_k_on/lam;
    [f, S] = simEHTTwitch(base0, t, tfine);
    sym(end+1) = struct('name', sprintf('Ca levels x%g, k_on /%g', lam, lam), 'dF', max(abs(f - fC))/max(fC), ...
        'dGate', max(abs(S.gate - SC.gate)), 'dBound', max(abs(S.bound - SC.bound)), 'dATP', abs(trapz(tfine, S.hyd)/trapz(tfine, SC.hyd) - 1)); %#ok<SAGROW>
end
t = thC;  t.tau1 = thC.tau2;  t.tau2 = thC.tau1;
[f, S] = simEHTTwitch(base0, t, tfine);
sym(end+1) = struct('name', 'tau1 <-> tau2', 'dF', max(abs(f - fC))/max(fC), 'dGate', max(abs(S.gate - SC.gate)), ...
    'dBound', max(abs(S.bound - SC.bound)), 'dATP', abs(trapz(tfine, S.hyd)/trapz(tfine, SC.hyd) - 1));
for k = 1:numel(sym)
    fprintf('N1 %-22s max rel |dF| %.1e  |d gate| %.1e  |d bound| %.1e  |d ATP| %.1e\n', sym(k).name, sym(k).dF, sym(k).dGate, sym(k).dBound, sym(k).dATP);
end

%% N2 -- equivalent refits under different fixed assumptions ---------------------------------------
% each alternative: label, fixed-parameter overrides (Ca / thin filament), base overrides (core),
% and a mask of which of the 6 parameters is held at a given value
A = struct('name', {}, 'fix', {}, 'core', {}, 'hold', {});
A(end+1) = struct('name', 'Cadia 0.05 uM', 'fix', struct('Cadia', 0.05, 'Camax', 1.0), 'core', struct(), 'hold', struct());
A(end+1) = struct('name', 'Cadia 0.20 uM', 'fix', struct('Cadia', 0.20, 'Camax', 1.0), 'core', struct(), 'hold', struct());
A(end+1) = struct('name', 'core rates x0.5', 'fix', fixC, 'core', struct('rateScale', 0.5), 'hold', struct());
A(end+1) = struct('name', 'core rates x2', 'fix', fixC, 'core', struct('rateScale', 2), 'hold', struct());
A(end+1) = struct('name', 'K_coop fixed 10', 'fix', fixC, 'core', struct(), 'hold', struct('i', 5, 'v', 10));
A(end+1) = struct('name', 'k_off fixed x0.5', 'fix', fixC, 'core', struct(), 'hold', struct('i', 6, 'v', xC(6) + log10(0.5)));
A(end+1) = struct('name', 'k_off fixed x2', 'fix', fixC, 'core', struct(), 'hold', struct('i', 6, 'v', xC(6) + log10(2)));
nA = numel(A);
cacheN = fullfile(res, 'eht_nonunique_fits.mat');
if exist(cacheN, 'file'), Cc = load(cacheN);  XA = Cc.XA;  SA = Cc.SA;     % rows with NaN SSR are (re)run
else, XA = nan(nA, 6);  SA = nan(nA, 1); end
todo = find(isnan(SA))';
if ~isempty(todo)
    lopt = optimoptions('lsqnonlin', 'Display', 'off', 'UseParallel', false, 'FiniteDifferenceStepSize', 1e-3, ...
        'StepTolerance', 1e-4, 'FunctionTolerance', 1e-6, 'MaxIterations', 60);
    xo = cell(numel(todo), 1);  so = zeros(numel(todo), 1);  At = A(todo);
    parfor q = 1:numel(todo)
        b = applyCore(base0, At(q).core);  m = true(1, 6);  x0 = xC; %#ok<PFBNS>
        if isfield(At(q).hold, 'i'), m(At(q).hold.i) = false;  x0(At(q).hold.i) = At(q).hold.v; end
        fr = @(z) resid6(expand(z, x0, m), b, At(q).fix, names6, tb, Fd);
        [z, rn] = lsqnonlin(fr, x0(m), lb6(m), ub6(m), lopt);
        xo{q} = expand(z, x0, m);  so(q) = rn;
    end
    for q = 1:numel(todo), XA(todo(q), :) = xo{q};  SA(todo(q)) = so(q); end
    save(cacheN, 'XA', 'SA', 'A');
end

%% diagnostics for fit C and every alternative ----------------------------------------------------
pc = [7 6.6 6.3 6.1 5.9 5.7 5.5 5.3 5 4.51];
alts = [struct('name', 'fit C', 'fix', fixC, 'core', struct(), 'hold', struct()), A];
Xall = [xC; XA];  Sall = [sum(FC.r0.^2); SA];
Rw = cell(numel(alts), 1);
parfor a = 1:numel(alts)
    b = applyCore(base0, alts(a).core);  th = toTh6(Xall(a, :), alts(a).fix, names6); %#ok<PFBNS>
    [fm, S] = simEHTTwitch(b, th, tfine);
    fmd = interp1(tfine, fm, tb);  sig = (fmd'*Fd)/(fmd'*fmd);
    Fss = arrayfun(@(p) steadyF(b, th, p), pc);
    [p50, nH] = hill(pc, Fss);
    [~, iCa] = max(S.Ca);  [~, iF] = max(fm);
    Rw{a} = struct('th', th, 'sigma', sig, 'rms', sqrt(mean((sig*fmd - Fd).^2))/max(Fd), 'fm', fm, 'S', S, ...
        'pCa50', p50, 'nH', nH, 'act', max(fm)/(Fss(end) - S.Fdia), 'gatePk', max(S.gate), 'boundPk', max(S.bound), ...
        'SRXmin', min(S.SRX), 'ATP', trapz(tfine, S.hyd), 'ATPsig', trapz(tfine, S.hyd)*sig, 'ATPdia', S.hyd(1), ...
        'lagCaF', tfine(iF) - tfine(iCa), 'koffEff', [min(th.tf_k_off*(1 + th.tf_K_coop*S.N)), max(th.tf_k_off*(1 + th.tf_K_coop*S.N))]);
end
R = [Rw{:}];
fprintf('N2 data SEM %.2f %% of peak\n', 100*sem);
for a = 1:numel(alts)
    t = R(a).th;
    fprintf('N2 %-17s rms %.2f %% | tau %3.0f/%3.0f ms  onset %+5.1f ms  k_on %6.2f  K_coop %5.2f  k_off %6.2f | pCa50 %.2f nH %.1f | sigma %.2f act %2.0f %% gate %.2f bound %.3f SRX %.3f | ATP %.2f (dia %.3f) ATPxsig %.1f | Ca->F lag %3.0f ms | k_off,eff %.0f-%.0f\n', ...
        alts(a).name, 100*R(a).rms, 1e3*t.tau1, 1e3*t.tau2, 1e3*t.shift, t.tf_k_on, t.tf_K_coop, t.tf_k_off, R(a).pCa50, R(a).nH, ...
        R(a).sigma, 100*R(a).act, R(a).gatePk, R(a).boundPk, R(a).SRXmin, R(a).ATP, R(a).ATPdia, R(a).ATPsig, 1e3*R(a).lagCaF, R(a).koffEff);
end
save(fullfile(res, 'eht_nonunique.mat'), 'sym', 'alts', 'Xall', 'Sall', 'R', 'tfine', 'tb', 'Fd', 'sem', 'pc');

%% N3 -- figure ----------------------------------------------------------------------------------------
nAl = numel(alts);  cm = lines(nAl);  cm(1, :) = [0 0 0];
fig = figure('Color', 'w', 'Position', [30 30 1650 900], 'Visible', 'off');
tl_ = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_, 'The fit is not unique: refits under different fixed assumptions (all 6 Ca/thin-filament parameters refitted each time)', 'FontWeight', 'normal');
ax = nexttile; hold(ax, 'on');
plot(ax, 1e3*tb, Fd, 'k.', 'MarkerSize', 10, 'DisplayName', 'data');
for a = 1:nAl
    plot(ax, 1e3*tfine, R(a).sigma*R(a).fm, '-', 'Color', cm(a, :), 'LineWidth', 1.2 + 0.8*(a == 1), ...
        'DisplayName', sprintf('%s (rms %.2f %%)', alts(a).name, 100*R(a).rms));
end
style(ax, 'A  Fitted twitches', 'active force (\muN)');  legend(ax, 'Box', 'off', 'FontSize', 7);  xlabel(ax, 'time from stimulus (ms)');
ax = nexttile; hold(ax, 'on');
for a = 1:nAl, plot(ax, 1e3*tfine, R(a).S.Ca, '-', 'Color', cm(a, :), 'LineWidth', 1.2 + 0.8*(a == 1)); end
style(ax, 'B  Ca transients behind them', '[Ca^{2+}]_i (\muM)');  xlabel(ax, 'time from stimulus (ms)');
ax = nexttile; hold(ax, 'on');
for a = 1:nAl, plot(ax, 1e3*tfine, R(a).S.gate, '-', 'Color', cm(a, :), 'LineWidth', 1.2 + 0.8*(a == 1)); end
style(ax, 'C  Attachment gate', 'gate');  xlabel(ax, 'time from stimulus (ms)');
ax = nexttile; hold(ax, 'on');                                                   % parameter ratios to fit C
pn = {'\tau (mean)', 'k_{on}', 'K_{coop}', 'k_{off}', 'k_{on}/k_{off}', '\sigma'};
val = @(r) [mean([r.th.tau1 r.th.tau2]), r.th.tf_k_on, r.th.tf_K_coop, r.th.tf_k_off, r.th.tf_k_on/r.th.tf_k_off, r.sigma];
v0 = val(R(1));
for a = 2:nAl, plot(ax, 1:numel(pn), val(R(a))./v0, 'o-', 'Color', cm(a, :), 'MarkerFaceColor', cm(a, :)); end
yline(ax, 1, 'k-');  set(ax, 'YScale', 'log', 'XTick', 1:numel(pn), 'XTickLabel', pn);  xlim(ax, [0.5 numel(pn) + 0.5]);
style(ax, 'D  Fitted values relative to fit C (log)', 'ratio to fit C');
ax = nexttile; hold(ax, 'on');                                                   % derived quantities
qn = {'pCa_{50}', 'n_H', 'activation', 'bound heads', 'ATP/head', 'ATP\times\sigma', 'diastolic ATP'};
qv = @(r) [10^(6 - r.pCa50), r.nH, r.act, r.boundPk, r.ATP, r.ATPsig, r.ATPdia];  % Ca50 in uM for the first
q0 = qv(R(1));
for a = 2:nAl, plot(ax, 1:numel(qn), qv(R(a))./q0, 'o-', 'Color', cm(a, :), 'MarkerFaceColor', cm(a, :)); end
yline(ax, 1, 'k-');  set(ax, 'YScale', 'log', 'XTick', 1:numel(qn), 'XTickLabel', [{'Ca_{50}'}, qn(2:end)]);  xlim(ax, [0.5 numel(qn) + 0.5]);
style(ax, 'E  Derived quantities relative to fit C (log)', 'ratio to fit C');
ax = nexttile; hold(ax, 'on');
b_ = barh(ax, 1:nAl, 100*[R.rms]);  b_.FaceColor = 'flat';  b_.CData = cm;
xline(ax, 100*sem, '--', 'data SEM', 'Color', [0.4 0.4 0.4]);
set(ax, 'YTick', 1:nAl, 'YTickLabel', {alts.name}, 'YDir', 'reverse');
style(ax, 'F  Fit error (RMS, % of peak)', '');  xlabel(ax, 'RMS error (% of peak)');
exportgraphics(fig, fullfile(res, 'eht_nonunique.png'), 'Resolution', 150);
fprintf('figure written\n');

%% ---- local functions ---------------------------------------------------------------------------------
function b = applyCore(b, core)
% rateScale multiplies the core's cycling rates (the set updateRates scales; getParams only applies
% params.xrate when updateModifiers is true, which the twitch wrapper does not use)
for f = string(fieldnames(core))'
    if f == "rateScale"
        for r = ["kah" "kamh" "ka" "kd" "k1" "k_1" "k2"], b.(r) = core.rateScale*b.(r); end
    else
        b.(f) = core.(f);
    end
end
end

function x = expand(z, x0, m)
x = x0;  x(m) = z;
end

function th = toTh6(x, fix, names)
% x(3) ('t0') shifts the whole model beat relative to the recorded stimulus (FitCEHTTwitch.m)
th = fix;
for i = 1:numel(names), th.(names{i}) = x(i); end
th.tf_k_on = 10^th.tf_k_on;  th.tf_k_off = 10^th.tf_k_off;
th.shift = th.t0;  th.t0 = 0;
end

function r = resid6(x, base, fix, names, tb, Fd)
fm  = simEHTTwitch(base, toTh6(x, fix, names), tb);
sig = max((fm' * Fd)/(fm' * fm), 0);
r   = (sig*fm - Fd)/max(Fd);
end

function F = steadyF(base, th, pc)
P = base;  P.CaModel = 'twoState';  P.CaInput = 'const';  P.pCa = pc;
for f = ["tf_k_on" "tf_k_off" "tf_K_coop"], if isfield(th, f), P.(f) = th.(f); end, end
P = getParams(P, [], true, false);  P.Vums = 0;  P.v = 0;
[~, y] = ode15s(@dPUdT_CombinedTransitions, [0 3], P.PU0(:), odeset('AbsTol', 1e-4, 'RelTol', 1e-3, 'MaxStep', 2e-3), P);
[~, o] = dPUdT_CombinedTransitions(3, y(end, :)', P);  F = o(1);
end

function [p50, nH] = hill(pc, F)
y = (F - F(1))/(F(end) - F(1));
q = lsqcurvefit(@(q, x) 1./(1 + 10.^(q(2)*(x - q(1)))), [5.8 3], pc, y, [4 0.5], [8 10], optimoptions('lsqcurvefit', 'Display', 'off'));
p50 = q(1);  nH = q(2);
end

function style(ax, ttl, yl)
title(ax, ttl, 'FontWeight', 'normal');  ylabel(ax, yl);  box(ax, 'off');  grid(ax, 'on');
set(ax, 'GridAlpha', 0.12, 'XColor', [0.32 0.32 0.31], 'YColor', [0.32 0.32 0.31], 'TickDir', 'out');
end
