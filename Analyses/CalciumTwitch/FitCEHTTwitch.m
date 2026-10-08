%% FitCEHTTwitch.m -- fit C (k_off free) in close-up: states, ATP, identifiability, other basins
%
% Fit C = the K3 refit of TestKoffEHTTwitch.m (results/eht_koff_bound.mat): 6 free parameters
% tau1, tau2, t0, tf_k_on, tf_K_coop, tf_k_off (+ sigma in closed form), Ca 0.10 -> 1.0 uM,
% ExampleHighATP core unchanged, SL 2.0 um, isometric.
%   C0  refit with t0 allowed negative (stimulus not synchronised with force); best = fit C
%   C1  twitch over 5 beats (periodicity), all internal states, net ATP turnover
%   C2  residual vs data noise, forward-difference Jacobian -> local CI and correlations
%   C3  steady force-pCa (3 s per level: the filament is slow) and the steady ATP rate at pCa 4.51
%   C4  one-at-a-time sensitivities at fit C (shape, force, ATP)
%   C5  other basins: refits with k_off free from low-activation starts (cached, ~30 min)
%   C6  close-up figure -> results/eht_fitC.png; numbers -> results/eht_fitC.mat
% Report: EHT_identification_report.md.

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
res  = fullfile(here, 'results');
L  = load(fullfile(res, 'eht_controls.mat'));  D = L.D;
KB = load(fullfile(res, 'eht_koff_bound.mat'));
PR = load(fullfile(res, 'eht_profile.mat'));
base = loadParams(fullfile(root, 'params', 'ExampleHighATP.m'));  base.SL0 = 2.0;
P0 = getParams(base, [], true, false);
sel = D.tb >= 0;  tb = D.tb(sel);  Fd = D.meanF(sel);  sdF = D.sdF(sel);
tfine = (0:0.0025:1)';
names6 = {'tau1', 'tau2', 't0', 'tf_k_on', 'tf_K_coop', 'tf_k_off'};
fix = struct('Cadia', 0.10, 'Camax', 1.0);
lb6 = [0.010 0.080 -0.150 log10(0.1) 1 log10(0.5)];   % t0 may be negative: the stimulus record is not synchronised
ub6 = [0.150 0.300 0.060 log10(3000) 20 log10(2000)];
xC = KB.x;
if isempty(gcp('nocreate')), refreshPool(5); end

%% C0 -- refit with the Ca onset free to precede the recorded stimulus ---------------------------
% Fit C (K3) had t0 on its bound of 0. The stimulus trace is not synchronised with force, so t0 is
% an unknown offset: it now shifts the whole model beat (simEHTTwitch 'shift'), -150..+60 ms.
cacheC0 = fullfile(res, 'eht_fitC_shift.mat');
x0s = [xC; xC; xC];  x0s(:, 3) = [0; -0.03; -0.08];
if exist(cacheC0, 'file'), Cc = load(cacheC0);  X0 = Cc.X0;  S0 = Cc.S0;
else
    lopt0 = optimoptions('lsqnonlin', 'Display', 'off', 'UseParallel', false, 'FiniteDifferenceStepSize', 1e-3, ...
        'StepTolerance', 1e-4, 'FunctionTolerance', 1e-6, 'MaxIterations', 60);
    X0 = zeros(size(x0s));  S0 = zeros(size(x0s, 1), 1);
    parfor q = 1:size(x0s, 1)
        [X0(q, :), S0(q)] = lsqnonlin(@(x) resid6(x, base, fix, names6, tb, Fd), x0s(q, :), lb6, ub6, lopt0); %#ok<PFBNS>
    end
    save(cacheC0, 'X0', 'S0', 'x0s', 'lb6', 'ub6');
end
for q = 1:size(X0, 1)
    t = toTh6(X0(q, :), fix, names6);
    fprintf('C0 start t0 %+.0f ms: SSR %.5f (fit C at t0 >= 0: %.5f) | t0 %+.1f ms tau %.0f/%.0f k_on %.2f K_coop %.2f k_off %.2f\n', ...
        1e3*x0s(q, 3), S0(q), KB.rn, 1e3*t.shift, 1e3*t.tau1, 1e3*t.tau2, t.tf_k_on, t.tf_K_coop, t.tf_k_off);
end
[sBest, qb] = min(S0);
if sBest < KB.rn, xC = X0(qb, :);  fprintf('C0: fit C updated to the shifted optimum\n'); end
thC = toTh6(xC, fix, names6);

%% C1 -- twitch, states, ATP -----------------------------------------------------------------
th5 = thC;  th5.nb = 5;  th5.pre = 2;
[fm, S] = simEHTTwitch(base, th5, tfine);
fmd = interp1(tfine, fm, tb);  sig = (fmd'*Fd)/(fmd'*fmd);
m = shape(tfine, fm);  mD = shape(tb, Fd);
ATP = trapz(tfine, S.hyd);  ATPdrx = trapz(tfine, S.hydDRX);  ATPsrx = trapz(tfine, S.hydSRX);
koffEff = thC.tf_k_off*(1 + thC.tf_K_coop*S.N);
fprintf('C1 fit C: sigma %.2f uN/kPa | Fdev %.2f kPa (x sigma = %.0f uN, data %.0f) | Fdia %.3f kPa (passive %.3f) | drift over beats %.2g kPa\n', ...
    sig, max(fm), sig*max(fm), max(Fd), S.Fdia, S.Fpas(1), S.drift);
fprintf('C1 metrics model/data: TTP %.0f/%.0f  10->pk %.0f/%.0f  RT50 %.0f/%.0f  RT90 %.0f/%.0f  W50 %.0f/%.0f ms\n', ...
    1e3*[m.TTP mD.TTP m.C10toPk mD.C10toPk m.RT50 mD.RT50 m.RT90 mD.RT90 m.W50 mD.W50]);
[~, iCa] = max(S.Ca);  [~, iG] = max(S.gate);  [~, iF] = max(fm);  [~, iB] = max(S.bound);
fprintf('C1 peaks (ms): Ca %.0f gate %.0f bound %.0f force %.0f | Ca 950 ms %.3f uM | 1-N %.4f -> %.3f | gate %.4f -> %.3f | bound %.4f -> %.3f\n', ...
    1e3*tfine([iCa iG iB iF]), interp1(tfine, S.Ca, 0.95), 1 - S.N(1), max(1 - S.N), S.gate(1), max(S.gate), S.bound(1), max(S.bound));
fprintf('C1 heads: PD %.3f -> min %.3f | PT %.3f -> max %.3f | SRX %.3f -> min %.3f | k_off,eff %.0f..%.0f /s\n', ...
    S.PD(1), min(S.PD), S.PT(1), max(S.PT), S.SRX(1), min(S.SRX), min(koffEff), max(koffEff));
fprintf('C1 ATP: %.3f /head/beat (DRX step %.3f, SRX step %.3f) | diastolic rate %.4f /s | peak %.2f /s | basal share %.0f %% | ATP x sigma %.2f | economy %.1f kPa/(ATP/head/beat)\n', ...
    ATP, ATPdrx, ATPsrx, S.hyd(1), max(S.hyd), 100*S.hyd(1)/ATP, ATP*sig, max(fm)/ATP);

%% C2 -- residual, noise, Jacobian -------------------------------------------------------------
r0 = resid6(xC, base, fix, names6, tb, Fd);
sem = sqrt(mean(D.sdFn(sel).^2))/sqrt(nnz(D.isCtrl));
ac1 = corr(r0(1:end-1), r0(2:end));
J = zeros(numel(r0), numel(xC));
parfor i = 1:numel(xC)
    xi = xC;  h = 1e-3*max(1, abs(xC(i)));  xi(i) = xi(i) + h;
    J(:, i) = (resid6(xi, base, fix, names6, tb, Fd) - r0)/h; %#ok<PFBNS>
end
s2 = sum(r0.^2)/(numel(r0) - numel(xC));
Cv = s2*pinv(J'*J);  sd = sqrt(diag(Cv))';  Rho = Cv./(sd'*sd);  sv = svd(J);
fprintf('C2 rms %.2f %% | template SEM %.2f %% (ratio %.1f) | lag-1 autocorrelation %.2f | cond(J) %.0f\n', ...
    100*sqrt(mean(r0.^2)), 100*sem, sqrt(mean(r0.^2))/sem, ac1, sv(1)/sv(end));
for i = 1:numel(xC)
    fprintf('C2 %-10s %9.4g  +-%.3g (95%%, %s)\n', names6{i}, xC(i), 1.96*sd(i), ternary(any(i == [4 6]), 'log10 units', 'linear'));
end
disp(array2table(round(Rho, 2), 'VariableNames', names6, 'RowNames', names6));

%% C3 -- steady force-pCa and steady ATP rate ------------------------------------------------------
pc = [7 6.6 6.3 6.1 5.9 5.7 5.5 5.3 5 4.51];  Fss = zeros(size(pc));  Hss = Fss;
parfor k = 1:numel(pc), [Fss(k), Hss(k)] = steadyState(base, thC, pc(k), 3); end %#ok<PFBNS>
[p50, nH] = hill(pc, Fss);
act = max(fm)/(Fss(end) - S.Fdia);
fprintf('C3 pCa50 %.2f (Ca50 %.2f uM) nH %.2f | Fsat %.1f kPa | twitch = %.0f %% of saturating developed force | ATP rate at pCa 4.51 %.2f /head/s (twitch peak %.2f, beat mean %.2f)\n', ...
    p50, 10^(6 - p50), nH, Fss(end), 100*act, Hss(end), max(S.hyd), ATP);

%% C4 -- sensitivities at fit C ----------------------------------------------------------------------
V = {'Cadia 0.05 uM',         struct('Cadia', 0.05);
     'Cadia 0.20 uM',         struct('Cadia', 0.20);
     'Ca amplitude x1.5',     struct('Camax', fix.Cadia + 1.5*(fix.Camax - fix.Cadia));
     'Ca amplitude /1.5',     struct('Camax', fix.Cadia + (fix.Camax - fix.Cadia)/1.5);
     'k_off x2 alone',        struct('tf_k_off', 2*thC.tf_k_off);
     'k_on,k_off x2 (same ss)', struct('tf_k_on', 2*thC.tf_k_on, 'tf_k_off', 2*thC.tf_k_off);
     'k_on,k_off /2 (same ss)', struct('tf_k_on', thC.tf_k_on/2, 'tf_k_off', thC.tf_k_off/2);
     'K_coop 10',             struct('tf_K_coop', 10);
     'kSE x0.1',              struct('kSE', 0.1*base.kSE);
     'ka x0.5',               struct('ka', 0.5*base.ka);
     'k2 x0.5',               struct('k2', 0.5*base.k2);
     'k1 x0.5',               struct('k1', 0.5*base.k1);
     'kmsrd x3',              struct('kmsrd', 3*base.kmsrd);
     'ksr0 x3',               struct('ksr0', 3*base.ksr0)};
SF = cell(size(V, 1), 1);  SH = zeros(size(V, 1), 1);  SS = cell(size(V, 1), 1);
parfor k = 1:size(V, 1)
    th = thC;  b = base;  ch = V{k, 2};
    for fn = string(fieldnames(ch))'
        if isfield(b, fn) && ~startsWith(fn, "tf_"), b.(fn) = ch.(fn); else, th.(fn) = ch.(fn); end
    end
    [SF{k}, Sk] = simEHTTwitch(b, th, tfine); %#ok<PFBNS>
    SH(k) = trapz(tfine, Sk.hyd);  SS{k} = struct('SRXmin', min(Sk.SRX), 'bound', max(Sk.bound));
end
sens = struct('name', V(:,1), 'F', SF, 'ATP', num2cell(SH), 'st', SS);
for k = 1:numel(sens)
    q = shape(tfine, sens(k).F);  fk = interp1(tfine, sens(k).F, tb);  sk = (fk'*Fd)/(fk'*fk);
    sens(k).m = q;  sens(k).rms = sqrt(mean((sk*fk - Fd).^2))/max(Fd);
    fprintf('C4 %-24s Fdev %+6.1f %%  ATP %+6.1f %%  TTP %3.0f  10->pk %3.0f  RT50 %3.0f  RT90 %3.0f ms  shape rms %4.1f %%  SRX min %.3f\n', ...
        sens(k).name, 100*(max(sens(k).F)/max(fm) - 1), 100*(sens(k).ATP/ATP - 1), 1e3*q.TTP, 1e3*q.C10toPk, 1e3*q.RT50, 1e3*q.RT90, ...
        100*sens(k).rms, sens(k).st.SRXmin);
end

%% C5 -- other basins: k_off free from low-activation starts -----------------------------------------
% starts: the profile's k_on 200 and 398 points with k_on, k_off both divided by 50 (same steady
% state, slow filament), and one generic start
xs = zeros(3, 6);  konStart = [200 398];
for q = 1:2
    j = find(abs(PR.X(:,4) - log10(konStart(q))) < 0.01, 1);  x = PR.X(j, :);
    xs(q, :) = [x(1:3), x(4) - log10(50), x(5), log10(P0.tf_k_off/50)];
end
xs(3, :) = [0.06 0.20 0.02 log10(3) 8 log10(3)];
cacheC5 = fullfile(res, 'eht_fitC_basins.mat');
if exist(cacheC5, 'file'), Cc = load(cacheC5);  XS = Cc.XS;  SSR = Cc.SSR;
else
    lopt = optimoptions('lsqnonlin', 'Display', 'off', 'UseParallel', false, 'FiniteDifferenceStepSize', 1e-3, ...
        'StepTolerance', 1e-4, 'FunctionTolerance', 1e-6, 'MaxIterations', 60);
    XS = zeros(size(xs));  SSR = zeros(size(xs, 1), 1);
    parfor q = 1:size(xs, 1)
        [XS(q, :), SSR(q)] = lsqnonlin(@(x) resid6(x, base, fix, names6, tb, Fd), min(max(xs(q, :), lb6), ub6), lb6, ub6, lopt); %#ok<PFBNS>
    end
    save(cacheC5, 'XS', 'SSR', 'xs');
end
basin = struct('x', num2cell(XS, 2), 'SSR', num2cell(SSR));
for q = 1:numel(basin)
    t = toTh6(XS(q, :), fix, names6);  [fq, Sq] = simEHTTwitch(base, t, tfine);
    fqd = interp1(tfine, fq, tb);  sq = (fqd'*Fd)/(fqd'*fqd);
    basin(q).sigma = sq;  basin(q).F = fq;  basin(q).gatePk = max(Sq.gate);  basin(q).ATP = trapz(tfine, Sq.hyd);
    fprintf('C5 start %d: SSR %.5f (fit C %.5f) | tau %.0f/%.0f t0 %.0f ms k_on %.2f K_coop %.2f k_off %.2f | sigma %.1f gate pk %.2f ATP %.2f ATPxsigma %.1f\n', ...
        q, SSR(q), KB.rn, 1e3*t.tau1, 1e3*t.tau2, 1e3*t.t0, t.tf_k_on, t.tf_K_coop, t.tf_k_off, sq, basin(q).gatePk, basin(q).ATP, basin(q).ATP*sq);
end
save(fullfile(res, 'eht_fitC.mat'), 'X0', 'S0', 'thC', 'xC', 'fm', 'S', 'sig', 'm', 'mD', 'ATP', 'ATPdrx', 'ATPsrx', 'koffEff', ...
    'r0', 'sem', 'ac1', 'J', 'Cv', 'Rho', 'pc', 'Fss', 'Hss', 'p50', 'nH', 'act', 'sens', 'basin', 'tfine', 'tb', 'Fd');

%% C6 -- close-up figure -------------------------------------------------------------------------------
Cb = [42 120 214]/255;  Co = [235 104 52]/255;  Cg = [27 175 122]/255;  Cy = [237 161 0]/255;  Cp = [150 90 200]/255;  G = [0.55 0.55 0.53];
fig = figure('Color', 'w', 'Position', [30 30 1650 1150], 'Visible', 'off');
tl_ = tiledlayout(fig, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_, sprintf(['Fit C in close-up: Plate 6 control EHTs, 1 Hz  (k_{off} %.1f s^{-1}, k_{on} %.2f \\muM^{-1}s^{-1}, ' ...
    'K_{coop} %.1f, Ca 0.10 \\rightarrow 1.0 \\muM, \\sigma %.1f \\muN/kPa)'], thC.tf_k_off, thC.tf_k_on, thC.tf_K_coop, sig), 'FontWeight', 'normal');

ax = nexttile([1 2]); hold(ax, 'on');                                    % A
fill(ax, 1e3*[tb; flipud(tb)], [Fd - sdF; flipud(Fd + sdF)], G, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'DisplayName', 'control wells, mean \pm SD (n = 15)');
plot(ax, 1e3*tb, Fd, 'k.', 'MarkerSize', 10, 'DisplayName', 'control mean');
plot(ax, 1e3*tfine, sig*fm, '-', 'Color', Co, 'LineWidth', 2, 'DisplayName', sprintf('fit C, rms %.2f %%', 100*sqrt(mean(r0.^2))));
style(ax, 'A  Twitch: data and fit C', 'active force (\muN)');  legend(ax, 'Box', 'off');  xlabel(ax, 'time from stimulus (ms)');

ax = nexttile; hold(ax, 'on');                                           % B residual
fill(ax, 1e3*[tb; flipud(tb)], 100*[-2*sem*ones(size(tb)); 2*sem*ones(size(tb))], G, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(ax, 1e3*tb, 100*r0, '-o', 'Color', Co, 'MarkerSize', 3, 'MarkerFaceColor', Co);  yline(ax, 0, '-', 'Color', G);
style(ax, sprintf('B  Residual (grey: \\pm2 SEM); lag-1 autocorr %.2f', ac1), 'model - data (% of peak)');  xlabel(ax, 'time from stimulus (ms)');

ax = nexttile; hold(ax, 'on');                                           % C log tail
plot(ax, 1e3*tb, max(Fd/max(Fd), 1e-4), 'k.', 'MarkerSize', 9);
plot(ax, 1e3*tfine, max(sig*fm/max(Fd), 1e-4), '-', 'Color', Co, 'LineWidth', 1.8);
set(ax, 'YScale', 'log');  ylim(ax, [3e-3 1.2]);
style(ax, 'C  Rise and relaxation on a log scale', 'F / F_{peak, data} (log)');  xlabel(ax, 'time from stimulus (ms)');

ax = nexttile; hold(ax, 'on');                                           % D Ca & thin filament
plot(ax, 1e3*tfine, S.Ca/max(S.Ca), '-', 'Color', Cb, 'LineWidth', 1.8, 'DisplayName', sprintf('Ca / Ca_{max} (peak %.0f ms)', 1e3*tfine(iCa)));
plot(ax, 1e3*tfine, 1 - S.N, '-', 'Color', Cp, 'LineWidth', 1.5, 'DisplayName', 'permissive sites 1-N');
plot(ax, 1e3*tfine, S.gate, '--', 'Color', Co, 'LineWidth', 1.5, 'DisplayName', 'attachment gate');
plot(ax, 1e3*tfine, fm/max(fm), ':', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'force / peak');
style(ax, 'D  Ca \rightarrow thin filament \rightarrow force', 'fraction');  legend(ax, 'Box', 'off', 'FontSize', 7, 'Location', 'east');  xlim(ax, [0 1000]);  xlabel(ax, 'time from stimulus (ms)');

ax = nexttile; hold(ax, 'on');                                           % E effective off rate
plot(ax, 1e3*tfine, koffEff, '-', 'Color', Cp, 'LineWidth', 1.8);
yyaxis(ax, 'right');  plot(ax, 1e3*tfine, 1e3./koffEff, '--', 'Color', G);  ylabel(ax, 'time constant (ms)');  ax.YColor = G;
yyaxis(ax, 'left');   ylabel(ax, 'k_{off}(1 + K_{coop}N) (s^{-1})');  ax.YColor = Cp;
title(ax, 'E  Effective site off-rate (solid) and its time constant', 'FontWeight', 'normal');  grid(ax, 'on');  box(ax, 'off');  xlabel(ax, 'time from stimulus (ms)');

ax = nexttile; hold(ax, 'on');                                           % F head pools
plot(ax, 1e3*tfine, S.PD, '-', 'Color', Cb, 'LineWidth', 1.6, 'DisplayName', 'primed PD (M.ADP.Pi)');
plot(ax, 1e3*tfine, S.bound, '-', 'Color', Co, 'LineWidth', 1.6, 'DisplayName', 'attached p1+p2');
plot(ax, 1e3*tfine, S.SRX, '-', 'Color', Cy, 'LineWidth', 1.6, 'DisplayName', 'SRX (P_{SR}+P_{SRD})');
plot(ax, 1e3*tfine, S.PT, '-', 'Color', Cg, 'LineWidth', 1.6, 'DisplayName', 'PT (M.ATP)');
style(ax, 'F  Myosin head pools', 'fraction of heads');  legend(ax, 'Box', 'off', 'FontSize', 8, 'Location', 'east');  xlabel(ax, 'time from stimulus (ms)');

ax = nexttile; hold(ax, 'on');                                           % G ATP
area(ax, 1e3*tfine, S.hyd, 'FaceColor', Cg, 'FaceAlpha', 0.25, 'EdgeColor', Cg, 'LineWidth', 1.6);
yline(ax, Hss(end), ':', 'Color', G);  ylim(ax, [0 1.7*Hss(end)]);
text(ax, 0.97, 0.97, sprintf(['AUC %.2f ATP/head/beat = mean %.2f s^{-1} at 1 Hz\\newlinediastolic rate %.3f s^{-1};  AUC \\times \\sigma = %.1f' ...
     '\\newlinedotted: steady rate at pCa 4.51, %.1f s^{-1}'], ATP, ATP, S.hyd(1), ATP*sig, Hss(end)), ...
     'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
style(ax, 'G  Net ATP turnover (forward - reverse hydrolysis)', 'ATP per head per s');  xlabel(ax, 'time from stimulus (ms)');

ax = nexttile; hold(ax, 'on');                                           % H force-pCa
pp = linspace(4.4, 7.2, 300);
plot(ax, pp, 1./(1 + 10.^(nH*(pp - p50))), '-', 'Color', Co);
plot(ax, pc, (Fss - Fss(1))/(Fss(end) - Fss(1)), 'o', 'Color', Co, 'MarkerFaceColor', Co);
xline(ax, 6 - log10(fix.Cadia), ':', 'Ca_{dia}', 'Color', G, 'LabelVerticalAlignment', 'middle');  xline(ax, 6 - log10(fix.Camax), ':', 'Ca_{max}', 'Color', G, 'LabelVerticalAlignment', 'middle');
yline(ax, act, '--', sprintf('twitch peak: %.0f %% of saturating', 100*act), 'Color', Cb, 'LabelHorizontalAlignment', 'right');
set(ax, 'XDir', 'reverse');  xlim(ax, [4.4 7.2]);
style(ax, sprintf('H  Force-pCa: pCa_{50} %.2f, n_H %.1f', p50, nH), 'normalised steady force');  xlabel(ax, 'pCa');

ax = nexttile; hold(ax, 'on');                                           % I force-Ca loop
can = (S.Ca - S.Ca(1))/(max(S.Ca) - S.Ca(1));
plot(ax, can, fm/max(fm), '-', 'Color', Co, 'LineWidth', 1.8);
plot(ax, can(1), 0, 'ko', 'MarkerFaceColor', 'k');  text(ax, 0.03, 0.05, 'stimulus', 'FontSize', 8);
style(ax, 'I  Force-Ca loop (normalised)', 'F / F_{peak}');  xlabel(ax, '(Ca - Ca_{dia}) / \DeltaCa');

ax = nexttile([1 2]); hold(ax, 'on');                                    % J sensitivities
cc = lines(numel(sens));
plot(ax, 1e3*tb, Fd/max(Fd), 'k.', 'MarkerSize', 8, 'DisplayName', 'data');
plot(ax, 1e3*tfine, fm/max(fm), '-', 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'fit C');
for k = 1:numel(sens)
    plot(ax, 1e3*tfine, sens(k).F/max(sens(k).F), '-', 'Color', cc(k,:), 'LineWidth', 1, 'DisplayName', ...
        sprintf('%s: F %+.0f %%, ATP %+.0f %%', sens(k).name, 100*(max(sens(k).F)/max(fm) - 1), 100*(sens(k).ATP/ATP - 1)));
end
style(ax, 'J  One-at-a-time sensitivities at fit C (each normalised to its peak)', 'F / F_{peak}');
legend(ax, 'Box', 'off', 'FontSize', 7, 'NumColumns', 1, 'Location', 'eastoutside');  xlabel(ax, 'time from stimulus (ms)');  xlim(ax, [0 800]);
exportgraphics(fig, fullfile(res, 'eht_fitC.png'), 'Resolution', 150);
fprintf('figure written\n');

%% ---- local functions ------------------------------------------------------------------------------
function th = toTh6(x, fix, names)
% x(3) ('t0') is applied as a shift of the whole model beat relative to the recorded stimulus;
% the transient itself starts at the model beat onset (equivalent to the old t0 for t0 >= 0)
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

function [F, H] = steadyState(base, th, pc, T)
P = base;  P.CaModel = 'twoState';  P.CaInput = 'const';  P.pCa = pc;
for f = ["tf_k_on" "tf_k_off" "tf_K_coop"], if isfield(th, f), P.(f) = th.(f); end, end
P = getParams(P, [], true, false);  P.Vums = 0;  P.v = 0;
[~, y] = ode15s(@dPUdT_CombinedTransitions, [0 T], P.PU0(:), odeset('AbsTol', 1e-4, 'RelTol', 1e-3, 'MaxStep', 2e-3), P);
[~, o, r] = dPUdT_CombinedTransitions(T, y(end, :)', P);
F = o(1);  H = (r(1) - r(2)) + (r(13) - r(14));
end

function [p50, nH] = hill(pc, F)
y = (F - F(1))/(F(end) - F(1));
q = lsqcurvefit(@(q, x) 1./(1 + 10.^(q(2)*(x - q(1)))), [5.8 3], pc, y, [4 0.5], [8 10], optimoptions('lsqcurvefit', 'Display', 'off'));
p50 = q(1);  nH = q(2);
end

function m = shape(tb, F)
[Fp, ip] = max(F);  up = F(1:ip);  dn = F(ip:end);  td = tb(ip:end);
c = @(tt, ff, lvl) interp1(ff(find(ff >= lvl, 1) + [-1 0]), tt(find(ff >= lvl, 1) + [-1 0]), lvl);
t10 = c(tb(1:ip), up, 0.1*Fp);  t50 = c(tb(1:ip), up, 0.5*Fp);
d = @(lvl) interp1(dn(find(dn <= lvl, 1) + [-1 0]), td(find(dn <= lvl, 1) + [-1 0]), lvl);
m = struct('TTP', tb(ip), 'C10toPk', tb(ip) - t10, 'RT50', d(0.5*Fp) - tb(ip), 'RT90', d(0.1*Fp) - tb(ip), 'W50', d(0.5*Fp) - t50);
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end

function style(ax, ttl, yl)
title(ax, ttl, 'FontWeight', 'normal');  ylabel(ax, yl);  box(ax, 'off');  grid(ax, 'on');
set(ax, 'GridAlpha', 0.12, 'XColor', [0.32 0.32 0.31], 'YColor', [0.32 0.32 0.31], 'TickDir', 'out');
end
