% SrxRateBridge.m -- put the sarcomere model's SRX next to the mant-ATP chase.
%
% The detached pool of dPUdT_CombinedTransitions and the kinetic core of the
% mant-ATP model (../mantATP-SRX, Modelica/mantATP.mo; chase subsystem in
% identifiability-analysis/chase_identifiability.py) have the SAME topology:
%
%     sarcomere      mant       meaning
%     P_SR           S_T        SRX . ATP
%     P_SRD          S_D        SRX . ADP.Pi
%     PT (PuATP)     D_T        DRX . ATP
%     PD (PuR)       D_D        DRX . ADP.Pi   (the actin-ready head)
%
%     PT  -> P_SR    ksr0*e^(-F/sigma2)       D_T -> S_T  kDT   (recruitment into SRX)
%     P_SRD -> PD    kmsrd*e^(+F/sigma_srd1)  S_D -> D_D  kSD   (escape from SRX)
%     PD  -> P_SRD   ksrd*e^(-F/sigma_srd2)   (none)
%     P_SR <-> P_SRD ksr2srd / ksrd2sr        kH / kmH    (hydrolysis in SRX)
%     PT <-> PD      kah / kamh               kH / kmH    (hydrolysis in DRX)
%     PD -> actin    ka*N_ov (cycle)          D_D -> R  kleak (basal Pi/ADP release)
%
% So the two can be compared rate for rate. This script (no fibre battery, a
% few short isometric runs only):
%   1. rate table: model at F = 0, at rest (passive only) and isometric, vs mant;
%   2. SRX fraction of the model's detached pool as a function of load F;
%   3. a VIRTUAL mant chase run on the model's own rates (label released from
%      DRX.ADP.Pi at the mant kleak): what the assay would see if the model's
%      rates held at rest, and how far they must be slowed for the chase to
%      turn double-exponential;
%   4. the same chase in the ACTIVE state (label released by the cross-bridge
%      cycle), lumped from the isometric fluxes.
% Output: results/rate_bridge.{txt,json,png}
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rs = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
if ~isfolder(rs); mkdir(rs); end

p0 = getParams(loadParams('robust1_seed'), [], true, false);
mf = str2func(p0.modelFcn);

%% --- isometric (active) and rest reference states at SL 2.2 --------------
pa = p0; pa.Velocity = [0; 0]; pa.SL0 = 2.2; pa = getParams(pa, pa.g, true, false);
[~, oa] = evaluateModel(mf, [-3 0], pa);
pr = pa; pr.ka = 0; pr = getParams(pr, pr.g, true, false);   % no attachment = relaxed
[~, orr] = evaluateModel(mf, [-5 0], pr);
S = @(o) struct('F', o.Force(end), 'SR', o.SR(end), 'SRD', o.SRD(end), 'PT', o.PuATP(end), ...
    'PD', o.PuR(end), 'att', o.p1_0(end) + o.p2_0(end), 'RD1', o.RD1(end), 'R1D', o.R1D(end), ...
    'R2T', o.R2T(end), 'RTD', o.RTD(end), 'RPT2SR', o.RPT2SR(end), 'RSRD2PD', o.RSRD2PD(end));
A = S(oa); R = S(orr);
A.SRX = A.SR + A.SRD; R.SRX = R.SR + R.SRD;

%% --- 1. rate table ---------------------------------------------------------
mant = struct('kH', 18, 'kmH', 1.8, 'kSD', 0.0047, 'kDT', 0.0244, 'kleak', 0.027, 'kR', 16);
kr = @(F) struct( ...
    'rec',  p0.ksr0  * exp(-F / p0.sigma2), ...        % DRX.ATP   -> SRX.ATP
    'esc',  p0.kmsrd * exp( F / p0.sigma_srd1), ...    % SRX.ADPPi -> DRX.ADPPi
    'dsr',  p0.ksrd  * exp(-F / p0.sigma_srd2), ...    % DRX.ADPPi -> SRX.ADPPi
    'hS',   p0.ksr2srd, 'mhS', p0.ksrd2sr, 'hD', p0.kah, 'mhD', p0.kamh);
Fset = [0, R.F, A.F];  Fnm = {'F = 0', sprintf('rest (%.1f kPa passive)', R.F), sprintf('isometric (%.1f kPa)', A.F)};
fid = fopen(fullfile(rs, 'rate_bridge.txt'), 'w'); c = onCleanup(@() fclose(fid));
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('Rate bridge: robust1_seed detached pool vs the mant-ATP fit (Walklate Fig 1A)\n\n');
out('%-34s %10s %10s %10s | %10s | model(F=0)/mant\n', 'step', Fnm{:}, 'mant');
rows = {'DRX->SRX recruitment (ksr0 | kDT)', 'rec', mant.kDT
        'SRX->DRX escape (kmsrd | kSD)',     'esc', mant.kSD
        'DRX.D->SRX.D (ksrd | none)',        'dsr', 0
        'hydrolysis in SRX (ksr2srd | kH)',  'hS',  mant.kH
        'hydrolysis in DRX (kah | kH)',      'hD',  mant.kH };
for i = 1:size(rows, 1)
    v = arrayfun(@(F) kr(F).(rows{i,2}), Fset);
    out('%-34s %10.4g %10.4g %10.4g | %10.4g | %8.3g\n', rows{i,1}, v, rows{i,3}, v(1)/rows{i,3});
end
Kmod = @(F) (kr(F).rec / kr(F).esc);      % crude SRX/DRX drive, same for mant
out('\nrecruit/escape ratio: model F=0 %.3g, rest %.3g, iso %.3g;  mant %.3g\n', ...
    Kmod(0), Kmod(R.F), Kmod(A.F), mant.kDT / mant.kSD);
out('force swing of the model''s ratio from rest to isometric: x%.3g\n', Kmod(R.F)/Kmod(A.F));
out('load needed on escape alone to span the mant->model gap (x%.3g): %.0f kPa (sigma_srd1 %.1f)\n', ...
    kr(0).esc/mant.kSD, p0.sigma_srd1*log(kr(0).esc/mant.kSD), p0.sigma_srd1);
out('\nSRX (SR+SRD) of the model: rest %.3f (SR %.3f, SRD %.3f) | isometric %.3f (SR %.3f, SRD %.3f)\n', ...
    R.SRX, R.SR, R.SRD, A.SRX, A.SR, A.SRD);
out('mant relaxed SRX (aging state S_T+S_D, chase_identifiability.json): 0.345\n');
out('isometric flux through SRX %.2f /s per head, pool %.3f -> mean SRX dwell %.1f ms (mant 1/kSD = %.0f s)\n', ...
    A.RPT2SR, A.SRX, 1e3*A.SRX/A.RPT2SR, 1/mant.kSD);

%% --- 2. SRX(F) of the closed detached pool (no actin) ---------------------
Fg = linspace(0, 100, 201); srxF = zeros(size(Fg));
for i = 1:numel(Fg)
    Q = iDetachedGenerator(kr(Fg(i)), 0, 0);    % closed 4-state chain
    pi_ = iStationary(Q); srxF(i) = pi_(1) + pi_(2);
end

%% --- 3. virtual mant chase on the model's rest rates ----------------------
% labelled states [S_T S_D D_T D_D]; label leaves D_D at kleak. Initial condition:
% stationary distribution of the aging cycle (D_D -kleak-> R -kR-> D_T).
kRest = kr(R.F);
[lamM, ampM] = iChase(iMantRates(mant), mant.kleak, mant.kR);
[lamR, ampR] = iChase(kRest, mant.kleak, mant.kR);
out('\nVirtual chase (label out of DRX.ADP.Pi at kleak = %.3g /s)\n', mant.kleak);
out('  mant rates  : slow rates %s /s, amplitudes %s\n', mat2str(lamM(1:2)', 3), mat2str(ampM(1:2)', 3));
out('  model rest  : slow rates %s /s, amplitudes %s\n', mat2str(lamR(1:2)', 3), mat2str(ampR(1:2)', 3));
% slow the model's SRX exchange (rec, esc, dsr) by mu at fixed equilibrium
mu = logspace(0, -6, 121); fSlow = zeros(size(mu)); rr = zeros(size(mu)); lslow = zeros(size(mu));
for i = 1:numel(mu)
    k = kRest; k.rec = k.rec*mu(i); k.esc = k.esc*mu(i); k.dsr = k.dsr*mu(i);
    [lam, amp] = iChase(k, mant.kleak, mant.kR);
    fSlow(i) = amp(1) / sum(amp(1:2)); rr(i) = lam(2)/lam(1); lslow(i) = lam(1);
end
iX = find(rr >= 3 & fSlow >= 0.1 & fSlow <= 0.9, 1);    % first mu: both slow modes carry >=10 % and are >=3x apart
out('  model rest, exchange slowed by mu: two phases resolvable (rate ratio>=3, each phase>=10%%) from mu = %.2g\n', mu(iX));
out('    -> escape %.3g /s, recruitment %.3g /s at that mu (mant %.3g / %.3g)\n', kRest.esc*mu(iX), kRest.rec*mu(iX), mant.kSD, mant.kDT);
[~, iM] = min(abs(log(lslow ./ lamM(1))));
out('  slow rate matches mant (%.3g /s) at mu = %.2g: slow amplitude %.2f (mant %.2f)\n', lamM(1), mu(iM), fSlow(iM), ampM(1)/sum(ampM(1:2)));

%% --- 4. virtual chase in the ACTIVE isometric state -----------------------
% attached pool lumped to one state Att: PD -> Att at RD1/PD, Att -> PD (weak
% detachment, keeps ADP.Pi) at R1D/att, Att -> PT at R2T/att (ADP released,
% fresh ATP bound = label gone).
kA = kr(A.F); kon = A.RD1/A.PD; kback = A.R1D/A.att; koff = A.R2T/A.att;
[lamA, ampA] = iChaseActive(kA, kon, kback, koff);
out('\nActive isometric chase (label leaves with ADP at strong detachment):\n');
out('  kon %.1f, kback %.1f, koff %.1f /s;  chase rates %s /s, amplitudes %s\n', kon, kback, koff, ...
    mat2str(lamA', 3), mat2str(ampA', 2));
out('  -> in contraction the SRX is not a kinetically separate population: every mode is < %.0f ms\n', 1e3/min(lamA));

%% --- figure & json --------------------------------------------------------
fg = figure('Color', 'w', 'Position', [50 50 1500 450], 'Visible', 'off');
tiledlayout(fg, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; box on; grid on;
plot(Fg, srxF, 'k-', 'LineWidth', 2);
plot(R.F, R.SRX, 'bo', 'MarkerFaceColor', 'b'); plot(A.F, A.SRX, 'rs', 'MarkerFaceColor', 'r');
yline(0.345, 'b--', 'mant relaxed 0.345'); xlabel('load F_{SR} (kPa)'); ylabel('SRX fraction (SR + SRD)');
legend({'closed detached pool (no actin)', 'model at rest', 'model isometric'}, 'Location', 'northeast');
title('Model SRX vs load');
nexttile; hold on; box on; grid on;
tt = logspace(-3, 3.5, 400);
yM = iChaseTrace(iMantRates(mant), mant.kleak, mant.kR, tt); yR = iChaseTrace(kRest, mant.kleak, mant.kR, tt);
kS = kRest; kS.rec = kS.rec*mu(iM); kS.esc = kS.esc*mu(iM); kS.dsr = kS.dsr*mu(iM);
yS = iChaseTrace(kS, mant.kleak, mant.kR, tt);
semilogx(tt, yM, 'b-', tt, yR, 'r-', tt, yS, 'k--', 'LineWidth', 2); set(gca, 'XScale', 'log');
xlabel('chase time (s)'); ylabel('bound label'); title('Virtual mant chase at rest');
legend({'mant fit', 'model rest rates', sprintf('model, exchange \\times %.1g', mu(iM))}, 'Location', 'southwest');
nexttile; hold on; box on; grid on;
ok2 = rr >= 3 & fSlow >= 0.1 & fSlow <= 0.9;
if any(ok2); patch(mu([find(ok2,1,'last') find(ok2,1,'last') find(ok2,1) find(ok2,1)]), [0 1 1 0], [0.85 0.92 1], 'EdgeColor', 'none'); end
semilogx(mu, fSlow, 'k-', 'LineWidth', 2); set(gca, 'XScale', 'log', 'YLim', [0 1.05]);
ylabel('share of the slowest phase'); xlabel('\mu (scale on the model''s SRX exchange rates, rest)');
text(mu(find(ok2,1)), 0.15, '  two phases resolvable', 'Color', [0.2 0.3 0.8]);
xline(mu(iM), 'b--', 'slow rate = mant'); title('When does the chase become two-phase?');
exportgraphics(fg, fullfile(rs, 'rate_bridge.png'), 'Resolution', 110);

J = struct('rest', R, 'iso', A, 'mant', mant, 'chase_mant', struct('rates', lamM(1:2), 'amp', ampM(1:2)), ...
    'chase_model_rest', struct('rates', lamR(1:2), 'amp', ampR(1:2)), 'mu_two_phase', mu(iX), ...
    'mu_match_slow', mu(iM), 'chase_active', struct('rates', lamA, 'amp', ampA), ...
    'SRX_vs_F', struct('F', Fg, 'SRX', srxF));
fid2 = fopen(fullfile(rs, 'rate_bridge.json'), 'w'); fprintf(fid2, '%s', jsonencode(J, 'PrettyPrint', true)); fclose(fid2);
fprintf('DONE SrxRateBridge\n');

%% ===========================================================================
function k = iMantRates(m)
    k = struct('rec', m.kDT, 'esc', m.kSD, 'dsr', 0, 'hS', m.kH, 'mhS', m.kmH, 'hD', m.kH, 'mhD', m.kmH);
end

function Q = iDetachedGenerator(k, kleak, kR)
% generator (columns = from) over [S_T S_D D_T D_D (R)]; R only if kR > 0
    Q = [ -(k.hS + 0*k.rec),  k.mhS,           k.rec,              0
           k.hS,             -(k.mhS + k.esc), 0,                  k.dsr
           0,                 0,              -(k.hD + k.rec),     k.mhD
           0,                 k.esc,           k.hD,             -(k.mhD + k.dsr + kleak) ];
    if kR > 0
        Q(5, :) = [0 0 0 kleak]; Q(:, 5) = [0; 0; kR; 0; -kR];
    end
end

function p = iStationary(Q)
    [V, D] = eig(Q); [~, j] = min(abs(diag(D))); p = real(V(:, j)); p = p / sum(p);
end

function [lam, amp] = iChase(k, kleak, kR)
    x0 = iStationary(iDetachedGenerator(k, kleak, kR)); x0 = x0(1:4) / sum(x0(1:4));
    Q = iDetachedGenerator(k, kleak, 0);
    [V, D] = eig(Q); lam = -real(diag(D)); c = V \ x0; amp = real(sum(V, 1)' .* c);
    [lam, o] = sort(lam); amp = amp(o);
end

function y = iChaseTrace(k, kleak, kR, t)
    x0 = iStationary(iDetachedGenerator(k, kleak, kR)); x0 = x0(1:4) / sum(x0(1:4));
    Q = iDetachedGenerator(k, kleak, 0); y = zeros(size(t));
    for i = 1:numel(t); y(i) = sum(expm(Q * t(i)) * x0); end
end

function [lam, amp] = iChaseActive(k, kon, kback, koff)
% states [S_T S_D D_T D_D Att]; label leaves from Att at koff
    Q = iDetachedGenerator(k, 0, 0);
    Q(5, 5) = 0; Q(4, 4) = Q(4, 4) - kon; Q(5, :) = [0 0 0 kon -(kback + koff)]; Q(:, 5) = [0; 0; 0; kback; -(kback + koff)];
    Qc = Q; Qc(3, 5) = koff;                          % closed cycle for the steady state
    x0 = iStationary(Qc);
    [V, D] = eig(Q); lam = -real(diag(D)); c = V \ x0; amp = real(sum(V, 1)' .* c);
    [lam, o] = sort(lam); amp = amp(o);
end
