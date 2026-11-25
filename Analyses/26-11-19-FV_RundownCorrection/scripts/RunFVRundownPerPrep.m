%% RunFVRundownPerPrep  Review round on #92: rundown per preparation, and why one factor fails above 3 ML/s.
%
% Data only. Run RunFVRundown.m first (it writes results/fv_rundown.mat).
%
% A. Rundown per preparation. Each preparation has one record per [ATP], so its
%    second/first ratio is rundown x ATP effect. It becomes an estimate of that
%    preparation's rundown once the ATP effect is taken as common to the
%    preparations of a pair (the order-balanced mean of section 4 of
%    RunFVRundown):   rho_p(v) = log(F2/F1) - s_p * A(v),
%    s_p = +1 when the low [ATP] ran second, -1 when it ran first.
%    Three readings: 0.5-1 ML/s (the ATP effect is smallest there), 0-3 ML/s,
%    and the first 1.10 Lo hold of each record.
% B. Above 3 ML/s. Order-balanced second - first in kPa at all eight speeds,
%    and an offset model in kPa on all 192 points:
%        F = s_p * R^[second] * g(v, ATP) + b * [second]
%    against the single factor (b = 0) and a free factor per speed.
% C. What the points above 3 ML/s can carry: the candidate fresh curves and
%    Vmax from a Hill fit, with and without those points.
% Output: results/fv_rundown_perprep.log, .mat, figures/fv_rundown_perprep.png, figures/fv_rundown_highspeed.png

here = fileparts(mfilename('fullpath'));  ana = fileparts(here);  root = fileparts(fileparts(ana));
addpath(fullfile(root, 'Analyses', 'FV_Beard2022'));
S0 = load(fullfile(ana, 'results', 'fv_rundown.mat'));
T = S0.T;  first = S0.first;  low = S0.low;  M = S0.M;  nP = numel(first);
D = loadBeard2022FV();
logf = fullfile(ana, 'results', 'fv_rundown_perprep.log');  if isfile(logf), delete(logf); end
diary(logf);
vAll = [0 0.5 1 2 3 4 5 6];  nV = numel(vAll);
sgn = 2*(first == 8) - 1;                       % +1: 8 mM first, low ATP second
F1 = nan(nP, nV);  F2 = F1;  H1 = nan(nP, 1);  H2 = H1;
for p = 1:nP
    for j = 1:nV
        F1(p, j) = T.F(T.prep == p & T.slot == 1 & T.v == vAll(j));
        F2(p, j) = T.F(T.prep == p & T.slot == 2 & T.v == vAll(j));
    end
    H1(p) = T.H(T.prep == p & T.slot == 1 & T.cyc == 1);  H2(p) = T.H(T.prep == p & T.slot == 2 & T.cyc == 1);
end

%% A. rundown per preparation
iv = 1:5;                                       % 0-3 ML/s, where the log is safe
dl = log(F2(:, iv)./F1(:, iv));  dh = log(H2./H1);
A = nan(nP, numel(iv));  Ah = nan(nP, 1);  rb = A;   % common ATP effect and common rundown of the pair
for a = [4 2]
    i8 = low == a & first == 8;  iL = low == a & first == a;
    A(low == a, :) = repmat((mean(dl(i8, :), 1) - mean(dl(iL, :), 1))/2, sum(low == a), 1);
    rb(low == a, :) = repmat((mean(dl(i8, :), 1) + mean(dl(iL, :), 1))/2, sum(low == a), 1);
    Ah(low == a) = (mean(dh(i8)) - mean(dh(iL)))/2;
    fprintf('%d mM pairs: ATP effect low/8 at v = %s: %s ; on the 1.10 Lo hold %.3f\n', a, mat2str(vAll(iv)), mat2str(round(exp(A(find(low == a, 1), :)), 3)), exp(Ah(find(low == a, 1))));
end
rv = dl - sgn.*A;                               % rundown per preparation and speed (log)
rLow = mean(rv(:, 2:3), 2);  rAll = mean(rv, 2);  rHi = mean(rv(:, 4:5), 2);  rH = dh - sgn.*Ah;
seAll = std(rv, 0, 2)/sqrt(numel(iv));
preps = unique({D.prep}, 'stable');
fprintf('\n== Rundown per preparation (second record / first record), the ATP effect of the pair taken as common\n');
fprintf('%-3s %-7s %9s %9s %14s %9s %9s\n', 'p', 'order', '0.5-1', '2-3', '0-3 (se)', 'hold', 'date');
for p = 1:nP
    fprintf('%-3d %d -> %d   %9.3f %9.3f %9.3f (%.3f) %7.3f   %s\n', p, first(p), setdiff([8 low(p)], first(p)), exp(rLow(p)), exp(rHi(p)), exp(rAll(p)), seAll(p)*exp(rAll(p)), exp(rH(p)), preps{p}(1:min(11, end)));
end
for a = [4 2]
    m = low == a;
    fprintf('%d mM pairs: 0.5-1 ML/s mean %.3f, range %.3f-%.3f, sd of log %.3f | 0-3 ML/s mean %.3f, range %.3f-%.3f, sd of log %.3f\n', a, ...
        exp(mean(rLow(m))), exp(min(rLow(m))), exp(max(rLow(m))), std(rLow(m)), exp(mean(rAll(m))), exp(min(rAll(m))), exp(max(rAll(m))), std(rAll(m)));
end
fprintf('all 12: 0.5-1 ML/s mean %.3f, sd of log %.3f | 0-3 ML/s mean %.3f, sd of log %.3f | hold mean %.3f, sd of log %.3f\n', exp(mean(rLow)), std(rLow), exp(mean(rAll)), std(rAll), exp(mean(rH)), std(rH));
% is the spread between preparations more than the scatter of the points? F test M1 against M1p
m1 = M(strcmp({M.name}, 'M1'));  m1p = M(strcmp({M.name}, 'M1p'));  n = 120;
Fst = ((m1.rss - m1p.rss)/(m1p.k - m1.k))/(m1p.rss/(n - m1p.k));
fprintf('F test, one factor against one per preparation: F(%d, %d) = %.2f, p = %.4f\n', m1p.k - m1.k, n - m1p.k, Fst, 1 - fcdf(Fst, m1p.k - m1.k, n - m1p.k));
fprintf('expected sd of a per-preparation estimate from point scatter alone (0-3 ML/s, 5 speeds): %.3f\n', m1p.sd*sqrt(2/numel(iv)));
c1 = corr(rLow, rHi);  c2 = corr(rLow, rH);  c3 = corr(rAll, rH);
fprintf('correlation between preparations: 0.5-1 against 2-3 ML/s r = %.2f; 0.5-1 ML/s against hold r = %.2f; 0-3 ML/s against hold r = %.2f (n = 12)\n', c1, c2, c3);
% without preparation 9 (2 -> 8, the outlier of the report)
k9 = (1:nP)' ~= 9;
fprintf('without preparation 9: sd of log 0-3 ML/s %.3f, range %.3f-%.3f; r(0.5-1, 2-3) = %.2f\n', std(rAll(k9)), exp(min(rAll(k9))), exp(max(rAll(k9))), corr(rLow(k9), rHi(k9)));
% between preparations in kPa: records run first against records run second, same [ATP]
fprintf('\n== Between preparations, kPa: records run second over records run first, same [ATP] (different preparations)\n');
for a = [8 4 2]
    for v = [0 1]
        x1 = T.F(T.atp == a & T.slot == 1 & T.v == v);  x2 = T.F(T.atp == a & T.slot == 2 & T.v == v);
        q = mean(x2)/mean(x1);  se = q*sqrt(var(x1)/numel(x1)/mean(x1)^2 + var(x2)/numel(x2)/mean(x2)^2);
        fprintf('%d mM, v = %g: first %.1f (n=%d), second %.1f (n=%d), ratio %.2f +- %.2f\n', a, v, mean(x1), numel(x1), mean(x2), numel(x2), q, se);
    end
end
% curves with a factor per preparation (0-3 ML/s reading)
f0 = S0.f0;
T.Fp = T.F./exp(f0(T.prep))./exp(rAll(T.prep).*(T.slot == 2));
fprintf('\n== Corrected with a factor per preparation against the single factor: mean and CV (%%) by [ATP]\n');
for a = [8 4 2]
    fprintf('%d mM  v: one factor mean (CV) | per preparation mean (CV)\n', a);
    for v = [0 0.5 1 2 3]
        m = T.atp == a & T.v == v;  x = T.Fc(m);  z = T.Fp(m);
        fprintf('   %3.1f: %.3f (%4.1f) | %.3f (%4.1f)\n', v, mean(x), 100*std(x)/mean(x), mean(z), 100*std(z)/mean(z));
    end
end

%% B. above 3 ML/s
fprintf('\n== Order-balanced second - first record within a preparation, kPa, v = %s\n', mat2str(vAll));
d = F2 - F1;  db = nan(2, nV);  seb = db;  ex = db;  jj = 0;
for a = [4 2]
    jj = jj + 1;  i8 = low == a & first == 8;  iL = low == a & first == a;
    db(jj, :) = (mean(d(i8, :), 1) + mean(d(iL, :), 1))/2;
    seb(jj, :) = sqrt(var(d(i8, :), 0, 1)/sum(i8) + var(d(iL, :), 0, 1)/sum(iL))/2;
    ex(jj, :) = (exp(S0.rho) - 1)*(mean(F1(i8, :), 1) + mean(F1(iL, :), 1))/2;
    fprintf('%d mM pairs: %s (se %s)\n   expected from the factor %.3f: %s\n', a, mat2str(round(db(jj, :), 2)), mat2str(round(seb(jj, :), 2)), exp(S0.rho), mat2str(round(ex(jj, :), 2)));
end
dAll = mean(db, 1);  seAllv = sqrt(sum(seb.^2, 1))/2;  exAll = mean(ex, 1);
fprintf('all pairs:   %s (se %s)\n   expected: %s\n   z of (observed - expected): %s\n', mat2str(round(dAll, 2)), mat2str(round(seAllv, 2)), mat2str(round(exAll, 2)), mat2str(round((dAll - exAll)./seAllv, 1)));
F1m = (mean(F1(first == 8, :), 1) + mean(F1(first ~= 8, :), 1))/2;
fprintf('as a ratio (1 + difference / mean first-record force): %s\n', mat2str(round(1 + dAll./F1m, 3)));
% the 8 mM records alone, between preparations: second against first, group means
g1 = arrayfun(@(v) mean(T.F(T.atp == 8 & T.slot == 1 & T.v == v)), vAll);
g2 = arrayfun(@(v) mean(T.F(T.atp == 8 & T.slot == 2 & T.v == v)), vAll);
pl = polyfit(g1, g2, 1);
fprintf('8 mM group means, kPa: first %s\n                       second %s\n   straight line second = %.3f * first + %.2f kPa, largest residual %.2f kPa\n', mat2str(round(g1, 2)), mat2str(round(g2, 2)), pl(1), pl(2), max(abs(g2 - polyval(pl, g1))));

% offset model on all 192 points, kPa, weighted
Y = T.F;  pI = T.prep;  s2 = double(T.slot == 2);
[~, vI] = ismember(T.v, vAll);  aI = (T.atp == 8)*1 + (T.atp == 4)*2 + (T.atp == 2)*3;  gI = (aI - 1)*nV + vI;
gm = accumarray(gI, Y, [3*nV 1], @mean);
w = 1./(1^2 + (0.10*gm(gI)).^2);                 % 1 kPa floor plus 10 % of the cell mean
fitm = struct('name', {}, 'k', {}, 'chi2', {}, 'aic', {}, 'par', {});
[c, ~] = bilinear(Y, w, pI, gI, nP, 3*nV);
fitm(end+1) = struct('name', 'no rundown', 'k', nP + 3*nV - 1, 'chi2', c, 'aic', c + 2*(nP + 3*nV - 1), 'par', []);
obj1 = @(R) bilinear(Y./R.^s2, w.*R.^(2*s2), pI, gI, nP, 3*nV);
R1 = fminsearch(obj1, 0.88);  c = obj1(R1);
fitm(end+1) = struct('name', 'one factor', 'k', nP + 3*nV, 'chi2', c, 'aic', c + 2*(nP + 3*nV), 'par', R1);
obj2 = @(q) bilinear((Y - q(2)*s2)./q(1).^s2, w.*q(1).^(2*s2), pI, gI, nP, 3*nV);
q2 = fminsearch(obj2, [0.8 1.5], optimset('TolX', 1e-5, 'TolFun', 1e-7));  c = obj2(q2);
fitm(end+1) = struct('name', 'factor + offset', 'k', nP + 3*nV + 1, 'chi2', c, 'aic', c + 2*(nP + 3*nV + 1), 'par', q2);
s8 = s2.*(T.atp == 8);  sL = s2.*(T.atp < 8);
obj3 = @(q) bilinear((Y - q(2)*s8 - q(3)*sL)./q(1).^s2, w.*q(1).^(2*s2), pI, gI, nP, 3*nV);
q3 = fminsearch(obj3, [q2 q2(2)], optimset('TolX', 1e-5, 'TolFun', 1e-7));  c = obj3(q3);
fitm(end+1) = struct('name', 'factor + offset by [ATP] of the second record', 'k', nP + 3*nV + 2, 'chi2', c, 'aic', c + 2*(nP + 3*nV + 2), 'par', q3);
obj4 = @(q) bilinear(Y./exp(q(vI)').^s2, w.*exp(q(vI)').^(2*s2), pI, gI, nP, 3*nV);
q4 = fminsearch(obj4, log(R1)*ones(1, nV), optimset('MaxFunEvals', 4000, 'MaxIter', 4000, 'TolX', 1e-5, 'TolFun', 1e-7));  c = obj4(q4);
fitm(end+1) = struct('name', 'factor per speed', 'k', nP + 3*nV - 1 + nV, 'chi2', c, 'aic', c + 2*(nP + 3*nV - 1 + nV), 'par', exp(q4));
fprintf('\n== Models in kPa on all %d points (0-6 ML/s), weights 1/(1 kPa^2 + (10 %% of the cell mean)^2)\n', numel(Y));
for i = 1:numel(fitm), fprintf('%-48s k = %2d  chi2 = %6.1f  AIC = %6.1f  parameters %s\n', fitm(i).name, fitm(i).k, fitm(i).chi2, fitm(i).aic, mat2str(round(fitm(i).par, 3))); end
% leave one preparation out, factor + offset
qL = nan(nP, 2);
for p = 1:nP
    k = pI ~= p;  [~, ~, pk] = unique(pI(k));
    o = @(q) bilinear((Y(k) - q(2)*s2(k))./q(1).^s2(k), w(k).*q(1).^(2*s2(k)), pk, gI(k), nP - 1, 3*nV);
    qL(p, :) = fminsearch(o, q2, optimset('TolX', 1e-5, 'TolFun', 1e-7));
end
fprintf('factor + offset, leave one preparation out: R %.3f..%.3f, offset %.2f..%.2f kPa (se from the jackknife: R %.3f, offset %.2f)\n', min(qL(:, 1)), max(qL(:, 1)), min(qL(:, 2)), max(qL(:, 2)), ...
    sqrt((nP - 1)/nP*sum((qL(:, 1) - mean(qL(:, 1))).^2)), sqrt((nP - 1)/nP*sum((qL(:, 2) - mean(qL(:, 2))).^2)));
qL3 = nan(nP, 3);
for p = 1:nP
    k = pI ~= p;  [~, ~, pk] = unique(pI(k));
    o = @(q) bilinear((Y(k) - q(2)*s8(k) - q(3)*sL(k))./q(1).^s2(k), w(k).*q(1).^(2*s2(k)), pk, gI(k), nP - 1, 3*nV);
    qL3(p, :) = fminsearch(o, q3, optimset('TolX', 1e-5, 'TolFun', 1e-7));
end
jk = @(x) sqrt((nP - 1)/nP*sum((x - mean(x)).^2));
fprintf('offset by [ATP] of the second record, leave one preparation out: 8 mM second %.2f..%.2f kPa (se %.2f), low ATP second %.2f..%.2f kPa (se %.2f)\n', min(qL3(:, 2)), max(qL3(:, 2)), jk(qL3(:, 2)), min(qL3(:, 3)), max(qL3(:, 3)), jk(qL3(:, 3)));
% the same in the raw numbers: 8 mM minus low ATP within a preparation, kPa, by order
d8 = nan(nP, nV);
for p = 1:nP, for j = 1:nV, d8(p, j) = T.F(T.prep == p & T.atp == 8 & T.v == vAll(j)) - T.F(T.prep == p & T.atp < 8 & T.v == vAll(j)); end, end
fprintf('8 mM minus low ATP within a preparation, kPa, at 4, 5, 6 ML/s: 8 mM first %s (se %s) | 8 mM second %s (se %s)\n', mat2str(round(mean(d8(first == 8, 6:8), 1), 2)), mat2str(round(std(d8(first == 8, 6:8), 0, 1)/sqrt(6), 2)), mat2str(round(mean(d8(first ~= 8, 6:8), 1), 2)), mat2str(round(std(d8(first ~= 8, 6:8), 0, 1)/sqrt(6), 2)));
% high-speed force is a property of the preparation: 6 ML/s force of the two records of a preparation
fprintf('force at 6 ML/s (kPa), first record against second record of the same preparation: r = %.2f; at 5 ML/s r = %.2f\n', corr(F1(:, 8), F2(:, 8)), corr(F1(:, 7), F2(:, 7)));
% corrected 8 mM curve, first against second records, single factor and factor + offset
[~, sO, gO] = obj2(q2);  [~, s1] = obj1(R1);
[~, ~, g1f] = obj1(R1);                              % unit: the fitted fresh 8 mM isometric force
T.Fo = (T.F - q2(2)*s2)./q2(1).^s2./sO(T.prep)/gO(1);  T.F1 = T.F./R1.^s2./s1(T.prep)/g1f(1);
fprintf('8 mM, corrected, records run first | run second, v = %s\n', mat2str(vAll));
fprintf('   one factor (kPa fit):   %s | %s\n', mat2str(round(arrayfun(@(v) mean(T.F1(T.atp == 8 & T.slot == 1 & T.v == v)), vAll), 3)), mat2str(round(arrayfun(@(v) mean(T.F1(T.atp == 8 & T.slot == 2 & T.v == v)), vAll), 3)));
fprintf('   factor + offset:        %s | %s\n', mat2str(round(arrayfun(@(v) mean(T.Fo(T.atp == 8 & T.slot == 1 & T.v == v)), vAll), 3)), mat2str(round(arrayfun(@(v) mean(T.Fo(T.atp == 8 & T.slot == 2 & T.v == v)), vAll), 3)));
fprintf('ATP ratio low/8 by order, factor + offset correction\n');
for a = [4 2], for v = [2 3 4 5]
    q = arrayfun(@(p) T.Fo(T.prep == p & T.atp == a & T.v == v)/T.Fo(T.prep == p & T.atp == 8 & T.v == v), find(low == a));
    f8 = first(low == a) == 8;
    q1 = arrayfun(@(p) T.F1(T.prep == p & T.atp == a & T.v == v)/T.F1(T.prep == p & T.atp == 8 & T.v == v), find(low == a));
    fprintf('  %d/8 at %g ML/s: one factor 8-first %.2f, low-first %.2f | factor + offset 8-first %.2f, low-first %.2f, all %.3f +- %.3f\n', a, v, mean(q1(f8)), mean(q1(~f8)), mean(q(f8)), mean(q(~f8)), mean(q), std(q)/sqrt(numel(q)));
end, end
fprintf('corrected curves, factor + offset, mean +- sem (fresh 8 mM isometric = 1), v = %s\n', mat2str(vAll));
for a = [8 4 2]
    fprintf('  %d mM: %s\n     sem %s\n', a, mat2str(round(arrayfun(@(v) mean(T.Fo(T.atp == a & T.v == v)), vAll), 3)), mat2str(round(arrayfun(@(v) std(T.Fo(T.atp == a & T.v == v))/sqrt(sum(T.atp == a & T.v == v)), vAll), 3)));
end
% what the traces say: plateau of the ramp, baseline, speed, 8 mM first against second
fprintf('\n== Traces, 8 mM records: force during the ramp (kPa above the slack baseline), run first | run second\n');
tg = [2 3 4 5 6 7];  cyc = [1 4 6];  vNom = [6 1 3 5 0.5 4 2];  TR = struct();
tt = -3:0.05:22;
for c = cyc
    X = nan(numel(D), numel(tt));  B = nan(numel(D), 1);  V = B;
    for k = 1:numel(D)
        if D(k).atp ~= 8, continue; end
        R = readIsoV(D(k).t, D(k).L, D(k).F);
        X(k, :) = interp1(D(k).t - R.tStart(c), movmean(D(k).F, 5), tt) - R.F0(c);  B(k) = R.F0(c);  V(k) = R.vMeas(c);
    end
    i1 = [D.atp] == 8 & [D.slot] == 1;  i2 = [D.atp] == 8 & [D.slot] == 2;
    TR.(sprintf('c%d', c)) = struct('t', tt, 'm1', mean(X(i1, :), 1), 'm2', mean(X(i2, :), 1), 'e1', std(X(i1, :), 0, 1)/sqrt(sum(i1)), 'e2', std(X(i2, :), 0, 1)/sqrt(sum(i2)));
    tEnd = (1.10 - 1.00)*0 + interp1([6 5 4], [7.3 8.6 13.1], vNom(c));
    wv = tt >= 2 & tt <= tEnd;
    fprintf('%g ML/s: plateau from 2 ms to the crossing %.2f | %.2f kPa; first to second half of the plateau %+.2f | %+.2f kPa; slack baseline %.2f | %.2f kPa; ramp speed %.3f | %.3f ML/s\n', vNom(c), ...
        mean(TR.(sprintf('c%d', c)).m1(wv)), mean(TR.(sprintf('c%d', c)).m2(wv)), halfdiff(TR.(sprintf('c%d', c)).m1(wv)), halfdiff(TR.(sprintf('c%d', c)).m2(wv)), mean(B(i1)), mean(B(i2)), mean(V(i1)), mean(V(i2)));
end
%% C. what the points above 3 ML/s can carry: candidate fresh curves, and Vmax from a Hill fit
fprintf('\n== Candidate fresh curves (fresh 8 mM isometric = 1) at v = %s\n', mat2str(vAll));
cn = {'first records only', 'factor + offset', 'single factor', 'uncorrected (#75 scale)'};
CC = struct();  rng(1);
for a = [8 4 2]
    ma = T.atp == a;
    cv = {T.Fo(ma & T.slot == 1), T.Fo(ma), T.F1(ma), T.Fn(ma)};  vv = {T.v(ma & T.slot == 1), T.v(ma), T.v(ma), T.v(ma)};
    cur = nan(4, nV);
    for i = 1:4, cur(i, :) = arrayfun(@(v) mean(cv{i}(vv{i} == v)), vAll); end
    sem = arrayfun(@(v) std(T.Fo(ma & T.v == v))/sqrt(sum(ma & T.v == v)), vAll);
    sys = (max(cur(1:3, :), [], 1) - min(cur(1:3, :), [], 1))/2;          % half the range of the three corrected candidates
    fprintf('%d mM\n', a);
    for i = 1:4, fprintf('  %-26s %s\n', cn{i}, mat2str(round(cur(i, :), 3))); end
    fprintf('  %-26s %s\n  %-26s %s\n  %-26s %s\n  %-26s %s\n', 'sem (factor + offset)', mat2str(round(sem, 3)), 'half range of candidates', mat2str(round(sys, 3)), ...
        'half range / sem', mat2str(round(sys./sem, 1)), 'half range, % of the force', mat2str(round(100*sys./cur(2, :), 0)));
    % Hill fits: F = F0 (1 - v/vmax)/(1 + v/(a vmax)); weights 1/(0.01^2 + (0.05 F)^2)
    fprintf('  Hill fit, Vmax (ML/s) and a/F0:   0-6 ML/s          0-3 ML/s\n');
    hv = nan(4, 2);
    for i = 1:4
        h6 = hillfit(vAll, cur(i, :));  h3 = hillfit(vAll(1:5), cur(i, 1:5));  hv(i, :) = [h6(3) h3(3)];
        fprintf('  %-26s        %5.2f (%.2f)      %5.2f (%.2f)\n', cn{i}, h6(3), h6(2), h3(3), h3(2));
    end
    % bootstrap over records, factor + offset curve
    recs = unique(T.rec(ma));  nB = 500;  vb = nan(nB, 2);
    Fm = nan(numel(recs), nV);
    for r = 1:numel(recs), for j = 1:nV, Fm(r, j) = T.Fo(T.rec == recs(r) & T.v == vAll(j)); end, end
    for b = 1:nB
        c = mean(Fm(randi(numel(recs), numel(recs), 1), :), 1);
        h6 = hillfit(vAll, c);  h3 = hillfit(vAll(1:5), c(1:5));  vb(b, :) = [h6(3) h3(3)];
    end
    q = prctile(vb, [2.5 50 97.5]);
    fprintf('  bootstrap over the %d records (factor + offset): Vmax 0-6 ML/s %.2f [%.2f, %.2f]; 0-3 ML/s only %.2f [%.2f, %.2f]\n', numel(recs), q(2, 1), q(1, 1), q(3, 1), q(2, 2), q(1, 2), q(3, 2));
    fprintf('  Vmax across the three corrected candidates, 0-6 ML/s: %.2f-%.2f ML/s\n', min(hv(1:3, 1)), max(hv(1:3, 1)));
    CC.(sprintf('a%d', a)) = struct('cur', cur, 'sem', sem, 'sys', sys, 'hv', hv, 'boot', q);
end
% how large the force is against what limits it
bsd = nan(numel(D), 1);
for k = 1:numel(D), R = readIsoV(D(k).t, D(k).L, D(k).F);  bsd(k) = std(R.F0(1:7)); end
fprintf('slack baseline, sd over the seven cycles of a record: mean %.2f kPa (range %.2f-%.2f)\n', mean(bsd), min(bsd), max(bsd));
for a = [8 4 2]
    fprintf('%d mM, mean force (kPa) at 3, 4, 5, 6 ML/s: %s ; between-preparation sd %s ; records with force below 2 kPa: %s of %d\n', a, ...
        mat2str(round(arrayfun(@(v) mean(T.F(T.atp == a & T.v == v)), [3 4 5 6]), 1)), mat2str(round(arrayfun(@(v) std(T.F(T.atp == a & T.v == v)), [3 4 5 6]), 1)), ...
        mat2str(arrayfun(@(v) sum(T.F(T.atp == a & T.v == v) < 2), [3 4 5 6])), sum(T.atp == a & T.v == 0));
end
diary off;
save(fullfile(ana, 'results', 'fv_rundown_perprep.mat'), 'CC', 'rv', 'rLow', 'rHi', 'rAll', 'rH', 'seAll', 'A', 'dAll', 'seAllv', 'exAll', 'fitm', 'q2', 'qL', 'g1', 'g2', 'pl', 'TR', 'T');

%% figure
col = containers.Map({8, 4, 2}, {[0 0 0], [0.85 0.33 0.10], [0 0.45 0.74]});
fig = figure('Visible', 'off', 'Position', [50 50 1500 820], 'Color', 'w');
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on;
[~, o] = sortrows([low == 2, first ~= 8]);
for i = 1:nP, p = o(i);  mk = 'o'; if first(p) ~= 8, mk = 's'; end
    errorbar(i, exp(rAll(p)), seAll(p)*exp(rAll(p)), mk, 'Color', col(low(p)), 'MarkerFaceColor', col(low(p)), 'MarkerSize', 7);
    plot(i + 0.22, exp(rLow(p)), mk, 'Color', col(low(p)), 'MarkerSize', 5);
    plot(i - 0.22, exp(rH(p)), 'x', 'Color', [0.5 0.5 0.5], 'MarkerSize', 6);
end
yline(exp(S0.rho), 'k--'); yline(1, ':'); set(gca, 'XTick', 1:nP, 'XTickLabel', arrayfun(@(p) sprintf('%d', p), o, 'UniformOutput', false)); xlim([0.3 nP + 0.7]);
xlabel('preparation (circle: 8 mM first, square: low ATP first)'); ylabel('second / first record, ATP effect removed');
title('Rundown per preparation (filled 0-3 ML/s, open 0.5-1 ML/s, x hold)'); grid on;
nexttile; hold on;
for p = 1:nP, mk = 'o'; if first(p) ~= 8, mk = 's'; end
    plot(exp(rLow(p)), exp(rHi(p)), mk, 'Color', col(low(p)), 'MarkerFaceColor', col(low(p)), 'MarkerSize', 7); text(exp(rLow(p)) + 0.008, exp(rHi(p)), sprintf('%d', p), 'FontSize', 8);
end
plot([0.7 1.3], [0.7 1.3], 'k:'); xlabel('rundown read at 0.5-1 ML/s'); ylabel('rundown read at 2-3 ML/s'); title(sprintf('Same preparation, two speed bands (r = %.2f)', c1)); grid on; axis tight;
nexttile; hold on;
errorbar(vAll, dAll, seAllv, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1.2, 'DisplayName', 'observed, order-balanced');
plot(vAll, exAll, 'r--', 'LineWidth', 1.2, 'DisplayName', sprintf('one factor %.2f', exp(S0.rho)));
plot(vAll, (q2(1) - 1)*F1m + q2(2), 'b-', 'LineWidth', 1.2, 'DisplayName', sprintf('factor %.2f + %.1f kPa', q2(1), q2(2)));
yline(0, ':', 'HandleVisibility', 'off'); xlabel('speed (ML/s)'); ylabel('second - first record (kPa)'); title('Rundown in kPa, within preparations'); legend('Box', 'off', 'Location', 'southeast'); grid on;
nexttile; hold on;
plot(g1, g2, 'ko', 'MarkerFaceColor', 'k'); xx = logspace(log10(1.5), log10(80), 60); plot(xx, polyval(pl, xx), 'b-'); plot(xx, exp(S0.rho)*xx, 'r--'); plot(xx, xx, 'k:');
for j = 1:nV, text(g1(j) + 1, g2(j) - 1.2, sprintf('%g', vAll(j)), 'FontSize', 8); end
xlabel('8 mM run first, group mean (kPa)'); ylabel('8 mM run second, group mean (kPa)'); title(sprintf('8 mM: second = %.2f first + %.1f kPa (blue); factor 0.88 (red)', pl(1), pl(2))); grid on;
set(gca, 'XScale', 'log', 'YScale', 'log'); xlim([1.5 80]); ylim([1.5 80]);
for c = [1 6]
    nexttile; hold on;  r = TR.(sprintf('c%d', c));
    fill([r.t fliplr(r.t)], [r.m1 - r.e1, fliplr(r.m1 + r.e1)], [0 0 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    fill([r.t fliplr(r.t)], [r.m2 - r.e2, fliplr(r.m2 + r.e2)], [0.8 0 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(r.t, r.m1, 'k-', 'LineWidth', 1.2); plot(r.t, r.m2, '-', 'Color', [0.8 0 0], 'LineWidth', 1.2);
    ylim([-4 14]); xlim([-1 20]); yline(0, ':'); xlabel('time from the step (ms)'); ylabel('force above slack baseline (kPa)');
    title(sprintf('%g ML/s ramp, 8 mM: run first black, second red', vNom(c))); grid on;
end
title(tl, 'Beard2022 isovelocity records: rundown per preparation, and the high-speed offset (black 8 mM, orange 4 mM pairs, blue 2 mM pairs)');
exportgraphics(fig, fullfile(ana, 'figures', 'fv_rundown_perprep.png'), 'Resolution', 150);
close(fig);

%% figure for part C: the tail above 3 ML/s
fig = figure('Visible', 'off', 'Position', [50 50 1500 430], 'Color', 'w');
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
hf = @(h, v) h(1)*(1 - v/h(3))./(1 + v/(h(2)*h(3)));  vf = linspace(0, 8, 200);
mkc = {'v', 'o', 's'};  cc = {[0 0.6 0.3], [0 0 0], [0.8 0 0]};
for a = [8 4 2]
    nexttile; hold on;  c = CC.(sprintf('a%d', a));
    h6 = hillfit(vAll, c.cur(2, :));  h3 = hillfit(vAll(1:5), c.cur(2, 1:5));
    plot(vf, max(hf(h3, vf), 1e-4), 'b--', 'LineWidth', 1.2, 'DisplayName', sprintf('Hill on 0-3 ML/s: V_{max} %.1f', h3(3)));
    plot(vf, max(hf(h6, vf), 1e-4), 'b-', 'LineWidth', 1.2, 'DisplayName', sprintf('Hill on 0-6 ML/s: V_{max} %.1f', h6(3)));
    for i = [1 3], plot(vAll + 0.06*(i - 2), c.cur(i, :), mkc{i}, 'Color', cc{i}, 'MarkerSize', 6, 'DisplayName', cn{i}); end
    errorbar(vAll, c.cur(2, :), sqrt(c.sem.^2 + c.sys.^2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 1, 'DisplayName', 'factor + offset, sem and candidate range');
    xline(3, ':', 'HandleVisibility', 'off'); set(gca, 'YScale', 'log'); ylim([0.01 1.3]); xlim([0 8]); grid on;
    xlabel('speed (ML/s)'); ylabel('force (fresh 8 mM isometric = 1)'); title(sprintf('%d mM', a)); legend('Box', 'off', 'Location', 'southwest', 'FontSize', 8);
end
title(tl, 'Above 3 ML/s the estimates of the fresh curve disagree, and V_{max} is set by those points');
exportgraphics(fig, fullfile(ana, 'figures', 'fv_rundown_highspeed.png'), 'Resolution', 150);
close(fig);

function [chi2, s, g] = bilinear(y, w, pI, gI, nP, nG)
% weighted least squares of y = s(pI) * g(gI) by alternating updates
s = ones(nP, 1);  g = accumarray(gI, w.*y, [nG 1])./accumarray(gI, w, [nG 1]);
for it = 1:200
    s = accumarray(pI, w.*y.*g(gI), [nP 1])./accumarray(pI, w.*g(gI).^2, [nP 1]);
    gn = accumarray(gI, w.*y.*s(pI), [nG 1])./accumarray(gI, w.*s(pI).^2, [nG 1]);
    if max(abs(gn - g)) < 1e-9, g = gn; break; end
    g = gn;
end
k = mean(s);  s = s/k;  g = g*k;
chi2 = sum(w.*(y - s(pI).*g(gI)).^2);
end

function d = halfdiff(x)
h = floor(numel(x)/2);  d = mean(x(h+1:end)) - mean(x(1:h));
end

function h = hillfit(v, F)
% Hill hyperbola F = F0 (1 - v/vmax)/(1 + v/(a vmax)); returns [F0, a/F0, vmax]
w = 1./(0.01^2 + (0.05*F).^2);
f = @(q, v) q(1)*(1 - v/exp(q(3)))./(1 + v/(exp(q(2))*exp(q(3))));
o = @(q) sum(w.*(F - f(q, v)).^2);
best = inf;
for v0 = [5 8 12]
    q = fminsearch(o, [F(1) log(0.3) log(v0)], optimset('MaxFunEvals', 4000, 'MaxIter', 4000, 'TolX', 1e-6, 'TolFun', 1e-10));
    if o(q) < best, best = o(q); h = [q(1) exp(q(2)) exp(q(3))]; end
end
end
