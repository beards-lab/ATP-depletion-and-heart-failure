%% RunFVRundown  Separate rundown from velocity in the Beard2022 isovelocity records (#92).
%
% Data only. 24 records = 12 preparations x 2 records (one at 8 mM, one at 4 or
% 2 mM; 3 preparations per running order). Per record: 7 isovelocity points
% (recorded order 6, 1, 3, 5, 0.5, 4, 2 ML/s), the isometric point at 1.00 Lo,
% and the 1.10 Lo hold force before each of the 7 cycles. All forces are
% baseline-subtracted (readIsoV F0).
%
% Factor model, log space, speeds 0-3 ML/s (forces well above the noise):
%   log F(prep p, record k, speed v) = f0(p) + g(v, ATP) + r + e
% with rundown forms
%   M0   r = 0
%   M1   r = rho * [second record]                       one factor
%   M1a  r = rho(pair) * [second record]                 one per ATP pair (8/4, 8/2)
%   M1p  r = rho(p) * [second record]                    one per preparation
%   M2   r = rho(v) * [second record]                    speed-dependent
%   M3   M1  + beta * log(H_c / H_1)                     within-record part from the holds
%   M3f  M1  + 1    * log(H_c / H_1)                     hold-normalised (beta fixed at 1)
%   M4   M1p + beta * log(H_c / H_1)
% Output: results/fv_rundown.log, results/fv_rundown.mat, figures/*.png.

here = fileparts(mfilename('fullpath'));  ana = fileparts(here);  root = fileparts(fileparts(ana));
addpath(fullfile(root, 'Analyses', 'FV_Beard2022'));
D = loadBeard2022FV();
logf = fullfile(ana, 'results', 'fv_rundown.log');  if isfile(logf), delete(logf); end
diary(logf);

%% 1. read every record
vNom = [6 1 3 5 0.5 4 2 0];  nR = numel(D);
preps = unique({D.prep}, 'stable');  nP = numel(preps);
T = table();
for k = 1:nR
    t = D(k).t;  F = D(k).F;
    R = readIsoV(t, D(k).L, F);
    assert(height(R) == 8, 'record %d: %d rows', k, height(R));
    H = nan(8, 1);
    for c = 1:7
        w = t >= R.tStart(c) - 20 & t < R.tStart(c) - 1;
        H(c) = mean(F(w)) - R.F0(c);
    end
    H(8) = H(7);
    tc = [R.tStart(1:7); t(end)]/1000;
    p = find(strcmp(preps, D(k).prep));
    T = [T; table(repmat(k, 8, 1), repmat(p, 8, 1), repmat(D(k).atp, 8, 1), repmat(D(k).slot, 8, 1), (1:8)', vNom', ...
        R.Fx - R.F0, H, H/H(1), tc, 'VariableNames', {'rec', 'prep', 'atp', 'slot', 'cyc', 'v', 'F', 'H', 'Hrel', 't'})]; %#ok<AGROW>
end
low = arrayfun(@(p) max(T.atp(T.prep == p & T.atp < 8)), (1:nP)');          % the low ATP of each preparation
first = arrayfun(@(p) T.atp(find(T.prep == p & T.slot == 1, 1)), (1:nP)');  % ATP of the first record
T.pair = low(T.prep);
% #75's normalisation: one scale per preparation, its 8 mM isometric force
Fref = arrayfun(@(p) T.F(T.prep == p & T.atp == 8 & T.v == 0), (1:nP)');
T.Fn = T.F./Fref(T.prep);
fprintf('Records %d, preparations %d. Orders: %s\n', nR, nP, strjoin(arrayfun(@(p) sprintf('%d->%d', first(p), setdiff([8 low(p)], first(p))), 1:nP, 'UniformOutput', false), ' '));
fprintf('Lowest baseline-subtracted force by speed (kPa), all records: %s (v = %s)\n', mat2str(round(arrayfun(@(v) min(T.F(T.v == v)), sort(vNom)), 2)), mat2str(sort(vNom)));

%% 2. within-record drift of the 1.10 Lo hold
fprintf('\n== Within-record drift: hold before cycle 7 over hold before cycle 1 (%%), by [ATP] and by record position\n');
dr = arrayfun(@(k) 100*(T.Hrel(T.rec == k & T.cyc == 7) - 1), (1:nR)');
ra = arrayfun(@(k) T.atp(find(T.rec == k, 1)), (1:nR)');  rs = arrayfun(@(k) T.slot(find(T.rec == k, 1)), (1:nR)');
rp = arrayfun(@(k) T.prep(find(T.rec == k, 1)), (1:nR)');
for a = [8 4 2]
    fprintf('%d mM: all %+5.1f +- %.1f (n=%d) | run first %+5.1f +- %.1f (n=%d) | run second %+5.1f +- %.1f (n=%d)\n', a, mean(dr(ra == a)), std(dr(ra == a)), sum(ra == a), ...
        mean(dr(ra == a & rs == 1)), std(dr(ra == a & rs == 1)), sum(ra == a & rs == 1), mean(dr(ra == a & rs == 2)), std(dr(ra == a & rs == 2)), sum(ra == a & rs == 2));
end
fprintf('by position: first records %+5.1f +- %.1f, second records %+5.1f +- %.1f\n', mean(dr(rs == 1)), std(dr(rs == 1)), mean(dr(rs == 2)), std(dr(rs == 2)));
% hold of the second record over the hold of the first record (first hold of each)
H1 = arrayfun(@(k) T.H(T.rec == k & T.cyc == 1), (1:nR)');
hr = arrayfun(@(p) H1(rp == p & rs == 2)/H1(rp == p & rs == 1), (1:nP)');
fprintf('first hold of the second record over first hold of the first record, by order:\n');
ord = {[8 4], [4 8], [8 2], [2 8]};
for i = 1:4
    m = first == ord{i}(1) & low == min(ord{i});
    fprintf('  %d -> %d: %s  mean %.3f\n', ord{i}(1), ord{i}(2), mat2str(round(hr(m)', 3)), mean(hr(m)));
end

%% 3. factor model
use = T.v <= 3;  S = T(use, :);  n = height(S);  y = log(S.F);
vs = unique(S.v);  as = [8 4 2];
Xp = double(S.prep == 1:nP);                                         % f0(p)
G = zeros(n, numel(vs)*3);  gi = 0;  gname = {};
for a = as, for v = vs', gi = gi + 1; G(:, gi) = S.atp == a & S.v == v; gname{gi} = sprintf('%g@%d', v, a); end, end %#ok<SAGROW>
G(:, 1) = [];  gname(1) = [];                                         % reference: v = 0 at 8 mM
s2 = double(S.slot == 2);  lh = log(S.Hrel);
M = struct('name', {}, 'X', {}, 'off', {}, 'k', {}, 'b', {}, 'rss', {}, 'sd', {}, 'aic', {}, 'bic', {}, 'cov', {});
add = @(M, name, X, off) [M, fitOne(name, [Xp G X], off, y)];
M = add(M, 'M0', zeros(n, 0), 0);
M = add(M, 'M1', s2, 0);
M = add(M, 'M1a', [s2.*(S.pair == 4), s2.*(S.pair == 2)], 0);
M = add(M, 'M1p', s2.*Xp, 0);
M = add(M, 'M2', s2.*double(S.v == vs'), 0);
M = add(M, 'M3', [s2, lh], 0);
M = add(M, 'M3f', s2, lh);
M = add(M, 'M4', [s2.*Xp, lh], 0);
fprintf('\n== Factor model on %d points (v <= 3 ML/s). k = parameters, sd = residual sd of log force\n', n);
fprintf('%-5s %3s %8s %8s %8s %8s\n', 'model', 'k', 'sd', 'AIC', 'BIC', 'dAIC');
for i = 1:numel(M), fprintf('%-5s %3d %8.4f %8.1f %8.1f %8.1f\n', M(i).name, M(i).k, M(i).sd, M(i).aic, M(i).bic, M(i).aic - min([M.aic])); end
nb = nP + size(G, 2);
m1 = M(2);   fprintf('M1  rho = %.3f +- %.3f  -> second record = %.3f of the first\n', m1.b(nb+1), sqrt(m1.cov(nb+1, nb+1)), exp(m1.b(nb+1)));
m1a = M(3);  fprintf('M1a rho(8/4 pairs) = %.3f +- %.3f (x%.3f), rho(8/2 pairs) = %.3f +- %.3f (x%.3f)\n', m1a.b(nb+1), sqrt(m1a.cov(nb+1, nb+1)), exp(m1a.b(nb+1)), m1a.b(nb+2), sqrt(m1a.cov(nb+2, nb+2)), exp(m1a.b(nb+2)));
m2 = M(5);   fprintf('M2  second/first by speed %s: %s (se %s)\n', mat2str(vs'), mat2str(round(exp(m2.b(nb+1:nb+numel(vs)))', 3)), mat2str(round(sqrt(diag(m2.cov(nb+1:nb+numel(vs), nb+1:nb+numel(vs))))', 3)));
m3 = M(6);   fprintf('M3  rho = %.3f +- %.3f, beta (hold exponent) = %.2f +- %.2f\n', m3.b(nb+1), sqrt(m3.cov(nb+1, nb+1)), m3.b(nb+2), sqrt(m3.cov(nb+2, nb+2)));
fprintf('M1p and M4 are rank deficient by 2: with a factor per preparation the mean rundown cannot be told from the mean ATP effect of each pair.\n');

%% 4. what the design separates: the raw second-minus-first difference by order
fprintf('\n== Raw log ratio second/first record of a preparation, mean by order, at v = %s\n', mat2str(vs'));
dl = nan(nP, numel(vs));
for p = 1:nP, for j = 1:numel(vs), dl(p, j) = log(S.F(S.prep == p & S.slot == 2 & S.v == vs(j))/S.F(S.prep == p & S.slot == 1 & S.v == vs(j))); end, end
for i = 1:4, m = first == ord{i}(1) & low == min(ord{i}); fprintf('  %d -> %d: %s\n', ord{i}(1), ord{i}(2), mat2str(round(mean(dl(m, :), 1), 3))); end
rhoP = mean(dl, 2);                     % per preparation, mean over the five speeds: rundown x ATP effect
fprintf('per preparation, second/first (mean over the speeds), by order:\n');
for i = 1:4, m = first == ord{i}(1) & low == min(ord{i}); fprintf('  %d -> %d: %s  sd of log %.3f\n', ord{i}(1), ord{i}(2), mat2str(round(exp(rhoP(m))', 3)), std(rhoP(m))); end
for a = [4 2]
    A = mean(dl(first == 8 & low == a, :), 1);  B = mean(dl(first == a & low == a, :), 1);
    fprintf('  %d mM pairs: rundown = half the sum %s ; ATP effect log(F_%d/F_8) = half the difference %s -> ratio %s\n', a, mat2str(round((A + B)/2, 3)), a, mat2str(round((A - B)/2, 3)), mat2str(round(exp((A - B)/2), 3)));
end
% speed dependence of rundown from the 8 mM records alone: first (6) against second (6), relative to the isometric point
r8 = nan(nP, numel(vs));
for p = 1:nP, for j = 1:numel(vs), r8(p, j) = log(S.F(S.prep == p & S.atp == 8 & S.v == vs(j))/S.F(S.prep == p & S.atp == 8 & S.v == 0)); end, end
d8 = mean(r8(first ~= 8, :), 1) - mean(r8(first == 8, :), 1);
se8 = sqrt(var(r8(first ~= 8, :), 0, 1)/sum(first ~= 8) + var(r8(first == 8, :), 0, 1)/sum(first == 8));
fprintf('8 mM records, F(v)/F(0): run second relative to run first, v = %s: %s (se %s)\n', mat2str(vs'), mat2str(round(exp(d8), 3)), mat2str(round(se8, 3)));

%% 5. corrected data and curves (correction of M1: divide second records by exp(rho); scale by exp(f0))
rho = m1.b(nb+1);  f0 = m1.b(1:nP);
T.Fc = T.F./exp(f0(T.prep))./exp(rho*(T.slot == 2));                 % fresh-state units: 1 = fresh 8 mM isometric
vAll = sort(vNom);
fprintf('\n== F-V per [ATP]: uncorrected (per-prep 8 mM isometric scale, #75) and corrected (M1), mean +- sem over records, and between-preparation CV (%%)\n');
C = struct();
for a = as
    fprintf('%d mM\n   v      uncorrected        CV   | corrected          CV   | P-V uncorr  corr\n', a);
    for j = 1:numel(vAll)
        m = T.atp == a & T.v == vAll(j);  x = T.Fn(m);  z = T.Fc(m);
        C.(sprintf('a%d', a))(j, :) = [vAll(j) mean(x) std(x)/sqrt(numel(x)) mean(z) std(z)/sqrt(numel(z)) std(x)/mean(x) std(z)/mean(z)];
        fprintf(' %4.1f   %.3f +- %.3f   %5.1f  | %.3f +- %.3f   %5.1f  | %.3f      %.3f\n', vAll(j), mean(x), std(x)/sqrt(numel(x)), 100*std(x)/mean(x), mean(z), std(z)/sqrt(numel(z)), 100*std(z)/mean(z), vAll(j)*mean(x), vAll(j)*mean(z));
    end
end
fprintf('8 mM corrected: records run first against run second (mean), v = %s\n   first  %s\n   second %s\n', mat2str(vAll), ...
    mat2str(round(arrayfun(@(v) mean(T.Fc(T.atp == 8 & T.slot == 1 & T.v == v)), vAll), 3)), mat2str(round(arrayfun(@(v) mean(T.Fc(T.atp == 8 & T.slot == 2 & T.v == v)), vAll), 3)));
fprintf('8 mM uncorrected (#75 scale): first %s\n                             second %s\n', ...
    mat2str(round(arrayfun(@(v) mean(T.Fn(T.atp == 8 & T.slot == 1 & T.v == v)), vAll), 3)), mat2str(round(arrayfun(@(v) mean(T.Fn(T.atp == 8 & T.slot == 2 & T.v == v)), vAll), 3)));

%% 6. order effect of the ATP ratios
fprintf('\n== ATP ratio low/8 within a preparation, by order, before and after the correction\n');
for a = [4 2]
    for v = [0 0.5 1 2 3 4]
        q = nan(nP, 2);
        for p = find(low == a)'
            q(p, 1) = T.F(T.prep == p & T.atp == a & T.v == v)/T.F(T.prep == p & T.atp == 8 & T.v == v);
            q(p, 2) = T.Fc(T.prep == p & T.atp == a & T.v == v)/T.Fc(T.prep == p & T.atp == 8 & T.v == v);
        end
        i8 = low == a & first == 8;  iL = low == a & first == a;
        fprintf('%d/8 at %3.1f ML/s: raw 8-first %.3f, low-first %.3f (gap %+.3f) | corrected %.3f, %.3f (gap %+.3f) | corrected all %.3f +- %.3f\n', a, v, ...
            mean(q(i8, 1)), mean(q(iL, 1)), mean(q(iL, 1)) - mean(q(i8, 1)), mean(q(i8, 2)), mean(q(iL, 2)), mean(q(iL, 2)) - mean(q(i8, 2)), mean(q(i8 | iL, 2)), std(q(i8 | iL, 2))/sqrt(sum(i8 | iL)));
    end
end

%% 7. checks: leave one preparation out; residuals against order
gM = exp(m1.b(nP+1:nb));  dmax = 0;  rhoL = nan(nP, 1);
for p = 1:nP
    k = S.prep ~= p;  Xk = [Xp(k, setdiff(1:nP, p)) G(k, :) s2(k)];
    b = Xk\y(k);  rhoL(p) = b(end);
    dmax = max(dmax, max(abs(exp(b(nP:end-1))./gM - 1)));
end
fprintf('\n== Leave one preparation out (M1): rho %.3f..%.3f (all: %.3f); largest change of any g(v, ATP): %.1f %%\n', min(rhoL), max(rhoL), rho, 100*dmax);
res0 = y - [Xp G]*M(1).b;  res1 = y - [Xp G s2]*m1.b;
fprintf('residual mean by record position: M0 first %+.4f second %+.4f | M1 first %+.4f second %+.4f\n', mean(res0(s2 == 0)), mean(res0(s2 == 1)), mean(res1(s2 == 0)), mean(res1(s2 == 1)));
fprintf('M1 residual sd by speed %s: %s\n', mat2str(vs'), mat2str(round(arrayfun(@(v) std(res1(S.v == v)), vs)', 3)));
fprintf('between-preparation sd of f0 (log fresh force scale): %.3f -> CV %.0f %%; mean fresh 8 mM isometric force %.1f kPa\n', std(f0), 100*std(f0), mean(exp(f0)));
diary off;
save(fullfile(ana, 'results', 'fv_rundown.mat'), 'T', 'M', 'C', 'rho', 'f0', 'rhoP', 'first', 'low', 'vs', 'dr', 'hr', 'dl');

%% figure
col = containers.Map({8, 4, 2}, {[0 0 0], [0.85 0.33 0.10], [0 0.45 0.74]});
fig = figure('Visible', 'off', 'Position', [50 50 1500 820], 'Color', 'w');
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on;
for a = as, c = C.(sprintf('a%d', a));
    errorbar(c(:, 1) - 0.04, c(:, 2), c(:, 3), 'o--', 'Color', [col(a) 0.5], 'MarkerSize', 4, 'DisplayName', sprintf('%d mM uncorrected', a));
    errorbar(c(:, 1) + 0.04, c(:, 4), c(:, 5), 's-', 'Color', col(a), 'LineWidth', 1.3, 'MarkerFaceColor', col(a), 'MarkerSize', 4, 'DisplayName', sprintf('%d mM corrected', a));
end
xlabel('speed (ML/s)'); ylabel('force (fresh 8 mM isometric = 1)'); title('F-V, mean \pm sem'); legend('Box', 'off'); grid on;
nexttile; hold on;
for a = as, c = C.(sprintf('a%d', a));
    plot(c(:, 1), c(:, 1).*c(:, 2), 'o--', 'Color', [col(a) 0.5], 'MarkerSize', 4);
    plot(c(:, 1), c(:, 1).*c(:, 4), 's-', 'Color', col(a), 'LineWidth', 1.3, 'MarkerFaceColor', col(a), 'MarkerSize', 4);
end
xlabel('speed (ML/s)'); ylabel('power (F_{ref} ML/s)'); title('P-V (dashed uncorrected, solid corrected)'); grid on;
nexttile; hold on;
for a = as, c = C.(sprintf('a%d', a)); k = c(:, 1) <= 4;
    plot(c(k, 1), 100*c(k, 6), 'o--', 'Color', [col(a) 0.5]); plot(c(k, 1), 100*c(k, 7), 's-', 'Color', col(a), 'LineWidth', 1.3, 'MarkerFaceColor', col(a));
end
xlabel('speed (ML/s)'); ylabel('between-preparation CV (%)'); title('Spread before (dashed) and after (solid)'); grid on;
nexttile; hold on;
mk = {'o', 's', 'd', '^'};  lab = {'8 \rightarrow 4', '4 \rightarrow 8', '8 \rightarrow 2', '2 \rightarrow 8'};
for i = 1:4, m = find(first == ord{i}(1) & low == min(ord{i}));
    plot(i + 0.12*(-1:1), exp(rhoP(m)), mk{i}, 'Color', col(min(ord{i})), 'MarkerFaceColor', col(min(ord{i})), 'MarkerSize', 7);
end
yline(exp(rho), 'k--'); yline(1, ':'); set(gca, 'XTick', 1:4, 'XTickLabel', lab); xlim([0.5 4.5]);
ylabel('second / first record, mean over 0-3 ML/s'); title('Per preparation: rundown \times ATP effect (dashed: fitted rundown)'); grid on;
nexttile; hold on;
for i = 1:4, m = first == ord{i}(1) & low == min(ord{i});
    ls = '-'; if ord{i}(1) ~= 8, ls = '--'; end
    plot(vs, exp(mean(dl(m, :), 1)), [mk{i} ls], 'Color', col(min(ord{i})), 'LineWidth', 1.2, 'DisplayName', lab{i});
end
yline(1, ':', 'HandleVisibility', 'off'); xlabel('speed (ML/s)'); ylabel('second record / first record, raw'); legend('Box', 'off', 'Location', 'best'); grid on;
title('What the balanced order separates');
nexttile; hold on;
for k = 1:nR
    m = T.rec == k & T.cyc <= 7;  ls = '-'; if rs(k) == 2, ls = '--'; end
    plot(T.t(m), T.Hrel(m), ['.' ls], 'Color', [col(ra(k)) 0.7]);
end
xlabel('time in the record (s)'); ylabel('1.10 L_o hold / first hold'); title('Within-record drift (solid: run first)'); grid on;
title(tl, 'Beard2022 isovelocity records: rundown separated from velocity (black 8 mM, orange 4 mM, blue 2 mM)');
exportgraphics(fig, fullfile(ana, 'figures', 'fv_rundown.png'), 'Resolution', 150);
close(fig);

function m = fitOne(name, X, off, y)
keep = any(X ~= 0, 1);  X = X(:, keep);
yy = y - off;  b = X\yy;  r = yy - X*b;  n = numel(y);  k = rank(X);
m.name = name;  m.X = [];  m.off = [];  m.k = k;  m.b = b;  m.rss = r'*r;  m.sd = sqrt(m.rss/(n - k));
m.aic = n*log(m.rss/n) + 2*k;  m.bic = n*log(m.rss/n) + k*log(n);
m.cov = m.sd^2*pinv(X'*X);
end
