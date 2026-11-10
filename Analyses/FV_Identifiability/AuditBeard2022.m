% AuditBeard2022.m -- is the Beard/Baker 2022 isovelocity dataset sound enough to
% carry a paper? Data only (no model run); seconds.
%
% Prints, and saves to results/audit.mat:
%   1. design: preparations, running order, n per [ATP]
%   2. per record: reference force, noise, drift of the 1.10 Lo hold force across
%      the seven cycles, whether the hold has settled before the next cycle
%   3. rundown between the two records of a preparation, from the balanced design
%   4. the ATP effect at each speed bracketed by running order
%   5. between-preparation spread against the size of the ATP effect
%   6. outlying preparations
% Figure: figures/audit_rundown.png
%
%   matlab -batch "run('Analyses/FV_Identifiability/AuditBeard2022.m')"

here = fileparts(mfilename('fullpath')); addpath(here);
addpath(fullfile(fileparts(here), 'FV_Beard2022'));
figP = fullfile(here, 'figures'); if ~isfolder(figP), mkdir(figP); end
D = loadBeard2022FV(); nS = numel(D);
ATP = [8 4 2];

%% 1. design
preps = unique({D.prep}, 'stable');
fprintf('\n1. DESIGN: %d preparations, %d records\n', numel(preps), nS);
ord = cellfun(@(p) D(find(strcmp({D.prep}, p), 1)).order, preps, 'uni', 0);
[uo, ~, io] = unique(ord); for i = 1:numel(uo), fprintf('   order %s: %d preparations\n', uo{i}, sum(io == i)); end
for a = ATP, fprintf('   %d mM: n = %d (first %d, second %d)\n', a, sum([D.atp] == a), sum([D.atp] == a & [D.slot] == 1), sum([D.atp] == a & [D.slot] == 2)); end

%% 2. per record
Q = table();
for k = 1:nS
    t = D(k).t; F = D(k).F; L = D(k).L;
    R = readIsoV(t, L, F); z = median(R.F0(1:7));
    h = nan(1, 7); unsettled = nan(1, 7); slackMs = nan(1, 7);
    for i = 1:7
        h(i) = mean(F(t >= R.tStart(i) - 60 & t <= R.tStart(i) - 5)) - z;
        late = mean(F(t >= R.tStart(i) - 20 & t <= R.tStart(i) - 5));
        early = mean(F(t >= R.tStart(i) - 60 & t <= R.tStart(i) - 45));
        unsettled(i) = (late - early) / h(i);
        seg = t > R.tSlack(i) + 2 & t < R.tSlack(i) + 200;
        j = find(seg & L > min(L(seg)) + 0.01, 1);
        if ~isempty(j), slackMs(i) = t(j) - R.tSlack(i); end
    end
    w = t < 400; noise = std(F(w) - movmean(F(w), 201));
    Fiso = R.Fx(end) - R.F0(end);
    Q = [Q; table(string(D(k).date), D(k).atp, D(k).slot, Fiso, h(1), h(7) / h(1) - 1, max(abs(unsettled(2:7))), ...
        noise, median(slackMs), z, 'VariableNames', ...
        {'prep', 'ATP', 'slot', 'Fiso100', 'F110', 'holdDrift', 'unsettled', 'noise', 'slackMs', 'zero'})]; %#ok<AGROW>
    D(k).R = R; D(k).h = h;
end
fprintf('\n2. PER RECORD (kPa; holdDrift = change of the 1.10 Lo hold force from cycle 1 to cycle 7, 3.1 s;\n');
fprintf('   unsettled = largest change of the hold force over its last 40 ms before a cycle, fraction of the hold)\n');
disp(Q);
for a = ATP
    q = Q(Q.ATP == a, :);
    fprintf('   %d mM: hold drift %+.1f %% +- %.1f (range %+.1f to %+.1f); unsettled <= %.2f %%; noise %.2f kPa; F110/Fiso100 %.2f +- %.2f\n', ...
        a, 100 * mean(q.holdDrift), 100 * std(q.holdDrift), 100 * min(q.holdDrift), 100 * max(q.holdDrift), ...
        100 * max(q.unsettled), mean(q.noise), mean(q.F110 ./ q.Fiso100), std(q.F110 ./ q.Fiso100));
end
for s = 1:2
    q = Q(Q.slot == s, :);
    fprintf('   slot %d (all [ATP]): hold drift %+.1f %% +- %.1f\n', s, 100 * mean(q.holdDrift), 100 * std(q.holdDrift));
end

%% 3 + 4. rundown and the ATP effect, bracketed by running order
% ratio r = F([ATP]) / F(8 mM) within a preparation. If the second record of a
% preparation is depressed by a factor rho (rundown) and the ATP effect is A:
%   8 mM first:  r = A * rho      low [ATP] first:  r = A / rho
% so A = sqrt(r_8first * r_lowfirst) and rho = sqrt(r_8first / r_lowfirst) if
% rho is the same in both orders. If it is not, A still lies between the two
% group means, because rundown only depresses.
v0 = D(1).summary.v; [vS, iv] = sort(v0); nV = numel(vS);
Fz = @(k) D(k).R.Fx(iv)' - D(k).R.F0(iv)';
B = table(); G = struct();
for a = [4 2]
    pr = unique({D([D.atp] == a).prep}); r = nan(numel(pr), nV + 1); first8 = false(numel(pr), 1);
    for ip = 1:numel(pr)
        k8 = find(strcmp({D.prep}, pr{ip}) & [D.atp] == 8); ka = find(strcmp({D.prep}, pr{ip}) & [D.atp] == a);
        r(ip, 1:nV) = Fz(ka) ./ Fz(k8); r(ip, nV + 1) = D(ka).h(1) / D(k8).h(1);   % last column: the 1.10 Lo hold
        first8(ip) = D(k8).slot == 1;
    end
    G.(sprintf('a%d', a)) = struct('r', r, 'first8', first8);
    lab = [arrayfun(@(v) sprintf('%g', v), vS', 'uni', 0), {'hold110'}];
    for j = 1:nV + 1
        a8 = r(first8, j); aL = r(~first8, j);
        sw = sqrt((var(a8) * (numel(a8) - 1) + var(aL) * (numel(aL) - 1)) / (numel(a8) + numel(aL) - 2));   % within-order sd
        B = [B; table(a, string(lab{j}), mean(a8), std(a8) / sqrt(numel(a8)), mean(aL), std(aL) / sqrt(numel(aL)), ...
            sqrt(mean(a8) * mean(aL)), sqrt(mean(a8) / mean(aL)), mean(r(:, j)), std(r(:, j)), sw, sw / sqrt(size(r, 1)), ...
            'VariableNames', {'ATP', 'v', 'r8first', 'sem8', 'rLowFirst', 'semLow', 'A', 'rho', 'mean', 'sdAll', 'sdWithinOrder', 'semBalanced'})]; %#ok<AGROW>
    end
end
fprintf('\n3+4. ATP EFFECT BY RUNNING ORDER (F([ATP])/F(8 mM) within a preparation, baseline-subtracted)\n');
fprintf('   A = sqrt(r8first*rLowFirst) is the order-balanced ATP effect, rho the rundown of the second record;\n');
fprintf('   semBalanced = within-order sd / sqrt(6), the uncertainty of the balanced mean once order is a fixed effect\n');
disp(B);

%% 5. between-preparation spread vs effect size
fprintf('\n5. SPREAD vs EFFECT (F-V points normalised by the preparation''s 8 mM isometric force at 1.00 Lo)\n');
SP = table();
for a = ATP
    ks = find([D.atp] == a); P = nan(numel(ks), nV);
    for c = 1:numel(ks)
        k = ks(c); k8 = find(strcmp({D.prep}, D(k).prep) & [D.atp] == 8);
        ref = D(k8).R.Fx(end) - D(k8).R.F0(end);
        P(c, :) = Fz(k) / ref;
    end
    for j = 1:nV
        SP = [SP; table(a, vS(j), numel(ks), mean(P(:, j)), std(P(:, j)), std(P(:, j)) / mean(P(:, j)), std(P(:, j)) / sqrt(numel(ks)) / mean(P(:, j)), ...
            'VariableNames', {'ATP', 'v', 'n', 'mean', 'sd', 'cv', 'relSem'})]; %#ok<AGROW>
    end
    Pn.(sprintf('a%d', a)) = P;
end
disp(SP);

%% 6. outlying preparations (8 mM F-V shape, leave-one-out z)
P8 = Pn.a8; k8s = find([D.atp] == 8); zo = nan(size(P8));
for c = 1:size(P8, 1)
    o = setdiff(1:size(P8, 1), c); zo(c, :) = (P8(c, :) - mean(P8(o, :))) ./ std(P8(o, :));
end
fprintf('\n6. OUTLIERS at 8 mM: leave-one-out z of each preparation''s F-V points (v = %s ML/s)\n', mat2str(vS(2:end)'));
for c = 1:numel(k8s)
    fprintf('   %-12s slot %d  Fiso %.1f kPa | %s | max |z| %.1f\n', D(k8s(c)).date, D(k8s(c)).slot, Q.Fiso100(k8s(c)), ...
        sprintf('%+5.1f ', zo(c, 2:end)), max(abs(zo(c, 2:end))));
end

%% figure
fig = figure('Visible', 'off', 'Position', [50 50 1250 420]); tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
col = [0 0 0; 0.85 0.45 0; 0.75 0 0];
nexttile(tl); hold on;                                % (a) hold force across the record
for ia = 1:3
    ks = find([D.atp] == ATP(ia)); H = cell2mat(arrayfun(@(k) D(k).h / D(k).h(1), ks', 'uni', 0));
    errorbar((1:7) + 0.08 * (ia - 2), mean(H), std(H), 'o-', 'Color', col(ia, :), 'MarkerFaceColor', 'w', 'LineWidth', 1.3, ...
        'DisplayName', sprintf('%d mM (n=%d)', ATP(ia), numel(ks)));
end
yline(1, 'k:', 'HandleVisibility', 'off'); xlabel('cycle within the record (0.5 s apart)'); ylabel('1.10 Lo hold force / first hold');
title('(a) Drift within a record'); legend('Location', 'southwest'); box on; grid on; xlim([0.5 7.5]);
for ia = 2:3                                          % (b, c) ATP ratio by running order
    a = ATP(ia); g = G.(sprintf('a%d', a)); nexttile(tl); hold on;
    b = B(B.ATP == a & B.v ~= "hold110", :); x = str2double(b.v);
    errorbar(x - 0.05, b.r8first, b.sem8, 's-', 'Color', [0.2 0.4 0.8], 'MarkerFaceColor', 'w', 'LineWidth', 1.3, 'DisplayName', '8 mM ran first (n=3)');
    errorbar(x + 0.05, b.rLowFirst, b.semLow, 'd-', 'Color', [0.1 0.6 0.3], 'MarkerFaceColor', 'w', 'LineWidth', 1.3, 'DisplayName', sprintf('%d mM ran first (n=3)', a));
    plot(x, b.A, 'k-', 'LineWidth', 2, 'DisplayName', 'order-balanced A');
    mf = fullfile(fileparts(here), 'FV_Beard2022', 'results', 'metrics.mat');
    if isfile(mf)
        M = load(mf, 'RT'); m = M.RT(M.RT.ATP == a, :);
        plot(m.v, m.model, 'x--', 'Color', col(ia, :), 'LineWidth', 1.5, 'MarkerSize', 9, 'DisplayName', 'model (gaps = slack)');
    end
    yline(1, 'k:', 'HandleVisibility', 'off'); xlabel('shortening velocity (ML/s)'); ylabel(sprintf('F(%d mM) / F(8 mM)', a));
    title(sprintf('(%s) %d/8 mM by running order (mean \\pm sem)', char('a' + ia - 1), a)); legend('Location', 'northeast', 'FontSize', 8);
    box on; grid on; xlim([-0.3 6.5]); ylim([0 1.6]);
end
exportgraphics(fig, fullfile(figP, 'audit_rundown.png'), 'Resolution', 150); close(fig);
if ~isfolder(fullfile(here, 'results')), mkdir(fullfile(here, 'results')); end
save(fullfile(here, 'results', 'audit.mat'), 'Q', 'B', 'SP', 'zo', 'G');
fprintf('\nfigure -> %s\n', fullfile(figP, 'audit_rundown.png'));
