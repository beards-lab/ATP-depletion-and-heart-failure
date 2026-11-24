% AnalyzeFvid.m -- which parameters do the Beard2022 isovelocity records constrain,
% compared with the slack battery? Local identifiability at the incumbent.
%
% Needs results/jac_beard and results/jac_battery (RunFvidJacobian.m). Minutes.
%
% NOTE (2026-11-10): MATLAB could not start when the Jacobian runs finished, so the
% numbers and figures in conclusions.md come from analyzeFvid.py, a Python twin
% that reads this script's data features (results/fvid*.mat) and reproduces its
% model feature vectors to 4e-13. This script has run on PARTIAL results only;
% it lacks the twin's later additions (tiers Bnr/Bok/Cok, single-parameter
% recovery) and writes the same figure names. RunFvidProfile.m takes its
% inputs from this script (results/fvid.mat, profile_plan.mat) or, when those
% are absent, from the twin (results/fvid_py.mat).
%
% Method
%   Feature vector f(theta) per data tier, each element with the standard error
%   of its mean, sem_i (between-preparation sd / sqrt(n)):
%     A   Beard2022 F-V points, 8/4/2 mM           (23 elements)
%     B8  Beard2022 full records, 8 mM only
%     B   Beard2022 full records, 8/4/2 mM         (points + trace samples)
%     C   slack battery (slack + ktr + passive + F-V), 8 mM, objective floors
%     D   B and C together
%   Jacobian in log-parameters by a one-sided x1.10 step:
%     J_ij = (f_i(theta_j x 1.10) - f_i(theta)) / (sem_i ln 1.10)
%   Fisher information F = J'J. For parameter j:
%     sFix   = 1/sqrt(F_jj)                    its log-uncertainty with all others known
%     sFree  = sqrt([(F + I/ln(10)^2)^-1]_jj)  with all others free, under a weak prior
%                                              of one decade on every parameter
%   Verdict (ln 2 = known to within a factor of two):
%     identified    sFree <= ln 2
%     sloppy        sFix <= ln 2 < sFree   (constrained only in combination)
%     unidentified  sFix > ln 2            (the data barely move when it changes)
%
% Output -> results/fvid.mat, figures/verdicts.csv, figures/fvid_*.png, tables below.

here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
cd(root); addpath(genpath(root));
wt = strsplit(genpath(fullfile(root, '.claude')), pathsep); wt = wt(~cellfun(@isempty, wt));
if ~isempty(wt), rmpath(wt{:}); end
figP = fullfile(here, 'figures'); if ~isfolder(figP), mkdir(figP); end
LN = log(1.10); SIG_PRIOR = log(10); THR = log(2); SD_FLOOR = 0.01;
[pBat, ~, ~, recs] = fvidBase();

%% 1. Beard2022: data features and their standard errors -------------------------
D = loadBeard2022FV(); nS = numel(D);
for k = 1:nS, R = readIsoV(D(k).t, D(k).L, D(k).F); D(k).iso = R.Fx(end) - R.F0(end); end
Fb = cell(1, 3);
for a = 1:3, r = load(fullfile(here, 'results', 'jac_beard', sprintf('e001_r%d.mat', a))); Fb{a} = r.r.F; end
[~, ~, S] = beardVector(recs, Fb, []);                  % the sample plan
dPts = nan(1, height(S.pt)); sPts = dPts; nPts = dPts; dTrc = nan(1, height(S.tr)); sTrc = dTrc; nTrc = dTrc;
for a = 1:3
    ks = find([D.atp] == recs(a).atp); P = []; T = [];
    for k = ks
        k8 = find(strcmp({D.prep}, D(k).prep) & [D.atp] == 8);
        rk = recs; rk(a).t = D(k).t; rk(a).L = D(k).L; Fk = {D(k).F, D(k).F, D(k).F};
        [p, t] = beardVector(rk, Fk, D(k8).iso, S);
        P(end+1, :) = p(S.pt.ATP == recs(a).atp); T(end+1, :) = t(S.tr.ATP == recs(a).atp); %#ok<SAGROW>
    end
    m = S.pt.ATP == recs(a).atp; dPts(m) = mean(P); sPts(m) = std(P); nPts(m) = size(P, 1);
    m = S.tr.ATP == recs(a).atp; dTrc(m) = mean(T); sTrc(m) = std(T); nTrc(m) = size(T, 1);
end
semPts = max(sPts, SD_FLOOR) ./ sqrt(nPts); semTrc = max(sTrc, SD_FLOOR) ./ sqrt(nTrc);

%% 2. model vectors for every evaluation -----------------------------------------
Eb = load(fullfile(here, 'results', 'jac_beard', 'evals.mat')); names = Eb.pool;
if any(arrayfun(@(e) ischar(e.name) && strcmp(e.name, 'pAtp_k2'), Eb.E)) && ~any(strcmp(names, 'pAtp_k2')), names{end+1} = 'pAtp_k2'; end  % the ATP gate is only in the Beard evaluations
P = numel(names);
Vp = nan(numel(Eb.E), height(S.pt)); Vt = nan(numel(Eb.E), height(S.tr));
% A record that is missing leaves NaN in its own elements only (the 8 mM record
% sets the force scale, so without it the whole evaluation is NaN).
for e = 1:numel(Eb.E)
    Fe = Fb; av = false(1, 3);
    for a = 1:3
        f = fullfile(here, 'results', 'jac_beard', sprintf('e%03d_r%d.mat', e, a));
        if ~isfile(f), continue; end
        try r = load(f); catch, continue; end
        if r.r.ok, Fe{a} = r.r.F; av(a) = true; end
    end
    if ~av(1), continue; end
    [vp, vt] = beardVector(recs, Fe, [], S);
    for a = find(~av), vp(S.pt.ATP == recs(a).atp) = NaN; vt(S.tr.ATP == recs(a).atp) = NaN; end
    Vp(e, :) = vp; Vt(e, :) = vt;
end
Ec = load(fullfile(here, 'results', 'jac_battery', 'evals.mat'));
r = load(fullfile(here, 'results', 'jac_battery', 'e001_r1.mat'));
[v0, dBat, semBat, labBat, skipped] = batteryVector(r.r.fm, r.r.fd, pBat.fn);
Vc = nan(numel(Ec.E), numel(v0)); Vc(1, :) = v0;
for e = 2:numel(Ec.E)
    f = fullfile(here, 'results', 'jac_battery', sprintf('e%03d_r1.mat', e));
    if ~isfile(f), continue; end
    try r = load(f); catch, continue; end
    if ~r.r.ok, continue; end
    v = batteryVector(r.r.fm, r.r.fd, pBat.fn); if numel(v) == numel(v0), Vc(e, :) = v; end
end
okB = isfinite(v0) & isfinite(dBat) & isfinite(semBat);
fprintf('Beard evaluations: %d of %d complete. Battery: %d of %d.\n', sum(all(isfinite(Vp), 2)), numel(Eb.E), sum(all(isfinite(Vc(:, okB)), 2)), numel(Ec.E));
fprintf('Battery entries left out of the information (priors, not data):\n'); fprintf('   %s\n', skipped{:});

%% 3. Jacobians --------------------------------------------------------------------
kindB = arrayfun(@(e) string(e.kind), Eb.E); kindC = arrayfun(@(e) string(e.kind), Ec.E);
col = @(V, E, kind, nm, sem, k, h) (V(find(kind == k & arrayfun(@(e) ischar(e.name) && strcmp(e.name, nm), E), 1), :) - V(1, :))' ./ sem(:) / h;
Jp = nan(height(S.pt), P); Jt = nan(height(S.tr), P); Jc = zeros(sum(okB), P); labOk = labBat(okB);
% A battery feature can vanish under a perturbation (the re-stretch peak detector
% finds no peak and returns NaN). That is a large change, not a small one, but it
% has no derivative; its sensitivity is set to 0, which under-states the battery.
for j = 1:P
    Jp(:, j) = col(Vp, Eb.E, kindB, names{j}, semPts, "jac", LN);
    Jt(:, j) = col(Vt, Eb.E, kindB, names{j}, semTrc, "jac", LN);
    if any(strcmp(Ec.pool, names{j}))
        c = col(Vc, Ec.E, kindC, names{j}, semBat, "jac", LN); c = c(okB); lost = ~isfinite(c);
        if any(lost) && ~all(lost)
            fprintf('  battery, %s x1.10: %d features vanish (%s); sensitivity set to 0\n', names{j}, sum(lost), strjoin(labOk(lost), ', '));
            c(lost) = 0;
        end
        Jc(:, j) = c;
    end
end
m8p = S.pt.ATP == 8; m8t = S.tr.ATP == 8;
% B8n3 and Bn3: the same records with every standard error rescaled to n = 3
% preparations (the slack battery's n), to separate what the protocol carries
% from what the larger sample carries.
n3 = @(J, n) J ./ sqrt(n(:) / 3);
tiers = struct('key', {'A', 'B8', 'B8n3', 'B', 'Bn3', 'C', 'D'}, ...
    'label', {'Beard F-V points, 8/4/2 mM', 'Beard records, 8 mM only', 'Beard records, 8 mM, as if n = 3', 'Beard records, 8/4/2 mM', ...
    'Beard records, 8/4/2 mM, as if n = 3', 'slack battery, 8 mM', 'Beard records + slack battery'}, ...
    'J', {Jp, [Jp(m8p, :); Jt(m8t, :)], [n3(Jp(m8p, :), nPts(m8p)); n3(Jt(m8t, :), nTrc(m8t))], [Jp; Jt], ...
    [n3(Jp, nPts); n3(Jt, nTrc)], Jc, [Jp; Jt; Jc]});
% A tier is analysed when at least 10 of its Jacobian columns are complete; all
% analysed tiers share one parameter set (the columns complete in every one of
% them), because the uncertainty of a parameter depends on which others are free.
fin = cell2mat(cellfun(@(J) all(isfinite(J), 1), {tiers.J}', 'uni', 0));
runT = sum(fin, 2) >= 10; use = all(fin(runT, :), 1); complete = all(runT) && all(use);
if ~all(runT), fprintf('PARTIAL RESULTS. Tiers without enough evaluations yet, left out: %s\n', strjoin({tiers(~runT).key}, ', ')); end
if ~all(use), fprintf('PARTIAL RESULTS. %d of %d parameters analysed; evaluations missing for: %s\n', sum(use), P, strjoin(names(~use), ', ')); end
tiers = tiers(runT); keys = {tiers.key}; has = @(k) any(strcmp(keys, k));
nmU = names(use); PU = sum(use);

%% 4. Fisher information and verdicts ----------------------------------------------
V = table(nmU(:), cellfun(@paramGroup, nmU(:), 'uni', 0), 'VariableNames', {'param', 'group'});
SP = struct();
for t = 1:numel(tiers)
    J = tiers(t).J(:, use); F = J' * J;
    ev = sort(eig((F + F') / 2), 'descend');
    sFix = 1 ./ sqrt(diag(F)); C = inv(F + eye(PU) / SIG_PRIOR^2); sFree = sqrt(diag(C));
    cls = repmat("unidentified", PU, 1); cls(sFix <= THR) = "sloppy"; cls(sFree <= THR) = "identified";
    if strcmp(tiers(t).key, 'C'), cls(strcmp(nmU, 'pAtp_k2')) = "not in tier"; end
    V.(['sFix_' tiers(t).key]) = sFix; V.(['sFree_' tiers(t).key]) = sFree; V.(['class_' tiers(t).key]) = cls;
    SP.(tiers(t).key) = struct('ev', ev, 'F', F, 'J', J, 'C', C, 'nEl', size(J, 1));
    fprintf('\n%-3s %-32s %4d elements | identified %2d, sloppy %2d, unidentified %2d | stiff directions (lambda > 1/ln2^2): %d | log10 lambda: max %.1f, median %.1f\n', ...
        tiers(t).key, tiers(t).label, size(J, 1), sum(cls == "identified"), sum(cls == "sloppy"), sum(cls == "unidentified"), ...
        sum(ev > 1 / THR^2), log10(ev(1)), log10(max(median(ev), 1e-12)));
end
[~, og] = sort(V.group); V = V(og, :);
fprintf('\nPer-parameter verdicts (sFree in log units: 0.69 = factor of two; 2.30 = the prior, i.e. no information):\n');
disp(V(:, [{'param', 'group'}, reshape([strcat('sFree_', keys); strcat('class_', keys)], 1, [])]));
fprintf('Identified parameters by mechanism group (count of group size):\n');
G = groupsummary(V, 'group', @(x) sum(x == "identified"), strcat('class_', keys)); disp(G);
only = @(a, b) strjoin(V.param(V.(['class_' a]) == "identified" & V.(['class_' b]) ~= "identified")', ', ');
bk = 'B'; if ~has('B'), bk = 'B8'; end                   % the Beard-record tier available
if has(bk) && has('C')
    fprintf('Identified by the slack battery but NOT by the Beard records (C not %s):\n   %s\n', bk, only('C', bk));
    fprintf('Identified by the Beard records but NOT by the slack battery (%s not C):\n   %s\n', bk, only(bk, 'C'));
    for c = ["sloppy", "unidentified"]
        fprintf('%s under both %s and C:\n   %s\n', c, bk, strjoin(V.param(V.(['class_' bk]) == c & V.class_C == c)', ', '));
    end
end
if has('D') && has('B') && has('C')
    fprintf('Identified only by the two together (D, neither B nor C):\n   %s\n', ...
        strjoin(V.param(V.class_D == "identified" & V.class_B ~= "identified" & V.class_C ~= "identified")', ', '));
end

%% 5. step and side check ------------------------------------------------------------
SC = table();
if has('B'), sets = {'B', [Vp, Vt], [semPts, semTrc], Eb.E, kindB};
else, sets = {'B8', [Vp(:, m8p), Vt(:, m8t)], [semPts(m8p), semTrc(m8t)], Eb.E, kindB}; end
sets = [sets; {'C', Vc(:, okB), semBat(okB), Ec.E, kindC}];
for s = 1:size(sets, 1)
    Vs = sets{s, 2}; se = sets{s, 3}; E = sets{s, 4}; kd = sets{s, 5};
    for e = find(kd == "half")
        nm = E(e).name; j = col(Vs, E, kd, nm, se, "jac", LN);
        jh = col(Vs, E, kd, nm, se, "half", log(1.05)); jm = col(Vs, E, kd, nm, se, "minus", -LN);
        ok3 = isfinite(j) & isfinite(jh) & isfinite(jm);    % features present at all three steps
        if mean(ok3) < 0.5, continue; end
        j = j(ok3); jh = jh(ok3); jm = jm(ok3);
        cs = @(a, b) (a' * b) / max(norm(a) * norm(b), eps);
        SC = [SC; table(string(sets{s, 1}), string(nm), norm(j), norm(jh) / norm(j), cs(j, jh), norm(jm) / norm(j), cs(j, jm), sum(~ok3), ...
            'VariableNames', {'tier', 'param', 'normJ', 'ratioHalf', 'cosHalf', 'ratioMinus', 'cosMinus', 'featuresLost'})]; %#ok<AGROW>
    end
end
fprintf('\nStep check: the Jacobian column from x1.05 and from x1/1.10 against the one used (x1.10).\n');
fprintf('ratio = norm of the alternative / norm used (1 = linear), cos = direction agreement (1 = same).\n');
disp(SC);
for tr = unique(SC.tier)'
    q = SC(SC.tier == tr, :); w = q.normJ > 1;         % columns that carry information
    fprintf('  tier %s: %d parameters (%d with |J| > 1). Half step: ratio median %.2f (%.2f-%.2f), cos median %.3f. Other side: ratio median %.2f (%.2f-%.2f), cos median %.3f\n', ...
        tr, height(q), sum(w), median(q.ratioHalf(w)), min(q.ratioHalf(w)), max(q.ratioHalf(w)), median(q.cosHalf(w)), ...
        median(q.ratioMinus(w)), min(q.ratioMinus(w)), max(q.ratioMinus(w)), median(q.cosMinus(w)));
end

%% 6. synthetic recovery ---------------------------------------------------------------
% Truth = the incumbent with every parameter moved by exp(0.10 N(0,1)); its features
% come from a real model run. The estimate is one regularised Gauss-Newton step
% from the incumbent, without and with measurement noise at the data's own sem.
RC = table(); rng(7); recov = struct();
synT = ["A", "B", "C"]; if ~complete, synT = strings(1, 0); fprintf('\n(synthetic recovery waits for the complete set)\n'); end
for tk = synT
    Fi = SP.(tk).F; J = SP.(tk).J; Cc = SP.(tk).C;
    switch tk
        case "A", Vs = Vp; se = semPts; E = Eb.E; kd = kindB;
        case "B", Vs = [Vp, Vt]; se = [semPts, semTrc]; E = Eb.E; kd = kindB;
        case "C", Vs = Vc(:, okB); se = semBat(okB); E = Ec.E; kd = kindC;
    end
    dT = []; dH = []; dN = [];
    for e = find(kd == "synth")
        fe = isfinite(Vs(e, :))';                          % features that exist at this truth
        fprintf('  tier %s, truth %d: %d of %d features exist\n', tk, e, sum(fe), numel(fe));
        if mean(fe) < 0.3, continue; end
        Jr = J(fe, :); Cr = inv(Jr' * Jr + eye(PU) / SIG_PRIOR^2);
        dtrue = zeros(PU, 1); [tf, loc] = ismember(nmU, E(e).name); dtrue(tf) = log(E(e).f(loc(tf)));
        rr = (Vs(e, fe) - Vs(1, fe))' ./ reshape(se(fe), [], 1);
        dH = [dH, Cr * (Jr' * rr)]; dT = [dT, dtrue]; %#ok<AGROW>
        en = zeros(PU, 1); for q = 1:200, en = en + (Cr * (Jr' * (rr + randn(size(rr)))) - dtrue).^2; end
        dN = [dN, sqrt(en / 200)]; %#ok<AGROW>
    end
    recov.(tk) = struct('dT', dT, 'dH', dH);
    cl = V.(['class_' char(tk)]); [~, ix] = ismember(V.param, nmU);        % V is sorted by group
    for c = ["identified", "sloppy", "unidentified"]
        m = ix(cl == c); if isempty(m) || isempty(dT), continue; end
        x = dT(m, :); y = dH(m, :);
        RC = [RC; table(tk, c, numel(m), rms(x(:)), rms(y(:) - x(:)), (x(:)' * y(:)) / (x(:)' * x(:)), rms(reshape(dN(m, :), [], 1)), ...
            sqrt(mean(diag(Cc(m, m)))), 'VariableNames', {'tier', 'class', 'nParam', 'rmsTruth', 'rmsError', 'slope', 'rmsErrorNoisy', 'predictedSd'})]; %#ok<AGROW>
    end
end
if complete
fprintf('\nSynthetic recovery (3 truths; all parameters moved at once, sd 0.10 in log):\n');
fprintf('rmsError = error of the noise-free estimate; slope = fraction of the true move recovered (1 = all, 0 = none);\n');
fprintf('rmsErrorNoisy = error with noise at the data sem added (200 draws); predictedSd = what the Fisher information predicts\n');
disp(RC);
end

%% 7. how well does the incumbent fit each tier? ------------------------------------------
fprintf('\nFit of the incumbent, chi2 per element against the standard error of the mean:\n');
zA = (Vp(1, :) - dPts) ./ semPts; zT = (Vt(1, :) - dTrc) ./ semTrc; zC = (v0(okB) - dBat(okB)) ./ semBat(okB);
fprintf('   A  F-V points        chi2/N = %6.1f (N = %d)\n', mean(zA.^2), numel(zA));
for a = [8 4 2], fprintf('        %d mM points     chi2/N = %6.1f\n', a, mean(zA(S.pt.ATP == a).^2)); end
fprintf('   B  trace samples     chi2/N = %6.1f (N = %d)\n', mean(zT.^2), numel(zT));
for kd = ["ramp", "redev", "hold110"]
    for a = [8 4 2], m = S.tr.kind == kd & S.tr.ATP == a; fprintf('        %-8s %d mM    chi2/N = %6.1f (N = %d)\n', kd, a, mean(zT(m).^2), sum(m)); end
end
fprintf('   C  slack battery     chi2/N = %6.1f (N = %d)\n', mean(zC.^2), numel(zC));

%% 8. the profile plan for stage 2 ---------------------------------------------------------
% The three parameters the slack battery identifies and the Beard records constrain
% worst. For each, move it by x2 and x1/2 and let every other parameter compensate
% as the Beard-record Fisher information says it can; RunFvidProfile simulates those
% six points on both data sets.
[~, ixV] = ismember(V.param, nmU); plan = struct([]);
if complete
cand = find(V.class_C == "identified" & V.class_B ~= "identified");
if numel(cand) < 3, cand = find(V.class_C == "identified"); end
[~, o] = sort(V.sFree_B(cand), 'descend'); worst = cand(o(1:min(3, numel(cand))));
FB = SP.B.F + eye(PU) / SIG_PRIOR^2; plan = struct('name', {}, 'f', {}, 'kind', {}, 'target', {}, 'delta', {}, 'predChi2B', {}, 'predChi2C', {});
for w = worst(:)'
    j = ixV(w); o = setdiff(1:PU, j);
    for dl = [-THR, THR]
        d = zeros(PU, 1); d(j) = dl; d(o) = -FB(o, o) \ (FB(o, j) * dl);
        d = max(min(d, log(3)), -log(3));                % keep the compensation inside x3
        plan(end+1) = struct('name', {nmU}, 'f', {exp(d')}, 'kind', 'profile', 'target', V.param{w}, 'delta', dl, ...
            'predChi2B', d' * SP.B.F * d, 'predChi2C', d' * SP.C.F * d); %#ok<SAGROW>
    end
end
fprintf('\nProfile plan (stage 2): %s\n', strjoin(unique({plan.target}, 'stable'), ', '));
for i = 1:numel(plan), fprintf('   %-28s x%.2f: predicted chi2 change, Beard records %.2f, slack battery %.1f\n', plan(i).target, exp(plan(i).delta), plan(i).predChi2B, plan(i).predChi2C); end
end

%% 9. figures ---------------------------------------------------------------------------------
set(0, 'DefaultAxesFontSize', 11);
cols = [0.85 0.45 0; 0.55 0.55 0.55; 0.8 0.8 0.8; 0.75 0 0; 0.95 0.6 0.6; 0.1 0.35 0.75; 0 0 0]; cols = cols(runT, :);
sfx = ''; if ~complete, sfx = '_partial'; end
fig = figure('Visible', 'off', 'Position', [50 50 720 430]); hold on;
for t = 1:numel(tiers)
    ev = SP.(tiers(t).key).ev; semilogy(1:PU, max(ev, 1e-6), '-o', 'Color', cols(t, :), 'MarkerSize', 3, 'LineWidth', 1.3, 'DisplayName', sprintf('%s  %s', tiers(t).key, tiers(t).label));
end
set(gca, 'YScale', 'log'); yline(1 / THR^2, 'k--', 'known to a factor of 2', 'HandleVisibility', 'off', 'LabelHorizontalAlignment', 'right');
yline(1 / SIG_PRIOR^2, 'k:', 'the prior (one decade)', 'HandleVisibility', 'off', 'LabelHorizontalAlignment', 'right');
xlabel('parameter combination (eigen-direction), sorted'); ylabel('Fisher information eigenvalue'); ylim([1e-4 1e8]); xlim([1 PU]);
title('How many parameter combinations each data set constrains'); legend('Location', 'northeast', 'FontSize', 8); box on; grid on;
exportgraphics(fig, fullfile(figP, ['fvid_spectrum' sfx '.png']), 'Resolution', 150); close(fig);

fig = figure('Visible', 'off', 'Position', [50 50 760 1250]);
M = cell2mat(cellfun(@(k) V.(['sFree_' k]), keys, 'uni', 0)); M(V.param == "pAtp_k2", strcmp(keys, 'C')) = NaN;
imagesc(log10(M), [log10(0.02) log10(SIG_PRIOR)]); colormap(gca, flipud(parula)); cb = colorbar; cb.Label.String = 'log_{10} uncertainty with all others free (log units)';
set(gca, 'XTick', 1:numel(keys), 'XTickLabel', keys, 'YTick', 1:PU, 'YTickLabel', strcat(V.param, "  [", V.group, "]"), 'TickLabelInterpreter', 'none', 'FontSize', 7);
hold on; cl = cellfun(@(k) V.(['class_' k]), keys, 'uni', 0); cl = [cl{:}];
[yy, xx] = find(cl == "identified"); plot(xx, yy, 'k.', 'MarkerSize', 9);
title({'Per-parameter uncertainty by data tier', 'dot = identified (within a factor of 2)'}, 'FontSize', 10);
exportgraphics(fig, fullfile(figP, ['fvid_parameters' sfx '.png']), 'Resolution', 150); close(fig);

if complete
fig = figure('Visible', 'off', 'Position', [50 50 1000 400]); tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for tk = ["B", "C"]
    nexttile(tl); hold on; R = recov.(tk); cl = V.(['class_' char(tk)]);
    cc = struct('identified', [0 0.5 0.2], 'sloppy', [0.9 0.6 0], 'unidentified', [0.7 0.7 0.7]);
    for c = ["unidentified", "sloppy", "identified"]
        m = ixV(cl == c); if isempty(m) || isempty(R.dT), continue; end
        x = R.dT(m, :); y = R.dH(m, :); plot(x(:), y(:), 'o', 'Color', cc.(c), 'MarkerFaceColor', cc.(c), 'MarkerSize', 4, 'DisplayName', char(c));
    end
    plot([-0.3 0.3], [-0.3 0.3], 'k-', 'HandleVisibility', 'off'); yline(0, 'k:', 'HandleVisibility', 'off');
    xlabel('true change (log)'); ylabel('recovered change (log)'); axis([-0.3 0.3 -0.3 0.3]); box on; grid on;
    title(sprintf('Synthetic recovery, tier %s', tk)); legend('Location', 'northwest', 'FontSize', 8);
end
exportgraphics(fig, fullfile(figP, 'fvid_recovery.png'), 'Resolution', 150); close(fig);
writetable(V, fullfile(figP, 'verdicts.csv'));
save(fullfile(here, 'results', 'profile_plan.mat'), 'plan');
end

save(fullfile(here, 'results', ['fvid' sfx '.mat']), 'V', 'SP', 'SC', 'RC', 'S', 'dPts', 'semPts', 'dTrc', 'semTrc', 'Vp', 'Vt', 'Vc', 'names', 'nmU', 'tiers', 'plan', 'labBat', 'okB', 'semBat', 'dBat');
fprintf('\nfigures -> %s\n', figP);

function g = paramGroup(nm)
% mechanism group of a parameter, for reading the table
if startsWith(nm, 'PieceWise'), g = 'strain-dependence shape';
elseif ismember(nm, {'ka', 'kd', 'k1', 'k_1', 'k2', 'k_2', 'kah', 'kamh'}), g = 'cycle rates';
elseif ismember(nm, {'kstiff1', 'kstiff2', 'kstiff1_n', 'kstiff2_n', 'estiff', 'dr', 'dr1', 'dr2', 'mu', 'mu_neg'}), g = 'cross-bridge mechanics';
elseif ismember(nm, {'kSE', 'ekSE', 'kSE_M', 'eta_M', 'x_M_slack'}), g = 'series elasticity';
elseif ismember(nm, {'ksr0', 'kmsrd', 'ksr2srd', 'ksrd2sr', 'sigma2', 'sigma_srd1', 'ksrd', 'sigma_srd2', 'SRXFromR2HighStrain', 'sSRXrip'}), g = 'SRX';
elseif ismember(nm, {'A0', 'v_ref_reg', 'tau_reg', 'P_bound_max'}), g = 'regulation';
elseif ismember(nm, {'k_cr', 'tau_cr', 'F_cr'}), g = 'recoil';
elseif ismember(nm, {'L_thick', 'L_thin', 'd_optimal', 'sigma_lattice', 'Lsc0'}), g = 'geometry';
elseif ismember(nm, {'gamma', 'k_pas'}), g = 'passive';
elseif strcmp(nm, 'pAtp_k2'), g = 'ATP gate';
else, g = 'other';
end
end
