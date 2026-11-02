% SummarizeSrxProtocols.m -- discriminability of the candidate protocols (RunSrxProtocols).
% Each read-out of each variant is compared with the incumbent (ref) and divided by an
% ASSUMED within-preparation measurement CV:
%   rates (k63)               8 %   (slack ktr ratio CV, SlackDataAnalysis)
%   force ratios              5 %   (normalised amplitudes 4-7 %, SlackDataAnalysis)
%   activation rise t50/t90  10 %   (no data here; a guess for a solution-switch / caged-Ca rig)
% |d|/CV >= 2 is called detectable in one preparation.
% Output: results/protocols.txt, results/protocols_summary.json, results/protocols.png
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsQ = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
S = load(fullfile(rsQ, 'protocols.mat'));
VP = S.VP; J = S.J; RES = S.RES; nV = numel(VP); nm = {VP.name};
CV = struct('rate', 0.08, 'force', 0.05, 'rise', 0.10);
jx  = @(v, kind, x) find([J.v] == v & strcmp({J.kind}, kind) & abs([J.x] - x) < 1e-12, 1);
val = @(i, f) iField(RES{i}, f);

fid = fopen(fullfile(rsQ, 'protocols.txt'), 'w'); c = onCleanup(@() fclose(fid));
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('Candidate protocols, ka-compensated variants (8 mM). [n] = |variant - ref| / CV; >= 2 detectable\n');
out('CV assumed: rates %.0f %%, force ratios %.0f %%, rise times %.0f %%\n\n', 100*CV.rate, 100*CV.force, 100*CV.rise);
for v = 1:nV; out('%-7s ka x%.3f  SRX_iso %.4f\n', nm{v}, S.kaV(v)/S.kaV(1), sum(S.icV(v, :))); end
nErr = sum(cellfun(@(r) isfield(r, 'err'), RES)); out('failed jobs: %d of %d\n', nErr, numel(RES));

% ---- A: activation step from rest
out('\nA  activation step from rest (relaxed steady state, then attachment on)\n');
out('%-7s %8s %8s %8s %8s\n', 'variant', 'SRXrest', 'Frest', 't50 ms', 't90 ms');
A = nan(nV, 3);
for v = 1:nV
    i = jx(v, 'A', 0); A(v, :) = [val(i, 'SRXrest'), val(i, 't50'), val(i, 't90')];
    out('%-7s %8.3f %8.2f %8.2f %8.2f', nm{v}, A(v, 1), val(i, 'Frest'), 1e3*A(v, 2), 1e3*A(v, 3));
    if v > 1; out('   [t50 %.1f, t90 %.1f]', abs(log(A(v, 2)/A(1, 2)))/CV.rise, abs(log(A(v, 3)/A(1, 3)))/CV.rise); end
    out('\n');
end

% ---- B: release-hold-restretch series
TB = S.TB; B = nan(nV, numel(TB)); Bs = B; Bm = B;
for v = 1:nV
    for k = 1:numel(TB)
        i = jx(v, 'B', TB(k)); B(v, k) = val(i, 'k63'); Bs(v, k) = val(i, 'SRXatRestretch'); Bm(v, k) = val(i, 'SRXmax');
    end
end
out('\nB  release 0.15 ML, hold T at short length, restretch in 2 ms; k63 after restretch (1/s)\n');
out('%-7s %s | [abs, shape k(T)/k(5 ms)]\n', 'variant', sprintf('%8.0fms', 1e3*TB));
for v = 1:nV
    out('%-7s %s', nm{v}, sprintf('%10.1f', B(v, :)));
    if v > 1
        sh = B(v, :)/B(v, 1); shr = B(1, :)/B(1, 1);
        out(' | [%.1f, %.1f]', max(abs(log(B(v, :)./B(1, :))))/CV.rate, max(abs(log(sh./shr)))/(sqrt(2)*CV.rate));
    end
    out('\n');
end
out('SRX peak during the hold / SRX at restretch:\n');
for v = 1:nV; out('%-7s %s\n', nm{v}, strjoin(arrayfun(@(k) sprintf('%.2f/%.2f', Bm(v, k), Bs(v, k)), 1:numel(TB), 'uni', 0), '  ')); end

% ---- C: paired release-restretch
DC = S.DC; C = nan(nV, numel(DC)); Cs = C;
for v = 1:nV
    for k = 1:numel(DC); i = jx(v, 'C', DC(k)); C(v, k) = val(i, 'ratio'); Cs(v, k) = val(i, 'SRXpre2'); end
end
out('\nC  two ktr manoeuvres (T = 20 ms), second/first k63 vs interval\n');
out('%-7s %s\n', 'variant', sprintf('%8.0fms', 1e3*DC));
for v = 1:nV
    out('%-7s %s', nm{v}, sprintf('%10.3f', C(v, :)));
    if v > 1; out(' | [%.1f]', max(abs(C(v, :) - C(1, :)))/(sqrt(2)*CV.rate)); end
    out('\n');
end

% ---- D: isometric SL series
SLs = S.SLs; D = nan(nV, numel(SLs)); Ds = D; i22 = find(abs(SLs - 2.2) < 1e-9);
for v = 1:nV
    for k = 1:numel(SLs); i = jx(v, 'D', SLs(k)); D(v, k) = val(i, 'F'); Ds(v, k) = val(i, 'SRX'); end
end
out('\nD  isometric SL series: F(SL)/F(2.2) [SRX]\n');
out('%-7s %s\n', 'variant', sprintf('%15.1fum', SLs));
for v = 1:nV
    out('%-7s %s', nm{v}, strjoin(arrayfun(@(k) sprintf('%8.3f[%.3f]', D(v, k)/D(v, i22), Ds(v, k)), 1:numel(SLs), 'uni', 0), ' '));
    if v > 1; out(' | [%.1f]', max(abs(D(v, :)/D(v, i22) - D(1, :)/D(1, i22)))/CV.force); end
    out('\n');
end

%% ---- figure ---------------------------------------------------------------------
cc = lines(nV);
fg = figure('Color', 'w', 'Position', [30 30 1500 900], 'Visible', 'off');
tiledlayout(fg, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for col = [2 3]
    nexttile; hold on; box on; grid on;
    for v = 1:nV
        r = RES{jx(v, 'A', 0)}; if isfield(r, 'tr'); plot(1e3*r.tr(:, 1), r.tr(:, col), 'Color', cc(v, :), 'LineWidth', 1.5); end
    end
    xlim([0 150]); xlabel('t after activation (ms)');
    if col == 2; ylabel('Force (kPa)'); title('A  activation step from rest'); legend(nm, 'Location', 'southeast');
    else; ylabel('SRX (SR + SRD)'); title('A  SRX drains on activation'); end
end
nexttile; hold on; box on; grid on;
for v = 1:nV; plot(1e3*TB, B(v, :), '-o', 'Color', cc(v, :), 'MarkerFaceColor', cc(v, :), 'LineWidth', 1.5); end
set(gca, 'XScale', 'log'); xlabel('hold at short length T (ms)'); ylabel('k_{63} after restretch (s^{-1})'); title('B  release-hold-restretch');
nexttile; hold on; box on; grid on;
for v = 1:nV
    r = RES{jx(v, 'B', 0.02)}; if isfield(r, 'tr'); plot(1e3*r.tr(:, 1), r.tr(:, 3), 'Color', cc(v, :), 'LineWidth', 1.5); end
end
xlim([-5 120]); xlabel('t (ms): release at 0, restretch 21-23 ms'); ylabel('SRX'); title('B  SRX through a 20 ms hold');
nexttile; hold on; box on; grid on;
for v = 1:nV; plot(1e3*DC, C(v, :), '-o', 'Color', cc(v, :), 'MarkerFaceColor', cc(v, :), 'LineWidth', 1.5); end
set(gca, 'XScale', 'log'); xlabel('interval between manoeuvres (ms)'); ylabel('k_2 / k_1'); title('C  paired release-restretch');
nexttile; hold on; box on; grid on;
for v = 1:nV; plot(SLs, D(v, :)/D(v, i22), '-o', 'Color', cc(v, :), 'MarkerFaceColor', cc(v, :), 'LineWidth', 1.5); end
xlabel('SL (\mum)'); ylabel('F / F(2.2)'); title('D  isometric length dependence');
exportgraphics(fg, fullfile(rsQ, 'protocols.png'), 'Resolution', 110);

Jn = struct('variants', {nm}, 'A_SRXrest', A(:, 1)', 'A_t50', A(:, 2)', 'A_t90', A(:, 3)', 'B_T', TB, 'B_k63', B, ...
    'B_SRXrestretch', Bs, 'C_D', DC, 'C_ratio', C, 'D_SL', SLs, 'D_F', D, 'D_SRX', Ds, 'CV', CV);
fj = fopen(fullfile(rsQ, 'protocols_summary.json'), 'w'); fprintf(fj, '%s', jsonencode(Jn, 'PrettyPrint', true)); fclose(fj);
fprintf('DONE SummarizeSrxProtocols\n');

function v = iField(r, f)
    v = NaN; if isstruct(r) && isfield(r, f) && ~isempty(r.(f)); v = r.(f); end
end
