% AnalyzeLevers.m — rank the lever scan by what actually matters:
% how much does a lever SLOW the post-restretch recovery, and what does that
% cost everywhere else?
%
% Needed slowing: ln(43.7/118.2) = -0.995 in dlnRsK. A lever is "useful" only
% if it delivers a decent fraction of that at a small dFeat. Merges
% results/lever_scan.mat (group 1) and results/lever_scan2.mat (group 2).

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

rows = struct([]);
S1 = fullfile(here,'results','lever_scan.mat');
S2 = fullfile(here,'results','lever_scan2.mat');
if isfile(S1); a = load(S1); rows = mergeRows(rows, a.rows); R0 = a.R0; end
if isfile(S2); b = load(S2); rows = mergeRows(rows, b.rows2); if ~exist('R0','var'); R0 = b.R0; end; end

ok = arrayfun(@(r) ~isempty(r.rsK) && isfinite(r.rsK) && isempty(r.err), rows);
r  = rows(ok);
fprintf('%d usable evaluations (of %d)\n', numel(r), numel(rows));

NEEDED = log(43.738 / mean(R0.rsK_m,'omitnan'));
fprintf('baseline rsK %.1f, data %.1f -> needed dlnRsK = %.3f\n\n', ...
        mean(R0.rsK_m,'omitnan'), 43.738, NEEDED);

%% ---- 1. everything that slows it, ranked by amount ----------------------
sl = r([r.dlnRsK] < -0.02);
[~, o] = sort([sl.dlnRsK]);
fprintf('LEVERS THAT SLOW THE POST-RESTRETCH RECOVERY\n');
fprintf('%-32s %6s %8s %7s %9s %8s %8s\n', 'param','fact','dlnRsK','%need','dFeat','ktr','steady');
for i = o
    fprintf('%-32s %6.2f %8.3f %6.0f%% %9.2f %8.1f %8.1f\n', ...
        sl(i).param, sl(i).factor, sl(i).dlnRsK, 100*sl(i).dlnRsK/NEEDED, ...
        sl(i).dFeat, sl(i).ktr, sl(i).steady);
end

%% ---- 2. efficiency: slowing per unit of damage --------------------------
fprintf('\nEFFICIENCY (dlnRsK per unit dFeat; only levers with dFeat > 0.05)\n');
eff = sl([sl.dFeat] > 0.05);
e   = arrayfun(@(x) -x.dlnRsK / x.dFeat, eff);
[~, o2] = sort(e, 'descend');
fprintf('%-32s %6s %8s %9s %10s\n','param','fact','dlnRsK','dFeat','slow/cost');
for i = o2(1:min(12,numel(o2)))
    fprintf('%-32s %6.2f %8.3f %9.2f %10.4f\n', eff(i).param, eff(i).factor, ...
        eff(i).dlnRsK, eff(i).dFeat, e(i));
end

%% ---- 3. the free ones: slow it AND improve the rest ---------------------
free = sl([sl.dFeat] < 0);
fprintf('\nFREE LEVERS (slow rsK *and* lower the total feature cost): %d\n', numel(free));
for i = 1:numel(free)
    fprintf('  %-30s x%-5.2g  dlnRsK %+6.3f  dFeat %+7.3f\n', ...
        free(i).param, free(i).factor, free(i).dlnRsK, free(i).dFeat);
end

%% ---- 4. how far can the best single lever get? -------------------------
[best, ib] = min([sl.dlnRsK]);
fprintf('\nBEST SINGLE LEVER: %s x%.2g -> dlnRsK %.3f = %.0f%% of what is needed\n', ...
    sl(ib).param, sl(ib).factor, best, 100*best/NEEDED);
fprintf('   at a cost of dFeat %+.2f (ktr %.1f, steady %.1f vs baseline %.1f)\n', ...
    sl(ib).dFeat, sl(ib).ktr, sl(ib).steady, mean(R0.steady_m,'omitnan'));

%% ---- figure ------------------------------------------------------------
fig = figure(705); clf; set(fig,'Position',[60 60 980 560],'Color','w'); hold on; box on;
x = [r.dlnRsK]; y = [r.dFeat];
% signed log on the cost axis: dFeat spans -0.2 .. +550, and a few levers are
% genuinely free (negative), so a plain log scale would drop them.
ys = sign(y) .* log10(1 + abs(y));
scatter(x, ys, 46, 'filled', 'MarkerFaceAlpha', 0.65);
big = abs(x) > 0.10 | y > 40;
for i = find(big)
    text(x(i), ys(i), ['  ' r(i).param], 'FontSize', 7.5, 'Interpreter','none');
end
xline(NEEDED, 'r--', 'LineWidth', 1.6, 'Label','needed slowing');
xline(0, 'k-'); yline(0, 'k-');
tk = [-0.2 0 0.5 2 10 50 200 550];
set(gca, 'YTick', sign(tk).*log10(1+abs(tk)), ...
         'YTickLabel', arrayfun(@(v) sprintf('%g',v), tk, 'UniformOutput', false));
xlabel('d ln(rsK)   (negative = slower recovery = better)');
ylabel('\Delta total feature cost  (signed log; positive = worse elsewhere)');
title('One-at-a-time lever scan: the useful corner (bottom-left) is empty');
exportgraphics(fig, fullfile(here,'results','lever_scan.png'), 'Resolution',140);
fprintf('\nwrote results/lever_scan.png\n');

function A = mergeRows(A, B)
    f = {'param','factor','rsK','rsKratio','ktr','steady','rsA','featTotal','dFeat','dlnRsK','err'};
    for i = 1:numel(B)
        s = struct();
        for k = 1:numel(f)
            if isfield(B, f{k}); s.(f{k}) = B(i).(f{k}); else; s.(f{k}) = []; end
        end
        if isempty(A); A = s; else; A(end+1) = s; end %#ok<AGROW>
    end
end
