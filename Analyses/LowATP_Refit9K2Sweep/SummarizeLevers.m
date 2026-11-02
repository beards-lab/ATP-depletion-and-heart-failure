% SummarizeLevers.m -- table, JSON and figures from results/levers.mat.
rootS = 'C:\home\git\ATP-depletion-and-heart-failure';
resS  = fullfile(rootS,'Analyses','LowATP_Refit9K2Sweep','results');
if ~exist('RL','var'); load(fullfile(resS,'levers.mat')); end            %#ok<LOAD>
if ~exist('fmB','var'); load(fullfile(resS,'levers_ref.mat')); end       %#ok<LOAD>

KEYL = {'steady','peak1_y','ktrfit','ktrP_k80','rsK80','t0','rsTfitE','vall_tq','FV2'};
PT8 = pooledTarget8mM(); KT8 = pooledKtrTarget(8); KT2 = pooledKtrTarget(2);
dataR = nan(1,numel(KEYL));
for i = 1:numel(KEYL)
    j = find(strcmp(PT8.name, KEYL{i}), 1);
    if ~isempty(j); dataR(i) = PT8.r(j); end
    j = find(strcmp(KT8.name, KEYL{i}), 1);
    if ~isempty(j); dataR(i) = KT2.val(j)/KT8.val(j); end
end
dataR(strcmp(KEYL,'FV2')) = FV2L.r;
mRatio = @(fm, fv) arrayfun(@(i) iRatio(fm, fm8L, KEYL{i}, fv, fv8L), 1:numel(KEYL));
rB = mRatio(fmB, fvB);
% "misfit" per feature: |log(model ratio / data ratio)|; negative change = better
misB = abs(log(rB ./ dataR));

levs = unique({RL.lever}, 'stable'); nLv = numel(levs);
T = struct(); rows = {};
fid = fopen(fullfile(resS,'levers.txt'),'w');
fprintf(fid, 'Lever sensitivity on %s, k2 gate g = %.2f (R2 at 2 mM); each lever x0.9 / x1.1 at 2 mM only\n', CORE, G_K2);
fprintf(fid, 'reference: score_rel %.4f  score_plain %.4f\n\n', sRB, sPB);
fprintf(fid, '%-8s %-3s %10s %10s %10s %10s | features moved toward data at the better side\n', ...
        'lever','tr','d_rel x0.9','d_rel x1.1','d_pl x0.9','d_pl x1.1');
Dm = nan(nLv, numel(KEYL));   % change of misfit per feature, at the better direction
best = nan(nLv,1); bestF = nan(nLv,1);
for k = 1:nLv
    iLo = find(strcmp({RL.lever}, levs{k}) & [RL.f] < 1, 1);
    iHi = find(strcmp({RL.lever}, levs{k}) & [RL.f] > 1, 1);
    dLo = RL(iLo).score_rel - sRB; dHi = RL(iHi).score_rel - sRB;
    [best(k), s] = min([dLo dHi]); ib = [iLo iHi]; ib = ib(s); bestF(k) = RL(ib).f;
    rK = mRatio(RL(ib).fm, RL(ib).fv2);
    Dm(k,:) = abs(log(rK ./ dataR)) - misB;
    good = KEYL(Dm(k,:) < -0.005); bad = KEYL(Dm(k,:) > 0.005);
    fprintf(fid, '%-8s %-3s %+10.4f %+10.4f %+10.4f %+10.4f | x%.1f better: %s | worse: %s\n', levs{k}, RL(iLo).tier, ...
            dLo, dHi, RL(iLo).score_plain-sPB, RL(iHi).score_plain-sPB, bestF(k), strjoin(good,','), strjoin(bad,','));
    T.(levs{k}) = struct('tier', RL(iLo).tier, 'd_rel', [dLo dHi], 'd_plain', [RL(iLo).score_plain-sPB RL(iHi).score_plain-sPB], ...
                         'best_f', bestF(k), 'better', {good}, 'worse', {bad}, 'ratios_best', rK);
end
fprintf(fid, '\nreference 2/8 ratios (model | data):\n');
for i = 1:numel(KEYL); fprintf(fid, '  %-9s %6.3f | %6.3f\n', KEYL{i}, rB(i), dataR(i)); end
fclose(fid);
J = struct('core', CORE, 'g_k2', G_K2, 'ref_rel', sRB, 'ref_plain', sPB, 'features', {KEYL}, ...
           'data_ratio', dataR, 'ref_ratio', rB, 'levers', T);
fw = fopen(fullfile(resS,'levers_summary.json'),'w'); fwrite(fw, jsonencode(J, 'PrettyPrint', true)); fclose(fw);
type(fullfile(resS,'levers.txt'));

%% ---- figure: tornado of d(score_rel) + feature heat map -------------------------
navy = [0 39 76]/255; red = [0.75 0.2 0.12];
lab = cellfun(@(l) sprintf('%s (%s)', l, T.(l).tier), levs, 'uni', 0);
f = figure('Color','w','Position',[50 50 1700 760],'Visible','off');
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
nexttile; hold on; box on; grid on;
D = cell2mat(cellfun(@(l) T.(l).d_rel, levs(:), 'uni', 0));
barh(1:nLv, D(:,1), 0.4, 'FaceColor', [0.55 0.65 0.8]); barh((1:nLv)+0.4, D(:,2), 0.4, 'FaceColor', navy);
set(gca, 'YTick', (1:nLv)+0.2, 'YTickLabel', lab, 'YDir', 'reverse', 'FontSize', 13, 'TickLabelInterpreter', 'none');
xline(0, 'k-'); xlabel('\Delta relevance-weighted 2 mM \chi^2 (negative = better)');
legend({'lever \times0.9 at 2 mM','lever \times1.1 at 2 mM'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
title(sprintf('One lever at a time, on top of R_2 \\times%.2f', G_K2), 'FontSize', 15);
nexttile;
imagesc(-Dm); colormap(gca, iRedBlue()); cl = max(abs(Dm(:))); if cl > 0; caxis([-cl cl]); end
set(gca, 'XTick', 1:numel(KEYL), 'XTickLabel', KEYL, 'YTick', 1:nLv, 'YTickLabel', ...
    cellfun(@(l) sprintf('%s x%.1f', l, T.(l).best_f), levs, 'uni', 0), 'FontSize', 13, 'TickLabelInterpreter', 'none');
xtickangle(35); cb = colorbar; cb.Label.String = 'blue = moves toward data';
title('Per feature, at each lever''s better direction', 'FontSize', 15);
exportgraphics(f, fullfile(resS,'levers_tornado.png'), 'Resolution', 200); close(f);

function r = iRatio(fm, fm8, nm, fv, fv8)
    if strcmp(nm, 'FV2'); r = fv/fv8; return; end
    a = iM(fm, nm); b = iM(fm8, nm); r = a/b;
end
function v = iM(st, nm)
    v = NaN;
    if isfield(st, nm); x = double(st.(nm)); x = x(isfinite(x)); if ~isempty(x); v = mean(x); end; end
end
function c = iRedBlue()
    % red (worse) -> white -> navy (better)
    n = 64;
    c = [linspace(0.75,1,n)' linspace(0.2,1,n)' linspace(0.12,1,n)'; ...
         linspace(1,0,n)'   linspace(1,0.15,n)' linspace(1,0.3,n)'];
end
