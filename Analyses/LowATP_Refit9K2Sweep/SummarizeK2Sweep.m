% SummarizeK2Sweep.m -- table, JSON and figures from results/k2sweep.mat.
% Called at the end of RunRefit9K2Sweep; can be re-run on its own:
%   matlab -batch "run('Analyses/LowATP_Refit9K2Sweep/SummarizeK2Sweep.m')"
rootS = 'C:\home\git\ATP-depletion-and-heart-failure';
resS  = fullfile(rootS,'Analyses','LowATP_Refit9K2Sweep','results');
if ~exist('RK2','var'); load(fullfile(resS,'k2sweep.mat')); end %#ok<LOAD>

gS  = [RK2.g]; nG = numel(gS);
sR  = [RK2.score_rel]; sP = [RK2.score_plain];
[~, ibS] = min(sR);

% key ATP-carrying features (+ force shape as a control that should NOT move)
KEY = { 'steady','Steady force (kPa)'
        'peak1_y','Re-stretch peak (kPa)'
        'ktrfit','Slack k_{tr} (s^{-1})'
        'ktrP_k80','k_{tr} protocol k80 (s^{-1})'
        'rsK80','Post-restretch rate rsK80 (s^{-1})'
        't0','Redevelopment delay t0 (s)'
        'rsTfitE','Recovery time rsTfitE (s)'
        'vall_tq','Valley time (s)' };
PT8S = pooledTarget8mM(); PT2S = pooledTarget2mM(); KT8S = pooledKtrTarget(8); KT2S = pooledKtrTarget(2);
mv = @(st, nm) iMean(st, nm);

J = struct('g', gS, 'p', [RK2.p], 'score_rel', sR, 'score_plain', sP, 'best_g', gS(ibS), ...
           'core', CORE, 'fv2_model', [RK2.fv2], 'fv2_data2', FV2K2.m2, 'fv2_data2_sd', FV2K2.s2, ...
           'fv2_data8', FV2K2.m8, 'fv2_model8', fv8K2, 'cost8', sum(c8K2));
J.ratios = struct(); J.values = struct();
fid = fopen(fullfile(resS,'k2sweep.txt'),'w');
fprintf(fid, 'k2 sweep on %s (only A2 detachment gated; g = R2(2 mM)/R2(8 mM))\n', CORE);
fprintf(fid, '8 mM cost on this objective (PS_ terms dropped): %.4f\n\n', sum(c8K2));
fprintf(fid, '%6s %7s %12s %12s %8s\n', 'g', 'p', 'score_rel', 'score_plain', 'FV2');
for k = 1:nG
    fprintf(fid, '%6.2f %7.3f %12.4f %12.4f %8.3f\n', gS(k), RK2(k).p, sR(k), sP(k), RK2(k).fv2);
end
fprintf(fid, '\nbest (relevance-weighted): g = %.2f\n\n', gS(ibS));
fprintf(fid, '%-12s %9s %9s %9s | %8s %8s | model 2/8 ratio per g\n', 'feature','data8','data2','model8','r data','sd');
for i = 1:size(KEY,1)
    nm = KEY{i,1};
    d8 = mv(fd8K2, nm); m8 = mv(fm8K2, nm); d2 = mv(RK2(1).fd, nm);
    m2 = arrayfun(@(R) mv(R.fm, nm), RK2);
    [rD, sdD] = iDataRatio(nm, PT8S, KT8S, KT2S);
    if ~isfinite(rD) && isfinite(d8) && d8 ~= 0; rD = d2/d8; end
    J.values.(nm) = struct('data8', d8, 'data2', d2, 'model8', m8, 'model2', m2);
    J.ratios.(nm) = struct('model', m2/m8, 'data', rD, 'data_sd', sdD, 'label', KEY{i,2});
    fprintf(fid, '%-12s %9.4g %9.4g %9.4g | %8.3f %8.3f | %s\n', nm, d8, d2, m8, rD, sdD, ...
            strtrim(sprintf('%6.3f ', m2/m8)));
end
J.ratios.FV2 = struct('model', [RK2.fv2]/fv8K2, 'data', FV2K2.r, 'data_sd', FV2K2.sr, 'label', 'FV2 P/P_0 at 2 ML/s');
fprintf(fid, '%-12s %9.4g %9.4g %9.4g | %8.3f %8.3f | %s\n', 'FV2', FV2K2.m8, FV2K2.m2, fv8K2, ...
        FV2K2.r, FV2K2.sr, strtrim(sprintf('%6.3f ', [RK2.fv2]/fv8K2)));

% largest remaining 2 mM terms at the best g
[cs, is] = sort(RK2(ibS).cost(:)', 'descend');
nmT = cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), fn2K2, 'uni', 0);
fprintf(fid, '\nlargest 2 mM terms at g = %.2f (plain chi-2, relevance a):\n', gS(ibS));
J.topTerms = {};
for k = 1:min(8, numel(cs))
    fprintf(fid, '  %-22s %8.4f   a = %.2f\n', nmT{is(k)}, cs(k), aK2(is(k)));
    J.topTerms{end+1} = sprintf('%s %.3f', nmT{is(k)}, cs(k));
end
fclose(fid);
J.messages = {};
fw = fopen(fullfile(resS,'k2sweep_summary.json'),'w'); fwrite(fw, jsonencode(J, 'PrettyPrint', true)); fclose(fw);
type(fullfile(resS,'k2sweep.txt'));

%% ---- figures ------------------------------------------------------------
navy = [0 39 76]/255; red = [0.75 0.2 0.12]; grey = [0.5 0.5 0.5];
f1 = figure('Color','w','Position',[50 50 1100 620],'Visible','off');
yyaxis left;  plot(gS, sR, '-o', 'Color', navy, 'LineWidth', 2.5, 'MarkerFaceColor', navy); ylabel('relevance-weighted 2 mM \chi^2');
yyaxis right; plot(gS, sP, '--s', 'Color', red, 'LineWidth', 2, 'MarkerFaceColor', red); ylabel('plain 2 mM \chi^2 (all terms)');
ax = gca; ax.YAxis(1).Color = navy; ax.YAxis(2).Color = red; set(ax, 'XDir', 'reverse', 'FontSize', 16, 'LineWidth', 1);
xlabel('g = R_2(2 mM) / R_2(8 mM)   (1 = no ATP effect)'); grid on; box on;
xline(gS(ibS), ':', sprintf('best g = %.2f', gS(ibS)), 'FontSize', 14, 'LabelVerticalAlignment', 'bottom');
title(sprintf('Low ATP on %s: only A2 detachment feels ATP', strrep(CORE,'_','\_')), 'FontSize', 18);
exportgraphics(f1, fullfile(resS,'k2sweep_cost.png'), 'Resolution', 200); close(f1);

f2 = figure('Color','w','Position',[50 50 1600 820],'Visible','off');
tl = tiledlayout(f2, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:size(KEY,1)
    nm = KEY{i,1}; V = J.values.(nm); nexttile; hold on; box on; grid on;
    sd2 = iDataSd(nm, PT2S, KT2S);
    if isfinite(V.data2)
        if isfinite(sd2); patch([min(gS) max(gS) max(gS) min(gS)], V.data2*[1-sd2 1-sd2 1+sd2 1+sd2], red, 'FaceAlpha', 0.15, 'EdgeColor', 'none'); end
        yline(V.data2, '-', 'Color', red, 'LineWidth', 1.5);
    end
    if isfinite(V.data8); yline(V.data8, '-', 'Color', navy, 'LineWidth', 1.5); end
    plot(gS, V.model2, '-o', 'Color', red, 'LineWidth', 2.5, 'MarkerFaceColor', red);
    plot(1, V.model8, 'o', 'Color', navy, 'MarkerSize', 10, 'LineWidth', 2.5, 'MarkerFaceColor', navy);
    set(gca, 'XDir', 'reverse', 'FontSize', 13); title(KEY{i,2}, 'FontSize', 14);
    if i > 4; xlabel('g'); end
end
lg = legend({'2 mM data \pm sd','2 mM data','8 mM data','2 mM model','8 mM model'}, 'Orientation', 'horizontal', 'FontSize', 13);
lg.Layout.Tile = 'south';
title(tl, 'Key features vs the A2-detachment factor g (lines: pooled data; markers: model)', 'FontSize', 17);
exportgraphics(f2, fullfile(resS,'k2sweep_features.png'), 'Resolution', 200); close(f2);

% trace figure: the POOLED data band in kPa (what the objective scores), not a single prep
PlotPooledTrace;

%% ---- helpers --------------------------------------------------------------
function v = iMean(st, nm)
    v = NaN;
    if isstruct(st) && isfield(st, nm)
        x = double(st.(nm)); x = x(isfinite(x));
        if ~isempty(x); v = mean(x); end
    end
end
function [r, sd] = iDataRatio(nm, PT8, KT8, KT2)
    r = NaN; sd = NaN;
    j = find(strcmp(PT8.name, nm), 1);
    if ~isempty(j); r = PT8.r(j); sd = PT8.sr(j); return; end
    j = find(strcmp(KT8.name, nm), 1);
    if ~isempty(j)
        r = KT2.val(j)/KT8.val(j); sd = r*sqrt(KT2.s(j)^2 + KT8.s(j)^2);
    end
end
function sd = iDataSd(nm, PT2, KT2)
    sd = NaN;
    j = find(strcmp(PT2.name, nm), 1); if ~isempty(j); sd = PT2.s(j); return; end
    j = find(strcmp(KT2.name, nm), 1); if ~isempty(j); sd = KT2.s(j); end
end
