% PlotPooledTrace.m -- the slack ladder against the POOLED data (the thing the
% objective scores), in absolute kPa, 8 and 2 mM.
%
% Data: each preparation's trace is rundown-lifted and divided by ONE divisor,
% its own rundown-lifted 8 mM pre-slack isometric (buildPooledTarget steps 1-2),
% then multiplied by the pooled level Fref -- exactly how the pooled feature
% target is built, so the band and the scored features share one scale.
% Model: absolute force (kPa), 8 mM incumbent and 2 mM at R2 x0.54.
% Replaces the earlier single-prep (03/27, raw kPa) trace, whose prep is 15 %
% stronger than the pool and made the model plateau look low.
rootS = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(rootS); addpath(genpath(rootS));
resS = fullfile(rootS,'Analyses','LowATP_Refit9K2Sweep','results');
evalc('T = buildPooledTarget(struct());');
Fref = T.scale.mean; lift = T.scale.lift; alpha = T.scale.perPrep(:);
DAYS = {'03_27_2026','04_03_2026','04_10_2026'}; CONDS = {'8mM','2mM'}; nS = 5;
TG = (-0.02 : 2e-4 : 0.57)';
BAND = nan(nS, numel(TG), 3, 2);
for c = 1:2
    for g = 1:nS
        V = nan(numel(DAYS), numel(TG));
        for p = 1:numel(DAYS)
            S = load(fullfile('data', sprintf('protocol_%s_%s_slack.mat', DAYS{p}, CONDS{c})));
            rl = S.velocitytable(S.velocitytable(:,2) < -20, 1);
            V(p,:) = interp1(S.datatable(:,1) - rl(g), S.datatable(:,3)*lift(p,c)/alpha(p)*Fref, TG, 'linear', NaN);
        end
        BAND(g,:,1,c) = mean(V,1,'omitnan');
        BAND(g,:,2,c) = BAND(g,:,1,c) - std(V,0,1,'omitnan');
        BAND(g,:,3,c) = BAND(g,:,1,c) + std(V,0,1,'omitnan');
    end
end
% model traces (simulated on the 03/27 protocol table -> same release times)
M = load(fullfile(resS,'k2sweep.mat'), 'RK2', 'tr8K2'); ig = find(abs([M.RK2.g]-0.54) < 1e-9, 1);
TR = {M.tr8K2, M.RK2(ig).tr};
REL = cell(1,2);
for c = 1:2
    S = load(fullfile('data', sprintf('protocol_%s_%s_slack.mat', DAYS{1}, CONDS{c})));
    REL{c} = S.velocitytable(S.velocitytable(:,2) < -20, 1);
end
navy = [0 39 76]/255; red = [0.75 0.2 0.12]; CC = {navy, red};
f = figure('Color','w','Position',[30 30 1800 600],'Visible','off');
tl = tiledlayout(f, 1, nS, 'Padding', 'compact', 'TileSpacing', 'compact');
for g = 1:nS
    nexttile; hold on; box on; grid on;
    for c = 1:2
        fill(1000*[TG; flipud(TG)], [BAND(g,:,2,c)'; flipud(BAND(g,:,3,c)')], CC{c}, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(1000*TG, BAND(g,:,1,c), '-', 'Color', CC{c}, 'LineWidth', 1.6, 'DisplayName', sprintf('%s pooled data \\pm SD (n=3)', CONDS{c}));
        % from the release on: the segments run in parallel from precomputed initial
        % states, so the stored trace has no samples before each release (a gap
        % that interpolation would bridge with a straight line)
        tr = TR{c}; r = REL{c}(g); m = tr.t - r >= 0 & tr.t - r <= TG(end);
        plot(1000*(tr.t(m) - r), tr.F(m), '--', 'Color', CC{c}, 'LineWidth', 2, 'DisplayName', sprintf('%s model', CONDS{c}));
    end
    ylim([-5 130]); xlim(1000*[TG(1) TG(end)]); set(gca, 'FontSize', 12);
    xlabel('time from release (ms)'); title(sprintf('slack cycle %d', g), 'FontSize', 13);
    if g == 1; ylabel('force (kPa, pooled level)'); legend('Location', 'south', 'FontSize', 10); end
end
title(tl, sprintf('Pooled data (each prep / own 8 mM isometric \\times F_{ref} = %.1f kPa) vs model: 8 mM incumbent, 2 mM R_2\\times0.54', Fref), 'FontSize', 14);
exportgraphics(f, fullfile(resS,'k2sweep_trace.png'), 'Resolution', 200); close(f);

% plateau check: mean force over the last 50 ms before each re-stretch-free end (steady window)
w = TG >= 0.50 & TG <= 0.56;
for c = 1:2
    tr = TR{c}; pm = nan(1,nS); pd = nan(1,nS);
    for g = 1:nS
        r = REL{c}(g); m = tr.t - r >= 0.50 & tr.t - r <= 0.56;
        pm(g) = mean(tr.F(m)); pd(g) = mean(BAND(g,w,1,c), 'omitnan');
    end
    fprintf('%s late plateau: data %s | model %s kPa\n', CONDS{c}, mat2str(pd,3), mat2str(pm,3));
end
