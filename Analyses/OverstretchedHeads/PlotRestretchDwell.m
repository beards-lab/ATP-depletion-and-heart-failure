%% PlotRestretchDwell.m -- the re-stretch, seen by the heads that carry it.
%
% For each configuration in CFG (name, override struct OR snapshot name), one serial slack run;
% for re-stretch SEG, aligned at the lengthening onset:
%   row 1  force, model vs pooled-prep data trace (the 03/27 protocol file)
%   row 2  attached mass beyond +15 nm, p1 (weak) and p2 (strong)
%   row 3  share of cross-bridge force carried by heads beyond +15 nm
% Writes results/restretch_dwell_<TAG>.png.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
if ~exist('CFG', 'var'), CFG = {'as fitted', struct()}; end
if ~exist('SEG', 'var'), SEG = [1 4]; end
if ~exist('TAG', 'var'), TAG = PARAMSET; end
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..')));
res = fullfile(here, 'results');

C = struct();
for c = 1:size(CFG, 1)
    ov = CFG{c, 2};
    if ischar(ov)                         % a snapshot name instead of overrides
        P = getParams(loadParams(ov), [], true, false);
    else
        P = getParams(loadParams(PARAMSET), [], true, false);
        f = fieldnames(ov); for k = 1:numel(f), P.(f{k}) = ov.(f{k}); end
    end
    P.RunSlackSegments = 'All'; P.MaxRunTime = 900; P.PlotEachSeparately = 0; P.PlotFeatureFitting = 0;
    P.ghostSave = ''; P.ghostLoad = ''; P.SaveBest = 0; P.mods = {}; P.g = ones(1,0);
    P = getParams(P, P.g, true, false);
    evalc('[~, o] = runSlackExperiment(P);');
    ss = P.ss; n = numel(o.t); [m1, m2, fsh] = deal(zeros(n, 1));
    for i = 1:n
        s  = P.s(:) - (-(o.SL(i) - o.LSE(i)) + o.LXBPivot(i))/2;
        p1 = max(0, o.PU(i, 1:ss)'); p2 = max(0, o.PU(i, ss+1:2*ss)');
        s1 = s + P.dr1; x = s + P.dr; e = sign(x).*abs(x).^P.estiff;
        fb = (P.kstiff1*(s1 >= 0) + P.kstiff1_n*(s1 < 0)).*s1.*p1 + (P.kstiff2*(s >= -P.dr) + P.kstiff2_n*(s < -P.dr)).*e.*p2;
        hi = s > 0.015;
        m1(i) = P.dS*sum(p1(hi)); m2(i) = P.dS*sum(p2(hi)); fsh(i) = sum(fb(hi)) / max(abs(sum(fb)), eps);
    end
    C(c).name = CFG{c,1}; C(c).t = o.t(:); C(c).SL = o.SL(:); C(c).F = o.Force(:);
    C(c).m1 = m1; C(c).m2 = m2; C(c).fsh = fsh; C(c).data = o.datatable;
end

dSL = [0; diff(C(1).SL) ./ max(diff(C(1).t), eps)];
on = find(dSL(2:end) > 0.5 & dSL(1:end-1) <= 0.5) + 1; on = on(C(1).t(on) > 70);
fig = figure('Color', 'w', 'Position', [80 80 520*numel(SEG) 760], 'Visible', 'off');
tl = tiledlayout(3, numel(SEG), 'TileSpacing', 'compact', 'Padding', 'compact');
col = lines(max(3, numel(C)));
for j = 1:numel(SEG)
    t0 = C(1).t(on(SEG(j))); win = [-0.005 0.045];
    D = C(1).data; wd = D(:,1) >= t0 + win(1) & D(:,1) <= t0 + win(2);
    ax = nexttile(tl, j); hold(ax, 'on');
    plot(ax, 1e3*(D(wd,1) - t0), D(wd,3), 'Color', [0.6 0.6 0.6], 'DisplayName', 'data 03/27');
    for c = 1:numel(C)
        w = C(c).t >= t0 + win(1) & C(c).t <= t0 + win(2);
        plot(ax, 1e3*(C(c).t(w) - t0), C(c).F(w), 'Color', col(c,:), 'LineWidth', 1.4, 'DisplayName', C(c).name);
    end
    title(ax, sprintf('re-stretch %d', SEG(j))); ylabel(ax, 'force (kPa)'); if j == 1, legend(ax, 'Location', 'southeast', 'Box', 'off'); end
    ax = nexttile(tl, numel(SEG) + j); hold(ax, 'on');
    for c = 1:numel(C)
        w = C(c).t >= t0 + win(1) & C(c).t <= t0 + win(2); tt = 1e3*(C(c).t(w) - t0);
        plot(ax, tt, C(c).m1(w), '-', 'Color', col(c,:), 'LineWidth', 1.4, 'DisplayName', [C(c).name ' p1']);
        plot(ax, tt, C(c).m2(w), ':', 'Color', col(c,:), 'LineWidth', 1.4, 'DisplayName', [C(c).name ' p2']);
    end
    ylabel(ax, 'attached beyond +15 nm'); if j == 1, legend(ax, 'Location', 'northeast', 'Box', 'off'); end
    ax = nexttile(tl, 2*numel(SEG) + j); hold(ax, 'on');
    for c = 1:numel(C)
        w = C(c).t >= t0 + win(1) & C(c).t <= t0 + win(2);
        plot(ax, 1e3*(C(c).t(w) - t0), 100*C(c).fsh(w), 'Color', col(c,:), 'LineWidth', 1.4);
    end
    ylabel(ax, 'XB force from >15 nm (%)'); xlabel(ax, 'time from re-stretch onset (ms)'); ylim(ax, [0 100]);
end
exportgraphics(fig, fullfile(res, ['restretch_dwell_' TAG '.png']), 'Resolution', 130);
save(fullfile(res, ['restretch_dwell_' TAG '.mat']), 'C', 'CFG', 'SEG');
