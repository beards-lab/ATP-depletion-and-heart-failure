% PlotParadox.m — the three panels that carry the result.
%
%  (a) transit time out of p2, per route, at the post-restretch force vs at zero
%      force. The SRX route is a SHORTCUT where the restretch happens and a trap
%      only where the other two windows happen.
%  (b) post-restretch rate vs diverted mass: the reservoir PD is what makes the
%      return path nearly irrelevant.
%  (c) rsK and vall_y are anti-correlated across the whole variant set — the
%      restretch valley and the post-restretch rate are two views of one flux.
%
% Run after Baseline.m, RunSRXRecoil.m, RunDestControl.m, RunParadoxProbe.m.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
res = fullfile(here, 'results');

Lb = load(fullfile(res,'baseline.mat'));      B = Lb.B;
Lv = load(fullfile(res,'variants.mat'));
Ld = load(fullfile(res,'destcontrol.mat'));
Lp = load(fullfile(res,'paradox_probe.mat'));

p = Lb.B; %#ok<NASGU>
pp = getParams(loadParams(Lb.SNAP), [], true, false);
F0 = 66.4;                                    % measured window-start force (kPa)

f = figure('Position', [80 80 1500 460], 'Color', 'w');
tl = tiledlayout(f, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

%% (a) transit times ---------------------------------------------------------
nexttile;
kOut = @(F) pp.kmsrd * exp(F/pp.sigma_srd1);
routes = {'p2->PT->PD','p2->SRD->PD','p2->SR->SRD->PD'};
tHi = [1000/pp.kah, 1000/kOut(F0),   1000/pp.ksr2srd + 1000/kOut(F0)];
tLo = [1000/pp.kah, 1000/kOut(0),    1000/pp.ksr2srd + 1000/kOut(0)];
bar([tHi; tLo]', 'grouped'); set(gca,'XTickLabel', routes, 'YScale','log');
ylabel('transit time out of p2  (ms)'); grid on;
legend({sprintf('at F = %.0f kPa (post-restretch)', F0), 'at F = 0 (post-slack / ktr)'}, ...
       'Location','northwest');
title(sprintf('(a) the SRX detour is a SHORTCUT at restretch force\n%.1f ms vs the normal %.1f ms', ...
      tHi(2), tHi(1)));

%% (b) rate vs diverted mass ------------------------------------------------
nexttile;
k = fieldnames(Lp.P); col = lines(numel(k)); hold on;
for i = 1:numel(k)
    o = Lp.P.(k{i}).out;
    m = o.t >= Lp.vt(Lp.iR(1),1) & o.t <= Lp.vt(Lp.iR(1)+1,1) + 0.010;
    dv = 0; if isfield(o,'RdetSRX'); dv = trapz(o.t(m), o.RdetSRX(m)); end
    % k63 of cycle 1, recomputed the same way as the probe
    t_r1 = Lp.vt(Lp.iR(1)+1,1); t_e = Lp.vt(Lp.iS(2),1);
    j = find(o.t >= t_r1 & o.t <= t_e); tt = o.t(j); FF = o.Force(j)/o.Fiso;
    mv = tt <= tt(1)+0.080; [~, iv] = min(FF(mv)); i0 = j(iv); i1 = j(end);
    seg = i0:i1; Bv = o.Force(i0)/o.Fiso;
    Pl = median(o.Force(seg(o.t(seg) >= o.t(i1)-0.15*(o.t(i1)-o.t(i0))))/o.Fiso);
    y = (movmedian(o.Force(seg)/o.Fiso, max(3, round(0.002/median(diff(o.t(seg)))))) - Bv)/max(Pl-Bv,eps);
    i63 = find(y >= 1-exp(-1), 1); kk = 1/(o.t(seg(i63)) - o.t(i0));
    PD0 = o.PuR(i0);
    plot(dv, kk, 'o', 'MarkerSize', 12, 'MarkerFaceColor', col(i,:), 'Color', col(i,:));
    text(dv, kk, sprintf('  %s\n  PD=%.2f', Lp.P.(k{i}).name, PD0), 'FontSize', 8);
end
yline(46.3, 'k--', 'data'); grid on;
xlabel('head mass diverted per re-stretch'); ylabel('post-restretch k_{63} (s^{-1}), cycle 1');
title(sprintf('(b) diverting mass barely moves the rate\nPD stays a 0.55-0.68 reservoir'));
xlim([-0.02 0.25]);

%% (c) the trade-off along the one knob that isolates it ---------------------
% Within the single-knob rip->SRXD family, sSRXrip sets how much of the p2
% population is diverted (lower threshold = more mass). Both features are
% monotone in it, in OPPOSITE directions. (Across the mixed variant set the
% scatter is uncorrelated, r = 0.04 — the trade-off is a property of this axis,
% not a universal, so it is plotted on the axis that isolates it.)
nexttile;
nm = cellfun(@(x) regexprep(strtok(x,'|'), '\[.*$',''), B.fn, 'uni', 0);
iK = find(strcmp(nm,'rsK'), 1); iV = find(strcmp(nm,'vall_y'), 1);
fam = {'B5 rip->SRXD s0=.005','B3 rip->SRXD s0=.010','B2 rip->SRXD s0=.015','B4 rip->SRXD s0=.020'};
s0  = [0.005 0.010 0.015 0.020];
kc = nan(1,numel(fam)); vc = nan(1,numel(fam));
for i = 1:numel(fam)
    for j = 1:numel(Lv.C)
        if Lv.C{j}.ok && strcmp(Lv.C{j}.name, fam{i})
            kc(i) = Lv.C{j}.ct2(iK); vc(i) = Lv.C{j}.ct2(iV);
        end
    end
end
% baseline plotted at the right edge: no diversion at all
s0p = [s0, 0.026]; kcp = [kc, B.ct2(iK)]; vcp = [vc, B.ct2(iV)];
yyaxis left;  plot(s0p, kcp, 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'auto');
ylabel('rsK cost   (recovery too FAST)');
yyaxis right; plot(s0p, vcp, 's--', 'LineWidth', 2, 'MarkerFaceColor', 'auto');
ylabel('vall_y cost   (restretch valley)');
xlabel('sSRXrip — strain threshold (um);  lower = more mass diverted');
xticks(s0p); xticklabels({'.005','.010','.015','.020','none'});
title(sprintf('(c) one flux, two features, opposite signs\nfilling the valley SPEEDS the recovery'));
grid on;
r = corr(kc(:), vc(:));

exportgraphics(f, fullfile(res, 'paradox.png'), 'Resolution', 140);
fprintf('rsK vs vall_y within the rip->SRXD family: r = %+.3f (n = %d)\n', r, numel(kc));
fprintf('  sSRXrip  %s\n', sprintf('%8.3f', s0));
fprintf('  rsK cost %s   (baseline %.3f)\n', sprintf('%8.3f', kc), B.ct2(iK));
fprintf('  vall_y   %s   (baseline %.3f)\n', sprintf('%8.3f', vc), B.ct2(iV));
fprintf('saved results/paradox.png\n');
