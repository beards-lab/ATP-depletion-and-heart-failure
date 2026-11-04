% PlotMomentsSummary.m
% Figure for conclusions.md: PDE vs Gaussian-moments force on four protocols,
% and the p2(s) distribution where the closure works (mid-ejection) and
% where it does not (just after the fast ktr restretch).

here = fileparts(mfilename('fullpath'));
params0 = momentsSetup('fastgate_best');
A = load(fullfile(here, 'results', 'aposteriori.mat'));
names = {'ktr', [0.98 1.3]; 'slack', [0.98 1.4]; 'beat', [0.95 2.3]; 'fv4', [0.99 1.08]};
fig = figure('Color', 'w', 'Position', [60 60 1250 720]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:4
    S = load(fullfile(here, 'results', ['ref_' names{k,1} '.mat']));
    o = A.R.(matlab.lang.makeValidName(names{k,1})).gauss3.o;
    ax = nexttile(tl);
    plot(ax, S.ref.t, S.ref.Force, 'k-', 'LineWidth', 1.6); hold(ax, 'on');
    plot(ax, o.t, o.Force, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
    xlim(ax, names{k,2}); grid(ax, 'on'); xlabel(ax, 't (s)'); ylabel(ax, 'Force (kPa)');
    r = A.R.(matlab.lang.makeValidName(names{k,1})).gauss3;
    title(ax, sprintf('%s: rms %.1f%%, max %.1f%% of F_{iso}', names{k,1}, r.rms, r.max));
    if k == 1, legend(ax, 'PDE (60 bins)', 'moments (3/state, Gaussian)', 'Location', 'southeast'); end
end
% distribution snapshots
snap = {'beat', 1.10, 'mid-ejection'; 'ktr', 1.0305, 'ktr, 0.5 ms after restretch'};
for k = 1:2
    S = load(fullfile(here, 'results', ['ref_' snap{k,1} '.mat'])); ref = S.ref;
    i = find(ref.t >= snap{k,2}, 1);
    s = params0.s(:) - (-(ref.SL(i) - ref.LSE(i)) + ref.LXBPivot(i))/2;
    p2 = max(0, ref.p2(i,:)');
    m0 = sum(p2)*params0.dS; mu = sum(s.*p2)*params0.dS/m0; sd = sqrt(sum((s-mu).^2.*p2)*params0.dS/m0);
    ax = nexttile(tl);
    bar(ax, 1e3*s, p2*params0.dS, 1, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none'); hold(ax, 'on');
    xs = linspace(mu - 4*sd, mu + 4*sd, 400);
    plot(ax, 1e3*xs, m0*params0.dS*exp(-(xs-mu).^2/(2*sd^2))/(sqrt(2*pi)*sd), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
    xlim(ax, 1e3*[mu - 5*sd, mu + 5*sd]); grid(ax, 'on');
    xlabel(ax, 'strain s (nm)'); ylabel(ax, 'p2 mass per bin');
    title(ax, sprintf('p2(s), %s (sd %.1f nm)', snap{k,3}, 1e3*sd));
end
exportgraphics(fig, fullfile(here, 'results', 'moments_vs_pde.png'), 'Resolution', 130);
