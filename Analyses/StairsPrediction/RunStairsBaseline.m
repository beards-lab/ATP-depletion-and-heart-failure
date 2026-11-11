% RunStairsBaseline.m — how well does the CURRENT best fit PREDICT the staircase?
%
% The staircase (protocol_03_27_2026_velocitytable_stairs) was never in the
% objective: params/spreadw_opt.m carries RunStairs = 0. So this is a genuine
% out-of-sample prediction test, and it doubles as the yardstick for the
% estiff experiment that follows.
%
% Outputs -> Analyses/StairsPrediction/results/baseline.mat + baseline.png

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resdir = fullfile(here, 'results');
if ~isfolder(resdir); mkdir(resdir); end

if ~exist('SNAP','var'); SNAP = 'spreadw_opt'; end

params0 = getParams(loadParams(SNAP), [], true, false);
params0.RunStairs = 1;
params0.PlotEachSeparately = false;
params0.ghostSave = '';  params0.ghostLoad = '';

fprintf('=== stairs prediction from %s ===\n', SNAP);
fprintf('estiff = %.6g   kstiff1=%.4g kstiff2=%.4g  dr=%.4g\n', ...
    params0.estiff, params0.kstiff1, params0.kstiff2, params0.dr);

tic; [E_stairs, out] = runStairsExperiment(params0); tRun = toc;
fprintf('E_stairs (MSE*10) = %.5f   [%.1f s, %d ODE points]\n', E_stairs, tRun, numel(out.t));

% ---------------- data & model on one instrument -------------------------
d  = load('data/protocol_03_27_2026_velocitytable_stairs.mat');
dt = d.datatable; vt = d.velocitytable;

% model force normalised exactly as runStairsExperiment does
i_move = find(vt(:,2) ~= 0, 1);
Fm_ss  = interp1(out.t, out.Force, vt(i_move,1), 'linear', 'extrap');
Fn     = out.Force(:) / Fm_ss;

fD = extractStairsAttributes(dt(:,1), dt(:,3), vt, 'st_', struct(), true);
fM = extractStairsAttributes(out.t(:), Fn,     vt, 'st_', struct(), false);

% ---------------- noise floor of the raw least-squares -------------------
% Residual scatter inside each hold, after removing a local linear trend.
res = [];
for i = 1:size(vt,1)-1
    if vt(i,2) ~= 0; continue; end
    m = dt(:,1) >= vt(i,1)+0.005 & dt(:,1) <= vt(i+1,1)-0.002;
    if nnz(m) < 20; continue; end
    p = polyfit(dt(m,1), dt(m,3), 1);
    res = [res; dt(m,3) - polyval(p, dt(m,1))]; %#ok<AGROW>
end
sigma = std(res);
fprintf('\n-- least-squares diagnostics --\n');
fprintf('per-sample noise sigma      = %.4f F0  (n=%d in-hold samples)\n', sigma, numel(res));
fprintf('irreducible MSE*10 floor    = %.5f   (a perfect model cannot go below)\n', 10*sigma^2);
fprintf('E_stairs / floor            = %.2f x\n', E_stairs/(10*sigma^2));

% ---------------- feature comparison -------------------------------------
K   = numel(fD.st_peak);
nEq = nnz(abs(fD.st_dSL - median(fD.st_dSL)) < 0.25*abs(median(fD.st_dSL)));
fprintf('\n-- per-step features (%d steps, %d equal-sized) --\n', K, nEq);
fprintf('%-4s %8s | %8s %8s | %8s %8s | %8s %8s | %8s %8s\n', ...
    'step','dSL', 'peakD','peakM', 'nadD','nadM', 'platD','platM', 'recD','recM');
for k = 1:K
    fprintf('%-4d %8.4f | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f\n', ...
        k, fD.st_dSL(k), fD.st_peak(k), fM.st_peak(k), fD.st_nadir(k), fM.st_nadir(k), ...
        fD.st_plat(k), fM.st_plat(k), fD.st_rec(k), fM.st_rec(k));
end

sc = {'st_peakFirst','st_peakMean','st_peakSlope','st_nadirMean', ...
      'st_platEnd','st_platSlope','st_platCurv','st_peakBig'};
fprintf('\n-- scalar summary --\n%-16s %10s %10s %10s\n','feature','data','model','rel.err');
for i = 1:numel(sc)
    a = fD.(sc{i}); b = fM.(sc{i});
    fprintf('%-16s %10.5f %10.5f %9.1f%%\n', sc{i}, a, b, 100*(b-a)/max(abs(a),eps));
end

% ---------------- figure --------------------------------------------------
fig = figure('Position',[80 80 1500 780],'Color','w');
tl = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

nexttile([1 2]); hold on; box on;
plot(dt(:,1), dt(:,3), '-', 'Color',[.75 .75 .75]);
plot(dt(:,1), sgolayfilt(dt(:,3),2,101), 'k-', 'LineWidth',1.4);
plot(out.t, Fn, 'b-', 'LineWidth',1.6);
xlim([dt(1,1) dt(end,1)]);
xlabel('t (s)'); ylabel('F / F_0');
legend({'data (raw)','data (smoothed)','model'},'Location','northwest');
title(sprintf('%s prediction of the staircase   |   E_{stairs}=%.4f (floor %.4f)', ...
    SNAP, E_stairs, 10*sigma^2),'Interpreter','none');

nexttile; hold on; box on;
plot(1:K, fD.st_peak, 'ks-','LineWidth',1.4,'MarkerFaceColor','k');
plot(1:K, fM.st_peak, 'bo-','LineWidth',1.4);
xlabel('step'); ylabel('peak - base  (F_0)');
title('elastic (T1-like) amplitude'); legend({'data','model'},'Location','best');

nexttile; hold on; box on;
plot(1:K, fD.st_plat, 'ks-','LineWidth',1.4,'MarkerFaceColor','k');
plot(1:K, fM.st_plat, 'bo-','LineWidth',1.4);
xlabel('step'); ylabel('plateau - F_0  (F_0)');
title('staircase envelope'); legend({'data','model'},'Location','best');

nexttile; hold on; box on;
plot(1:K, fD.st_nadir, 'ks-','LineWidth',1.4,'MarkerFaceColor','k');
plot(1:K, fM.st_nadir, 'bo-','LineWidth',1.4);
yline(0,'k:');
xlabel('step'); ylabel('nadir - base  (F_0)');
title('phase-2 undershoot'); legend({'data','model'},'Location','best');

nexttile; hold on; box on;
plot(1:K, fD.st_rec, 'ks-','LineWidth',1.4,'MarkerFaceColor','k');
plot(1:K, fM.st_rec, 'bo-','LineWidth',1.4);
xlabel('step'); ylabel('plateau - nadir  (F_0)');
title('delayed recovery'); legend({'data','model'},'Location','best');

exportgraphics(fig, fullfile(resdir,'baseline.png'), 'Resolution',140);
save(fullfile(resdir,'baseline.mat'), 'E_stairs','fD','fM','sigma','SNAP','Fn','out');
fprintf('\nsaved -> %s\n', fullfile(resdir,'baseline.png'));
