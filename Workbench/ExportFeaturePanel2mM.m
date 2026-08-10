% ExportFeaturePanel2mM.m — 2 mM feature panels, FULL params0.fn (FV + slack + PS_),
% with GHOST comparisons.
%
% Two gotchas this script exists to handle:
%  1. plotFeatures.m:53 does `figure(80085); clf;` -- it makes its OWN figure at the
%     default size, so pre-creating a large figure does nothing. Resize 80085 AFTER
%     the call and let the 'flow' tiledlayout reflow into it.
%  2. RunBakersExp is a SCRIPT and clobbers common workspace names (F, E, i, ...).
%     Everything persistent here is prefixed pnl_ so it survives the call.
%
% Panel A: ghost = previous best low-ATP (Michaelis gates)  vs  sim = current best
% Panel B: ghost = current best                             vs  sim = current + passive x1.04
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root)); LoadData;
pnl_outdir=fullfile(root,'Analyses','LowATP_AbsoluteFit','results','panels_2mM');
if ~exist(pnl_outdir,'dir'); mkdir(pnl_outdir); end

pnl_PREV = 'lowATP_only_opt';        % previous best: frozen 8 mM + Michaelis gates
pnl_CUR  = 'lowATP_powergate2_opt';  % current best: 21 ATP power gates

pnl_tmp  = getParams(loadParams(pnl_CUR), [], true, false);
pnl_FN   = pnl_tmp.fn; [~,pnl_iu]=unique(pnl_FN,'stable'); pnl_FN=pnl_FN(pnl_iu);
pnl_KPAS = pnl_tmp.k_pas;
fprintf('using the snapshot fn: %d terms (FV and PS_ included)\n', numel(pnl_FN));

pnl_cfg = { pnl_PREV, 1.00; pnl_CUR, 1.00; pnl_CUR, 1.04 };
pnl_feats = cell(1,3); pnl_data = []; pnl_cost = nan(1,3);

for pnl_k = 1:3
    fprintf('\n===== %s (k_pas x%.2f) =====\n', pnl_cfg{pnl_k,1}, pnl_cfg{pnl_k,2});
    params0 = getParams(loadParams(pnl_cfg{pnl_k,1}), [], true, false);
    params0.MgATP=2; params0.MgADP=0; params0.Pi=0.8;
    params0.k_pas = pnl_KPAS*pnl_cfg{pnl_k,2};
    params0.RunForceVelocity=true; params0.RunForceVelocityTime=false;
    params0.RunSlack=true; params0.RunSlackPassive=true;
    params0.RunKtr=false; params0.RunStairs=false;
    params0.EvalFeatures=true; params0.BreakOnODEUnstable=false; params0.OptimizeOn='Feats';
    params0.RunSlackSegments='AllPar'; params0.MaxRunTime=400;
    params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0;
    params0.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
    params0.fn = pnl_FN;
    RunBakersExp;                       % clobbers plain names -- pnl_* survive
    pnl_feats{pnl_k} = features_model;
    if isempty(pnl_data); pnl_data = features_data; end
    pnl_cost(pnl_k) = sum(evalFeatureCost(features_data, features_model, pnl_FN, 2));
    fprintf('  cost = %.4f\n', pnl_cost(pnl_k));
end

fprintf('\ncosts: prev %.4f | current %.4f | current+passive %.4f\n', pnl_cost);

bigPanel(pnl_data, pnl_feats{2}, pnl_feats{1}, pnl_FN, pnl_outdir, ...
    sprintf('2 mM ATP — data (o), current best 21 power gates (x) [%.2f], GHOST = previous best Michaelis gates (+) [%.2f]', ...
            pnl_cost(2), pnl_cost(1)), ...
    '2mM_plotFeatures_panel.png');

bigPanel(pnl_data, pnl_feats{3}, pnl_feats{2}, pnl_FN, pnl_outdir, ...
    sprintf('2 mM ATP — data (o), current + passive \\times1.04 (x) [%.2f], GHOST = current best (+) [%.2f]', ...
            pnl_cost(3), pnl_cost(2)), ...
    '2mM_plotFeatures_passive.png');
disp('PANEL EXPORT DONE');

% ================= local functions =================
function bigPanel(dat, sim, ghost, FN, outdir, ttl, fname)
    close all force;
    plotFeatures(dat, sim, ghost, FN);          % -> figure 80085
    f = figure(80085);
    set(f,'Units','pixels','Position',[10 10 2600 1460],'Color','w');
    set(findall(f,'-property','FontSize'),'FontSize',13);
    sgtitle(f, ttl, 'FontSize',16,'FontWeight','bold');
    drawnow; pause(1); drawnow;
    dst = fullfile(outdir, fname);
    exportgraphics(f, dst, 'Resolution', 150);
    ii = imfinfo(dst);
    fprintf('saved %s  (%d x %d px)\n', fname, ii.Width, ii.Height);
end
