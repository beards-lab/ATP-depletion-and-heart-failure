% Export8mMPanels.m — 8 mM counterpart of Export2mMPanels.m.
% Produces, into Analyses/SummaryFigure/results/panels_8mM/:
%   8mM_plotFeatures_panel.png   all scored slack features, model vs data
%   8mM_Slack_full.png           the whole slack protocol, data vs model
%   8mM_Slack_zoom1.png          zoom on the FIRST slack segment, features marked
% and caches the traces to panels_8mM/slack8_cache.mat so the plots can be
% re-tuned without re-simulating.

root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root)); LoadData;
outdir = fullfile(root,'Analyses','SummaryFigure','results','panels_8mM');
if ~exist(outdir,'dir'); mkdir(outdir); end

BASE = 'rskR2_w025_opt';                 % the frozen high-ATP baseline
p0 = getParams(loadParams(BASE), [], true, false);
p0.EvalFeatures = true; p0.BreakOnODEUnstable = false; p0.OptimizeOn = 'Feats';
p0.RunSlackSegments = 'All'; p0.MaxRunTime = 400;
p0.PlotEachSeparately = 0; p0.PlotFeatureFitting = 0; p0.PlotFullscreen = 0;
p0.RunForceVelocity = false; p0.RunForceVelocityTime = false;
p0.RunStairs = false; p0.RunKtr = false; p0.RunSlackPassive = false;
p0.RunSlack = true;

fprintf('\n===== slack experiment, 8 mM =====\n');
[E_slack, out_slack, fm, fd] = runSlackExperiment(p0);
fprintf('E_slack = %s\n', mat2str(round(E_slack,4)));
fprintf('out_slack fields: %s\n', strjoin(fieldnames(out_slack)', ', '));

save(fullfile(outdir,'slack8_cache.mat'), 'out_slack','fm','fd','E_slack','-v7.3');

%% ---- 1. scored-feature panel -------------------------------------------
fn = p0.fn; [~,iu] = unique(fn,'stable'); fn = fn(iu);
fnS = fn(~(startsWith(fn,'FV') | startsWith(fn,'PS_')));
figure('Units','pixels','Position',[20 20 2560 1400],'Color','w');
plotFeatures(fd, fm, [], fnS);
sgtitle('8 mM ATP — model vs data, all scored slack features');
drawnow;
exportgraphics(gcf, fullfile(outdir,'8mM_plotFeatures_panel.png'), 'Resolution', 200);
fprintf('  saved 8mM_plotFeatures_panel.png\n');

disp('EXPORT 8mM PART 1 DONE');
