% ExportFeaturePanel2mM.m — export the 2 mM plotFeatures panel BIG.
% plotFeatures.m:53 does `figure(80085); clf;` -- it makes its OWN figure at the
% default size, so pre-creating a large figure does nothing. The fix is to resize
% figure 80085 AFTER the call and let the 'flow' tiledlayout reflow into it.
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root)); LoadData;
outdir=fullfile(root,'Analyses','LowATP_AbsoluteFit','results','panels_2mM');
if ~exist(outdir,'dir'); mkdir(outdir); end

q = getParams(loadParams('lowATP_powergate2_opt'), [], true, false);
q.MgATP=2; q.MgADP=0; q.Pi=0.8;
q.RunForceVelocity=false; q.RunForceVelocityTime=false; q.RunKtr=false;
q.RunStairs=false; q.RunSlackPassive=false; q.RunSlack=true;
q.EvalFeatures=true; q.BreakOnODEUnstable=false;
q.RunSlackSegments='All'; q.MaxRunTime=400;
q.PlotEachSeparately=0; q.PlotFeatureFitting=0;
q.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
[~,~,fm,fd] = runSlackExperiment(q);

fn = q.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
fnS = fn(~(startsWith(fn,'FV')|startsWith(fn,'PS_')));
fprintf('plotting %d feature tiles\n', numel(fnS));

close all force;
plotFeatures(fd, fm, [], fnS);          % creates/uses figure 80085

f = figure(80085);
set(f,'Units','pixels','Position',[10 10 2560 1440],'Color','w');
% enlarge text now that the tiles are big
set(findall(f,'-property','FontSize'),'FontSize',13);
tl = findobj(f,'Type','tiledlayout');
if ~isempty(tl); tl.TileSpacing='compact'; tl.Padding='compact'; end
sgtitle(f, '2 mM ATP — model (x) vs data (o), all scored slack features', ...
        'FontSize',20,'FontWeight','bold');
drawnow; pause(1); drawnow;                       % let the flow layout settle

dst = fullfile(outdir,'2mM_plotFeatures_panel.png');
exportgraphics(f, dst, 'Resolution', 150);
i = imfinfo(dst);
fprintf('saved %s  (%d x %d px)\n', dst, i.Width, i.Height);
disp('PANEL EXPORT DONE');
