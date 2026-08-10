% Export2mMPanels.m — export each experiment panel at 2 mM as its own high-res PNG.
% RunBakersExp is run once per experiment with a single RunXXX flag set, as suggested.
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root)); LoadData;
outdir=fullfile(root,'Analyses','LowATP_AbsoluteFit','results','panels_2mM');
if ~exist(outdir,'dir'); mkdir(outdir); end

BASE = 'lowATP_powergate2_opt';        % best low-ATP parametrisation
p0 = getParams(loadParams(BASE), [], true, false);
p0.MgATP=2; p0.MgADP=0; p0.Pi=0.8;
p0.EvalFeatures=true; p0.BreakOnODEUnstable=false; p0.OptimizeOn='Feats';
p0.RunSlackSegments='All'; p0.MaxRunTime=400;
p0.PlotEachSeparately=1; p0.PlotFeatureFitting=0; p0.PlotFullscreen=0;
% 2 mM protocol tables where they exist
p0.velocitytableonfile      = 'protocol_03_27_2026_2mM_slack.mat';
p0.ktr_velocitytableonfile  = 'protocol_03_27_2026_velocitytable_ktr_2mM.mat';
p0.stairs_velocitytableonfile='protocol_03_27_2026_velocitytable_stairs.mat';

RUNS = { 'ForceVelocity','RunForceVelocity'; ...
         'Slack',        'RunSlack'; ...
         'Stairs',       'RunStairs'; ...
         'Ktr',          'RunKtr' };

for i = 1:size(RUNS,1)
    nm = RUNS{i,1}; flag = RUNS{i,2};
    fprintf('\n===== %s (2 mM) =====\n', nm);
    close all force;
    params0 = p0;
    params0.RunForceVelocity=false; params0.RunForceVelocityTime=false;
    params0.RunSlack=false; params0.RunStairs=false; params0.RunKtr=false;
    params0.RunSlackPassive=false;
    params0.(flag) = true;
    try
        RunBakersExp;
        figs = findobj('Type','figure');
        figs = flipud(figs(:));
        k = 0;
        for f = figs'
            if isempty(findobj(f,'Type','axes')); continue; end
            k = k + 1;
            fn = fullfile(outdir, sprintf('2mM_%s_%02d.png', nm, k));
            try
                % blow the figure up to a large canvas BEFORE export, otherwise the
                % axes keep their small on-screen layout and detail stays unreadable
                set(f,'Units','pixels','Position',[30 30 2400 1350],'Color','w');
                drawnow;
                exportgraphics(f, fn, 'Resolution', 200);
                fprintf('  saved %s\n', fn);
            catch ME
                fprintf('  [skip fig] %s\n', ME.message);
            end
        end
        if k==0; fprintf('  (no figure produced)\n'); end
    catch ME
        fprintf('  FAILED: %s\n', ME.message);
    end
end

% ---- plotFeatures panel (2 mM, slack features) at high res ----
fprintf('\n===== plotFeatures panel (2 mM) =====\n');
close all force;
q = p0; q.RunForceVelocity=false; q.RunForceVelocityTime=false; q.RunStairs=false;
q.RunKtr=false; q.RunSlackPassive=false; q.RunSlack=true;
q.PlotEachSeparately=0;
[~,~,fm,fd] = runSlackExperiment(q);
fn = q.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
fnS = fn(~(startsWith(fn,'FV')|startsWith(fn,'PS_')));
figure('Units','pixels','Position',[20 20 2560 1400],'Color','w');
plotFeatures(fd, fm, [], fnS);
sgtitle('2 mM ATP — model vs data, all scored slack features');
drawnow;
exportgraphics(gcf, fullfile(outdir,'2mM_plotFeatures_panel.png'), 'Resolution', 200);
fprintf('  saved 2mM_plotFeatures_panel.png\n');
disp('EXPORT DONE');
