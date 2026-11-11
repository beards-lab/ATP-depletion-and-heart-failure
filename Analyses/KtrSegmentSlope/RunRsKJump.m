% RunRsKJump.m -- where does rsK80 jump between two parameter sets that
% differ by 1e-6 (the campaign's in-memory best vs its %g snapshot on disk)?
% Instruments fitRestretchRecovery through the RSK80_DIAG global.
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunRsKJump.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','KtrSegmentSlope','results');
cl=parcluster('local'); jsl=fullfile(tempdir,sprintf('mljob_rj_%d',feature('getpid'))); if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation=jsl; if isempty(gcp('nocreate')); parpool(cl,5); end
S=load(fullfile(root,'params','refit7_state.mat')); st=S.state;
pm = st.best_params;
pd = getParams(loadParams('refit7_opt'),[],true,false);
for f = {'RunSlackSegments','MaxRunTime','PlotEachSeparately','PlotFeatureFitting','EvalFeatures','ghostSave','ghostLoad','SaveBest'}
    pd.(f{1}) = pm.(f{1});
end
P = {pm, pd}; NAME = {'memory','disk'};
global RSK80_DIAG
D = cell(1,2); FM = cell(1,2);
for k = 1:2
    RSK80_DIAG = struct([]);
    [~,~,fm] = runSlackExperiment(P{k});
    D{k} = RSK80_DIAG; FM{k} = fm;
    fprintf('%-7s %d estimator calls | rsK80 %s | rsT20 %s\n', NAME{k}, numel(D{k}), ...
            mat2str(round(fm.rsK80,2)), mat2str(round(1e3*fm.rsT20,2)));
end
RSK80_DIAG = [];
n = min(numel(D{1}), numel(D{2}));
fprintf('\n%4s | %-7s %8s %8s %8s %8s %8s\n','call','set','tB(ms)','B','t20(ms)','t80(ms)','k80');
jump = zeros(1,n);
for i = 1:n
    for k = 1:2
        d = D{k}(i);
        fprintf('%4d | %-7s %8.2f %8.3f %8.2f %8.2f %8.2f\n', i, NAME{k}, 1e3*d.tB, d.B, 1e3*d.t20, 1e3*d.t80, d.k80);
    end
    jump(i) = abs(D{1}(i).k80 - D{2}(i).k80);
end
[~,iw] = max(jump);
fprintf('\nlargest jump: call %d, |dk80| = %.2f\n', iw, jump(iw));
f = figure('Position',[40 40 1100 450],'Color','w'); hold on; box on; grid on
C = {[0 0.45 0.74],[0.85 0.33 0.1]};
for k = 1:2
    d = D{k}(iw);
    plot(1e3*d.tt, d.zx, '-', 'Color', C{k}, 'LineWidth', 1.5);
    plot(1e3*[d.t20 d.t80], [0.2 0.8], 'o', 'Color', C{k}, 'MarkerFaceColor', C{k});
end
yline(0.2,':'); yline(0.8,':'); xlim([0 80]);
xlabel('ms after the recovery minimum'); ylabel('z = (F - B)/(steady - B)');
legend({'memory z','memory t20/t80','disk z','disk t20/t80'},'Location','southeast','Box','off');
title(sprintf('rsK80 estimator, call %d: k80 %.1f (memory) vs %.1f (disk)', iw, D{1}(iw).k80, D{2}(iw).k80));
exportgraphics(f, fullfile(resdir,'rsk80_jump.png'),'Resolution',110);
save(fullfile(resdir,'rsk80_jump.mat'),'D','FM','jump');
fprintf('DONE\n');
