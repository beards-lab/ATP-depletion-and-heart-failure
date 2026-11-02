% RegenPooledTargets.m -- after the floored-campaign merge: regenerate the 2 mM
% pooled target with the CURRENT feature set, and prove the 8 mM one would not move.
%
%   8 mM: written to a staging folder only and compared row by row with the live
%         Model/pooledTarget8mM.m (never installed from here).
%   2 mM: written straight to Model/pooledTarget2mM.m (the Sep 16 file lacks
%         ktrfit, ktr80, rsTfitE, tOn20 and the vertex features).
%   matlab -batch "run('Analyses/LowATP_Refit9K2Sweep/RegenPooledTargets.m')"
rootR = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(rootR); addpath(genpath(rootR));
diary(fullfile(rootR,'Analyses','LowATP_Refit9K2Sweep','results','regen_targets.log')); diary on;
stage = fullfile(tempdir, 'regen8_stage'); if ~isfolder(stage); mkdir(stage); end
iWrite(8, fullfile(stage, 'pooledTarget8mM.m'));
clear pooledTarget8mM; live = pooledTarget8mM();
addpath(stage, '-begin'); clear pooledTarget8mM; new = pooledTarget8mM(); rmpath(stage); clear pooledTarget8mM;
fprintf('8 mM staged vs live: name %d val %d s %d z %d Fref %d\n', isequal(new.name, live.name), ...
        isequaln(new.val, live.val), isequaln(new.s, live.s), isequaln(new.z, live.z), isequal(new.Fref, live.Fref));
if isequal(new.name, live.name)
    fprintf('   max |dval|/|val| = %.3g   max |ds| = %.3g\n', max(abs(new.val(:)-live.val(:))./max(abs(live.val(:)),eps)), max(abs(new.s-live.s)));
end
iWrite(2, '');
clear pooledTarget2mM; T2 = pooledTarget2mM();
fprintf('2 mM target: %d features, Fref %.2f, cond %d\n', numel(T2.name), T2.Fref, T2.cond);
for q = {'steady','ktrfit','ktr80','rsK80','rsTfitE','t0'}
    j = find(strcmp(T2.name, q{1}), 1);
    if ~isempty(j); fprintf('   %-8s val %s  s %.3f\n', q{1}, mat2str(round(T2.val(j,:),4)), T2.s(j)); end
end
K8 = pooledKtrTarget(8); K2 = pooledKtrTarget(2);
for j = 1:numel(K8.name)
    fprintf('ktr protocol %-9s 8 mM %.3f (s %.3f) | 2 mM %.3f (s %.3f) | ratio %.3f\n', K8.name{j}, K8.val(j), K8.s(j), K2.val(j), K2.s(j), K2.val(j)/K8.val(j));
end
F = fv2Target();
fprintf('FV2 P/P0: 8 mM %.3f +- %.3f | 2 mM %.3f +- %.3f | ratio %.3f +- %.3f\n', F.m8, F.s8, F.m2, F.s2, F.r, F.sr);
fprintf('DONE\n'); diary off;

function iWrite(COND, OUTF) %#ok<INUSD>
    % WritePooledTarget is a script that clears its workspace -> run it in here
    OPTS = struct(); %#ok<NASGU>
    run(fullfile('C:\home\git\ATP-depletion-and-heart-failure','Analyses','DataPooling','WritePooledTarget.m'));
end
