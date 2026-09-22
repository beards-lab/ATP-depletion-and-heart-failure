% RunVerify5.m -- closing verification of the refit5 campaign.
%
%   1. reload params/refit5_opt.m from disk and confirm it reproduces the
%      recorded 0.6885 (a snapshot that does not reproduce is not a result)
%   2. per-feature breakdown against refit5_seed
%   3. the ktr LADDER -- the point of the whole exercise: does the slope now
%      fall, and by how much, against the pooled -1.85/segment
%   4. held out: the staircase (RunStairs = 0 throughout) and the ktr protocol
%      recovery on all three days via validateKtrProtocol
%   5. the P_bound_max the search settled on, and whether the occupancy gate
%      is any steeper than it was
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunVerify5.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~isfolder(resdir); mkdir(resdir); end

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_v5_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 6); end

TAGS = {'refit5_seed','refit5_opt'};
PT = pooledTarget8mM();
jk = find(strcmp(PT.name,'ktr'),1);
ktr_pool = PT.val(jk,:);
b = polyfit(1:5, ktr_pool, 1); SLOPE_DATA = b(1);

R = struct();
for k = 1:2
    p = getParams(loadParams(TAGS{k}), [], true, false);
    p.RunForceVelocity=true; p.RunSlack=true; p.RunStairs=false; p.RunKtr=true;
    p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
    p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
    p.RunSlackSegments='AllPar'; p.MaxRunTime=600; p.OptimizeOn='Feats';
    p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0; p.mods={}; p.g=ones(1,0);
    [R(k).cost,~,R(k).fc,R(k).fn] = evaluateBakersExp(1, p, true, false);
    R(k).p = p;
    % repeat, same process: the objective must be deterministic
    R(k).cost2 = evaluateBakersExp(1, p, true, false);
    % the ktr ladder
    ps = p; ps.RunSlackSegments='All';
    [~,~,fm] = runSlackExperiment(ps);
    R(k).ktr = fm.ktr(:)';
    bb = polyfit(1:5, R(k).ktr, 1); R(k).slope = bb(1);
    R(k).fm = fm;
    % held out: staircase
    pst = p; pst.RunStairs = 1; pst = getParams(pst, pst.g, true);
    R(k).Estairs = runStairsExperiment(pst);
    fprintf('%-14s cost %.4f (repeat %.4f, delta %.2g) | ktr slope %+0.3f | stairs %.4f\n', ...
            TAGS{k}, R(k).cost, R(k).cost2, R(k).cost2-R(k).cost, R(k).slope, R(k).Estairs);
end

fprintf('\n=== reproduction ===\n');
fprintf('  refit5_opt on disk    %.4f\n', R(2).cost);
fprintf('  campaign log recorded 0.6885\n');
fprintf('  difference            %+.4f  (%.2f %%)\n', R(2).cost-0.6885, 100*(R(2).cost-0.6885)/0.6885);

%% ---- per-feature -------------------------------------------------------
c1 = R(1).fc(:); c2 = R(2).fc(:); d = c2-c1;
lab = cellfun(@(s) regexprep(s,'\|.*$',''), R(1).fn, 'UniformOutput', false);
[~,od] = sort(d,'ascend');
fprintf('\n=== per-feature (refit5_seed -> refit5_opt) ===\n');
fprintf('%-24s %10s %10s %11s\n','feature','seed','opt','delta');
for i = od(:)'
    if abs(d(i)) < 1e-6 && c1(i) < 1e-6; continue; end
    fprintf('%-24s %10.4f %10.4f %+11.4f\n', lab{i}, c1(i), c2(i), d(i));
end
fprintf('%-24s %10.4f %10.4f %+11.4f\n','TOTAL',sum(c1),sum(c2),sum(d));
fprintf('%-24s %10.4f %10.4f %+11.4f   <- HELD OUT\n','stairs', ...
        R(1).Estairs, R(2).Estairs, R(2).Estairs-R(1).Estairs);

%% ---- the ktr ladder ----------------------------------------------------
fprintf('\n=== the ktr ladder (SL 2.04 -> 1.88 um) ===\n');
fprintf('%-16s %7s %7s %7s %7s %7s   %10s %8s\n','source','seg1','seg2','seg3','seg4','seg5','slope','span');
fprintf('%-16s %7.2f %7.2f %7.2f %7.2f %7.2f   %+10.3f %7.1f%%\n','POOLED DATA', ...
        ktr_pool, SLOPE_DATA, 100*(ktr_pool(5)-ktr_pool(1))/ktr_pool(1));
for k = 1:2
    fprintf('%-16s %7.2f %7.2f %7.2f %7.2f %7.2f   %+10.3f %7.1f%%\n', TAGS{k}, ...
            R(k).ktr, R(k).slope, 100*(R(k).ktr(5)-R(k).ktr(1))/R(k).ktr(1));
end
fprintf('  slope closed: %.0f %% of the gap\n', ...
        100*(R(2).slope-R(1).slope)/(SLOPE_DATA-R(1).slope));

%% ---- the occupancy gate ------------------------------------------------
fprintf('\n=== occupancy gate ===\n');
for k = 1:2
    ps = R(k).p; ps.RunSlackSegments='All';
    [~,os] = runSlackExperiment(ps);
    Pb = os.p1_0(:) + os.p2_0(:);
    fprintf('%-14s P_bound_max %.4f | P_bound med %.4f max %.4f | f_sat %.3f / %.3f | dfsat/dPb %.2f\n', ...
            TAGS{k}, ps.P_bound_max, median(Pb), max(Pb), ...
            1-median(Pb)/ps.P_bound_max, 1-max(Pb)/ps.P_bound_max, -1/ps.P_bound_max);
end

%% ---- held out: the ktr protocol ----------------------------------------
fprintf('\n=== HELD OUT: ktr protocol ===\n');
try
    validateKtrProtocol('refit5_opt');
catch ME
    fprintf(2,'validateKtrProtocol failed: %s\n', ME.message);
end

save(fullfile(resdir,'verify5.mat'),'R','TAGS','ktr_pool','SLOPE_DATA','lab');
fprintf('\nDONE\n');
