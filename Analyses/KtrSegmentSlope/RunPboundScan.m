% RunPboundScan.m -- what does pulling P_bound_max back into the physiological
% range actually cost, on the corrected objective?
%
% P_bound_max is the binding-site occupancy ceiling in
%   f_sat = 1 - P_bound/P_bound_max                 (OccupancyForm 'linear')
% multiplying the attachment flux. Rat cardiac max-Ca isometric bound fraction is
% ~0.12-0.15 (Analyses/BindingSiteOccupancy). refit4_opt carries 1.208, i.e. the
% gate never bites, and that is exactly why the 40 % SL gradient in
% N_overlap*f_lattice sets only the TARGET of force redevelopment and never its
% RATE -- the model's ktr is flat across the slack ladder while every
% preparation's falls.
%
% Before spending a campaign on it: sweep the ceiling on the refit5 objective and
% record the full cost, the per-term breakdown, the ktr ladder slope, and the
% actual bound fraction the model runs at (so we know whether the ceiling is even
% reachable). Seeding the campaign at a good value is worth far more than hoping
% the search rediscovers it from 1.208.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunPboundScan.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~isfolder(resdir); mkdir(resdir); end

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_pbs_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 6); end

p0 = getParams(loadParams('refit5_seed'), [], true, false);
p0.RunForceVelocity=true; p0.RunSlack=true; p0.RunStairs=false; p0.RunKtr=true;
p0.RunForceVelocityTime=false; p0.RunSlackPassive=true; p0.EvalFeatures=true;
p0.BreakOnODEUnstable=false; p0.PlotEachSeparately=0; p0.PlotFeatureFitting=0;
p0.RunSlackSegments='AllPar'; p0.MaxRunTime=300; p0.OptimizeOn='Feats';
p0.ghostSave=''; p0.ghostLoad=''; p0.SaveBest=0; p0.mods={}; p0.g=ones(1,0);

PB = [0.12 0.15 0.20 0.30 0.50 0.80 p0.P_bound_max];
nP = numel(PB);
C = nan(nP,1); SL = nan(nP,1); KT = nan(nP,5); CT = cell(nP,1); NM = {};
PBMAXOBS = nan(nP,1);

for i = 1:nP
    p = p0; p.P_bound_max = PB(i);
    p = getParams(p, p.g, true, false);
    try
        [c,~,ct,nm] = evaluateBakersExp(1, p, true, false);
        C(i) = c; CT{i} = ct; NM = nm;
        % the ktr ladder, from the same parameters
        ps = p; ps.RunSlackSegments = 'All';
        [~, os, fm] = runSlackExperiment(ps);
        KT(i,:) = fm.ktr(:)';
        b = polyfit(1:5, KT(i,:), 1); SL(i) = b(1);
        % how close does the model actually get to the ceiling?
        if isfield(os,'outputs') && ~isempty(os.outputs)
            % columns: [Force F_active F_passive N_overlap p1_0 p2_0 ...]
            PBMAXOBS(i) = max(os.outputs(:,5) + os.outputs(:,6));
        end
    catch ME
        fprintf(2,'P_bound_max %.3f FAILED: %s\n', PB(i), ME.message);
    end
    fprintf('P_bound_max %6.3f   cost %8.4f   ktr slope %+7.3f   ktr %5.1f ... %5.1f   max P_bound %.3f\n', ...
            PB(i), C(i), SL(i), KT(i,1), KT(i,5), PBMAXOBS(i));
end

fprintf('\n=== per-term breakdown (rows = P_bound_max) ===\n');
keep = find(cellfun(@(x) ~isempty(x), CT));
big = [];
for i = keep'; big = union(big, find(CT{i} > 0.005)); end
fprintf('%8s %8s', 'P_bmax', 'TOTAL');
for j = big(:)'; fprintf(' %12s', NM{j}(1:min(12,end))); end
fprintf('\n');
for i = keep'
    fprintf('%8.3f %8.4f', PB(i), C(i));
    for j = big(:)'; fprintf(' %12.4f', CT{i}(j)); end
    fprintf('\n');
end

fprintf('\n=== data target: ktr slope -1.848, ladder 51.61 49.50 47.87 46.32 43.96 ===\n');
[~,ibest] = min(C);
fprintf('best cost at P_bound_max = %.3f (%.4f, slope %+.3f)\n', PB(ibest), C(ibest), SL(ibest));

save(fullfile(resdir,'pbound_scan.mat'),'PB','C','SL','KT','CT','NM','PBMAXOBS');
fprintf('\nDONE\n');
