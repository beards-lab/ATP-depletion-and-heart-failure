% RunFinalVerify.m -- close out the ktrfence campaign.
%   1. cost, reloaded from the snapshot (NOT the campaign header -- %g loss)
%   2. the fence: did it stay inert, i.e. did tau_reg hold?
%   3. the held-out ktr protocol, all three days
%   4. the ktr OVERSTRETCH PEAK -- not in the objective, and the thing that
%      falsified the ApoPool mechanism. Measured, not eyeballed.
%   5. what the refit moved, per feature, against the pooled target
%   6. which parameters actually moved
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RefitStrategy','results');
c = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_fv_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
c.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(c, 5); end

TAGS = {'ktrfence_seed','ktrfence_opt'};
kd  = load(fullfile(root,'data','protocol_03_27_2026_velocitytable_ktr.mat'));
KVT = kd.velocitytable; KDT = kd.datatable;

fprintf('\n================ 1-2. COST AND THE FENCE ================\n');
P = cell(1,2); C = zeros(1,2); FC = cell(1,2);
for k = 1:2
    p = getParams(loadParams(TAGS{k}), [], true, false);
    p.PlotEachSeparately=0; p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0;
    p.EvalFeatures=1; p.BreakOnODEUnstable=0; p.MaxRunTime=600;
    p = getParams(p, p.g, true);
    [C(k),~,FC{k},FN] = evaluateBakersExp([], p, false, false);
    P{k} = p;
    iF = find(strcmp(FN,'ktrP_k63'),1);
    fprintf('%-16s cost %.4f   tau_reg %8.4g ms   MaxwellTensionOnly %d   fence %.6f\n', ...
            TAGS{k}, C(k), p.tau_reg*1e3, p.UseMaxwellTensionOnly, FC{k}(iF));
end
fprintf('campaign header said 0.5271; reloaded %.4f  (%%g snapshot round-trip)\n', C(2));
fprintf('improvement: %.4f -> %.4f  (%.1f %%)\n', C(1), C(2), 100*(1-C(2)/C(1)));

fprintf('\n================ 3. HELD-OUT ktr PROTOCOL ================\n');
for k = 1:2; validateKtrProtocol(TAGS{k}); end

fprintf('\n================ 4. THE ktr OVERSTRETCH PEAK (not scored) ================\n');
% the 1.05 ML overstretch lives between the restretch ramp start and the end
% of the final release: KVT rows end-4 .. end-1
tA = KVT(end-4,1); tB = KVT(end-1,1);
Fiso_d = 1.0;                                  % the data trace is already normalised
md = KDT(:,1) >= tA & KDT(:,1) <= tB;
pk_d = max(KDT(md,3))/Fiso_d;
fprintf('window %.4f .. %.4f s  (the 0.799 -> 1.0509 -> 1.000 ML manoeuvre)\n', tA, tB);
fprintf('%-16s %12s %12s\n','','peak F/F_iso','vs data');
fprintf('%-16s %12.3f %12s\n','DATA 03/27', pk_d, '--');
for k = 1:2
    pk_ = P{k}; pk_.SL0 = 2*KVT(1,4); pk_ = getParams(pk_, pk_.g, true);
    pk_.Velocity = KVT(:,2);
    [~, o] = evaluateModel(str2func(pk_.modelFcn), KVT(:,1), pk_);
    Fi = interp1(o.t, o.Force, KVT(2,1), 'linear','extrap');
    m  = o.t >= tA & o.t <= tB;
    pk = max(o.Force(m))/Fi;
    fprintf('%-16s %12.3f %11.2fx\n', TAGS{k}, pk, pk/pk_d);
end

fprintf('\n================ 5. PER-FEATURE, vs THE POOLED TARGET ================\n');
lab = cellfun(@(s) regexprep(s,'\|.*$',''), FN, 'UniformOutput', false);
d = FC{2} - FC{1}; [~,o] = sort(d,'descend');
fprintf('%-24s %10s %10s %10s\n','feature','seed','opt','delta');
for i = o(:)'
    if abs(d(i)) < 1e-6 && FC{1}(i) < 1e-6; continue; end
    fprintf('%-24s %10.4f %10.4f %+10.4f\n', lab{i}, FC{1}(i), FC{2}(i), d(i));
end
fprintf('%-24s %10.4f %10.4f %+10.4f\n','TOTAL', C(1), C(2), C(2)-C(1));

fprintf('\n================ 6. WHAT MOVED ================\n');
a = P{1}; b = P{2}; rows = {};
for f = fieldnames(a)'
    va = a.(f{1}); vb = b.(f{1});
    if ~isnumeric(va) || ~isscalar(va) || ~isnumeric(vb) || ~isscalar(vb); continue; end
    if va == 0 || ~isfinite(va) || ~isfinite(vb); continue; end
    r = vb/va;
    if abs(r-1) > 0.02; rows(end+1,:) = {f{1}, va, vb, r}; end %#ok
end
[~,ix] = sort(abs(cell2mat(rows(:,4))-1),'descend');
fprintf('%-30s %14s %14s %8s\n','parameter','seed','opt','ratio');
for i = ix(:)'
    fprintf('%-30s %14.6g %14.6g %8.3f\n', rows{i,1}, rows{i,2}, rows{i,3}, rows{i,4});
end
fprintf('\n%d of %d pooled parameters moved by more than 2 %%\n', numel(ix), numel(P{1}.leanPool));
save(fullfile(resdir,'finalverify.mat'),'C','FC','FN','lab','rows');
fprintf('\nDONE\n');
