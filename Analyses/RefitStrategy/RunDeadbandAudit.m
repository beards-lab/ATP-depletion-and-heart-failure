% RunDeadbandAudit.m -- how far can the deadband be narrowed?
%
% Evaluated at the seed we would actually refit from: lean_opt + the two ktr
% fixes (= L3 of RunRefitStrategy.m, Et 0.5169 on the lean objective).
%
% THE QUESTION. Every data term carries fl = FLOOR_CHI2 * nEl * s^2 with
% FLOOR_CHI2 = 1, i.e. "being within one between-preparation SD is free". That
% was set when the target was ONE preparation. It no longer is: params0.
% UsePooledTarget = 1 replaces the slack features with the MEAN of three
% preparations (Model/pooledTargetFeatures), and the uncertainty on a mean of
% three is s/sqrt(3), not s. So for the pooled terms the defensible deadband is
% s^2/3, i.e. FLOOR_CHI2 = 1/3.
%
% The terms pooledTargetFeatures explicitly does NOT touch -- F-V, the passive
% PS_* -- still have a single-source target, so their deadband stays at 1.
%
% Run: matlab -batch "run('Analyses/RefitStrategy/RunDeadbandAudit.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RefitStrategy','results');
if ~exist(resdir,'dir'); mkdir(resdir); end
c = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_db_%d', feature('getpid')));
if ~exist(jsl,'dir'); mkdir(jsl); end
c.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(c, 5); end

p = getParams(loadParams('lean_opt'), [], true, false);
p.PlotEachSeparately=0; p.ghostSave=''; p.ghostLoad=''; p.BreakOnODEUnstable=0; p.SaveBest=0;
p.UseMaxwellTensionOnly = 1; p.x_M_slack = 0.0108565;   % fix 1 (x_M_slack MUST come too)
p.tau_reg = 2e-4;                                        % fix 2
p = getParams(p, p.g, true);

FN   = p.fn(:)';                       % the REAL entries, floors and all
FN0  = regexprep(FN, '\|fl:[0-9eE.+-]+', '');
lab  = cellfun(@(s) regexprep(s,'\|.*$',''), FN, 'UniformOutput', false);

q = p; q.fn = FN;  [Et_fl,~,fc_fl,~] = evaluateBakersExp([], q, false, false);
q = p; q.fn = FN0; [Et_nf,~,fc_nf,~] = evaluateBakersExp([], q, false, false);

% weights and floors, parsed off the entries
n = numel(FN); w = ones(1,n); fl = zeros(1,n);
for i = 1:n
    tk = split(FN{i},'|');
    for j = 2:numel(tk)
        if startsWith(tk{j},'fl:'); fl(i) = str2double(extractAfter(tk{j},'fl:')); end
        v = str2double(tk{j}); if ~isnan(v); w(i) = v; end
    end
end
raw    = fc_nf(:)';                    % unweighted chi-2, deadband removed
isData = fl > 0;

% which terms have a POOLED (mean-of-3) target. pooledTargetFeatures.m states
% it leaves F-V, the ktr protocol and the passive PS_* alone.
POOLED = {'ktr','A','t0_crossing','restretchSlopeStart','peak1_y','peak1_dSL', ...
          'vall_y','vall_t','peak2','steady','vall2_dy','rsK80'};
isPooled = ismember(lab, POOLED) & isData;

fprintf('\n===== SEED: lean_opt + both ktr fixes =====\n');
fprintf('Et(floored, as shipped) %.4f    Et(no deadband at all) %.4f\n', Et_fl, Et_nf);

fprintf('\n%-22s %6s %9s %9s %8s %10s %10s  %s\n', ...
        'feature','pooled','raw chi2','floor','prep-w','w*floor','w*post','status');
[~,ord] = sort(raw./max(fl,eps),'descend');
nIn = 0; head = 0;
for i = ord
    if ~isData(i); continue; end
    pw = sqrt(raw(i)/fl(i)); post = max(0, raw(i)-fl(i));
    if post <= 0; st='INSIDE (pays 0)'; nIn=nIn+1; head=head+w(i)*raw(i); else; st=''; end
    fprintf('%-22s %6s %9.4f %9.4f %8.2f %10.4f %10.4f  %s\n', lab{i}, ...
            string(isPooled(i)), raw(i), fl(i), pw, w(i)*fl(i), w(i)*post, st);
end
fprintf('\n%d of %d data terms sit INSIDE the 1-prep-width deadband and pay nothing.\n', nIn, nnz(isData));
fprintf('Their un-deadbanded weighted cost would be %.4f, against a floored total of %.4f.\n', head, Et_fl);

%% ---- the chi-square identity: is every term commensurate? ----
fprintf('\n===== w x floor: what one preparation-width COSTS, per term =====\n');
wf = w.*fl; med = median(wf(isData));
[~,o2] = sort(wf,'descend');
for i = o2
    if ~isData(i); continue; end
    fprintf('  %-22s %10.4f   %8.1fx median\n', lab{i}, wf(i), wf(i)/med);
end
fprintf('  median = %.4f  (= lambda*nEl for a term weighted lambda/s^2)\n', med);

%% ---- sweeps ----
iFV = find(strcmp(lab,'FV_normpowerAvg'),1);
sweep = @(fpool, fother) iSweep(raw, fl, w, isData, isPooled, fpool, fother, iFV);
fprintf('\n===== SWEEP 1: narrow EVERYTHING (the naive reading) =====\n');
fprintf('  %8s %12s %8s %12s %10s\n','FLOOR','Et','n live','FV term','FV share');
for f = [1 0.5 1/3 0.25 0.1 0]
    [Et2,live,fv] = sweep(f,f);
    fprintf('  %8.3f %12.4f %8d %12.4f %9.1f%%\n', f, Et2, live, fv, 100*fv/max(Et2,eps));
end
fprintf('\n===== SWEEP 2: narrow only the POOLED terms; others stay at 1 =====\n');
fprintf('  %8s %12s %8s %12s %10s\n','FLOOR','Et','n live','FV term','FV share');
for f = [1 0.5 1/3 0.25 0.1 0]
    [Et2,live,fv] = sweep(f,1);
    fprintf('  %8.3f %12.4f %8d %12.4f %9.1f%%\n', f, Et2, live, fv, 100*fv/max(Et2,eps));
end

fprintf(['\nFLOOR_CHI2 = 1/3 = %.4f is the SEM of a mean of three ' ...
         '(pooledTarget8mM: nP = 3, dof = 4.58).\n'], 1/3);
save(fullfile(resdir,'deadband.mat'),'FN','lab','raw','fl','w','isData','isPooled','Et_fl','Et_nf');
fprintf('\nDONE\n');

function [Et2, live, fv] = iSweep(raw, fl, w, isData, isPooled, fpool, fother, iFV)
    cpost = raw; live = 0;
    for i = 1:numel(raw)
        if ~isData(i); continue; end
        f_i = fother; if isPooled(i); f_i = fpool; end
        cpost(i) = max(0, raw(i) - f_i*fl(i));
        if cpost(i) > 0; live = live + 1; end
    end
    Et2 = sum(w.*cpost);
    fv  = w(iFV)*cpost(iFV);
end
