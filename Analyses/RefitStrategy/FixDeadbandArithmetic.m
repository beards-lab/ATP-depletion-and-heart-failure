% FixDeadbandArithmetic.m -- redo the deadband audit's arithmetic.
%
% evaluateBakersExp's 3rd output is the WEIGHTED per-feature cost, not the raw
% chi-2. The first pass multiplied by w a second time, which put F-V at 14.7
% preparation-widths when it is actually at 0.35. No re-simulation needed --
% results/deadband.mat holds the weighted costs.
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RefitStrategy','results');
L = load(fullfile(resdir,'deadband.mat'));
lab=L.lab; fl=L.fl; w=L.w; isData=L.isData; isPooled=L.isPooled;
wraw = L.raw;                 % WEIGHTED, no deadband
chi2 = wraw ./ w;             % <- the actual chi-2

fprintf('\n===== SEED: lean_opt + both ktr fixes =====\n');
fprintf('Et floored (as shipped) %.4f   Et with NO deadband %.4f\n', L.Et_fl, L.Et_nf);
iFV = find(strcmp(lab,'FV_normpowerAvg'),1);
fprintf('Without the deadband F-V alone is %.4f of %.4f = %.1f%% of the objective.\n', ...
        wraw(iFV), L.Et_nf, 100*wraw(iFV)/L.Et_nf);

fprintf('\n%-22s %6s %10s %9s %8s %12s\n','feature','pooled','chi2','floor','prep-w','w*(chi2-fl)');
pw = sqrt(chi2./max(fl,eps));
[~,ord] = sort(pw,'descend');
nIn=0;
for i = ord
    if ~isData(i); continue; end
    post = w(i)*max(0, chi2(i)-fl(i));
    if post<=0; nIn=nIn+1; st='  INSIDE'; else; st=''; end
    fprintf('%-22s %6s %10.5f %9.4f %8.2f %12.4f%s\n', lab{i}, string(isPooled(i)), ...
            chi2(i), fl(i), pw(i), post, st);
end
fprintf('\n%d of %d data terms are inside the 1-prep-width deadband.\n', nIn, nnz(isData));

%% ---- sweeps, done right ----
fprintf('\n===== SWEEP 1: one FLOOR_CHI2 for every data term =====\n');
fprintf('  %8s %10s %8s %10s %10s   %s\n','FLOOR','Et','n live','FV term','FV share','newly live');
prev = pw <= 1;                                  % inside at FLOOR = 1
for f = [1 0.5 1/3 0.25 0.1 0.05 0]
    [Et2,live,fvt,inside] = iEval(chi2, fl, w, isData, isPooled, f, f, iFV);
    newly = lab(isData & ~inside & prev);
    fprintf('  %8.3f %10.4f %8d %10.4f %9.1f%%   %s\n', f, Et2, live, fvt, ...
            100*fvt/max(Et2,eps), strjoin(newly,' '));
    prev = inside;
end

fprintf('\n===== SWEEP 2: floor = s^2/n, n = how many sources the target averages =====\n');
fprintf('pooled slack terms n=3 -> 1/3 ; F-V averages 2 lab curves -> 1/2 ; PS_* n=1 -> 1\n');
[Et2,live,fvt,inside] = iEval(chi2, fl, w, isData, isPooled, 1/3, 1, iFV);
% F-V gets 1/2 explicitly
cpost = chi2;
for i=1:numel(chi2)
    if ~isData(i); continue; end
    f_i = 1; if isPooled(i); f_i = 1/3; end
    if i==iFV; f_i = 0.5; end
    cpost(i) = max(0, chi2(i)-f_i*fl(i));
end
Et3 = sum(w.*cpost) + sum(wraw(~isData));
live3 = nnz(cpost(isData)>0);
fprintf('  Et = %.4f   n live = %d   FV term = %.4f\n', Et3, live3, w(iFV)*cpost(iFV));
fprintf('  live terms: %s\n', strjoin(lab(isData & cpost>0),' '));

fprintf('\n===== what the objective is MADE OF, per setting =====\n');
for f = [1 1/3]
    cpost = chi2;
    for i=1:numel(chi2)
        if ~isData(i); continue; end
        f_i = f; if ~isPooled(i); f_i = 1; end
        cpost(i) = max(0, chi2(i)-f_i*fl(i));
    end
    tot = sum(w.*cpost) + sum(wraw(~isData));
    fprintf('\n  FLOOR_CHI2 = %.3f on pooled terms (others 1):  Et = %.4f\n', f, tot);
    [~,o] = sort(w.*cpost,'descend');
    for i = o(1:min(8,numel(o)))
        if w(i)*cpost(i) <= 0; continue; end
        fprintf('    %-22s %8.4f  (%4.1f%%)\n', lab{i}, w(i)*cpost(i), 100*w(i)*cpost(i)/tot);
    end
    nd = wraw(~isData); ndl = lab(~isData);
    [~,o2] = sort(nd,'descend');
    for k = o2(1:min(4,numel(o2)))
        if nd(k) <= 0; continue; end
        fprintf('    %-22s %8.4f  (%4.1f%%)  [guard/boundary, no deadband]\n', ndl{k}, nd(k), 100*nd(k)/tot);
    end
end
save(fullfile(resdir,'deadband_fixed.mat'),'chi2','pw','fl','w','lab','isData','isPooled');
fprintf('\nDONE\n');

function [Et2, live, fvt, inside] = iEval(chi2, fl, w, isData, isPooled, fpool, fother, iFV)
    cpost = chi2; inside = true(size(chi2));
    for i = 1:numel(chi2)
        if ~isData(i); continue; end
        f_i = fother; if isPooled(i); f_i = fpool; end
        cpost(i) = max(0, chi2(i) - f_i*fl(i));
        inside(i) = cpost(i) <= 0;
    end
    Et2  = sum(w.*cpost);
    live = nnz(~inside & isData);
    fvt  = w(iFV)*cpost(iFV);
end
