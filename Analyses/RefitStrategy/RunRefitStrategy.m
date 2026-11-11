% RunRefitStrategy.m -- everything needed to choose the refit, in one batch job.
%
% (A) Can the LEAN lineage carry the two ktr fixes, and what do they cost on the
%     lean (floored, pooled) objective?
% (B) The deadband. FLOOR_CHI2 = 1 today. Which features are inside it and
%     therefore paying nothing, how far inside, and what narrowing would buy.
%
% Run:  matlab -batch "run('Analyses/RefitStrategy/RunRefitStrategy.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RefitStrategy','results');
if ~exist(resdir,'dir'); mkdir(resdir); end

% own job storage so this cannot collide with the desktop session's cluster
c = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_refit_%d', feature('getpid')));
if ~exist(jsl,'dir'); mkdir(jsl); end
c.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(c, 5); end

b = getParams(loadParams('lean_opt'), [], true, false);
b.PlotEachSeparately=0; b.ghostSave=''; b.ghostLoad=''; b.BreakOnODEUnstable=0; b.SaveBest=0;
kd = load(fullfile(root,'data','protocol_03_27_2026_velocitytable_ktr.mat'));
KVT = kd.velocitytable; KDT = kd.datatable; NOISE = 0.0222;

%% ================= (A) lean + the two ktr fixes =================
% NB x_M_slack MUST come back with the flag. BuildLeanSeed zeroed it when it
% DELETED the mechanism, and dPUdT_CombinedTransitions.m:322 divides by it, so
% flipping the flag alone gives Inf. XMS is fastgate_best's value.
XMS = 0.0108565;
mt  = @(p) setfield(setfield(p,'UseMaxwellTensionOnly',1),'x_M_slack',XMS);
fg  = @(p) setfield(p,'tau_reg',2e-4);
V = { 'L0  lean_opt as is',       @(p) p
      'L1  + MaxwellTensionOnly', mt
      'L2  + fast gate',          fg
      'L3  + both',               @(p) fg(mt(p)) };

fprintf('\n=============== (A) LEAN + KTR FIXES ===============\n');
fprintf('%-28s %10s %10s %10s %12s\n','variant','Et(lean)','dEt','ktr k63','rmse/noise');
A = struct('name',{},'Et',{},'fc',{},'k63',{},'rn',{});
for v = 1:size(V,1)
    p = V{v,2}(b); p = getParams(p, p.g, true);
    [Et,~,fc,FN] = evaluateBakersExp([], p, false, false);
    [k63, rn] = iKtr(p, KVT, KDT, NOISE);
    A(v) = struct('name',V{v,1},'Et',Et,'fc',fc,'k63',k63,'rn',rn);
    if v == 1; d = '--'; else; d = sprintf('%+.4f', Et-A(1).Et); end
    fprintf('%-28s %10.4f %10s %10.1f %11.2fx\n', V{v,1}, Et, d, k63, rn);
end
fprintf('%-28s %10s %10s %10.1f %11.2fx\n','DATA 03/27','--','--',49.4,1.00);
lab = cellfun(@(s) regexprep(s,'\|.*$',''), FN, 'UniformOutput', false);
fprintf('\n-- per-feature, L3 vs L0 --\n');
dd = A(4).fc - A(1).fc; [~,o] = sort(abs(dd),'descend');
for i = o(1:min(8,numel(o)))
    fprintf('%-24s %8.4f -> %8.4f  (%+.4f)\n', lab{i}, A(1).fc(i), A(4).fc(i), dd(i));
end

%% ================= (B) the deadband =================
% Recover each entry's floor, then re-evaluate with floors STRIPPED to get the
% raw chi-2. floor = FLOOR_CHI2 * nEl * s^2 with FLOOR_CHI2 = 1, so the floor IS
% the chi-2 of a model exactly one preparation-width away: raw/floor is the
% model's distance in SQUARED preparation-widths.
FN_nofloor = regexprep(FN, '\|fl:[0-9eE.+-]+', '');
fprintf('\n=============== (B) DEADBAND AUDIT ===============\n');
BASES = {'L0 lean_opt', b; 'L3 lean + both fixes', fg(mt(b))};
B = struct();
for s = 1:size(BASES,1)
    p = BASES{s,2}; p = getParams(p, p.g, true);
    p.fn = FN;         [Et_fl,~,fc_fl,~] = evaluateBakersExp([], p, false, false);
    p.fn = FN_nofloor; [Et_nf,~,fc_nf,~] = evaluateBakersExp([], p, false, false);

    % weights and floors, parsed straight off the fn entries
    w = zeros(1,numel(FN)); fl = zeros(1,numel(FN));
    for i = 1:numel(FN)
        tk = split(FN{i},'|');
        for j = 2:numel(tk)
            if startsWith(tk{j},'fl:'); fl(i) = str2double(extractAfter(tk{j},'fl:')); end
            n = str2double(tk{j});
            if ~isnan(n); w(i) = n; end
        end
        if w(i) == 0; w(i) = 1; end
    end
    raw    = fc_nf;        % unweighted chi-2, no deadband
    isData = fl > 0;       % only the floored entries are data terms

    B(s).name=BASES{s,1}; B(s).Et_fl=Et_fl; B(s).Et_nf=Et_nf;
    B(s).raw=raw; B(s).fl=fl; B(s).w=w; B(s).isData=isData; B(s).fc_fl=fc_fl;

    fprintf('\n--- %s   Et(floored) %.4f   Et(no floor) %.4f ---\n', BASES{s,1}, Et_fl, Et_nf);
    fprintf('%-22s %9s %9s %8s %10s %10s  %s\n', ...
            'feature','raw chi2','floor','prep-w','w*floor','w*post','status');
    [~,ord] = sort(raw./max(fl,eps),'descend');
    nInside = 0; headroom = 0;
    for i = ord
        if ~isData(i); continue; end
        pw   = sqrt(raw(i)/fl(i));
        post = max(0, raw(i)-fl(i));
        if post <= 0
            st = 'INSIDE (pays 0)'; nInside = nInside+1; headroom = headroom + w(i)*raw(i);
        else
            st = '';
        end
        fprintf('%-22s %9.4f %9.4f %8.2f %10.4f %10.4f  %s\n', ...
                lab{i}, raw(i), fl(i), pw, w(i)*fl(i), w(i)*post, st);
    end
    fprintf('  -> %d of %d data terms are INSIDE the 1-prep-width deadband and pay nothing;\n', ...
            nInside, nnz(isData));
    fprintf('     their un-floored weighted cost would be %.4f (vs the floored total %.4f)\n', ...
            headroom, Et_fl);

    fprintf('\n  FLOOR_CHI2 sweep (same weights, same model -- only the deadband moves):\n');
    fprintf('  %8s %12s %8s %12s %12s\n','FLOOR','Et','n live','FV term','FV share');
    iFV = find(strcmp(lab,'FV_normpowerAvg'), 1);
    for fc2 = [1 0.5 0.25 0.1 0]
        cpost = raw; live = 0;
        for i = 1:numel(FN)
            if isData(i)
                cpost(i) = max(0, raw(i) - fc2*fl(i));
                if cpost(i) > 0; live = live+1; end
            end
        end
        Et2   = sum(w.*cpost);
        fvterm = w(iFV)*cpost(iFV);
        fprintf('  %8.2f %12.4f %8d %12.4f %11.1f%%\n', fc2, Et2, live, fvterm, 100*fvterm/max(Et2,eps));
    end
    % the same sweep with FV held at one prep-width
    fprintf('\n  ... and with FV alone pinned at FLOOR_CHI2 = 1:\n');
    fprintf('  %8s %12s %8s %12s %12s\n','FLOOR','Et','n live','FV term','FV share');
    for fc2 = [1 0.5 0.25 0.1 0]
        cpost = raw; live = 0;
        for i = 1:numel(FN)
            if isData(i)
                f_i = fc2; if i == iFV; f_i = 1; end
                cpost(i) = max(0, raw(i) - f_i*fl(i));
                if cpost(i) > 0; live = live+1; end
            end
        end
        Et2    = sum(w.*cpost);
        fvterm = w(iFV)*cpost(iFV);
        fprintf('  %8.2f %12.4f %8d %12.4f %11.1f%%\n', fc2, Et2, live, fvterm, 100*fvterm/max(Et2,eps));
    end
end

%% ---- the chi-square identity, and where FV sits relative to it ----
fprintf('\n=============== (C) IS THE OBJECTIVE A CHI-SQUARE? ===============\n');
fprintf('For a term weighted w = lambda/s^2 with floor = nEl*s^2, w*floor = lambda*nEl:\n');
fprintf('every term costs the SAME to sit one preparation-width out. Deviations from\n');
fprintf('that are pins, not measurements.\n\n');
w = B(1).w; fl = B(1).fl; isData = B(1).isData;
wf = w.*fl;
med = median(wf(isData));
fprintf('%-22s %12s %12s\n','feature','w*floor','x median');
[~,ord] = sort(wf,'descend');
for i = ord
    if ~isData(i); continue; end
    fprintf('%-22s %12.4f %11.1fx\n', lab{i}, wf(i), wf(i)/med);
end
fprintf('\nmedian w*floor = %.4f  (= lambda*nEl for the chi-square-consistent terms)\n', med);

save(fullfile(resdir,'refitstrategy.mat'),'A','B','FN','FN_nofloor','lab');
fprintf('\nDONE\n');

function [k63, rn] = iKtr(p, KVT, KDT, NOISE)
    pk = p; pk.SL0 = 2*KVT(1,4); pk = getParams(pk, pk.g, true); pk.Velocity = KVT(:,2);
    [~, ok] = evaluateModel(str2func(pk.modelFcn), KVT(:,1), pk);
    FisoK = interp1(ok.t, ok.Force, KVT(2,1), 'linear', 'extrap');
    f  = extractKtrAttributes(ok.t, ok.Force/FisoK, KVT);
    Fi = interp1(ok.t, ok.Force/FisoK, KDT(:,1), 'linear', 'extrap');
    m  = KDT(:,1) >= KVT(end-1,1);
    k63 = f.ktrP_k63;
    rn  = sqrt(mean((Fi(m)-KDT(m,3)).^2))/NOISE;
end
