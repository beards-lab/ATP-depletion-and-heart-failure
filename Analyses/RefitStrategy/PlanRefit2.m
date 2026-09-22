% PlanRefit2.m -- pick the restart point and size the new floors, on numbers.
%
% Three questions:
%  (1) which paramset is actually the best fit on SLACK and FV? Score all three
%      candidates on the same (lean, pooled, floored) objective, broken into
%      blocks, so "best fit" is measured rather than asserted.
%  (2) how far can the slack and FV floors narrow? A floor only bites when the
%      feature is near it, so report each term's distance in preparation-widths
%      and the floor at which it would start paying.
%  (3) A/steady is what reads optically, and it is NOT a feature. Its spread
%      across preparations should be SMALLER than A's or steady's, because the
%      per-preparation strength is common-mode and cancels in the ratio. If it
%      is, the ratio earns a small floor honestly and is the right way to make
%      the optical offset matter.
%
%   matlab -batch "run('Analyses/RefitStrategy/PlanRefit2.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RefitStrategy','results');
c = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_p2_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
c.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(c, 5); end

REF  = getParams(loadParams('ktrfence_seed'), [], true, false);   % owns the objective
FN   = REF.fn(:)';
lab  = cellfun(@(s) regexprep(s,'\|.*$',''), FN, 'UniformOutput', false);
SLACKBLK = {'ktr','A','t0_crossing','restretchSlopeStart','peak1_y','peak1_dSL', ...
            'vall_y','vall_t','peak2','steady','vall2_dy','rsK80[1:4]'};
PASBLK   = {'PS_restretchPeak','PS_rampupRMSE','PS_holdDecayRMSE','PS_steady22','PS_steady20'};

%% ---- (1) which restart point ------------------------------------------
CAND = {'fastgate_best','ktrfence_seed','ktrfence_opt'};
fprintf('\n=============== (1) BEST FIT ON SLACK AND FV ===============\n');
fprintf('all scored on the SAME objective (lean fn, pooled target, current floors)\n\n');
fprintf('%-16s %9s %9s %9s %9s %9s\n','paramset','TOTAL','slack','FV','passive','other');
S = struct();
for i = 1:numel(CAND)
    p = getParams(loadParams(CAND{i}), [], true, false);
    p.fn = FN; p.UsePooledTarget = REF.UsePooledTarget; p.RunKtr = 1;
    p.PlotEachSeparately=0; p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0;
    p.EvalFeatures=1; p.BreakOnODEUnstable=0; p.MaxRunTime=600;
    p = getParams(p, p.g, true);
    [Et,~,fc,~] = evaluateBakersExp([], p, false, false);
    isS = ismember(lab, SLACKBLK); isP = ismember(lab, PASBLK);
    isF = strcmp(lab,'FV_normpowerAvg');
    S(i).tag=CAND{i}; S(i).Et=Et; S(i).fc=fc; S(i).p=p;
    S(i).slack=sum(fc(isS)); S(i).fv=sum(fc(isF)); S(i).pas=sum(fc(isP));
    S(i).other=Et - S(i).slack - S(i).fv - S(i).pas;
    fprintf('%-16s %9.4f %9.4f %9.4f %9.4f %9.4f\n', CAND{i}, Et, S(i).slack, S(i).fv, S(i).pas, S(i).other);
end
[~,best] = min([S.slack] + [S.fv]);
fprintf('\nbest on slack+FV combined: %s (%.4f)\n', S(best).tag, S(best).slack+S(best).fv);

%% ---- (2) how far can the floors narrow ---------------------------------
fprintf('\n=============== (2) FLOOR HEADROOM at %s ===============\n', S(best).tag);
p = S(best).p;
FN0 = regexprep(FN, '\|fl:[0-9eE.+-]+', '');
q = p; q.fn = FN0; [~,~,fc_nf,~] = evaluateBakersExp([], q, false, false);
w = ones(1,numel(FN)); fl = zeros(1,numel(FN));
for i=1:numel(FN)
    tk = split(FN{i},'|');
    for j=2:numel(tk)
        if startsWith(tk{j},'fl:'); fl(i)=str2double(extractAfter(tk{j},'fl:')); end
        vv = str2double(tk{j}); if ~isnan(vv); w(i)=vv; end
    end
end
chi2 = fc_nf(:)'./w; isData = fl>0;
fprintf('%-22s %10s %10s %8s %12s  %s\n','feature','chi2','floor now','prep-w','bites below','');
[~,ord] = sort(chi2./max(fl,eps),'descend');
for i = ord
    if ~isData(i); continue; end
    pw = sqrt(chi2(i)/fl(i));
    bite = chi2(i)/fl(i);                   % scale factor on the CURRENT floor
    note = '';
    if pw > 1; note = 'already paying'; end
    if strcmp(lab{i},'FV_normpowerAvg'); note = [note '   <-- w = 1804, handle with care']; end
    fprintf('%-22s %10.5f %10.5f %8.2f %12.2fx  %s\n', lab{i}, chi2(i), fl(i), pw, bite, note);
end

%% ---- (3) the A/steady ratio, and what spread it really has -------------
fprintf('\n=============== (3) A/steady ACROSS PREPARATIONS ===============\n');
DAYS = {'03_27','04_03','04_10'};
Aa = nan(3,5); Ss = nan(3,5);
for d = 1:3
    F = load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d}))).features_data;
    Aa(d,:) = F.A(:)'; Ss(d,:) = F.steady(:)';
end
Rr = Aa./Ss;
sA = std(Aa,0,1)./mean(Aa,1);
sS = std(Ss,0,1)./mean(Ss,1);
sR = std(Rr,0,1)./mean(Rr,1);
fprintf('%-14s %8s %8s %8s %8s %8s\n','segment','1','2','3','4','5');
fprintf('%-14s %8.3f %8.3f %8.3f %8.3f %8.3f\n','rel SD of A',    sA);
fprintf('%-14s %8.3f %8.3f %8.3f %8.3f %8.3f\n','rel SD steady',  sS);
fprintf('%-14s %8.3f %8.3f %8.3f %8.3f %8.3f\n','rel SD A/steady',sR);
fprintf('\nmean relative SD:  A %.3f   steady %.3f   A/steady %.3f  ->  ratio is %.1fx tighter\n', ...
        mean(sA), mean(sS), mean(sR), mean(sA)/mean(sR));
PT = pooledTarget8mM();
jA=find(strcmp(PT.name,'A'),1); jS=find(strcmp(PT.name,'steady'),1);
tgtR = PT.val(jA,:)./PT.val(jS,:);
fprintf('pooled A/steady target: %s\n', mat2str(round(tgtR,4)));
LAMBDA = 0.0106294; nEl = 5;
sRm = mean(sR);
flR = (1/3)*nEl*sRm^2;  wR = LAMBDA/sRm^2;
fprintf('proposed  A_over_steady|fl:%.4g|%.4g   (n=3 -> s^2/3 ; w*fl = %.6g = lambda*nEl)\n', ...
        flR, wR, wR*flR);
for i=1:numel(CAND)
    fm = []; %#ok
end

%% ---- (4) the ktr fence, retargeted on the 3-day mean -------------------
fprintf('\n=============== (4) ktr FENCE, RETARGETED ===============\n');
kd = zeros(1,3);
for d=1:3
    K = load(fullfile(root,'data',sprintf('protocol_%s_2026_velocitytable_ktr.mat',DAYS{d})));
    f = extractKtrAttributes(K.datatable(:,1), K.datatable(:,3), K.velocitytable);
    kd(d) = f.ktrP_k63;
end
sK = std(kd)/mean(kd);
fprintf('per day %s  mean %.2f  rel SD %.3f\n', mat2str(round(kd,1)), mean(kd), sK);
fprintf('the current fence targets 03/27 (%.1f) -- WRONG population for an s from 3 days.\n', kd(1));
fprintf('retargeted on the mean with n=3: fl = s^2/3 = %.4g, w = lambda/s^2 = %.4g\n', ...
        sK^2/3, LAMBDA/sK^2);
save(fullfile(resdir,'planrefit2.mat'),'S','chi2','fl','w','lab','sA','sS','sR','kd','sK');
fprintf('\nDONE\n');
