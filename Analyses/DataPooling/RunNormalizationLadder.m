% RunNormalizationLadder.m -- which normalisation makes the three preparations
% agree, and what does it cost the ATP signal?
%
% THE PROBLEM
% Absolute tension (kPa) is force divided by an ESTIMATED cross-sectional area.
% That estimate is the least reproducible number in the experiment, so every
% force-dimension observable of preparation p carries a common unknown
% multiplicative factor alpha_p. Fitting absolute kPa fits alpha, not the model.
% Dividing it out is normalisation -- but by WHICH divisor, at WHICH scope, and
% after which corrections?
%
% THE MODEL
%   y_{p,c,k} = alpha_p^{d_k} * ( mu_{c,k} + pi_{p,k} ) * exp(eps)
%   p = preparation (3)  c = condition (8 / 2 mM)  k = feature
%   d_k = 1 for force-dimension features, 0 for rates / times / ratios
%   pi  = the preparation's own PASSIVE (titin) contribution, measured
%         independently in its PNB + mavacamten run
% alpha_p is a property of the PREPARATION, not of the run, so both ATP
% conditions of a day share it. That is what makes the 2 mM / 8 mM ratio an
% alpha-free observable, and it is why the `scope` axis below is not a matter of
% taste. Scale-free features carry d_k = 0 and no scheme can touch them: they
% are the built-in control group and the floor on what normalisation can buy.
%
% FOUR CROSSED FACTORS
%   anchor    what estimates alpha_p        none | Fss22 | Fiso20a | Fiso20b
%                                           | Fpeak | Fsteady | MoR
%   scope     per run, or per preparation   run | prep
%   rundown   lift the run recorded second  off | on
%   passive   subtract the day's own PNB+Mava force  off | on
%
% THE CRITERION -- and why "lowest spread wins" is not enough on its own
% Any divisor drawn from the data can shrink the spread of the very feature it
% came from, to zero in the limit. Every scheme is therefore scored
% LEAVE-ONE-FEATURE-OUT: the scale is re-estimated with feature k held out
% before k's spread is measured. Reported per scheme:
%   sBar   LOFO between-preparation relative SD, RMS over force features
%   sMed   the same, median (RMS is dominated by the two loose features)
%   VE     fraction of between-preparation log-variance the scheme removes
%   zATP   discriminability of the 2 mM / 8 mM force effect AFTER normalising
% A scheme that wins on sBar by destroying zATP has removed the signal, not the
% nuisance.
%
% Output -> results/normalisation_ladder.mat + printed tables

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here, 'results'); if ~isfolder(resDir); mkdir(resDir); end

PHI = 0.452;      % permanent fraction of the within-run decline (RundownCorrection M.4)

R = poolLoadRuns();
days = unique({R.day}, 'stable'); nP = numel(days);

%% ---- feature registry --------------------------------------------------
% class F = force-dimension (d_k = 1, carries alpha)
%       0 = scale-free (rate / time / length / ratio) -- the control group
% passiveProxy: which PNB+Mava feature is this feature's passive contribution.
% `A` and `Am` are recovery amplitudes at the SLACK length and the passive run
% has no recovery to fit, so they borrow `vall_y`, the passive force at that
% same length.
FEAT = { 'A','F','vall_y'; 'Am','F','vall_y'; 'peak1_y','F','peak1_y'
         'vall_y','F','vall_y'; 'peak2','F','peak2'; 'steady','F','steady'
         'vall2_dy','F','vall2_dy'; 'rsA','F','rsA'
         'restretchSlopeStart','F','restretchSlopeStart'
         'ktr','0',''; 'rsK','0',''; 'rsK63','0',''; 'ktr2_k','0',''
         't0','0',''; 'peak1_t','0',''; 'vall_t','0',''; 'vall2_t','0',''
         'rsT0','0',''; 'peak1_dSL','0',''; 'rsR2','0','' };
isF = strcmp(FEAT(:,2),'F');
nF  = size(FEAT,1); kF = find(isF); nEl = 5;

%% ---- assemble Y (nP x nCond x nF x nEl), anchors, passive --------------
COND = {'8mM','2mM'};
ancNames = {'Fss22','Fiso20a','Fiso20b','Fpeak','Fsteady'};
Y   = nan(nP, numel(COND), nF, nEl);
PAS = nan(nP, nF, nEl);                        % the day's PNB+Mava features
ANC = struct(); for a = 1:numel(ancNames); ANC.(ancNames{a}) = nan(nP,numel(COND)); end
ANCP = nan(nP,1);                              % passive level of the anchor
[ORD, SLOPE, TACT] = deal(nan(nP,numel(COND)));
for p = 1:nP
    for c = 1:numel(COND)
        i = find(strcmp({R.day},days{p}) & strcmp({R.cond},COND{c}));
        for k = 1:nF
            v = double(R(i).feat.(FEAT{k,1})(:))'; m = min(nEl,numel(v));
            Y(p,c,k,1:m) = v(1:m);
        end
        for a = 1:numel(ancNames); ANC.(ancNames{a})(p,c) = R(i).A.(ancNames{a}); end
        ORD(p,c) = R(i).order; SLOPE(p,c) = R(i).slope; TACT(p,c) = R(i).Tact;
    end
    ip = find(strcmp({R.day},days{p}) & strcmp({R.cond},'PNB'));
    for k = kF(:)'
        v = double(R(ip).feat.(FEAT{k,3})(:))'; m = min(nEl,numel(v));
        PAS(p,k,1:m) = v(1:m);
    end
    ANCP(p) = R(ip).A.Fss22;                   % passive at SL 2.2, same window
end

%% ---- rundown lift of the run recorded second --------------------------
lift = ones(nP,numel(COND));
for p = 1:nP
    ce = ORD(p,:)==1; cl = ORD(p,:)==2;
    lift(p,cl) = 1 + PHI*abs(SLOPE(p,ce))*TACT(p,ce) / ANC.Fss22(p,cl);
end
fprintf('=== rundown lift applied to the run recorded second ===\n');
for p = 1:nP
    fprintf('  %-6s  %s -> %s   lift = %.4f (+%.1f%%)\n', days{p}, ...
        COND{ORD(p,:)==1}, COND{ORD(p,:)==2}, lift(p,ORD(p,:)==2), ...
        100*(lift(p,ORD(p,:)==2)-1));
end
fprintf('\n=== passive (PNB+Mava) share of each preparation, at SL 2.2 ===\n');
for p = 1:nP
    fprintf('  %-6s  passive %5.2f kPa  =  %4.1f%% of its 8 mM Fss22 (%.1f kPa)\n', ...
        days{p}, ANCP(p), 100*ANCP(p)/ANC.Fss22(p,1), ANC.Fss22(p,1));
end

%% ---- the ladder --------------------------------------------------------
SCHEME = {};
for a = [ancNames {'MoR'} {'none'}]
    for sc = {'run','prep'}
        for rd = [false true]
            for pa = [false true]
                if strcmp(a{1},'none') && strcmp(sc{1},'run'); continue; end
                SCHEME(end+1,:) = {a{1}, sc{1}, rd, pa}; %#ok<SAGROW>
            end
        end
    end
end
nS = size(SCHEME,1);

res = struct('name',{},'anchor',{},'scope',{},'rundown',{},'passive',{}, ...
             'sK',{},'sBar',{},'sMed',{},'VE',{},'rATP',{},'sATP',{},'zATP',{});

% baseline spread: raw kPa, nothing removed
s_raw = nan(1,nF); for k = 1:nF; s_raw(k) = relSD(squeeze(Y(:,1,k,:))); end

for s = 1:nS
    [anc, scope, rd, pa] = deal(SCHEME{s,1}, SCHEME{s,2}, SCHEME{s,3}, SCHEME{s,4});

    Yn = Y; ancV = ANC;
    if pa                                       % subtract the measured passive
        for c = 1:numel(COND); Yn(:,c,kF,:) = Yn(:,c,kF,:) - reshape(PAS(:,kF,:),nP,1,numel(kF),nEl); end
        for a = 1:numel(ancNames); ancV.(ancNames{a}) = ancV.(ancNames{a}) - ANCP; end
    end
    if rd                                       % lift the degraded run
        for c = 1:numel(COND); Yn(:,c,kF,:) = bsxfun(@times, Yn(:,c,kF,:), lift(:,c)); end
        for a = 1:numel(ancNames); ancV.(ancNames{a}) = ancV.(ancNames{a}) .* lift; end
    end

    [alpha, alphaLOO] = estimateScale(Yn, ancV, anc, scope, kF, FEAT);

    sK = nan(1,nF);
    for k = 1:nF
        V = squeeze(Yn(:,1,k,:));
        if isF(k); V = V ./ alphaLOO(:,1,k); end
        sK(k) = relSD(V);
    end
    sBar = sqrt(mean(sK(kF).^2,'omitnan'));
    sMed = median(sK(kF),'omitnan');
    VE   = 1 - sum(sK(kF).^2,'omitnan')/sum(s_raw(kF).^2,'omitnan');

    % ATP effect on the force LEVEL. Only meaningful under `prep` scope: a
    % per-RUN scale divides each condition by its own level, so whatever it
    % reports about the ATP level is an artefact of the gauge, not a
    % measurement. Marked NaN there rather than printed as a number.
    if strcmp(scope,'prep')
        rp = nan(nP,1);
        for p = 1:nP
            v8 = squeeze(Yn(p,1,kF,:))./alpha(p,1);
            v2 = squeeze(Yn(p,2,kF,:))./alpha(p,2);
            q  = v2(:)./v8(:); q = q(isfinite(q) & q > 0);
            rp(p) = exp(mean(log(q)));
        end
        rATP = mean(rp); sATP = std(rp);
        zATP = abs(rATP-1)/sqrt(sATP^2 + sBar^2);
    else
        [rATP, sATP, zATP] = deal(NaN);
    end

    res(end+1) = struct('name', sprintf('%s/%s%s%s', anc, scope, tern(rd,'+rd',''), tern(pa,'+pas','')), ...
        'anchor',anc,'scope',scope,'rundown',rd,'passive',pa, ...
        'sK',sK,'sBar',sBar,'sMed',sMed,'VE',VE,'rATP',rATP,'sATP',sATP,'zATP',zATP); %#ok<SAGROW>
end

%% ---- report: the ladder -----------------------------------------------
fprintf('\n=== normalisation ladder, sorted by LOFO spread (n = %d preparations) ===\n', nP);
fprintf('(ATP columns are blank for `run` scope: a per-run scale defines the ATP\n');
fprintf(' level ratio away, so nothing there is a measurement.)\n');
fprintf('%-24s %7s %7s %6s %8s %8s %7s\n','scheme','sBar','sMed','VE','ATPratio','sd(ATP)','z(ATP)');
[~,ord] = sort([res.sBar]);
for i = ord
    if isnan(res(i).rATP); atp = sprintf('%8s %8s %7s','-','-','-');
    else; atp = sprintf('%8.3f %8.3f %7.2f', res(i).rATP, res(i).sATP, res(i).zATP); end
    fprintf('%-24s %6.1f%% %6.1f%% %6.2f %s\n', res(i).name, ...
        100*res(i).sBar, 100*res(i).sMed, res(i).VE, atp);
end

%% ---- is the ranking of the anchors real, at n = 3? ---------------------
% The schemes are scored on the SAME nine force features, so compare them
% paired: bootstrap the feature set (and, inside it, the elements) and ask how
% often each scheme beats the leader. Anything that wins less than ~95 % of the
% time is not distinguishable from the leader on this dataset.
rng(7); NB = 4000;
cmp = {'Fss22/prep+rd','Fsteady/prep+rd','MoR/prep+rd','Fiso20a/prep+rd', ...
       'Fiso20b/prep+rd','Fpeak/prep+rd','none/prep+rd'};
idx = cellfun(@(nm) find(strcmp({res.name},nm),1), cmp);
SK  = cell2mat(arrayfun(@(r) r.sK(kF)', res(idx), 'uni', 0));   % nFF x nCmp
win = zeros(1,numel(idx));
for b = 1:NB
    j = randi(size(SK,1), size(SK,1), 1);
    [~, w] = min(sqrt(mean(SK(j,:).^2, 1, 'omitnan')));
    win(w) = win(w) + 1;
end
fprintf('\n=== paired feature bootstrap: which anchor is really best? (%d draws) ===\n', NB);
fprintf('%-20s %8s %10s\n','scheme','sBar','P(best)');
for i = 1:numel(idx)
    fprintf('%-20s %7.1f%% %9.1f%%\n', cmp{i}, 100*res(idx(i)).sBar, 100*win(i)/NB);
end

%% ---- report: marginal effect of each factor ---------------------------
fprintf('\n=== marginal effect of each factor (mean over all other levels) ===\n');
fprintf('%-14s %-10s %7s %7s %8s\n','factor','level','sBar','sMed','z(ATP)');
show = @(lab,lev,m) fprintf('%-14s %-10s %6.1f%% %6.1f%% %8.2f\n', lab, lev, ...
        100*mean([res(m).sBar]), 100*mean([res(m).sMed]), mean([res(m).zATP],'omitnan'));
for a = [ancNames {'MoR'} {'none'}]
    show('anchor', a{1}, strcmp({res.anchor},a{1}));
end
for sc = {'run','prep'}; show('scope', sc{1}, strcmp({res.scope},sc{1})); end
show('rundown','off',~[res.rundown]); show('rundown','on',[res.rundown]);
show('passive','off',~[res.passive]); show('passive','on',[res.passive]);

%% ---- report: per-feature, for the leading schemes ---------------------
pick = {'none/prep','Fss22/prep','Fss22/prep+rd','Fss22/prep+rd+pas','MoR/prep+rd+pas'};
fprintf('\n=== per-feature LOFO spread ===\n%-22s', 'feature');
for q = pick; fprintf(' %18s', q{1}); end; fprintf('\n');
for k = [kF(:)' find(~isF)']
    fprintf('%-22s', [FEAT{k,1} tern(isF(k),'','  (scale-free)')]);
    for q = pick
        j = find(strcmp({res.name},q{1}),1);
        fprintf(' %17.1f%%', 100*res(j).sK(k));
    end
    fprintf('\n');
end

%% ---- is ONE scale per preparation enough? ------------------------------
% Centre the log force features by feature (8 mM), then take the SVD. A pure
% common scale is a rank-1 matrix whose right singular vector is FLAT: every
% feature loads equally. Both must hold.
Lg = log(reshape(squeeze(Y(:,1,kF,:)), nP, []));
Lg = Lg(:, all(isfinite(Lg),1) & all(imag(Lg)==0,1)); Lg = Lg - mean(Lg,1);
[~,Sv,V] = svd(Lg,'econ'); sv = diag(Sv).^2; sv = sv/sum(sv);
cosFlat = abs(V(:,1)'*ones(size(V,1),1)/sqrt(size(V,1)));
fprintf('\n=== is a single scale per preparation enough? ===\n');
fprintf('  variance in PC1 of the centred log force features : %.1f%%\n', 100*sv(1));
fprintf('  |cos| between the PC1 loading and the flat direction : %.3f\n', cosFlat);
fprintf('  rank-1 AND flat  =>  the between-prep difference IS one common scale\n');

%% ---- how well is the spread itself determined? -------------------------
nu = nP-1; ci = sqrt(nu./[chi2q(0.975,nu) chi2q(0.025,nu)]);
fprintf('\n=== how many measurements do we actually have? ===\n');
fprintf('  n = %d preparations -> nu = %d dof for EVERY spread estimate\n', nP, nu);
fprintf('  95%% CI on any s : [%.2f x s, %.2f x s]\n', ci(1), ci(2));
fprintf('  a feature measured at s = 10%% is only known to lie in %.1f-%.1f%%\n', ...
        100*0.10*ci(1), 100*0.10*ci(2));

save(fullfile(resDir,'normalisation_ladder.mat'), 'res','FEAT','isF','kF','days', ...
     'COND','ANC','ANCP','PAS','ORD','SLOPE','TACT','lift','Y','s_raw','sv','cosFlat','ci');
fprintf('\nsaved -> %s\n', fullfile(resDir,'normalisation_ladder.mat'));

%% ======================= helpers ========================================
function [alpha, alphaLOO] = estimateScale(Yn, ancV, anc, scope, kF, FEAT)
% alpha    : nP x nCond
% alphaLOO : nP x nCond x nFeat -- the divisor for feature k, estimated WITHOUT
%            using feature k. Identical to alpha for anchors external to the
%            feature set; genuinely leave-one-out for MoR and Fsteady.
    [nP,nC,nFt,~] = size(Yn);
    alpha = ones(nP,nC); alphaLOO = ones(nP,nC,nFt);
    switch anc
        case 'none'
            return
        case 'MoR'
            for k = [0 kF(:)']
                use = setdiff(kF,k);
                a = mor(Yn,1,use); a2 = tern(strcmp(scope,'run'), mor(Yn,2,use), a);
                if k == 0; alpha = [a a2];
                else;      alphaLOO(:,1,k) = a; alphaLOO(:,2,k) = a2; end
            end
        otherwise
            v = ancV.(anc);
            switch scope
                case 'run',  alpha = v ./ mean(v(:,1));
                case 'prep', alpha = repmat(v(:,1)/mean(v(:,1)), 1, nC);
            end
            for k = 1:nFt; alphaLOO(:,:,k) = alpha; end
            if strcmp(anc,'Fsteady')       % the anchor IS the `steady` feature
                ks = find(strcmp(FEAT(:,1),'steady'),1);
                alphaLOO(:,:,ks) = NaN;    % nothing left to rebuild it from
            end
    end
end

function a = mor(Yn, c, use)
% Median-of-ratios scale: the single factor that best lines this preparation up
% with the others across ALL retained force features at once. Reference-free
% (DESeq size factors / the same idea in proteomics). Gauge: geometric mean 1.
    nP = size(Yn,1);
    Lg = reshape(log(squeeze(Yn(:,c,use,:))), nP, []);
    Lg(~isfinite(Lg)) = NaN;
    a = exp(median(Lg - mean(Lg,1,'omitnan'), 2, 'omitnan'));
    a = a / exp(mean(log(a)));
end

function s = relSD(V)
% RMS over elements of the between-preparation relative SD -- matches the
% sum-over-elements form of evalFeatureCost, so 1/s^2 is a chi-square weight.
    mu = mean(V,1,'omitnan'); sd = std(V,0,1,'omitnan');
    g = isfinite(mu) & abs(mu) > 0 & isfinite(sd);
    if ~any(g); s = NaN; return; end
    s = sqrt(mean((sd(g)./abs(mu(g))).^2));
end

function o = tern(c,a,b); if c; o = a; else; o = b; end; end

function x = chi2q(p,nu)
    if exist('chi2inv','file') == 2; x = chi2inv(p,nu); return; end
    T = struct('p0_025_2',0.050636, 'p0_975_2',7.377759);   % nu = 2 only
    x = T.(sprintf('p%s_%d', strrep(sprintf('%.3f',p),'.','_'), nu));
end
