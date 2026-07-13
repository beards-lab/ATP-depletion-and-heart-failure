% featureSensitivity.m
% =========================================================================
% FEATURE-RESOLVED sensitivity of the best 2-state fit (params_2state_a2hop.m,
% cost ~2.95). Unlike calcSensitivities (total-cost only), this perturbs every
% parameter by +25% and records the elasticity of EACH feature:
%     S(feature, param) = (dfeature/feature0) / (dparam/param0)
% and, using the data targets, whether moving that param TOWARD or AWAY from data.
% Output: for each residual feature, the ranked levers + whether up/down helps.
%
% Writes: Analyses/FeatureSensitivity/featureSensitivity_report.txt
%         Analyses/FeatureSensitivity/featureSensitivity.mat  (S matrix etc.)
% Run:  cd(root); addpath(genpath('.')); featureSensitivity
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));

SNAP   = 'params/params_2state_a2hop.m';
REPORT = fullfile('Analyses','FeatureSensitivity','featureSensitivity_report.txt');
MRT    = 90;   PERT = 0.25;   % +25% one-sided perturbation

params0 = getParams(); run(SNAP); pbase = params0;

% ---- parameters to probe (2-state comprehensive pool) ----
plist = {'ka','kd','k1','k_1','k2','k_2','kah', ...
    'kstiff1','kstiff2','kstiff1_n','kstiff2_n','kSE','kSE_M','eta_M','mu_neg', ...
    'dr','dr2','ksr0','kmsrd','ksr2srd','ksrd2sr','sigma1','sigma2','sigma_srd1','sigma_srd2', ...
    'A0','v_ref_reg','tau_reg','P_bound_max','k_pas','kA2hop','sA2hop', ...
    'PieceWiseStrainDepParams__2','PieceWiseStrainDepParams__3','PieceWiseStrainDepParams__4', ...
    'PieceWiseStrainDep2Params__2','PieceWiseStrainDep2Params__3','PieceWiseStrainDep2Params__4', ...
    'PieceWiseStrainDep2Params__5','PieceWiseStrainDep2Params__6','PieceWiseStrainDep2Params__7', ...
    'PieceWiseStrainDepR1DParams__2','PieceWiseStrainDepR1DParams__3','PieceWiseStrainDepR1DParams__4', ...
    'PieceWiseStrainDepR21Params__2','PieceWiseStrainDepR21Params__3'};

% ---- features to track (mean over segments) + FV points ----
featNames = {'ktr','A','peak1_y','peak1_dSL','vall_y','peak2','ovrsht_dy','steady', ...
    'restretchSlopeStart','t0_crossing','FV1','FV2','FV3','FV4'};

% ---- base evaluation ----
[~, cv0, fm0, fd] = costOfSnap(SNAP, pbase.fn, MRT);
f0 = featvec(fm0, featNames);
fdv = featvec(fd,  featNames);   % data targets (FV from fd; slack feats from fd)

np = numel(plist); nf = numel(featNames);
S = nan(nf, np);
for j = 1:np
    p = plist{j};
    bv = seedval(pbase, p);
    try
        [~,~,fmj,~] = costOfSnap(SNAP, pbase.fn, MRT, struct(p, bv*(1+PERT)));
        fj = featvec(fmj, featNames);
        S(:,j) = ((fj - f0) ./ f0) / PERT;   % elasticity
    catch e
        fprintf('  %s perturb ERROR: %s\n', p, e.message);
    end
    fprintf('done %d/%d: %s\n', j, np, p);
end

% ---- write report: per feature, ranked levers + helpful direction ----
fid = fopen(REPORT,'w');
fprintf(fid, '==== featureSensitivity on %s (base cost from cv) ====\n', SNAP);
fprintf(fid, 'S = elasticity (dfeat/feat0)/(dparam/param0), +25%% one-sided.\n');
fprintf(fid, 'For each feature: model0, data, gap(=data-model); then top levers.\n');
fprintf(fid, '  "UP helps" = increasing param moves feature toward data; "DOWN helps" = decreasing does.\n\n');

for i = 1:nf
    gap = fdv(i) - f0(i);
    fprintf(fid, '---- %s : model=%.4g  data=%.4g  gap=%+.4g ----\n', featNames{i}, f0(i), fdv(i), gap);
    [~, ord] = sort(abs(S(i,:)), 'descend', 'MissingPlacement','last');
    kk = 0;
    for t = 1:np
        j = ord(t); s = S(i,j);
        if ~isfinite(s) || abs(s) < 0.05; continue; end
        dir = 'UP helps'; if sign(s) ~= sign(gap); dir = 'DOWN helps'; end
        if gap == 0; dir = '(matched)'; end
        fprintf(fid, '    %-32s S=%+6.2f   %s\n', plist{j}, s, dir);
        kk = kk + 1; if kk >= 6; break; end
    end
    fprintf(fid, '\n');
end
fclose(fid);
save(fullfile('Analyses','FeatureSensitivity','featureSensitivity.mat'), 'S','plist','featNames','f0','fdv','cv0');
type(REPORT);
disp('DONE featureSensitivity');

% ---- helpers ----
function v = seedval(p, name)
    us = strfind(name, '__');
    if isempty(us); v = p.(name);
    else, f = name(1:us(1)-1); idx = str2double(name(us(1)+2:end)); a = p.(f); v = a(idx); end
end
function fv = featvec(fm, names)
    fv = nan(numel(names),1);
    for i = 1:numel(names)
        n = names{i};
        if strncmp(n,'FV',2)
            k = str2double(n(3)) + 1;  % FV1->index2 ... (drop v=0)
            if isfield(fm,'FV_fnorm') && numel(fm.FV_fnorm)>=k; fv(i) = fm.FV_fnorm(k); end
        elseif isfield(fm,n) && ~isempty(fm.(n))
            fv(i) = mean(fm.(n),'omitnan');
        end
    end
end
