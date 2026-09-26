% RunNumericalNoise.m -- how much does each scored feature move when the
% parameters move by a part in a million?
%
% WHY. refit7's round-2 best scored 0.3355 in memory and 0.3707 from its own
% snapshot on disk, whose parameters differ by at most 1.1e-6 (the %g
% rounding in writeParamsToMFile). RunRsKJump showed the estimator is not the
% cause -- the anchor minimum is identical -- the model TRACE moves: segment
% 4's re-stretch recovery crossed 80 % at 33.73 ms in one and 32.12 ms in the
% other, a 4 % change in rsK80, which the tight rsK80 term (s = 0.071) plus the
% floor subtraction turned into a 2x change in cost.
%
% A term that rewards differences smaller than the model's own
% reproducibility is optimising noise: the search keeps the lucky evaluation
% and the gain does not survive a reload. So measure that noise per feature,
% and fold it into the deadband:
%       fl = nEl * ( s_between^2 / 3  +  s_num^2 )
% i.e. the same chi-2 construction, with the numerical spread added in
% quadrature to the uncertainty of the 3-heart mean.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunNumericalNoise.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~exist('NREP','var'); NREP = 8; end
if ~exist('EPS_REL','var'); EPS_REL = 1e-6; end
cl=parcluster('local'); jsl=fullfile(tempdir,sprintf('mljob_nn_%d',feature('getpid'))); if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation=jsl; if isempty(gcp('nocreate')); parpool(cl,5); end

S = load(fullfile(root,'params','refit7_state.mat')); st = S.state;
base = st.best_params; base.mods = {}; base.g = ones(1,0);
pool = base.leanPool;
FEAT = {'ktr80','ktrfit','ktr80_shape','ktrfit_shape','rsK80','rsKfit','rsT20','rsTfit','steady','vall_t','peak1_dSL','t0'};
rng(7);
R = struct('fc',{},'feat',{},'kP',{});
for r = 0:NREP
    p = base;
    if r > 0                                  % r = 0 is the unperturbed point
        for i = 1:numel(pool)
            nm = pool{i}; us = strfind(nm,'__');
            if isempty(us)
                if isfield(p,nm) && isscalar(p.(nm)) && p.(nm) ~= 0
                    p.(nm) = p.(nm) * (1 + EPS_REL*randn);
                end
            else
                b = nm(1:us(1)-1); ix = str2double(nm(us(1)+2:end));
                if isfield(p,b) && numel(p.(b)) >= ix && p.(b)(ix) ~= 0
                    v = p.(b); v(ix) = v(ix) * (1 + EPS_REL*randn); p.(b) = v;
                end
            end
        end
    end
    p = getParams(p, p.g, true, false);
    params0 = p; LoadData; RunBakersExp;                          %#ok<NASGU>
    [fc,~,~] = evalFeatureCost(features_data, features_model, p.fn);
    fe = struct();
    for q = FEAT; fe.(q{1}) = features_model.(q{1})(:)'; end
    kP = NaN; if isfield(features_model,'ktrP_k80'); kP = features_model.ktrP_k80; end
    R(end+1) = struct('fc', fc(:)', 'feat', fe, 'kP', kP); %#ok<SAGROW>
    fprintf('rep %d  cost %.4f\n', r, sum(fc));
end
fd = features_data;
lab = cellfun(@(s) regexprep(regexprep(s,'\|.*$',''),'\[.*$',''), p.fn, 'uni', 0);

%% ---- the objective's own noise ------------------------------------------
C  = vertcat(R.fc);
tot = sum(C,2);
fprintf('\nTOTAL: mean %.4f  SD %.4f (%.1f %%)  range %.4f..%.4f\n', mean(tot), std(tot), ...
        100*std(tot)/mean(tot), min(tot), max(tot));
fprintf('%-22s %9s %9s %9s\n','term','mean','SD','SD/mean');
for j = 1:numel(lab)
    if mean(C(:,j)) < 1e-6 && std(C(:,j)) < 1e-6; continue; end
    fprintf('%-22s %9.4f %9.4f %8.1f%%\n', lab{j}, mean(C(:,j)), std(C(:,j)), 100*std(C(:,j))/max(mean(C(:,j)),eps));
end

%% ---- per-feature numerical spread, relative to the DATA mean -----------
% (evalFeatureCost normalises residuals by mean(data), so this is the unit
% the floor lives in)
fprintf('\n%-14s %10s %12s   per-element SD / mean(data)\n','feature','s_num','(elements)');
SNUM = struct();
for q = FEAT
    X = vertcat(arrayfun(@(r) r.feat.(q{1}), R, 'uni', 0)); X = vertcat(X{:});
    if isfield(fd, q{1}); dm = mean(double(fd.(q{1})(:)), 'omitnan'); else; dm = mean(X(:), 'omitnan'); end
    se = std(X, 0, 1) / abs(dm);
    SNUM.(q{1}) = sqrt(mean(se.^2));             % RMS over elements
    fprintf('%-14s %10.4f   %s\n', q{1}, SNUM.(q{1}), mat2str(round(100*se,2)));
end
kPs = [R.kP]; SNUM.ktrP_k80 = std(kPs) / pooledKtrTarget().val(2);
fprintf('%-14s %10.4f\n', 'ktrP_k80', SNUM.ktrP_k80);

save(fullfile(resdir,'numerical_noise.mat'),'R','SNUM','tot','lab','EPS_REL','NREP');
fprintf('\nDONE\n');
