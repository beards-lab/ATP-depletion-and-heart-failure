% RunLeverScan2.m — the group RunLeverScan.m did not reach (it died after
% `kSE x0.6` while three MATLAB jobs were competing for 12 cores).
%
% These are the MECHANICS levers: series elasticity, stiffness, the Maxwell
% dashpot and the strain offsets, plus the A2 re-stretch hop. They matter here
% because the post-restretch window is the one window where the series element
% is loaded throughout, so a compliance lever could in principle slow the
% force rise without touching the cross-bridge rates.
%
% Merges into the same results/lever_scan.mat (as `rows2`).

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if isempty(gcp('nocreate')); parpool('local', 6); end

SNAP    = 'params/optfull7_opt_mov.m';
FACTORS = [0.6 1.7];
OUT     = fullfile(here, 'results', 'lever_scan2.mat');

cand = {'kSE','ekSE','kstiff2','kstiff1','kstiff2_n','kstiff1_n', ...
        'mu_neg','mu','kSE_M','eta_M','estiff', ...
        'dr','dr1','dr2', ...
        'kA2hop','sA2hop', ...
        'A0','v_ref_reg','tau_reg','ka','P_bound_max', ...
        'LXBpivot','L_thick','L_thin'};

base = getParams(loadParams(SNAP), [], true, false);
base.MaxRunTime = 400;
opts = struct('slackOnly', true);

fprintf('=== reference run ===\n');
R0 = evalRs(base, opts);
fprintf('featTotal %.3f | rsK %.1f (x%.2f) | ktr %.1f | steady %.1f\n', ...
    R0.featTotal, mean(R0.rsK_m), R0.rsK_ratio, mean(R0.ktr_m), mean(R0.steady_m));

rows2 = struct('param',{},'factor',{},'rsK',{},'rsKratio',{},'ktr',{}, ...
               'steady',{},'rsA',{},'featTotal',{},'dFeat',{},'dlnRsK',{},'err',{});
n = 0; tAll = tic;
for ic = 1:numel(cand)
    nm = cand{ic};
    v0 = getp(base, nm);
    if isempty(v0) || ~isfinite(v0) || v0 == 0
        fprintf('%-30s SKIP (absent or zero)\n', nm); continue;
    end
    for f = FACTORS
        n = n + 1;
        p = setp(base, nm, v0*f);
        p = getParams(p, [], true, false);
        R = evalRs(p, opts);
        rows2(n).param = nm; rows2(n).factor = f;
        if ~R.ok
            rows2(n).err = R.err;
            fprintf('%-30s x%-4.2g  FAILED: %s\n', nm, f, R.err);
        else
            rows2(n).err = '';
            rows2(n).rsK = mean(R.rsK_m,'omitnan');  rows2(n).rsKratio = R.rsK_ratio;
            rows2(n).ktr = mean(R.ktr_m,'omitnan');  rows2(n).steady   = mean(R.steady_m,'omitnan');
            rows2(n).rsA = mean(R.rsA_m,'omitnan');  rows2(n).featTotal = R.featTotal;
            rows2(n).dFeat  = R.featTotal - R0.featTotal;
            rows2(n).dlnRsK = log(rows2(n).rsK / mean(R0.rsK_m,'omitnan'));
            fprintf(['%-30s x%-4.2g  rsK %6.1f (x%4.2f)  dlnRsK %+6.3f | ktr %5.1f | ' ...
                     'steady %5.1f | feat %7.3f (%+7.3f)\n'], nm, f, rows2(n).rsK, ...
                rows2(n).rsKratio, rows2(n).dlnRsK, rows2(n).ktr, rows2(n).steady, ...
                rows2(n).featTotal, rows2(n).dFeat);
        end
        save(OUT, 'rows2', 'R0', 'SNAP', 'FACTORS', 'cand');
    end
end
fprintf('\n=== group 2 done, %d evals in %.1f min ===\n', n, toc(tAll)/60);

function v = getp(p, name)
    us = strfind(name, '__');
    if isempty(us)
        if isfield(p,name) && isnumeric(p.(name)) && isscalar(p.(name)); v = p.(name); else; v = []; end
    else
        f = name(1:us(1)-1); i = str2double(name(us(1)+2:end));
        if isfield(p,f) && numel(p.(f)) >= i; v = p.(f)(i); else; v = []; end
    end
end
function p = setp(p, name, val)
    us = strfind(name, '__');
    if isempty(us); p.(name) = val;
    else; f = name(1:us(1)-1); i = str2double(name(us(1)+2:end)); p.(f)(i) = val; end
end
