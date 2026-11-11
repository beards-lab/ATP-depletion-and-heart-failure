% RunComboCost.m — what does the k2 + eta_M combination cost on the real
% protocol and the full battery?
%
% RunSlopeUnderLevers found that two levers together fix BOTH halves of the
% post-restretch defect on the synthetic amplitude sweep:
%   k2 x0.6 + eta_M x0.3  ->  rsK ~49 s^-1 at EVERY amplitude, slope +21
%                             (baseline: 80->137 s^-1, slope +714; data ~44-52,
%                              slope -126)
% k2 sets the LEVEL, eta_M removes the AMPLITUDE DEPENDENCE. That is a much
% better outcome than the one-at-a-time scan suggested, so it has to be checked
% where it counts: on the REAL protocol, with the full battery, including the
% passive experiment where a Maxwell-element change should send its bill.
%
% Also tests whether the known force/rate compensations (ka, kstiff2) can pay
% for k2's side effects, since that is what the optimiser will try to do.
%
% Writes results/combo_cost.mat

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if isempty(gcp('nocreate')); parpool('local', 5); end

SNAP = 'params/optfull7_opt_mov.m';
base = getParams(loadParams(SNAP), [], true, false);
base.MaxRunTime = 400;
opts = struct('slackOnly', false);

VAR = { 'baseline',                  {}; ...
        'eta_M x0.3',                {'eta_M',0.3}; ...
        'k2 x0.6',                   {'k2',0.6}; ...
        'k2 .6 + etaM .3',           {'k2',0.6,'eta_M',0.3}; ...
        'k2 .7 + etaM .3',           {'k2',0.7,'eta_M',0.3}; ...
        'k2 .6 + etaM .3 + ka 1.4',  {'k2',0.6,'eta_M',0.3,'ka',1.4}; ...
        'k2 .6 + etaM .3 + kst2 .85',{'k2',0.6,'eta_M',0.3,'kstiff2',0.85}; ...
        'k2 .7 + etaM .3 + ka 1.2',  {'k2',0.7,'eta_M',0.3,'ka',1.2} };

res = struct('name',{},'featTotal',{},'rsK',{},'rsKratio',{},'ktr',{}, ...
             'ktrratio',{},'steady',{},'byName',{},'err',{});
R0 = [];
for iv = 1:size(VAR,1)
    nm = VAR{iv,1}; spec = VAR{iv,2};
    p = base;
    for ip = 1:2:numel(spec); p.(spec{ip}) = p.(spec{ip}) * spec{ip+1}; end
    p = getParams(p, [], true, false);
    R = evalRs(p, opts);
    n = numel(res)+1; res(n).name = nm;
    if ~R.ok
        res(n).err = R.err;
        fprintf('%-28s FAILED: %s\n', nm, R.err); continue;
    end
    if isempty(R0); R0 = R; end
    res(n).err = ''; res(n).featTotal = R.featTotal;
    res(n).rsK = mean(R.rsK_m,'omitnan'); res(n).rsKratio = R.rsK_ratio;
    res(n).ktr = mean(R.ktr_m,'omitnan'); res(n).ktrratio = R.ktr_ratio;
    res(n).steady = mean(R.steady_m,'omitnan'); res(n).byName = R.byName;
    fprintf(['%-28s rsK %6.1f (x%4.2f) | ktr %5.1f (x%4.2f) | steady %5.1f | ' ...
             'feat %8.3f\n'], nm, res(n).rsK, res(n).rsKratio, res(n).ktr, ...
             res(n).ktrratio, res(n).steady, res(n).featTotal);
end

save(fullfile(here,'results','combo_cost.mat'), 'res', 'VAR', 'SNAP');

%% ---- which features pay? ------------------------------------------------
ok = cellfun(@isempty, {res.err});
r  = res(ok);
if numel(r) > 1
    fns = fieldnames(r(1).byName);
    fprintf('\nPER-FEATURE COST (weighted, L1) — baseline then deltas\n');
    fprintf('%-22s %8s', 'feature', 'base');
    for i = 2:numel(r); fprintf(' %11s', shortname(r(i).name)); end
    fprintf('\n');
    for k = 1:numel(fns)
        b = r(1).byName.(fns{k});
        row = sprintf('%-22s %8.3f', fns{k}, b);
        big = false;
        for i = 2:numel(r)
            d = r(i).byName.(fns{k}) - b;
            row = [row sprintf(' %+11.3f', d)]; %#ok<AGROW>
            if abs(d) > 0.5; big = true; end
        end
        if big || b > 0.5; fprintf('%s\n', row); end
    end
    fprintf('%-22s %8.3f', 'TOTAL', r(1).featTotal);
    for i = 2:numel(r); fprintf(' %+11.3f', r(i).featTotal - r(1).featTotal); end
    fprintf('\n');
end

function s = shortname(n)
    s = strrep(strrep(strrep(n,' ',''),'+','_'),'.','');
    if numel(s) > 11; s = s(1:11); end
end
