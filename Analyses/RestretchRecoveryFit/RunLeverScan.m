% RunLeverScan.m — one-at-a-time sweep: which parameters move the
% post-restretch redevelopment rate rsK, and what do they cost elsewhere?
%
% Base = params/optfull7_opt_mov.m, where rsK is 2.70x too fast while the
% recovery AMPLITUDE (rsA) and the post-slack rate (ktr) are already right.
% So we are looking for a lever with a large d(ln rsK) and a small
% d(featTotal) — i.e. one that is not welded to force level or to ktr.
%
% Slack-only evaluation (47 s vs 122 s); the top hits are re-checked on the
% full battery afterwards. Writes results/lever_scan.mat incrementally so a
% crash or a kill loses at most one evaluation.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if ~exist(fullfile(here, 'results'), 'dir'); mkdir(fullfile(here, 'results')); end

SNAP    = 'params/optfull7_opt_mov.m';
FACTORS = [0.6 1.7];
OUT     = fullfile(here, 'results', 'lever_scan.mat');

% ---- candidate levers, grouped by the mechanism they probe --------------
cand = { ...
  ... % A. strain-dependent detachment — the leading suspect. R2 reaches
  ... %    2100-4250/s at restretch strains, so every stretched bridge is gone
  ... %    in <0.5 ms and the recovery is pure fast re-attachment.
  'PieceWiseStrainDep2X_logOffset', 'PieceWiseStrainDep2Params__5', ...
  'PieceWiseStrainDep2Params__6', 'PieceWiseStrainDep2Params__7', ...
  'PieceWiseStrainDep2X__4', 'PieceWiseStrainDep2X__5', ...
  'PieceWiseStrainDepR1DX_logOffset', 'PieceWiseStrainDepR1DParams__3', ...
  'PieceWiseStrainDepR1DParams__4', ...
  'PieceWiseStrainDepR21X_logOffset', 'PieceWiseStrainDepX_logOffset', ...
  ... % B. attachment side — refill speed
  'ka', 'A0', 'v_ref_reg', 'tau_reg', 'P_bound_max', ...
  ... % C. core cycle rates
  'k2', 'kd', 'k1', 'k_1', 'kah', 'kamh', 'k_2', ...
  ... % D. SRX manifold (predicted NOT to work — included to re-test at this base)
  'ksr0', 'kmsrd', 'ksr2srd', 'ksrd2sr', 'sigma2', 'sigma_srd1', 'sigma1', ...
  ... % E. mechanics / series elasticity
  'kSE', 'ekSE', 'kstiff2', 'kstiff1', 'mu_neg', 'kSE_M', 'eta_M', 'dr2', 'dr', ...
  ... % F. A2 restretch hop
  'kA2hop', 'sA2hop'};

% A pool speeds the 5 slack cycles up ~4x. Only evaluateModel runs on the
% workers (extractSlackAttributes is client-side), but start a FRESH pool
% anyway so no worker holds a stale compiled copy of an edited model file.
% 4 workers, not 6: the overnight optimizer runs concurrently on the same box.
if isempty(gcp('nocreate'))
    parpool('local', 4);
end

base = getParams(loadParams(SNAP), [], true, false);
base.MaxRunTime = 300;
opts = struct('slackOnly', true);

fprintf('=== reference run ===\n');
R0 = evalRs(base, opts);
fprintf('featTotal %.3f | rsK %.1f (x%.2f) | ktr %.1f (x%.2f) | steady %.1f\n', ...
    R0.featTotal, mean(R0.rsK_m), R0.rsK_ratio, mean(R0.ktr_m), R0.ktr_ratio, mean(R0.steady_m));

rows = struct('param', {}, 'factor', {}, 'val', {}, 'rsK', {}, 'rsKratio', {}, ...
              'ktr', {}, 'ktrratio', {}, 'steady', {}, 'rsA', {}, ...
              'featTotal', {}, 'dFeat', {}, 'dlnRsK', {}, 'err', {}, 'time', {});

n = 0; tAll = tic;
for ic = 1:numel(cand)
    nm = cand{ic};
    v0 = getp(base, nm);
    if isempty(v0) || ~isfinite(v0) || v0 == 0
        fprintf('%-38s SKIP (absent or zero: cannot scale)\n', nm); continue;
    end
    for f = FACTORS
        n = n + 1;
        p = setp(base, nm, v0 * f);
        p = getParams(p, [], true, false);      % re-resolve derived fields
        R = evalRs(p, opts);
        rows(n).param  = nm;   rows(n).factor = f;   rows(n).val = v0*f;
        rows(n).time   = R.time;
        if ~R.ok
            rows(n).err = R.err;
            fprintf('%-30s x%-4.2g  FAILED: %s\n', nm, f, R.err);
        else
            rows(n).err       = '';
            rows(n).rsK       = mean(R.rsK_m, 'omitnan');
            rows(n).rsKratio  = R.rsK_ratio;
            rows(n).ktr       = mean(R.ktr_m, 'omitnan');
            rows(n).ktrratio  = R.ktr_ratio;
            rows(n).steady    = mean(R.steady_m, 'omitnan');
            rows(n).rsA       = mean(R.rsA_m, 'omitnan');
            rows(n).featTotal = R.featTotal;
            rows(n).dFeat     = R.featTotal - R0.featTotal;
            rows(n).dlnRsK    = log(rows(n).rsK / mean(R0.rsK_m, 'omitnan'));
            fprintf(['%-30s x%-4.2g  rsK %6.1f (x%4.2f)  dlnRsK %+6.3f | ' ...
                     'ktr %5.1f | steady %5.1f | feat %7.3f (%+7.3f)\n'], ...
                nm, f, rows(n).rsK, rows(n).rsKratio, rows(n).dlnRsK, ...
                rows(n).ktr, rows(n).steady, rows(n).featTotal, rows(n).dFeat);
        end
        save(OUT, 'rows', 'R0', 'SNAP', 'FACTORS', 'cand');
    end
end

fprintf('\n=== scan done, %d evals in %.1f min ===\n', n, toc(tAll)/60);

%% ---- ranking: efficiency = how much rsK slowing per unit cost elsewhere --
ok  = ~cellfun(@isempty, {rows.rsK});
r   = rows(ok);
slw = [r.dlnRsK] < 0;                       % only levers that SLOW the recovery
r   = r(slw);
[~, ord] = sort([r.dlnRsK]);                % most slowing first
fprintf('\nTop levers that SLOW the post-restretch recovery:\n');
fprintf('%-30s %6s %8s %9s %9s %8s\n', 'param', 'factor', 'dlnRsK', 'dFeat', 'ktrRatio', 'steady');
for i = ord(1:min(20, numel(ord)))
    fprintf('%-30s %6.2f %8.3f %9.3f %9.2f %8.1f\n', ...
        r(i).param, r(i).factor, r(i).dlnRsK, r(i).dFeat, r(i).ktrratio, r(i).steady);
end

% ---------------------------------------------------------------- helpers
function v = getp(p, name)
    us = strfind(name, '__');
    if isempty(us)
        if isfield(p, name) && isnumeric(p.(name)) && isscalar(p.(name))
            v = p.(name);
        else
            v = [];
        end
    else
        f = name(1:us(1)-1); i = str2double(name(us(1)+2:end));
        if isfield(p, f) && numel(p.(f)) >= i; v = p.(f)(i); else; v = []; end
    end
end

function p = setp(p, name, val)
    us = strfind(name, '__');
    if isempty(us)
        p.(name) = val;
    else
        f = name(1:us(1)-1); i = str2double(name(us(1)+2:end));
        p.(f)(i) = val;
    end
end
