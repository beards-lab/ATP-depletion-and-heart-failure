% RunAmpNoMaxwell.m — the residual after the Maxwell transient is removed:
% is it a LEVEL error (uniformly too fast -> rate constants can fix it) or an
% AMPLITUDE-DEPENDENCE error (needs new mechanism)?
%
% Synthetic slack/re-stretch protocols of varying depth, same estimator as
% rsK, run with and without the Maxwell dashpot.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if isempty(gcp('nocreate')); parpool('local', 4); end

SNAP   = 'params/rskR2_w025_opt.m';
DEPTHS = [0.080 0.060 0.040 0.030 0.020 0.015 0.010];
MIN_AMP = 3.0;

p0 = getParams(loadParams(SNAP), [], true, false);
p0.MaxRunTime = 300;
fcn = str2func(p0.modelFcn);

V_REL = -36; V_STR = 3.9; T_ACT = 2.0; T_HOLD = 0.298; T_REC = 0.400;

fprintf('%-9s | %-22s | %-22s\n', '', 'Maxwell ON', 'Maxwell OFF');
fprintf('%-9s %10s %10s %10s %10s\n', 'depth', 'rsK', 'A(kPa)', 'rsK', 'A(kPa)');
K = nan(numel(DEPTHS),2); A = nan(numel(DEPTHS),2);

for id = 1:numel(DEPTHS)
    d = DEPTHS(id);
    t_rel = T_ACT; dt_rel = d/abs(V_REL); t_h0 = t_rel+dt_rel;
    t_str = t_h0+T_HOLD; dt_str = d/V_STR; t_h1 = t_str+dt_str; t_end = t_h1+T_REC;
    vt = [-2.0,0,0,1.10; t_rel,V_REL,0,1.10; t_h0,0,0,1.10-d; ...
          t_str,V_STR,0,1.10-d; t_h1,0,0,1.10; t_end,0,0,1.10];

    fprintf('%-9.3f', d);
    for im = 1:2
        p = p0; p.UseMaxwellDashpot = (im==1);
        p.Velocity = vt(:,2); p = getParams(p, p.g, true);
        try
            [~, out] = evaluateModel(fcn, vt(:,1), p);
        catch
            fprintf('%10s%10s','FAIL',''); continue;
        end
        nr = makeMonotonous(out.t); T = out.t(nr); F = out.Force(nr);
        s = fitWin(T, F, t_h1, t_end, 0.080);
        K(id,im) = s.k; A(id,im) = s.A;
        fprintf('%10.1f%10.2f', s.k, s.A);
    end
    fprintf('\n');
end

fprintf('\n%-14s %10s %10s %10s\n','','slope','intercept','r');
for im = 1:2
    if im==1; lbl='Maxwell ON'; else; lbl='Maxwell OFF'; end
    ok = A(:,im) >= MIN_AMP & isfinite(K(:,im));
    if nnz(ok) < 3; continue; end
    pp = polyfit(DEPTHS(ok)', K(ok,im), 1); rr = corr(DEPTHS(ok)', K(ok,im));
    fprintf('%-14s %10.0f %10.1f %10.2f\n', lbl, pp(1), pp(2), rr);
end
fprintf('%-14s %10.0f %10.1f %10s\n','DATA (protocol)', -94, 53.2, '-');

fprintf(['\nREAD: the DATA slope is -94 s^-1/ML at a level of ~44 s^-1. If\n' ...
         '"Maxwell OFF" is FLAT and merely too high, the residual is a level\n' ...
         'error and rate constants can close it. If it still RISES steeply with\n' ...
         'amplitude, a stretch-engaged mechanism is genuinely missing.\n']);

save(fullfile(here,'results','amp_no_maxwell.mat'), 'K','A','DEPTHS','SNAP');

function s = fitWin(T, F, ta, tb, vallWin)
    m = T >= ta & T <= tb; t = T(m); y = F(m);
    s = struct('k',NaN,'A',NaN);
    if numel(t) < 20; return; end
    nW = max(3, nnz(t <= t(1)+vallWin)); [B,iB] = min(y(1:nW));
    tt = t(iB:end)-t(iB); yy = y(iB:end)-B;
    if numel(tt) < 20; return; end
    Aa = median(yy(tt >= tt(end)-0.15*(tt(end)-tt(1))));
    if ~isfinite(Aa) || Aa <= 0; return; end
    obj = @(q) sum((yy - Aa.*(1-exp(-abs(q(1)).*max(tt-max(q(2),0),0)))).^2);
    q = fminsearch(obj, [45,0.004], optimset('Display','off'));
    s.k = abs(q(1)); s.A = Aa;
end
