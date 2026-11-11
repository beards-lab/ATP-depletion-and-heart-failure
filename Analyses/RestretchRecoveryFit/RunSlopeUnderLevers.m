% RunSlopeUnderLevers.m — do the best rsK levers move the amplitude SLOPE, or
% only the intercept?
%
% The gap splits into ×1.8 "the model's small-signal relaxation is too fast"
% and ×1.7 "it speeds up further when stretched hard" — and only the second is
% structural (the data's amplitude slope is NEGATIVE: -126 s^-1/ML, while the
% model's is +714). A lever that lowers rsK uniformly buys the first half and
% leaves the second untouched.
%
% The prediction for the Maxwell levers (eta_M, kSE_M) is that they do exactly
% that, because a linear series filter delays force transmission by the same
% factor whatever the amplitude. Asserted in conclusions.md — verified here.
%
% Writes results/slope_under_levers.mat + .png

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if isempty(gcp('nocreate')); parpool('local', 4); end

SNAP   = 'params/optfull7_opt_mov.m';
DEPTHS = [0.080 0.060 0.040 0.030 0.020 0.015 0.010];
% Each row: {label, {param,factor; ...}}. The last rows are COMBINATIONS —
% the single-lever runs showed eta_M flattening the SLOPE while k2 lowers the
% LEVEL, i.e. they address the two halves of the gap separately, so the
% obvious question is what they do together.
VAR = { 'baseline',            {}; ...
        'eta_M x0.6',          {'eta_M',0.6}; ...
        'eta_M x0.3',          {'eta_M',0.3}; ...
        'kSE_M x0.6',          {'kSE_M',0.6}; ...
        'k2 x0.6',             {'k2',0.6}; ...
        'k2 .6 + eta_M .3',    {'k2',0.6,'eta_M',0.3}; ...
        'k2 .7 + eta_M .3',    {'k2',0.7,'eta_M',0.3}; ...
        'k2 .6+etaM .3+kSEM .6', {'k2',0.6,'eta_M',0.3,'kSE_M',0.6} };

p0 = getParams(loadParams(SNAP), [], true, false);
p0.MaxRunTime = 300;
fcn = str2func(p0.modelFcn);
V_REL = -36; V_STR = 3.9; T_ACT = 2.0; T_HOLD = 0.298; T_REC = 0.400;

S = struct('name',{},'depth',{},'k',{},'A',{});
for iv = 1:size(VAR,1)
    nm = VAR{iv,1}; spec = VAR{iv,2};
    fprintf('\n--- %s ---\n', nm);
    for d = DEPTHS
        t_rel=T_ACT; t_h0=t_rel+d/abs(V_REL); t_str=t_h0+T_HOLD;
        t_h1=t_str+d/V_STR; t_end=t_h1+T_REC;
        vt = [-2,0,0,1.10; t_rel,V_REL,0,1.10; t_h0,0,0,1.10-d; ...
              t_str,V_STR,0,1.10-d; t_h1,0,0,1.10; t_end,0,0,1.10];
        p = p0;
        for ip = 1:2:numel(spec); p.(spec{ip}) = p.(spec{ip}) * spec{ip+1}; end
        p.Velocity = vt(:,2); p = getParams(p, p.g, true);
        try
            [~, out] = evaluateModel(fcn, vt(:,1), p);
        catch ME
            fprintf('  %.3f FAILED: %s\n', d, ME.message); continue;
        end
        nr = makeMonotonous(out.t);
        s  = fitWin(out.t(nr), out.Force(nr), t_h1, t_end, 0.080);
        n = numel(S)+1; S(n).name=nm; S(n).depth=d; S(n).k=s.k; S(n).A=s.A;
        fprintf('  depth %.3f  k %6.1f  A %6.2f\n', d, s.k, s.A);
    end
end
save(fullfile(here,'results','slope_under_levers.mat'), 'S', 'VAR', 'DEPTHS', 'SNAP');

fprintf('\n%-14s %10s %8s %12s %12s\n','variant','slope','r','k@0.080','k@0.010');
fig = figure(706); clf; set(fig,'Position',[70 70 760 460],'Color','w'); hold on; box on;
co = lines(size(VAR,1));
for iv = 1:size(VAR,1)
    nm = VAR{iv,1};
    m = strcmp({S.name}, nm);
    a = [S(m).depth]'; k = [S(m).k]';
    ok = isfinite(a) & isfinite(k);
    if nnz(ok) < 3; continue; end
    pf = polyfit(a(ok),k(ok),1); rr = corr(a(ok),k(ok));
    kb = k(a==0.080); ks = k(a==0.010);
    fprintf('%-14s %10.0f %8.2f %12.1f %12.1f\n', nm, pf(1), rr, ...
            ternary(isempty(kb),NaN,kb), ternary(isempty(ks),NaN,ks));
    plot(a, k, 'o-', 'Color', co(iv,:), 'LineWidth',1.6, 'DisplayName', ...
         sprintf('%s (slope %+.0f)', nm, pf(1)));
end
fprintf('%-14s %10.0f %8.2f %12s %12s\n','DATA',-126,-0.31,'~52','n/a');
yline(43.7,'k--','data ~44 s^{-1}','LineWidth',1.2);
xlabel('re-stretch amplitude (ML)'); ylabel('rsK (s^{-1})');
legend('Location','northwest','FontSize',8);
title('Do the best rsK levers flatten the slope, or just shift it down?');
exportgraphics(fig, fullfile(here,'results','slope_under_levers.png'),'Resolution',140);
fprintf('\nwrote results/slope_under_levers.png\n');

function v = ternary(c,a,b); if c; v=a; else; v=b; end; end

function s = fitWin(T, F, ta, tb, vallWin)
    m = T >= ta & T <= tb; t = T(m); y = F(m);
    s = struct('k',NaN,'A',NaN);
    if numel(t) < 20; return; end
    nW = max(3, nnz(t <= t(1)+vallWin)); [B, iB] = min(y(1:nW));
    tt = t(iB:end)-t(iB); yy = y(iB:end)-B;
    if numel(tt) < 20; return; end
    A = median(yy(tt >= tt(end)-0.15*(tt(end)-tt(1))));
    if ~isfinite(A) || A <= 0; return; end
    obj = @(q) sum((yy - A.*(1-exp(-abs(q(1)).*max(tt-max(q(2),0),0)))).^2);
    q = fminsearch(obj, [45, 0.004], optimset('Display','off'));
    s.k = abs(q(1)); s.A = A;
end
