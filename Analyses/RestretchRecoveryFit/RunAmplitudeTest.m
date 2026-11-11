% RunAmplitudeTest.m — is the model's fast post-restretch recovery a
% LARGE-PERTURBATION (nonlinear) effect, or a STRUCTURAL one?
%
% THE ARGUMENT
% ------------
% Each of the three redevelopment windows is a relaxation back to the SAME
% isometric steady state, from a different displaced state. For a relaxation
% that is even approximately linear, the observed rate belongs to the SYSTEM
% (its slowest engaged eigenvalue); only the mode AMPLITUDES depend on where
% you start. The data behaves exactly like that: ktr, post-slack and
% post-restretch all give ~41-47 s^-1 despite amplitudes of 1.00, 0.74 and
% 0.22 F_iso. The model does not: it gives ~52, ~52 and ~118 s^-1.
%
% So either
%   (A) the re-stretch is a big enough perturbation to leave the linear
%       regime, and the model's fast rate is a nonlinear excursion  -- in
%       which case shrinking the perturbation must make rsK converge onto the
%       post-slack rate; or
%   (B) the two windows engage structurally different processes, and rsK
%       stays fast however small the re-stretch is.
%
% (A) and (B) call for completely different fixes, so this is worth one
% experiment. We build synthetic slack/re-stretch protocols of decreasing
% depth -- always releasing and re-stretching by the SAME amount, so the
% muscle returns to the same length -- and track rsK vs the post-slack rate
% as the excursion shrinks.
%
% Writes results/amplitude_test.mat + results/amplitude_test.png.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if ~exist(fullfile(here,'results'),'dir'); mkdir(fullfile(here,'results')); end

if isempty(gcp('nocreate')); parpool('local', 3); end

SNAP  = 'params/optfull7_opt_mov.m';
p0    = getParams(loadParams(SNAP), [], true, false);
p0.MaxRunTime = 300;
fcn   = str2func(p0.modelFcn);

% Depths as a fraction of the reference length (ML units). The real protocol's
% first cycle releases 0.080 ML; go down two decades from there.
DEPTHS = [0.080 0.060 0.040 0.030 0.020 0.015 0.010 0.007 0.005 0.003];

% A fit is only meaningful while the window still has an amplitude well above
% the noise/round-off floor. Below this the reported "rate" is estimator
% noise, not a property of the model, and must not be read as convergence.
MIN_AMP_KPA = 3.0;

V_REL  = -36;      % release velocity (ML/s), as in the real protocol cycle 1
V_STR  =  3.9;     % re-stretch velocity (ML/s)
T_ACT  = 2.0;      % settle time before the release
T_HOLD = 0.298;    % hold at short length (as in the protocol)
T_REC  = 0.400;    % recovery window after the re-stretch

res = struct('depth',{},'kPostSlack',{},'kPostRestretch',{},'ampSlack',{}, ...
             'ampRestretch',{},'F0',{});

for id = 1:numel(DEPTHS)
    d = DEPTHS(id);
    t_rel  = T_ACT;
    dt_rel = d / abs(V_REL);
    t_h0   = t_rel + dt_rel;
    t_str  = t_h0 + T_HOLD;
    dt_str = d / V_STR;
    t_h1   = t_str + dt_str;
    t_end  = t_h1 + T_REC;

    % velocitytable: [time, velocity, (unused), ML]
    vt = [ -2.0,   0,      0, 1.10; ...
           t_rel,  V_REL,  0, 1.10; ...
           t_h0,   0,      0, 1.10-d; ...
           t_str,  V_STR,  0, 1.10-d; ...
           t_h1,   0,      0, 1.10; ...
           t_end,  0,      0, 1.10];

    p = p0; p.Velocity = vt(:,2);
    p = getParams(p, p.g, true);
    fprintf('depth %.3f ML ... ', d);
    tt = tic;
    try
        [~, out] = evaluateModel(fcn, vt(:,1), p);
    catch ME
        fprintf('FAILED: %s\n', ME.message); continue;
    end

    nr = makeMonotonous(out.t);
    T  = out.t(nr); F = out.Force(nr);

    % post-slack window: from the end of the release ramp to the re-stretch
    kS = fitWin(T, F, t_h0, t_str, 0.0);
    % post-restretch window: from the end of the re-stretch ramp, anchored at
    % the minimum in the following 80 ms (same convention as rsK)
    kR = fitWin(T, F, t_h1, t_end, 0.080);

    n = numel(res)+1;
    res(n).depth          = d;
    res(n).kPostSlack     = kS.k;
    res(n).kPostRestretch = kR.k;
    res(n).ampSlack       = kS.A;
    res(n).ampRestretch   = kR.A;
    res(n).r2Slack        = kS.r2;
    res(n).r2Restretch    = kR.r2;
    res(n).F0             = interp1(T, F, t_rel-0.005);
    res(n).okSlack        = kS.A >= MIN_AMP_KPA;
    res(n).okRestretch    = kR.A >= MIN_AMP_KPA;

    flag = '';
    if ~res(n).okRestretch; flag = '  <-- restretch amp below floor'; end
    if ~res(n).okSlack;     flag = [flag '  <-- slack amp below floor']; end
    fprintf(['%.1f s | postSlack k %6.1f (A %6.2f) | postRestretch k %6.1f ' ...
             '(A %5.2f) | ratio %5.2f%s\n'], ...
        toc(tt), kS.k, kS.A, kR.k, kR.A, kR.k/kS.k, flag);
end

save(fullfile(here,'results','amplitude_test.mat'), 'res', 'DEPTHS', 'SNAP');

%% ---- report -------------------------------------------------------------
fprintf('\n%-9s %11s %8s %11s %8s %8s %7s\n', 'depth ML', 'k_postRstr', ...
        'A_restr', 'k_postSlack', 'A_slack', 'ratio', 'use?');
for n = 1:numel(res)
    use = 'yes'; if ~res(n).okRestretch; use = 'NOISE'; end
    fprintf('%-9.3f %11.1f %8.2f %11.1f %8.2f %8.2f %7s\n', res(n).depth, ...
        res(n).kPostRestretch, res(n).ampRestretch, ...
        res(n).kPostSlack, res(n).ampSlack, ...
        res(n).kPostRestretch/res(n).kPostSlack, use);
end

% Small-signal limit: the smallest depth whose amplitude is still above floor.
okr = [res.okRestretch];
if any(okr)
    idx  = find(okr, 1, 'last');
    kLin = res(idx).kPostRestretch;
    kBig = res(1).kPostRestretch;
    fprintf(['\nSMALL-SIGNAL LIMIT (smallest usable depth %.3f ML, A %.1f kPa):' ...
             ' k = %.1f s^-1\n'], res(idx).depth, res(idx).ampRestretch, kLin);
    fprintf('LARGE PERTURBATION (%.3f ML): k = %.1f s^-1\n', res(1).depth, kBig);
    fprintf('\nDECOMPOSITION of the %.2fx model/data gap at the real protocol depth:\n', kBig/43.7);
    fprintf('   linear   : model small-signal %.0f vs data %.0f          = x%.2f\n', kLin, 43.7, kLin/43.7);
    fprintf('   nonlinear: large-stretch speed-up %.0f -> %.0f            = x%.2f\n', kLin, kBig, kBig/kLin);
end

fprintf(['\nREADING: the post-restretch rate tends to the model''s SMALL-SIGNAL\n' ...
         'relaxation rate as the perturbation shrinks. Whatever that limit is,\n' ...
         'it is a property of the model at its operating point -- if it already\n' ...
         'exceeds the data''s ~44 s^-1, then part of the gap is a plain rate\n' ...
         'error (parameters can fix it) and the remainder is a large-stretch\n' ...
         'nonlinearity (they cannot).\n']);

fig = figure(702); clf; set(fig,'Position',[80 80 900 380],'Color','w');
subplot(1,2,1); hold on; box on;
plot([res.depth], [res.kPostSlack],     'o-', 'LineWidth',1.6, 'DisplayName','post-slack');
plot([res.depth], [res.kPostRestretch], 's-', 'LineWidth',1.6, 'DisplayName','post-restretch');
set(gca,'XScale','log'); xlabel('release/re-stretch depth (ML)'); ylabel('k (s^{-1})');
legend('Location','northwest'); title('rate vs perturbation size');
subplot(1,2,2); hold on; box on;
plot([res.depth], [res.kPostRestretch]./[res.kPostSlack], 'k^-','LineWidth',1.6);
yline(1,'r--','equal rates');
set(gca,'XScale','log'); xlabel('depth (ML)'); ylabel('k_{postRestretch} / k_{postSlack}');
title('do the two windows converge?');
exportgraphics(fig, fullfile(here,'results','amplitude_test.png'), 'Resolution',140);
fprintf('wrote results/amplitude_test.png\n');

% ---------------------------------------------------------------- helper
function s = fitWin(T, F, ta, tb, vallWin)
% First-order fit over [ta,tb] with the baseline anchored at the minimum in
% the first VALLWIN (0 = anchor at the window start) and A fixed at the
% window's plateau. Same estimator as Model/extractSlackAttributes' rsK.
    m = T >= ta & T <= tb;
    t = T(m); y = F(m);
    s = struct('k',NaN,'A',NaN,'t0',NaN,'r2',NaN);
    if numel(t) < 20; return; end
    if vallWin > 0
        nW = max(3, nnz(t <= t(1)+vallWin));
        [B, iB] = min(y(1:nW));
    else
        B = y(1); iB = 1;
    end
    tt = t(iB:end) - t(iB); yy = y(iB:end) - B;
    if numel(tt) < 20; return; end
    A = median(yy(tt >= tt(end) - 0.15*(tt(end)-tt(1))));
    if ~isfinite(A) || A <= 0; return; end
    obj = @(q) sum((yy - A.*(1-exp(-abs(q(1)).*max(tt-max(q(2),0),0)))).^2);
    q = fminsearch(obj, [45, 0.004], optimset('Display','off'));
    s.k = abs(q(1)); s.t0 = max(q(2),0); s.A = A;
    pred = A.*(1-exp(-s.k.*max(tt-s.t0,0)));
    s.r2 = 1 - sum((yy-pred).^2)/max(sum((yy-mean(yy)).^2),eps);
end
