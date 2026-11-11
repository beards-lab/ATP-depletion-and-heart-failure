%% ProbeOverstretch.m -- how many attached heads sit beyond the physiological strain range, and when?
%
% One serial slack run and one ktr-protocol run on PARAMSET (its own numerics),
% then, for every output time, the strain of every bin is rebuilt from
% (SL, LSE, LXBPivot) and the attached mass and active force are split by
% |strain|. Force per bin follows the RHS exactly (UseNegativeKstiff split,
% estiff exponent on p2). Reported: the share of attached heads and of active
% force beyond LIM (default 10/15/20 nm), overall and in each re-stretch window.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
if ~exist('LIM', 'var'), LIM = [0.010 0.015 0.020]; end
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..')));
res = fullfile(here, 'results');

P = getParams(loadParams(PARAMSET), [], true, false);
P.RunSlackSegments = 'All'; P.MaxRunTime = 900; P.PlotEachSeparately = 0; P.PlotFeatureFitting = 0;
P.ghostSave = ''; P.ghostLoad = ''; P.SaveBest = 0; P.mods = {}; P.g = ones(1,0);
P = getParams(P, P.g, true, false);

evalc('[~, oS, fmS] = runSlackExperiment(P);');
Ktr_mean = [37.7928 29.0 25.8033];
evalc('[~, oK, fmK] = runKtrExperiment(P, Ktr_mean);');

S = struct();
for pr = {'slack', 'ktr'}
    if strcmp(pr{1}, 'slack'), o = oS; else, o = oK; end
    [A, Fa, frac, ffrac, smax, smin] = splitByStrain(o, P, LIM);
    S.(pr{1}) = struct('t', o.t(:), 'SL', o.SL(:), 'A', A, 'F', Fa, 'frac', frac, 'ffrac', ffrac, 'smax', smax, 'smin', smin, 'Ftot', o.Force(:));
    fprintf('\n=== %s protocol (%d time points) ===\n', pr{1}, numel(o.t));
    for j = 1:numel(LIM)
        [mf, im] = max(frac(:, j)); [mF, iF] = max(abs(ffrac(:, j)));
        fprintf('  beyond +-%2.0f nm: peak share of attached heads %5.1f %% (t %.3f s, attached %.3f) | peak share of active force %5.1f %% (t %.3f s)\n', ...
            1e3*LIM(j), 100*mf, o.t(im), A(im), 100*mF, o.t(iF));
    end
    fprintf('  occupied strain range (density > 1e-3 of peak): %.0f .. %.0f nm\n', 1e3*min(smin), 1e3*max(smax));
end

% re-stretch windows of the slack protocol: from each lengthening onset to +150 ms
o = S.slack; dSLdt = [0; diff(o.SL) ./ max(diff(o.t), eps)];
on = find(dSLdt(2:end) > 0.5 & dSLdt(1:end-1) <= 0.5) + 1;      % lengthening onsets
on = on(o.t(on) > 70);
fprintf('\n=== slack re-stretches: share beyond +-15 nm (j = 2) ===\n');
fprintf('%4s %8s %10s %12s %12s %14s\n', 'seg', 't_on', 'peak F', 'heads @peakF', 'force @peakF', 'max heads 0-150ms');
W = [];
for k = 1:numel(on)
    w = o.t >= o.t(on(k)) & o.t <= o.t(on(k)) + 0.150;
    iw = find(w); [Fp, ip] = max(o.Ftot(iw)); ip = iw(ip);
    W(k, :) = [o.t(on(k)) Fp o.frac(ip, 2) o.ffrac(ip, 2) max(o.frac(iw, 2))]; %#ok<SAGROW>
    fprintf('%4d %8.3f %10.2f %11.1f%% %11.1f%% %13.1f%%\n', k, W(k,1), Fp, 100*W(k,3), 100*W(k,4), 100*W(k,5));
end
S.restretch = W; S.LIM = LIM; S.fmS = fmS; S.fmK = fmK;
save(fullfile(res, ['overstretch_' PARAMSET '.mat']), 'S');

function [A, Fa, frac, ffrac, smax, smin] = splitByStrain(o, P, LIM)
ss = P.ss; dS = P.dS; n = numel(o.t);
A = zeros(n,1); Fa = zeros(n,1); frac = zeros(n, numel(LIM)); ffrac = frac; smax = zeros(n,1); smin = smax;
for i = 1:n
    s  = P.s(:) - (-(o.SL(i) - o.LSE(i)) + o.LXBPivot(i))/2;
    p1 = max(0, o.PU(i, 1:ss)');  p2 = max(0, o.PU(i, ss+1:2*ss)');
    s1 = s + P.dr1;  x = s + P.dr;  e = sign(x).*abs(x).^P.estiff;
    f  = (P.kstiff1*(s1 >= 0) + P.kstiff1_n*(s1 < 0)).*s1.*p1 + (P.kstiff2*(s >= -P.dr) + P.kstiff2_n*(s < -P.dr)).*e.*p2;
    m  = p1 + p2;
    A(i) = dS*sum(m);  Fa(i) = dS*sum(f);
    for j = 1:numel(LIM)
        out = abs(s) > LIM(j);
        frac(i,j)  = dS*sum(m(out)) / max(A(i), eps);
        ffrac(i,j) = dS*sum(f(out)) / max(abs(Fa(i)), eps);
    end
    occ = m > 1e-3*max(m);
    if any(occ), smax(i) = max(s(occ)); smin(i) = min(s(occ)); end
end
end
