%% ProbeRestretchDwell.m -- WHERE and WHEN do heads detach during a re-stretch?
%
% ProbeOverstretch says how much attached mass sits beyond a strain limit.
% This asks the question the physiology is about: at what strain does a head
% LET GO, which state is it in when it does, and how long does a head that has
% been dragged past the limit survive there. Everything is rebuilt per bin from
% the saved state vector with the RHS's own rate functions (2-state tree,
% UsePieceWiseStrainDep, fast-ppval path), and checked against the logged
% XB force so a wrong strain reconstruction cannot go unnoticed.
%
% Per re-stretch window (lengthening onset -> +40 ms), for the slack protocol
% and the ktr overstretch, and -- as the isolation control -- the isometric
% plateau and the shortening phases:
%   * stretch (um SL) to the model's first peak, next to the data's peak1_dSL
%   * detachment-flux-weighted strain, p1 and p2 separately (the "strain at
%     which heads let go"); share of detachment happening beyond 10/15/20 nm
%   * lifetime of heads beyond 15 nm = mass / outflux (p1 and p2)
%   * share of active force carried by heads beyond 15 nm at the peak
%
%   PARAMSET (default refit9_opt), OVR (struct of overrides, default none),
%   TAG (results suffix). Needs no pool.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
if ~exist('OVR', 'var'), OVR = struct(); end
if ~exist('TAG', 'var'), TAG = PARAMSET; end
if ~exist('RUNKTR', 'var'), RUNKTR = true; end
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..')));
res = fullfile(here, 'results');

P = getParams(loadParams(PARAMSET), [], true, false);
f = fieldnames(OVR); for k = 1:numel(f), P.(f{k}) = OVR.(f{k}); end
P.RunSlackSegments = 'All'; P.MaxRunTime = 900; P.PlotEachSeparately = 0; P.PlotFeatureFitting = 0;
P.ghostSave = ''; P.ghostLoad = ''; P.SaveBest = 0; P.mods = {}; P.g = ones(1,0);
P = getParams(P, P.g, true, false);

evalc('[~, oS, fmS, fdS] = runSlackExperiment(P);');
Q = struct();
Q.slack = binFluxes(oS, P);
if RUNKTR
    evalc('[~, oK, fmK] = runKtrExperiment(P, [37.7928 29.0 25.8033]);');
    Q.ktr = binFluxes(oK, P);
end

%% slack: re-stretch windows vs data
o = Q.slack; dSLdt = [0; diff(o.SL) ./ max(diff(o.t), eps)];
on = find(dSLdt(2:end) > 0.5 & dSLdt(1:end-1) <= 0.5) + 1;
on = on(o.t(on) > 70);
fprintf('\n[%s] check: max |F_rebuilt - FXB| / max FXB = %.2e\n', TAG, o.fErr);
fprintf('\n%-4s %9s %9s | %8s %8s | %8s %8s | %6s %6s %6s | %7s %7s | %7s\n', 'seg', 'dSLpk_m', 'dSLpk_d', ...
    'sDet_p1', 'sDet_p2', 'sAtt_pk', 'p1:p2', '>10nm', '>15nm', '>20nm', 'tau1_15', 'tau2_15', 'F>15@pk', 'F>15max');
W = zeros(numel(on), 14);
for k = 1:numel(on)
    i0 = on(k); t0 = o.t(i0);
    w  = find(o.t >= t0 & o.t <= t0 + 0.040);
    ip = firstPeak(o.Ftot(w)); ip = w(ip);
    wd = w(o.t(w) <= o.t(ip) + 0.010);                       % onset -> peak + 10 ms
    dt = [diff(o.t(wd)); 0];
    J1 = o.J1(wd, :) .* dt;  J2 = o.J2(wd, :) .* dt;  S = o.S(wd, :);
    sd1 = sum(J1 .* S, 'all') / sum(J1, 'all');
    sd2 = sum(J2 .* S, 'all') / sum(J2, 'all');
    Jt = J1 + J2;
    b = @(L) sum(Jt(S > L), 'all') / sum(Jt, 'all');
    % lifetime of heads beyond 15 nm, averaged over the window: mass / outflux
    m1 = sum(o.M1(wd, :) .* (S > 0.015) .* dt, 'all'); m2 = sum(o.M2(wd, :) .* (S > 0.015) .* dt, 'all');
    q1 = sum((J1 + o.J1h(wd, :) .* dt) .* (S > 0.015), 'all');   % rebinding also ends a head's stay beyond 15 nm
    q2 = sum((J2 + o.J21(wd, :) .* dt) .* (S > 0.015), 'all');
    sAtt = sum((o.M1(ip,:) + o.M2(ip,:)) .* o.S(ip,:)) / sum(o.M1(ip,:) + o.M2(ip,:));
    r12 = sum(o.M1(ip,:)) / max(sum(o.M2(ip,:)), eps);
    Fpk = sum(o.Fb(ip, o.S(ip,:) > 0.015)) / max(abs(sum(o.Fb(ip,:))), eps);
    Fsh = sum(o.Fb(w,:) .* (o.S(w,:) > 0.015), 2) ./ max(abs(sum(o.Fb(w,:), 2)), eps);   % share of XB force from >15 nm heads
    Fmax15 = max(Fsh);
    dSLm = o.SL(ip) - o.SL(i0);
    dSLd = NaN; if isfield(fdS, 'peak1_dSL') && numel(fdS.peak1_dSL) >= k, dSLd = fdS.peak1_dSL(k); end
    W(k, :) = [k dSLm dSLd sd1 sd2 sAtt r12 b(0.010) b(0.015) b(0.020) m1/max(q1,eps) m2/max(q2,eps) Fpk Fmax15];
    fprintf('%-4d %9.4f %9.4f | %6.1fnm %6.1fnm | %6.1fnm %8.2f | %5.0f%% %5.0f%% %5.0f%% | %5.2fms %5.2fms | %6.0f%% %6.0f%%\n', ...
        k, dSLm, dSLd, 1e3*sd1, 1e3*sd2, 1e3*sAtt, r12, 100*W(k,8), 100*W(k,9), 100*W(k,10), 1e3*W(k,11), 1e3*W(k,12), 100*Fpk, 100*Fmax15);
end
Q.W = W;
Q.Wcols = {'seg','dSLpk_model','dSLpk_data','sDet_p1','sDet_p2','sAtt_pk','p1_over_p2','detFrac10','detFrac15','detFrac20','tau1_15','tau2_15','Fshare15_pk','Fshare15_max40ms'};
fprintf('model peak1_dSL (feature, vs SLslack): %s | data %s\n', mat2str(fmS.peak1_dSL(:)',3), mat2str(fdS.peak1_dSL(:)',3));

%% isolation: which protocol phases put heads beyond 10/15 nm at all?
ph = phaseMasks(o, on);
fprintf('\n%-26s %10s %10s %10s %10s\n', 'slack phase (s>0 side)', 'max >10nm', 'max >15nm', 'mean >15nm', 'max <-15nm');
for q = fieldnames(ph)'
    m = ph.(q{1});
    fprintf('%-26s %9.1f%% %9.1f%% %9.2f%% %9.1f%%\n', q{1}, 100*max(o.frac10(m)), 100*max(o.frac15(m)), 100*mean(o.frac15(m)), 100*max(o.fneg15(m)));
end
if RUNKTR
    k = Q.ktr; dK = [0; diff(k.SL) ./ max(diff(k.t), eps)];
    fprintf('%-26s %9.1f%% %9.1f%%\n', 'ktr: whole protocol', 100*max(k.frac10), 100*max(k.frac15));
    fprintf('%-26s %9.1f%% %9.1f%%\n', 'ktr: lengthening', 100*max(k.frac10(dK > 0.05)), 100*max(k.frac15(dK > 0.05)));
end
Q.fmS = fmS; Q.fdS = fdS; Q.P = P; Q.OVR = OVR;
Q.slack = rmfield(Q.slack, {'J1','J2','J21','J1h','M1','M2','S','Fb'});   % keep the file small
if RUNKTR, Q.ktr = rmfield(Q.ktr, {'J1','J2','J21','J1h','M1','M2','S','Fb'}); end
save(fullfile(res, ['dwell_' TAG '.mat']), 'Q');

function B = binFluxes(o, P)
% Per time x bin: strain, masses, detachment flux densities, force. 2-state,
% UsePieceWiseStrainDep tree of dPUdT_CombinedTransitions (fast-ppval path).
assert(P.NumberOfStates == 2 && P.UsePieceWiseStrainDep, 'probe covers the 2-state piecewise tree only');
ss = P.ss; dS = P.dS; n = numel(o.t); cl = P.ClampStrainDepExtrapolation;
kx = @(pp, x) ppval(pp, min(max(x, pp.breaks(1)), pp.breaks(end)) .* cl + x .* ~cl);
B = struct('t', o.t(:), 'SL', o.SL(:), 'Ftot', o.Force(:));
[B.S, B.M1, B.M2, B.J1, B.J2, B.J21, B.Fb, B.J1h] = deal(zeros(n, ss));
B.frac10 = zeros(n,1); B.frac15 = zeros(n,1); B.fneg15 = zeros(n,1); fx = zeros(n,1);
for i = 1:n
    s  = P.s(:) - (-(o.SL(i) - o.LSE(i)) + o.LXBPivot(i))/2;
    p1 = max(0, o.PU(i, 1:ss)');  p2 = max(0, o.PU(i, ss+1:2*ss)');
    r1D = P.kd  * kx(P.PieceWiseStrainDepR1D, s);
    r21 = P.k_1 * kx(P.PieceWiseStrainDepR21, s);
    r2  = P.k2  * kx(P.PieceWiseStrainDep2, s + P.dr2 - P.dr);
    rip = P.k2rip * min(1e9, max(0, exp(P.alphaRip*(s + P.dr3))));
    r1h = zeros(size(s));                                              % yield-gate rebinding (leaves the bin, not the filament)
    if isfield(P, 'UseYieldGate') && P.UseYieldGate
        g = 1 ./ (1 + exp(-(s - P.s_yield) / P.w_yield));
        hop = 0; if isfield(P, 'YieldHop1'), hop = min(1, max(0, P.YieldHop1)); end   % rebinding is not detachment
        r1D = r1D + (1 - hop) * P.k_yield1 * g;  rip = rip + P.k_yield2 * g;  r1h = hop * P.k_yield1 * g;
    end
    s1 = s + P.dr1;  x = s + P.dr;  e = sign(x).*abs(x).^P.estiff;
    fb = (P.kstiff1*(s1 >= 0) + P.kstiff1_n*(s1 < 0)).*s1.*p1 + (P.kstiff2*(s >= -P.dr) + P.kstiff2_n*(s < -P.dr)).*e.*p2;
    B.S(i,:) = s; B.M1(i,:) = p1*dS; B.M2(i,:) = p2*dS;
    B.J1(i,:) = r1D.*p1*dS;  B.J1h(i,:) = r1h.*p1*dS;  B.J2(i,:) = (r2 + rip).*p2*dS;  B.J21(i,:) = r21.*p2*dS;  B.Fb(i,:) = fb*dS;
    m = p1 + p2; A = max(sum(m), eps);
    B.frac10(i) = sum(m(s > 0.010))/A;  B.frac15(i) = sum(m(s > 0.015))/A;   % STRETCHED side only
    B.fneg15(i) = sum(m(s < -0.015))/A;
    fx(i) = sum(fb)*dS;
end
B.fErr = max(abs(fx(:) - o.FXB(:))) / max(abs(o.FXB));
end

function ip = firstPeak(y)
% first local maximum whose prominence clears 5 % of the window's range
% (the extractSlackAttributes rule for a smooth trace)
[~, loc] = findpeaks(y(:), 'MinPeakProminence', 0.05*range(y));
if isempty(loc), [~, ip] = max(y); else, ip = loc(1); end
end

function ph = phaseMasks(o, on)
% Isometric plateau (last 50 ms before each release), shortening (slack),
% the re-stretch itself, and the post-restretch hold.
dSL = [0; diff(o.SL) ./ max(diff(o.t), eps)];
ph.isometric_prerelease = false(size(o.t)); ph.shortening = dSL < -0.05 & o.t > 70;
ph.restretch_lengthening = dSL > 0.05 & o.t > 70; ph.post_restretch_hold = false(size(o.t));
rel = find(dSL(2:end) < -0.05 & dSL(1:end-1) >= -0.05) + 1; rel = rel(o.t(rel) > 70);
for k = 1:numel(rel), ph.isometric_prerelease = ph.isometric_prerelease | (o.t > o.t(rel(k)) - 0.050 & o.t < o.t(rel(k)) - 0.001); end
for k = 1:numel(on)
    tEnd = o.t(on(k)); j = find(o.t > tEnd & abs(dSL) < 0.05, 1);
    if ~isempty(j), ph.post_restretch_hold = ph.post_restretch_hold | (o.t >= o.t(j) & o.t < o.t(j) + 0.100); end
end
end
