function w = oscWindow(day, cond, type, cyc, t, L, F, ta, tb, vRamp, fs)
%OSCWINDOW  One recovery window: exponential fit, residual, residual spectrum.
%
%   The exponential is the SAME estimator the objective uses
%   (Model/extractSlackAttributes/fitRestretchRecovery): baseline B = the
%   minimum in the first 80 ms, amplitude A = end-of-window plateau minus B
%   held FIXED, only k and t0 free. What is left over is the ripple.
%
%   Returns [] when the window is too short to say anything.

w = [];
m = t >= ta & t <= tb;
if nnz(m) < 40; return; end
tt = t(m) - ta;  yy = F(m);  ll = L(m);

nMed = max(3, round(0.002*fs));
mv = tt <= 0.080;
if nnz(mv) < 3; mv = true(size(tt)); end
[B, iB] = min(movmedian(yy(mv), nMed));
tt = tt(iB:end) - tt(iB);  yy = yy(iB:end);  ll = ll(iB:end);
if numel(tt) < 40 || tt(end) <= 0; return; end

P = median(yy(tt >= tt(end) - 0.15*tt(end)));
A = P - B;
z = yy - B;

if isfinite(A) && A > 0
    obj  = @(q) sum((z - A.*(1-exp(-abs(q(1)).*max(tt-max(q(2),0),0)))).^2);
    q    = fminsearch(obj, [45, 0.004], optimset('Display','off','MaxFunEvals',4e3));
    k    = abs(q(1));  t0 = max(q(2),0);
    pred = A.*(1-exp(-k.*max(tt-t0,0)));
else
    % PNB / passive: no active rise. The "recovery" is a flat line, so the
    % ripple is measured against a slow trend instead of against a rise.
    A = NaN; k = NaN; t0 = NaN;
    pred = movmean(z, max(3, round(0.050*fs)));
end
res = z - pred;
r2  = 1 - sum(res.^2)/max(sum((z-mean(z)).^2), eps);

[f0, pk, frac, fax, pxx] = rippleSpectrum(res, fs);
resL = ll - movmean(ll, max(3, round(0.050*fs)));
[f0L, pkL, ~, ~, pxxL]   = rippleSpectrum(resL, fs);

% amplitude scale for the ripple: against the RISE when there is one, against
% the local force level when there is not (PNB).
scale = A; if ~isfinite(scale) || scale <= 0; scale = max(median(yy), eps); end

w = struct('day',day, 'cond',cond, 'type',type, 'cyc',cyc, ...
    't',tt, 'z',z, 'pred',pred, 'res',res, 'fs',fs, 'L',ll, ...
    'A',A, 'B',B, 'k',k, 't0',t0, 'r2',r2, 'vRamp',vRamp, 'Fmed',median(yy), ...
    'f0',f0, 'pkAmp',pk, 'bandFrac',frac, 'fax',fax, 'pxx',pxx, ...
    'f0L',f0L, 'pkAmpL',pkL, 'pxxL',pxxL, ...
    'resRMS',rms(res), 'relRMS',rms(res)/scale, 'span',tt(end));
end
