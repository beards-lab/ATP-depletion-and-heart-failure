function [f0, pkAmp, bandFrac, fax, pxx] = rippleSpectrum(res, fs, band)
%RIPPLESPECTRUM  Dominant frequency of a recovery residual.
%
%   [F0, PKAMP, BANDFRAC, FAX, PXX] = rippleSpectrum(RES, FS)
%
%   Welch PSD of RES (Hann, 50 % overlap, at most 4 segments so a 280 ms
%   window still resolves ~15 Hz). Searched above 8 Hz so the exponential's
%   own residual curvature -- which lives at 1/span -- cannot win.
%
%   F0        frequency of the largest peak (Hz)
%   PKAMP     equivalent sinusoid amplitude at F0, in the units of RES
%   BANDFRAC  fraction of the residual variance inside F0 +/- 25 %
%   FAX, PXX  the periodogram itself

if nargin < 3 || isempty(band); band = [8 Inf]; end
f0 = NaN; pkAmp = NaN; bandFrac = NaN; fax = []; pxx = [];
res = res(:); res = res(isfinite(res));
n = numel(res); if n < 32; return; end

nseg = max(64, 2^floor(log2(n/2)));           % >= 2 averages when possible
nseg = min(nseg, n);
[pxx, fax] = pwelch(res - mean(res), hann(nseg), round(nseg/2), max(256,nseg), fs);

hi = min(band(2), 0.45*fs);
sel = fax >= band(1) & fax <= hi;
if ~any(sel); return; end
[~, j] = max(pxx(sel));
fsel = fax(sel);  f0 = fsel(j);

df = median(diff(fax));
inb = fax >= 0.75*f0 & fax <= 1.25*f0;
Pin = sum(pxx(inb))*df;
bandFrac = Pin / max(sum(pxx)*df, eps);
pkAmp = sqrt(2*Pin);                          % sinusoid amplitude = sqrt(2*P)
end
