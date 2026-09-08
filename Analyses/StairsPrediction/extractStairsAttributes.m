function [features, seg] = extractStairsAttributes(t, F, velocitytable, prefix, features, isData)
%EXTRACTSTAIRSATTRIBUTES  Features of a staircase (step-stretch/hold) protocol.
%
%   [FEATURES, SEG] = extractStairsAttributes(T, F, VELOCITYTABLE, PREFIX,
%                                             FEATURES, ISDATA)
%
%   The staircase is a train of identical small ramp-stretches separated by
%   isometric holds, closed by one larger stretch. Each step produces the
%   classic three-phase response:
%
%       base --(ramp)--> PEAK --(fast decay)--> NADIR --(slow rise)--> PLATEAU
%
%   Every amplitude below is a WINDOW MEAN, never an argmax/argmin. On the
%   4 kHz experimental trace the point-to-point noise is sigma ~ 0.023 F0, so a
%   max() over a ramp window is biased upward by ~1 sigma (~15 % of the signal)
%   while a 25-sample mean has SE ~ 0.005 F0 and no bias. The same windows are
%   applied to the (smooth) model output, so data and model are measured with
%   one instrument.
%
%   Windows, all relative to the ramp of step k (Tr = ramp duration):
%     base    : [t_r0 - 20 ms, t_r0 -  2 ms]      pre-step isometric level
%     peak    : [t_r0 + 0.55*Tr, t_r0 + 0.85*Tr]  late-ramp elastic plateau
%     nadir   : [t_r1 +  8 ms,   t_r1 + 18 ms]    phase-2 undershoot
%     plateau : [t_h1 - 20 ms,   t_h1 -  2 ms]    end-of-hold level
%
%   Features written (PREFIX prepended, default 'st_'), all normalised to the
%   pre-staircase isometric force F0 = base of step 1:
%     st_peak   (1xK) peak_k  - base_k   elastic/T1 amplitude of each step
%     st_nadir  (1xK) nadir_k - base_k   phase-2 undershoot (negative = below base)
%     st_plat   (1xK) plat_k  - F0       cumulative staircase envelope vs F0
%     st_rec    (1xK) plat_k  - nadir_k  delayed force recovery within the hold
%     st_SL     (1xK) SL at the end of step k (the x-axis of the envelope)
%     st_peakSlope    slope of st_peak over the EQUAL-SIZED steps, per step
%                     (>0 = strain stiffening, <0 = strain softening)
%     st_peakFirst    st_peak(1)         amplitude scale, first step
%     st_platEnd      st_plat(end)       total force gain over the staircase
%     st_platCurv     curvature of the equal-step envelope (2nd difference mean)
%
%   ISDATA (default false) applies a Savitzky-Golay pre-filter; it affects only
%   the diagnostic peak location, not the window means.
%
%   See also: runStairsExperiment, evalFeatureCost, extractSlackAttributes

BASE_PRE_S   = 0.020;   % base window length before ramp onset
GUARD_S      = 0.002;   % guard band kept clear of every transition
PEAK_FRAC    = [0.55 0.85];  % fraction of the ramp used for the peak mean
NADIR_LAG_S  = [0.008 0.018];% lag after ramp end for the phase-2 mean
PLAT_PRE_S   = 0.020;   % plateau window length before the next ramp

if nargin < 4 || isempty(prefix);   prefix   = 'st_';  end
if nargin < 5;                      features = struct(); end
if nargin < 6 || isempty(isData);   isData   = false;  end

t = t(:); F = F(:);
if isData && numel(F) > 41
    F = sgolayfilt(F, 2, 21);   % ~5 ms at 4 kHz: kills point noise, keeps the ramp
end

% ---- segment the protocol from the velocity table -----------------------
% Rows with v ~= 0 open a ramp; the ramp ends at the next row; the hold ends
% at the row after that (or at the end of the record).
vt   = velocitytable;
iRmp = find(vt(:,2) ~= 0);
iRmp = iRmp(iRmp < size(vt,1));
K    = numel(iRmp);

seg = struct('t_r0',[],'t_r1',[],'t_h1',[],'dSL',[]);
for k = 1:K
    i = iRmp(k);
    seg(k).t_r0 = vt(i,1);
    seg(k).t_r1 = vt(i+1,1);
    if i+2 <= size(vt,1)
        seg(k).t_h1 = vt(i+2,1);
    else
        seg(k).t_h1 = t(end);
    end
    seg(k).dSL = vt(i,2) * (seg(k).t_r1 - seg(k).t_r0);   % ML units
end

wmean = @(a,b) iWinMean(t, F, a, b);

base = nan(1,K); peak = nan(1,K); nad = nan(1,K); plat = nan(1,K);
for k = 1:K
    Tr = seg(k).t_r1 - seg(k).t_r0;
    base(k) = wmean(seg(k).t_r0 - BASE_PRE_S, seg(k).t_r0 - GUARD_S);
    peak(k) = wmean(seg(k).t_r0 + PEAK_FRAC(1)*Tr, seg(k).t_r0 + PEAK_FRAC(2)*Tr);
    nad(k)  = wmean(seg(k).t_r1 + NADIR_LAG_S(1), ...
                    min(seg(k).t_r1 + NADIR_LAG_S(2), seg(k).t_h1 - GUARD_S));
    plat(k) = wmean(max(seg(k).t_h1 - PLAT_PRE_S, seg(k).t_r1 + GUARD_S), ...
                    seg(k).t_h1 - GUARD_S);
end

F0 = base(1);
if ~isfinite(F0) || F0 == 0; F0 = 1; end

st_peak  = (peak - base) / F0;
st_nadir = (nad  - base) / F0;
st_plat  = (plat - F0)   / F0;
st_rec   = (plat - nad)  / F0;

% ---- which steps are the equal-sized ones? ------------------------------
% The staircase closes with one larger stretch; keep it out of every trend.
dSL   = [seg.dSL];
isEq  = abs(dSL - median(dSL)) < 0.25*abs(median(dSL));
kEq   = find(isEq);

pk = polyfit(kEq(:), st_peak(kEq)', 1);
pl = polyfit(kEq(:), st_plat(kEq)', 1);

features.([prefix 'peak'])   = st_peak;
features.([prefix 'nadir'])  = st_nadir;
features.([prefix 'plat'])   = st_plat;
features.([prefix 'rec'])    = st_rec;
features.([prefix 'dSL'])    = dSL;
features.([prefix 'peakSlope'])  = pk(1);
features.([prefix 'peakFirst'])  = st_peak(kEq(1));
features.([prefix 'peakMean'])   = mean(st_peak(kEq));
features.([prefix 'nadirMean'])  = mean(st_nadir(kEq));
features.([prefix 'platEnd'])    = st_plat(end);
features.([prefix 'platSlope'])  = pl(1);
features.([prefix 'platCurv'])   = mean(diff(st_plat(kEq), 2));
features.([prefix 'peakBig'])    = st_peak(end);
end

% ------------------------------------------------------------------------
function m = iWinMean(t, F, a, b)
%IWINMEAN Mean of F over t in [a,b]; falls back to interpolation if the
%   window holds fewer than two samples (coarse ODE output).
m = NaN;
if ~(b > a); return; end
sel = t >= a & t <= b;
if nnz(sel) >= 2
    m = mean(F(sel));
elseif numel(t) >= 2
    tq = linspace(a, b, 25);
    m  = mean(interp1(t, F, tq, 'linear', 'extrap'));
end
end
