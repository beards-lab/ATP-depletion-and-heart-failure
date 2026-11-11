function F = robustSlackFeatures(t, y, SL, vt, useSmoothing, tL, LSE)
%ROBUSTSLACKFEATURES Current slack timing/rate features and robust replacements.
%
%   F = robustSlackFeatures(T, Y, SL, VT, USESMOOTHING)
%   F = robustSlackFeatures(T, Y, SL, VT, USESMOOTHING, TL, LSE)   % + model t0_crossing
%
% Same windows, anchors and smoothers as extractSlackAttributes (reproduced
% here so both versions read the identical trace), plus replacements for the
% features that read a single sample. The rule behind every replacement:
%
%   * a LEVEL CROSSING on a steep flank is well conditioned: its error is
%     noise / slope, and linear interpolation makes it continuous in the trace;
%   * an extremum VALUE is well conditioned (first-order blind to where the
%     extremum is);
%   * an extremum LOCATION is not: on a flat top/bottom its error is
%     sqrt(noise / curvature), and read off samples it is quantized to one
%     sample -- 1 ms on the model's 1 kHz read, 0.25 ms on the data.
%
% So each timing that was anchored on an argmin/argmax sample is re-anchored
% on a protocol time or on a flank crossing, and each peak location is read
% by a local quadratic fit instead of by the sample index. Every replacement
% measures from the PROTOCOL time (velocity-table row), never from the first
% sample of its window: that sample is itself quantized to the read grid.
%
% Fields (one column per slack cycle):
%   current (as scored now)      robust replacement
%   rsK80   2-crossing rate       rsKfit    band regression (already in the code)
%   rsT20   t20 - t(argmin)       rsT20e    t20 - ramp end   (protocol anchor)
%                                 rsTfitE   band-fit t20 - ramp end
%   peak1_t   argmax sample       peak1_tq  quadratic vertex (+-1 ms)
%   peak1_dSL SL at argmax sample peak1_dSLq SL interpolated at the vertex
%   vall_t    argmin sample       vall_tq   quadratic vertex (+-1.5 ms)
%                                 tFall50   mid-level crossing on the falling flank
%   ktr80   2-crossing hold rate  ktrfit    band regression (already in the code)
%   t0_crossing (model: first     tOn20     hold force crosses 20 % of its rise,
%     output sample with LSE>=0)            from the release start (model & data)
%                                 t0x       model only: LSE zero crossing, interpolated
%   restretchSlopeStart (>=3      rsSlope1ms fixed 1 ms window on a 0.25 ms grid
%     samples in first 1 ms)
%   rsK80                         rsKls/rsT0ls fixed-amplitude least squares, raw samples
%   peak1_y / vall_y (max/min     peak1_yq / vall_yq  the same parabola's vertex value
%     of raw samples)

if nargin < 5 || isempty(useSmoothing); useSmoothing = false; end
if nargin < 7; tL = []; LSE = []; end
t = t(:); y = y(:); SL = SL(:);
[t, iu] = unique(t); y = y(iu); SL = SL(iu);
tL = tL(:); LSE = LSE(:);

RESTRETCH_VALL_WIN_S = 0.080; FORCE_THRESHOLD_KPA = 10;
PEAK_PROM_K_SIGMA = 3; PEAK_PROM_MIN_FRAC = 0.05; PEAK_MIN_DISTANCE_S = 0.001;
hp = y - movmedian(y, 9); noise_sigma = 1.4826 * median(abs(hp - median(hp)));

segs = find(vt(:, 2) < -1)';
nS = numel(segs);
names = {'rsK80','rsKfit','rsT20','rsT20e','rsTfit','rsTfitE','peak1_t','peak1_tq', ...
         'peak1_dSL','peak1_dSLq','vall_t','vall_tq','tFall50','ktr80','ktrfit', ...
         'tOn20','t0x','restretchSlopeStart','rsSlope1ms','steady','peak1_y','vall_y', ...
         'peak1_yq','vall_yq','rsKls','rsT0ls'};
for k = 1:numel(names); F.(names{k}) = nan(1, nS); end

for i = 1:nS
    if segs(i) + 4 > size(vt, 1); break; end
    vs = vt(segs(i):segs(i)+4, 1);
    t_seg = vs(1);

    % ---- slack hold (vs2..vs3): SLslack, hold redevelopment --------------
    win = t > vs(2) & t < vs(3) & y > FORCE_THRESHOLD_KPA;
    SLslack = NaN; if any(win); s_ = SL(win); SLslack = s_(end); end
    winH = t > vs(2) & t < vs(3);
    [F.ktr80(i), F.ktrfit(i), t20h] = holdRate(t(winH), y(winH));
    F.tOn20(i) = t20h - t_seg;
    if ~isempty(LSE)
        m = tL >= vs(2) & tL <= vs(3);
        tl = tL(m); ls = LSE(m);
        j = find(ls >= 0, 1);
        if ~isempty(j)
            if j > 1 && ls(j) ~= ls(j-1)
                F.t0x(i) = tl(j-1) + (0 - ls(j-1)) / (ls(j) - ls(j-1)) * (tl(j) - tl(j-1)) - t_seg;
            else
                F.t0x(i) = tl(j) - t_seg;
            end
        end
    end

    % ---- re-stretch ramp (vs3..vs4): slope, peak, valley ------------------
    w = t >= vs(3) & t <= vs(4);
    tw = t(w); yw = y(w); sw = SL(w);
    if numel(tw) >= 4
        r0 = tw > tw(1) & tw < tw(1) + 1e-3; r0(1:3) = true;
        p = polyfit(sw(r0), yw(r0), 1); F.restretchSlopeStart(i) = p(1);
        tg = (vs(3):2.5e-4:vs(3) + 1e-3)';          % protocol-anchored, fixed 1 ms
        p = polyfit(interp1(t, SL, tg), interp1(t, y, tg), 1); F.rsSlope1ms(i) = p(1);

        prom = max(PEAK_PROM_K_SIGMA * noise_sigma, PEAK_PROM_MIN_FRAC * (max(yw) - min(yw)));
        [pk, pt] = findpeaks(yw, tw, 'MinPeakDistance', PEAK_MIN_DISTANCE_S, 'MinPeakProminence', prom);
        if ~isempty(pk)
            F.peak1_y(i) = pk(1);
            F.peak1_t(i) = pt(1) - tw(1);
            F.peak1_dSL(i) = sw(find(tw >= pt(1), 1)) - SLslack;
            [tq, F.peak1_yq(i)] = vertexTime(tw, yw, pt(1), 1e-3);
            F.peak1_tq(i) = tq - vs(3);                 % protocol anchor, not the first sample
            F.peak1_dSLq(i) = interp1(t, SL, tq) - SLslack;
            after = tw > pt(1);
            if any(after)
                ta = tw(after); ya = yw(after);
                [vv, iv] = min(ya);
                F.vall_y(i) = vv;
                F.vall_t(i) = ta(iv) - tw(1);
                [tv, F.vall_yq(i)] = vertexTime(tw, yw, ta(iv), 1.5e-3);
                F.vall_tq(i) = tv - vs(3);
                % falling flank between the peak and the valley, crossed at mid level
                nm = max(1, round(1e-3 / median(diff(tw))));
                yf = movmedian(yw, nm);
                lev = 0.5 * (pk(1) + vv);
                jf = find(tw > pt(1) & tw <= ta(iv));
                F.tFall50(i) = crossDown(tw(jf), yf(jf), lev) - vs(3);
            end
        end
    end

    % ---- post-restretch hold (vs4..vs5): steady, recovery rate/delay ------
    ws = t >= vs(5) - 0.02 & t <= vs(5);
    steady = median(y(ws)); F.steady(i) = steady;
    w = t >= vs(4) & t <= vs(5);
    tw = t(w) - t_seg; yw = y(w);
    if useSmoothing
        ya = smoothdata(yw, 'movmedian', max(3, round(0.005 / median(diff(tw)))));
    else
        ya = yw;
    end
    [F.rsK80(i), F.rsKfit(i), F.rsT20(i), F.rsT20e(i), F.rsTfit(i), F.rsTfitE(i), F.rsKls(i), F.rsT0ls(i)] = ...
        recovery(tw, yw, ya, steady, RESTRETCH_VALL_WIN_S, vs(4) - t_seg);
end
end

% =========================================================================
function [k80, kfit, t20, t20e, tfit, tfitE, kls, t0ls] = recovery(t, y, ya, steady, vallWin, tEnd)
% fitRestretchRecovery's rsK80/rsT20/rsKfit/rsTfit, plus the ramp-end anchored delays
    k80 = NaN; kfit = NaN; t20 = NaN; t20e = NaN; tfit = NaN; tfitE = NaN; kls = NaN; t0ls = NaN;
    ok = isfinite(t) & isfinite(y) & isfinite(ya); t = t(ok); y = y(ok); ya = ya(ok);
    if numel(t) < 20; return; end
    nWin = nnz(t <= t(1) + vallWin); if nWin < 3; nWin = numel(t); end
    [B, iB] = min(ya(1:nWin));
    tt = t(iB:end) - t(iB); yy = y(iB:end) - B;
    Aext = steady - B;
    if ~(isfinite(Aext) && Aext > 0) || numel(tt) < 20; return; end
    nx = max(3, round(0.005 / median(diff(tt))));
    zx = movmedian(yy, nx) / Aext;
    a = firstCrossing(tt, zx, 0.20); b = firstCrossing(tt, zx, 0.80);
    if isfinite(a) && isfinite(b) && b > a; k80 = log(4) / (b - a); end
    t20 = a;
    t20e = a + t(iB) - tEnd;                     % anchored on the ramp END (protocol time), not the argmin
    [kfit, tfit] = bandRateFit(tt, zx, 0.20, 0.80);
    tfitE = tfit + t(iB) - tEnd;
    % Least squares on the RAW samples, amplitude fixed at Aext (as for rsK80),
    % rate and onset free, over the first passage up to z = 0.9 of the smoothed
    % shape -- every sample of the rise, none of the later overshoot. For white
    % noise this is the efficient estimator; the two crossings use two points.
    i9 = find(zx >= 0.9, 1); if isempty(i9), i9 = numel(tt); end
    if i9 >= 10 && isfinite(a)
        xx = tt(1:i9); zv = yy(1:i9) / Aext;
        k0 = kfit; if ~isfinite(k0), k0 = 45; end
        obj = @(q) sum((zv - (1 - exp(-abs(q(1)) .* max(xx - q(2), 0)))).^2);
        q = fminsearch(obj, [k0, a - log(1/0.8)/k0], ...
                       optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-12, 'MaxFunEvals', 4e3, 'MaxIter', 4e3));
        kls = abs(q(1)); t0ls = q(2) + t(iB) - tEnd;
    end
end

function [k80, kfit, t20] = holdRate(t, y)
% extractSlackAttributes' holdRate80, also returning the absolute 20 % crossing
    k80 = NaN; kfit = NaN; t20 = NaN;
    ok = isfinite(t) & isfinite(y); t = t(ok); y = y(ok);
    if numel(t) < 40; return; end
    nx = max(3, round(0.005 / median(diff(t))));
    ys = movmedian(y, nx);
    half = find(t <= t(1) + 0.5 * (t(end) - t(1)));
    [B, iB] = min(ys(half));
    P = median(ys(t >= t(end) - 0.040));
    if ~isfinite(B) || ~isfinite(P) || P - B <= 0; return; end
    z = (ys(iB:end) - B) / (P - B); tt = t(iB:end);
    t20 = firstCrossing(tt, z, 0.20); t80 = firstCrossing(tt, z, 0.80);
    if isfinite(t20) && isfinite(t80) && t80 > t20; k80 = log(4) / (t80 - t20); end
    kfit = bandRateFit(tt, z, 0.20, 0.80);
end

function [tq, yq] = vertexTime(t, y, t0, hw)
% vertex of a least-squares parabola through the samples within +-hw of t0;
% falls back to t0 when the fit is not a clean extremum inside the window
    m = abs(t - t0) <= hw;
    if nnz(m) < 3
        [~, o] = sort(abs(t - t0)); m = false(size(t)); m(o(1:min(3, numel(t)))) = true;
    end
    x = t(m) - t0; p = polyfit(x, y(m), 2);
    tq = t0; yq = polyval(p, 0);
    if p(1) ~= 0
        xv = -p(2) / (2 * p(1));
        if abs(xv) <= hw; tq = t0 + xv; yq = polyval(p, xv); end
    end
end

function tx = crossDown(tt, z, lev)
    tx = NaN; i = find(z <= lev, 1);
    if isempty(i); return; end
    if i == 1; tx = tt(1); return; end
    z1 = z(i-1); z2 = z(i);
    if z2 == z1; tx = tt(i); else; tx = tt(i-1) + (lev - z1) / (z2 - z1) * (tt(i) - tt(i-1)); end
end

function tx = firstCrossing(tt, z, lev)
    tx = NaN; i = find(z >= lev, 1);
    if isempty(i); return; end
    if i == 1; tx = tt(1); return; end
    z1 = z(i-1); z2 = z(i);
    if z2 == z1; tx = tt(i); else; tx = tt(i-1) + (lev - z1) / (z2 - z1) * (tt(i) - tt(i-1)); end
end

function [k, tlo] = bandRateFit(tt, z, lo, hi)
% verbatim from extractSlackAttributes
    k = NaN; tlo = NaN;
    i1 = find(z >= lo, 1); i2 = find(z >= hi, 1);
    if isempty(i1) || isempty(i2) || i2 - i1 < 3; return; end
    idx = i1:i2; zz = min(z(idx), 1 - 1e-6);
    x = tt(idx); yv = log(1 - zz); w = (1 - zz).^2;
    W = sum(w); xm = sum(w.*x)/W; ym = sum(w.*yv)/W;
    sxx = sum(w.*(x - xm).^2); if ~(sxx > 0); return; end
    slope = sum(w.*(x - xm).*(yv - ym)) / sxx; if ~(slope < 0); return; end
    k = -slope; a = ym - slope*xm; tlo = (log(1 - lo) - a) / slope;
end
