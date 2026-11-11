function R = readIsoV(t, L, F, opts)
% readIsoV  Read the isovelocity points from ONE continuous Beard2022 record.
%
%   R = readIsoV(t, L, F)          t (ms), L (Lo), F (kPa) column vectors
%   R = readIsoV(t, L, F, opts)    opts.halfWin (ms, default 0.5), opts.Lx (1.00)
%
% The same reader is applied to data and to model output, so any reader bias is
% common to both. The protocol per cycle: isometric hold at ~1.10 Lo -> fast
% step -> constant-velocity shortening through 1.00 Lo -> slack -> re-stretch to
% 1.10. The last cycle is the v = 0 point: a step to 1.00 Lo and an isometric
% hold to the end of the record.
%
% R is a table with one row per cycle, in recorded order:
%   tStart  (ms)   first sample with L below the 1.10 hold (cycle onset)
%   tCross  (ms)   first sample with L <= Lx after onset (NaN for v = 0)
%   vMeas   (ML/s) shortening speed fitted over the ramp around tCross
%   Fx      (kPa)  mean force over tCross +- halfWin  (v = 0: mean over the
%                  last 200 ms of the record)
%   Lx      (Lo)   mean length over the same window
%   tSlack  (ms)   onset of the post-ramp slack (end of the ramp)
%   F0      (kPa)  force baseline: mean force while the preparation hangs slack
%                  after the ramp (L within 0.01 Lo of the cycle's shortest
%                  length, from 3 ms after the release). The workbook's own
%                  `zero` column is the same quantity. The v = 0 row takes
%                  the baseline of the last cycle.

if nargin < 4, opts = struct(); end
if ~isfield(opts, 'halfWin'), opts.halfWin = 0.5; end
if ~isfield(opts, 'Lx'), opts.Lx = 1.00; end
t = t(:); L = L(:); F = F(:);

Lhold = median(L(t < 400));                        % the 1.10 Lo hold
below = L < Lhold - 0.01;
on = find(below & ~[false; below(1:end-1)]);       % rising edges of "below hold"
% keep onsets that follow at least 50 ms at the hold (rejects ringing)
keep = true(size(on));
for i = 1:numel(on)
    prevHold = find(~below(1:on(i)-1), 1, 'last');
    if isempty(prevHold), continue; end
    j0 = find(t >= t(on(i)) - 50, 1);
    keep(i) = all(~below(j0:on(i)-1));
end
on = on(keep);

n = numel(on);
R = table(nan(n,1), nan(n,1), nan(n,1), nan(n,1), nan(n,1), nan(n,1), nan(n,1), ...
    'VariableNames', {'tStart', 'tCross', 'vMeas', 'Fx', 'Lx', 'tSlack', 'F0'});
% smoothed speed (1 ms moving mean) for the slack detector, ML/s
dt = median(diff(t));
vS = -movmean([0; diff(L)] / dt * 1000, max(1, round(1/dt)));
for i = 1:n
    i0 = on(i);
    R.tStart(i) = t(i0);
    iEnd = numel(t); if i < n, iEnd = on(i+1) - 1; end
    seg = i0:iEnd;
    ic = seg(find(L(seg) <= opts.Lx, 1));
    if isempty(ic), continue; end
    R.tCross(i) = t(ic);
    w = abs(t - t(ic)) <= opts.halfWin;
    R.Fx(i) = mean(F(w)); R.Lx(i) = mean(L(w));
    % slack onset: first sample after the crossing shortening faster than 30 ML/s
    % (the fastest ramp is 6 ML/s, the slack release is ~400 ML/s)
    is = find(t(seg) > t(ic) & vS(seg) > 30, 1);
    if ~isempty(is), R.tSlack(i) = t(seg(is)); else, R.tSlack(i) = t(iEnd); end
    % ramp speed: straight-line fit to L over the ramp before the slack, from
    % 2 ms after the onset step to 0.5 ms before the slack
    wr = seg(t(seg) >= t(i0) + 2 & t(seg) <= R.tSlack(i) - 1 & L(seg) > opts.Lx - 0.03);
    wr = wr(t(wr) >= t(ic) - max(10, 0.5 * (t(ic) - t(i0))));   % ramp part only
    p = polyfit(t(wr)/1000, L(wr), 1); R.vMeas(i) = -p(1);
    % force baseline while slack (skip 3 ms of ringing after the release)
    sl = seg(t(seg) > R.tSlack(i) + 3 & L(seg) <= min(L(seg)) + 0.01);
    if ~isempty(sl), R.F0(i) = mean(F(sl)); end
end
% v = 0 point: the last cycle re-stretches to ~1.00 Lo and holds there to the end
w = t > t(end) - 200;
if abs(median(L(w)) - opts.Lx) < 0.01
    R(end+1, :) = {NaN, NaN, 0, mean(F(w)), mean(L(w)), NaN, R.F0(n)};
end
end
