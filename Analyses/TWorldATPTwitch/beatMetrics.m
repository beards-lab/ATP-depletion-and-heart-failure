function m = beatMetrics(t, y)
%BEATMETRICS  Transient metrics of one beat. T in s (from the stimulus), Y the trace.
%   dia = value at t(1); peak, amp = peak - dia; ttp (ms); d50/d90 = time from the peak until
%   the trace has fallen 50 % / 90 % of the way back to dia (ms); tau = exponential decay time
%   constant fitted between 80 % and 20 % of amp on the way down (ms).
t = t(:);  y = y(:);
dia = y(1);  [pk, ip] = max(y);  amp = pk - dia;
td = t(ip:end);  yd = y(ip:end);
dn = @(f) 1e3*(td(find(yd <= dia + (1 - f)*amp, 1)) - t(ip));
m.dia = dia;  m.peak = pk;  m.amp = amp;  m.ttp = 1e3*t(ip);
m.d50 = val(dn(0.5));  m.d90 = val(dn(0.9));
w = yd <= dia + 0.8*amp & yd >= dia + 0.2*amp & td <= td(find(yd <= dia + 0.2*amp, 1) + 0);
if nnz(w) > 3
    p = polyfit(td(w), log(yd(w) - dia), 1);  m.tau = -1e3/p(1);
else
    m.tau = NaN;
end
end
function v = val(x), if isempty(x), v = NaN; else, v = x; end, end
