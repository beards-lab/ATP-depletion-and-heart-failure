function m = caTransientMetrics(kind, opt, T)
%CATRANSIENTMETRICS  Shape metrics of a periodic caTransient KIND with options OPT.
%   M = CATRANSIENTMETRICS(KIND, OPT, T) returns, for one beat of period T (s):
%     dia, peak (uM), amp; ttp (s from beat onset, i.e. from the stimulus), ttp_onset
%     (from the transient's own onset t0 if present), t10_90 rise, dec50, dec90 (peak ->
%     50 %/90 % decay, s), ctd50, ctd90 (durations at 50 %/10 % of amplitude), tau_d_fit
%     (mono-exponential tau fitted on the decay from 80 % to 20 %), and the fraction of
%     the amplitude still present at the next stimulus (carry).
if nargin < 3, T = 1; end
opt.T = T;
t  = linspace(0, T, 20001)';  t = t(1:end-1);
Ca = caTransient(kind, t, opt);
[pk, ip] = max(Ca);  dia = min(Ca);  amp = pk - dia;
up = Ca(1:ip);  tu = t(1:ip);  dn = Ca(ip:end);  td = t(ip:end);
lvl = @(f) dia + f*amp;
cu = @(f) tu(find(up >= lvl(f), 1));
cd_ = @(f) td(find(dn <= lvl(f), 1));
t0 = 0;  if isfield(opt, 't0'), t0 = opt.t0; end
m.dia = dia;  m.peak = pk;  m.amp = amp;  m.ttp = t(ip);  m.ttp_onset = t(ip) - t0;
m.t10_90 = cu(0.9) - cu(0.1);
m.dec50 = cd_(0.5) - t(ip);  m.dec90 = cd_(0.1) - t(ip);  m.dec80 = cd_(0.2) - t(ip);
m.t80rise = cu(0.8) - t0;  m.dt90_10 = cd_(0.1) - cd_(0.9);
m.ctd50 = cd_(0.5) - cu(0.5);  m.ctd90 = cd_(0.1) - cu(0.1);
w = td >= cd_(0.8) & td <= cd_(0.2);
pfit = polyfit(td(w), log(max(dn(w) - dia, 1e-12)), 1);  m.tau_d_fit = -1/pfit(1);
m.carry = (Ca(end) - dia)/amp;
end
