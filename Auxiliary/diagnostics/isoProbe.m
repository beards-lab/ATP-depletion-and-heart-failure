function R = isoProbe(p0)
% ISOPROBE  Fast isometric-only measurement of the attached/ATPase/ktr triad.
%   Used for parameter sweeps (no force-velocity runs). See mechProbe for Vmax.
p = p0; p.Velocity = 0; p.MgATP = 8; p.Ca = 1000; p.ValuesInTime = 1;
p.UseSlack = 0; p.UseKtrProtocol = 0;
p.RunSlack = 0; p.RunForceVelocity = 0; p.RunKtr = 0; p.RunStairs = 0;
p = getParams(p);
[~, o] = evaluateModel(@dPUdT_CombinedTransitions, [0 0.6], p);
t = o.t(:); F = o.Force(:); Fss = F(end);
R.att   = o.p1_0(end) + o.p2_0(end);
R.XTOR  = o.RTD(end);
R.PT    = o.PuATP(end);
R.PD    = o.PuR(end);
R.SRX   = o.SR(end) + o.SRD(end);
R.k2eff = (o.R2T(end) + o.R2D(end)) / max(o.p2_0(end), eps);
R.F     = Fss;
w = t > 0.002 & F < 0.92*Fss & F > 0.08*Fss;
if nnz(w) >= 3
    pf = polyfit(t(w), log(1 - F(w)/Fss), 1); R.ktr = -pf(1);
else
    R.ktr = NaN;
end
end
