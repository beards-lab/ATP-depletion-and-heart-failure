function P = fv2Read(t, F, vt)
%FV2READ P/P0 of the FV2 ramp: mean force over [t1+5 ms, t2] / isometric
%   (median over 2-40 ms). t1/t2 = start/end of the -2 ML/s ramp, taken from
%   the protocol's velocity table (the row whose speed is between -5 and -0.5 ML/s).
i  = find(vt(:,2) < -0.5 & vt(:,2) > -5, 1);
t1 = vt(i,1); t2 = vt(i+1,1);
Fss = median(F(t >= 0.002 & t < 0.040));
w   = t >= t1 + 0.005 & t <= t2;
P   = mean(F(w)) / Fss;
end
