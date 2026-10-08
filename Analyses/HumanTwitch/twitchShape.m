function m = twitchShape(tb, F)
%TWITCHSHAPE  TTP (from stimulus), 10 %->peak, RT50, RT90 (from peak), width at 50 %; NaN if undefined.
[Fp, ip] = max(F);  up = F(1:ip);  dn = F(ip:end);  tu = tb(1:ip);  td = tb(ip:end);
m.Fpk = Fp;  m.TTP = tb(ip);
m.t10 = crossUp(tu, up, 0.1*Fp);  m.t50 = crossUp(tu, up, 0.5*Fp);
m.C10toPk = tb(ip) - m.t10;
m.RT50 = crossDn(td, dn, 0.5*Fp) - tb(ip);  m.RT90 = crossDn(td, dn, 0.1*Fp) - tb(ip);
m.W50 = crossDn(td, dn, 0.5*Fp) - m.t50;
end
function t = crossUp(tt, ff, lvl)
k = find(ff >= lvl, 1);  if isempty(k) || k == 1, t = NaN; return; end
t = interp1(ff(k-1:k), tt(k-1:k), lvl);
end
function t = crossDn(tt, ff, lvl)
k = find(ff <= lvl, 1);  if isempty(k) || k == 1, t = NaN; return; end
t = interp1(ff(k-1:k), tt(k-1:k), lvl);
end
