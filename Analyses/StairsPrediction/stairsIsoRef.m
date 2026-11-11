function ref = stairsIsoRef(base, vt)
%STAIRSISOREF  Isometric operating point of the staircase run, for compensation.
%
%   Returns the frozen isometric p2 distribution, the strain axis the ODE
%   actually uses at that instant (which needs the DYNAMIC out.LXBPivot, not
%   params.LXBpivot), and a handle pw(e) = sum sign(x)|x|^e p2 dS, so a caller
%   can rescale kstiff2 to hold the isometric active force fixed while changing
%   the exponent.

p = base; p.RunStairs = 1; p.PlotEachSeparately = false;
p.ghostSave = ''; p.ghostLoad = '';
[~, o] = runStairsExperiment(p);

ss   = numel(o.params.s);
tMov = vt(find(vt(:,2)~=0,1),1);
[~,i0] = min(abs(o.t - (tMov - 0.005)));

p2   = o.PU(i0, ss+1:2*ss)';
cSh  = (-(o.SL(i0)-o.LSE(i0)) + o.LXBPivot(i0))/2;
x    = o.params.s(:) - cSh + base.dr;
dS   = o.params.s(2) - o.params.s(1);

ref.p2   = p2;
ref.x    = x;
ref.dS   = dS;
ref.pw   = @(e) dS*sum(sign(x).*abs(x).^e.*p2);
ref.Fiso = o.Force(i0);
ref.Fp1  = o.Force(i0) - base.kstiff2*ref.pw(base.estiff);   % everything not p2
ref.out  = o;
end
