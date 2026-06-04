function ktr = ktrProbe(p0)
% KTRPROBE  Brenner-style force-redevelopment rate.
%   Activates to isometric steady state, then instantaneously detaches the
%   bound bridges (p1=p2=0) while holding the upstream pools (PT,PD,SRX) and
%   the length at their steady-state values, and fits the force redevelopment
%   F(t)=Fss*(1-exp(-ktr*t)). This isolates attachment/power-stroke kinetics
%   from the slow hydrolysis/SRX filling, matching the experimental ktr.
p = p0; p.Velocity = 0; p.MgATP = 8; p.Ca = 1000; p.ValuesInTime = 1;
p.UseSlack = 0; p.UseKtrProtocol = 0;
p.RunSlack = 0; p.RunForceVelocity = 0; p.RunKtr = 0; p.RunStairs = 0;
p = getParams(p);
[~, o] = evaluateModel(@dPUdT_CombinedTransitions, [0 0.6], p);

ss = p.ss;
% steady-state pools
SR = o.SR(end); SRD = o.SRD(end); PD = o.PuR(end);
att = o.p1_0(end) + o.p2_0(end);
SLss = o.SL(end); LSEss = o.LSE(end);

% Build a fresh detached IC (correct length, no restart/regrid issues).
% Layout (2-state): [p1(1:ss) p2(1:ss) P_SR NP SL LSE PD P_SRD x_dash]
PU0 = zeros(1, 2*ss + 7);
PU0(2*ss+1) = SR;            % P_SR at SS
PU0(2*ss+2) = 0;            % NP
PU0(2*ss+3) = SLss;        % SL held
PU0(2*ss+4) = LSEss;       % LSE at SS
PU0(2*ss+5) = PD + att;    % freed bridges return to PD (post-hydrolysis, ready)
PU0(2*ss+6) = SRD;          % P_SRD at SS
PU0(2*ss+7) = 0;            % x_dash

p.PU0 = PU0; p.UseSpaceExtension = 0;
[~, o2] = evaluateModel(@dPUdT_CombinedTransitions, [0 0.15], p);
t = o2.t(:); F = o2.Force(:); Fss = F(end);
w = t > 0.001 & F < 0.92*Fss & F > 0.08*Fss;
if nnz(w) >= 3
    pf = polyfit(t(w), log(1 - F(w)/Fss), 1); ktr = -pf(1);
else
    ktr = NaN;
end
end
