function R = mechProbe(params0, fcn)
% MECHPROBE  Measure the attached pool / ATPase / ktr / Vmax tradeoff.
%   R = mechProbe(params0[, fcn]) runs activated isometric contraction and a
%   small force-velocity sweep, and returns the four coupled quantities the
%   modeller must reconcile, using the model's own feature conventions
%   (extractSlackAttributes): attached_ss = p1_0+p2_0, XTOR = RTD (hydrolysis
%   flux = ATP turnover per head /s).
%
%   Fields:
%     attached      p1_0 + p2_0 at isometric steady state  (target > 0.30)
%     p1 / p2       attached sub-pools
%     XTOR          ATP turnover (RTD, gross hydrolysis flux) [1/s]  (target < 5)
%     XTOR_net      net hydrolysis (RTD - RDT) [1/s]
%     PT / PD       free pools (PuATP pre-hydrolysis / PuR post-hydrolysis)
%     SRX / NP      parked / non-permissive
%     F             isometric force
%     k2eff         effective isometric detachment rate from p2 [1/s]
%     dwell_ms      mean attached lifetime = attached/XTOR [ms]
%     ktr           force-redevelopment rate from activation rise [1/s] (~50)
%     Vmax          force-velocity zero-force intercept [ML/s]
if nargin < 2, fcn = @dPUdT_CombinedTransitions; end

% ---- isometric activation to steady state ----
p = params0; p.Velocity = 0; p.MgATP = 8; p.Ca = 1000; p.ValuesInTime = 1;
p.UseSlack = 0; p.UseKtrProtocol = 0;
p.RunSlack = 0; p.RunForceVelocity = 0; p.RunKtr = 0; p.RunStairs = 0;
p = getParams(p);
[~, o] = evaluateModel(fcn, [0 0.6], p);
t = o.t(:); F = o.Force(:);

R.attached = o.p1_0(end) + o.p2_0(end);
R.p1 = o.p1_0(end); R.p2 = o.p2_0(end);
R.XTOR = o.RTD(end);
R.XTOR_net = o.RTD(end) - o.RDT(end);
R.PT = o.PuATP(end); R.PD = o.PuR(end);
R.SRX = o.SR(end) + o.SRD(end); R.NP = o.NP(end);
R.F = F(end);
% effective isometric detachment from p2 (R2T+R2D are strain-integrated fluxes)
detflux = o.R2T(end) + o.R2D(end);
R.k2eff = detflux / max(o.p2_0(end), eps);
R.dwell_ms = 1000 * R.attached / max(R.XTOR, eps);

% ---- ktr from the force-redevelopment (activation) rise ----
Fss = F(end);
win = t > 0.002 & F < 0.92*Fss & F > 0.08*Fss;
if nnz(win) >= 3
    pf = polyfit(t(win), log(1 - F(win)/Fss), 1);
    R.ktr = -pf(1);
else
    R.ktr = NaN;
end

% ---- Vmax from a constant-velocity force-velocity sweep ----
vels = [0 -0.5 -1 -1.5 -2 -3];
Fv = nan(size(vels));
for i = 1:numel(vels)
    pv = params0; pv.Velocity = vels(i); pv.MgATP = 8; pv.Ca = 1000; pv.ValuesInTime = 1;
    pv.UseSlack = 0; pv.UseKtrProtocol = 0;
    pv.RunSlack = 0; pv.RunForceVelocity = 0; pv.RunKtr = 0; pv.RunStairs = 0;
    pv = getParams(pv);
    try
        [~, ov] = evaluateModel(fcn, [0 0.25], pv);
        Fv(i) = ov.FXB(end) + ov.FXBPassive(end);
    catch
        Fv(i) = NaN;
    end
end
R.FV_v = vels; R.FV_f = Fv;
% zero-force intercept of FXB (active) -> Vmax; interpolate on active force
Fa = Fv; % includes small passive; subtract passive proxy at v=0? use active only:
good = ~isnan(Fa);
vg = vels(good); fg = Fa(good);
% find first sign change to F<=0 (active force ~ 0)
R.Vmax = NaN;
for i = 1:numel(vg)-1
    if fg(i) > 0 && fg(i+1) <= 0 || (fg(i) > 0.5 && fg(i+1) < 0.5)
        % linear interp to F=0 (use 0.5 as near-zero floor due to passive)
        R.Vmax = interp1([fg(i) fg(i+1)], [vg(i) vg(i+1)], 0, 'linear', 'extrap');
        break;
    end
end
end
