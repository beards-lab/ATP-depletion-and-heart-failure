function R = srxProbe(params0, fcn)
% SRXPROBE  Characterize SRX (super-relaxed) dynamics for a given paramset.
%   R = srxProbe(params0[, fcn]) runs activated isometric simulations and
%   returns metrics describing the DRX<->SRX behaviour: steady-state pool
%   size at the active working point, the force-sensitivity swing, and the
%   intrinsic equilibration time.
%
%   Fields of R:
%     F_active     active isometric steady-state force (the F_SR driver)
%     SR_active / SRD_active / SRX_active   SRX pool at active force
%     SRX_rest / F_rest   SRX pool at a short-SL (low passive force) probe
%     swing        SRX_rest / SRX_active  (force-sensitivity dynamic range)
%     tau_eig_active_ms / tau_eig_rest_ms   slow-mode relaxation time (ms)
%
%   The relaxation time is the slowest eigenvalue of the [P_SR,P_SRD]
%   exchange matrix (PT,PD buffered). For the symmetric design the slow
%   mode equals the mobilization rate, so tau_rest ~ 1000/kmsr [ms]; this
%   was confirmed against direct simulation (kmsr=60->16.7ms, 100->10.0ms).
%
%   See also: getParams, evaluateModel, dPUdT_CombinedTransitions
if nargin < 2, fcn = @dPUdT_CombinedTransitions; end

R.F_active = runIso(params0, fcn, params0.SL0, @(o) o.FXB(end)+o.FXBPassive(end));
[R.SRX_active, R.SR_active, R.SRD_active] = isoSRX(params0, fcn, params0.SL0);

% Resting probe: short SL -> low passive force -> upper end of the SRX swing
R.F_rest = runIso(params0, fcn, 1.85, @(o) o.FXB(end)+o.FXBPassive(end));
R.SRX_rest = isoSRX(params0, fcn, 1.85);
R.swing = R.SRX_rest / max(R.SRX_active, eps);

R.tau_eig_active_ms = srxEigTau(params0, R.F_active);
R.tau_eig_rest_ms   = srxEigTau(params0, 0);
end

function v = runIso(p0, fcn, SL, pick)
o = isoRun(p0, fcn, SL); v = pick(o);
end

function [SRX, SR, SRD] = isoSRX(p0, fcn, SL)
o = isoRun(p0, fcn, SL);
SR = o.SR(end); SRD = o.SRD(end); SRX = SR + SRD;
end

function o = isoRun(p0, fcn, SL)
p = p0; p.Velocity = 0; p.MgATP = 8; p.Ca = 1000; p.ValuesInTime = 1;
p.UseSlack = 0; p.UseKtrProtocol = 0;
p.RunSlack = 0; p.RunForceVelocity = 0; p.RunKtr = 0; p.RunStairs = 0;
p.SL0 = SL; p = getParams(p);
[~, o] = evaluateModel(fcn, [0 0.6], p);
end

function tau_ms = srxEigTau(p, F)
% Slow-mode relaxation time of the [P_SR,P_SRD] subsystem with PT,PD buffered.
a_mob = p.kmsr  * exp(F/p.sigma1);       % SR  -> PT  (mobilize)
d_mob = p.kmsrd * exp(F/p.sigma_srd1);   % SRD -> PD  (mobilize)
c1 = p.ksr2srd; c2 = p.ksrd2sr;          % SR <-> SRD coupling
M = [ -(a_mob+c1), c2 ; c1, -(d_mob+c2) ];
tau_ms = 1000 / abs(max(real(eig(M))));
end
