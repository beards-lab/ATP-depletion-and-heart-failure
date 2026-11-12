function [f, outputs] = dPUdT_HostWall(t, PU, params)
%DPUDT_HOSTWALL  Probe wrapper: one TriSeg wall's sarcomere, driven open-loop by the host.
%
%   Wraps the UNCHANGED dPUdT_CombinedTransitions. Two things are added, both outside
%   the model file:
%   1. Ca activation: one appended state N (non-permissive thin-filament fraction),
%      exactly the factorized gate of Analyses/CalciumTwitch/dPUdT_CaGate.m
%      (gate scales params.ka; gate = 1 reproduces the calibrated model).
%   2. Host kinematics: the sarcomere follows the host's Ls(t), dSL/dt = v_host(t).
%      The series element and the contractile element are the model's own
%      (dLSE/dt = v_host - velHS). The ONLY thing that changes between modes is the
%      velocity the velocity-driven mechanisms see (A_reg drive, A2Reattaching gate,
%      mu2*vel, dynamic passive, catch bond, stretch activation):
%        params.host.velMode = 'imposed'  -> v_host(t)      (what evaluateModel does)
%                              'velHS'    -> contractile-element velocity (two passes)
%                              'zero'     -> 0 (velocity-driven mechanisms switched off)
%      SL and LSE are integrated with v_host in every mode, so the kinematics are
%      identical across modes.
%
%   params.host: Ls(t)/v(t) handles of one wall (hostSignals), velMode.
%   params.tf:   k_on, k_off, K_coop, CaFcn(t) (uM).
%   OUTPUTS = [core outputs, gate, N, Ca, v_host, v_mech, velHS].

n    = numel(PU);
N    = min(1, max(0, PU(n)));
core = PU(1:n-1);
ss   = params.ss;  Ns = params.NumberOfStates;
iSL  = Ns*ss + 3;  iLSE = Ns*ss + 4;

% --- Ca gate (as dPUdT_CaGate) ---
tf    = params.tf;
bound = params.dS*sum(max(0, core(1:Ns*ss)));
Ca    = tf.CaFcn(t);
Pu    = max(0, 1 - N - bound);
Jon   = tf.k_on*Ca*N*(1 + tf.K_coop*(1 - N));
Joff  = tf.k_off*Pu*(1 + tf.K_coop*N);
dN    = Joff - Jon;
gate  = Pu/max(1 - bound, eps);

p = params;
p.ka = params.ka*gate;
vh = params.host.v(t);

switch params.host.velMode
    case 'imposed'
        vm = vh;
    case 'zero'
        vm = 0;
    case 'velHS'
        p.Vums = 0;
        f0 = dPUdT_CombinedTransitions(t, core, p);
        vm = 0 - f0(iLSE);             % dLSE = vel - velHS  ->  velHS = -dLSE at vel = 0
    otherwise
        error('velMode %s', params.host.velMode);
end
p.Vums = vm;
[fc, outputs] = dPUdT_CombinedTransitions(t, core, p);
velHS = vm - fc(iLSE);
% kinematics: the sarcomere follows the host whatever the mechanisms saw
fc(iSL)  = vh;
fc(iLSE) = vh - velHS;
f = [fc; dN];
outputs = [outputs, gate, N, Ca, vh, vm, velHS];
end
