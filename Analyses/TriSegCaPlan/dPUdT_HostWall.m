function [f, outputs] = dPUdT_HostWall(t, PU, params)
%DPUDT_HOSTWALL  Probe wrapper: one TriSeg wall's sarcomere, driven open-loop by the host.
%
%   Wraps the UNCHANGED dPUdT_CombinedTransitions. State = [core; N; vf]. Added,
%   all outside the model file:
%   1. Ca activation: N (non-permissive thin-filament fraction), exactly the
%      factorized gate of Analyses/CalciumTwitch/dPUdT_CaGate.m (gate scales
%      params.ka; gate = 1 reproduces the calibrated model).
%   2. Host kinematics: the sarcomere follows the host's Ls(t), dSL/dt = v_host(t).
%      The series element and the contractile element are the model's own
%      (dLSE/dt = v_host - velHS). The ONLY thing that changes between modes is the
%      velocity the velocity-driven mechanisms see (A_reg drive, A2Reattaching gate,
%      mu2*vel, dynamic passive, catch bond, stretch activation):
%        params.host.velMode = 'imposed'   v_host(t)   (what evaluateModel does)
%                              'zero'      0 (velocity-driven mechanisms off)
%                              'velHS'     contractile-element velocity, two passes
%                                          (an algebraic loop: velHS depends on vel)
%                              'filtHS'    vf, a first-order low-pass of velHS
%                              'filtHost'  vf, a first-order low-pass of v_host
%      SL and LSE are integrated with v_host in every mode, so the kinematics are
%      identical across modes.
%   3. vf: dvf/dt = (x - vf)/params.host.tauV, x = velHS ('filtHS') or v_host (else).
%
%   params.host: Ls(t)/v(t) handles of one wall (hostSignals), velMode, tauV (s).
%   params.tf:   k_on, k_off, K_coop, CaFcn(t) (uM).
%   OUTPUTS = [core outputs, gate, N, Ca, v_host, v_mech, velHS].

n    = numel(PU);
vf   = PU(n);
N    = min(1, max(0, PU(n-1)));
core = PU(1:n-2);
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
    case {'filtHS', 'filtHost'}
        vm = vf;
    otherwise
        error('velMode %s', params.host.velMode);
end
p.Vums = vm;
[fc, outputs] = dPUdT_CombinedTransitions(t, core, p);
velHS = vm - fc(iLSE);
% kinematics: the sarcomere follows the host whatever the mechanisms saw
fc(iSL)  = vh;
fc(iLSE) = vh - velHS;
if strcmp(params.host.velMode, 'filtHS'), dvf = (velHS - vf)/params.host.tauV;
else,                                     dvf = (vh - vf)/params.host.tauV;
end
f = [fc; dN; dvf];
outputs = [outputs, gate, N, Ca, vh, vm, velHS];
end
