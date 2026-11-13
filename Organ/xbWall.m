function [dw, sigma, aux] = xbWall(t, w, Ls, Ca, W)
%XBWALL  One wall's sarcomere: the cross-bridge model, given its length from outside.
%
% PROVENANCE, block by block (decision 0003, point 5):
%   [XB]    Model/dPUdT_CombinedTransitions.m is CALLED, never copied. Whatever that
%           file becomes as the cross-bridge model advances is what runs here.
%   [GLUE]  everything else in this file:
%           - state packing: the model's state holds SL and LSE; in a wall the
%             sarcomere length Ls comes from the geometry, so the wall carries the
%             contractile length Lsc = SL - LSE instead and hands the model
%             SL = Ls, LSE = Ls - Lsc;
%           - the Ca gate: the factorised thin-filament gate of
%             Analyses/CalciumTwitch/dPUdT_CaGate.m (N = non-permissive fraction;
%             the gate scales the attachment rate ka; gate = 1 is the calibrated
%             model). Its constants k_on, k_off, K_coop are the HOST's
%             (model.m l.110-113, 412-413), with one calibration factor on k_on;
%           - the velocity the model's velocity-driven mechanisms see: a first-order
%             low-pass of dLs/dt with time constant tauV, obtained without
%             differentiating the geometry from one extra state Lf (low-passed Ls):
%             dLf/dt = (Ls - Lf)/tauV,  v = (Ls - Lf)/tauV  (#84 replaces this).
%
%   w = [core (the model's state vector; its SL and LSE slots are unused here and
%        kept at their initial values); Lsc; Lf; N]
%   Ls (um), Ca (uM) from the host;  W from setupTrisegXB (W.P = model parameters)
%   dw     time derivative of w
%   sigma  total wall stress (kPa): the model's series-element force, its output 1
%   aux    [F_active, F_passive, gate, N, v, velHS, LSE, bound]

nC  = W.nCore;
Lsc = w(nC + 1);  Lf = w(nC + 2);  N = min(1, max(0, w(nC + 3)));

% [GLUE] kinematics handed to the model
v    = (Ls - Lf)/W.tauV;
core = w(1:nC);
core(W.iSL)  = Ls;
core(W.iLSE) = Ls - Lsc;

% [GLUE] Ca gate
P     = W.P;
bound = P.dS*sum(max(0, core(1:W.nBins)));
Pu    = max(0, 1 - N - bound);
Jon   = W.tf.k_on*Ca*N*(1 + W.tf.K_coop*(1 - N));
Joff  = W.tf.k_off*Pu*(1 + W.tf.K_coop*N);
gate  = Pu/max(1 - bound, eps);
P.ka   = P.ka*gate;
P.Vums = v;

% [XB] the cross-bridge model, unchanged
[fc, o] = dPUdT_CombinedTransitions(t, core, P);

% [GLUE] back to the wall's state: the model returns dLSE/dt = v - velHS
velHS = v - fc(W.iLSE);
fc(W.iSL)  = 0;                       % unused slots
fc(W.iLSE) = 0;
dw    = [fc; velHS; v; Joff - Jon];   % dLsc/dt = velHS; dLf/dt = (Ls - Lf)/tauV = v
sigma = o(1);
if nargout > 2
    aux = [o(2), o(3), gate, N, v, velHS, Ls - Lsc, bound];
end
end
