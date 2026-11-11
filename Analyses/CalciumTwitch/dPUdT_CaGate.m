function [f, outputs] = dPUdT_CaGate(t, PU, params)
%DPUDT_CAGATE  Calcium-gated probe wrapper around dPUdT_CombinedTransitions.
%
%   [F, OUTPUTS] = DPUDT_CAGATE(T, PU, PARAMS) appends ONE state to the core
%   state vector -- N, the non-permissive fraction of thin-filament regulatory
%   units -- and gates attachment by the permissive fraction of FREE sites.
%   The core model is not modified: the gate enters by scaling params.ka,
%   which is exact because every other attachment modifier in the core
%   multiplies ka_eff (stretch activation, A_reg, A_def, velocity Gaussian).
%
%   Thin-filament kinetics (Tewari 2016 / Marzban 2020 / Campbell 2018 form):
%       Jon  = k_on*Ca*N*(1 + K_coop*(1-N))           Ca in uM
%       Joff = k_off*Pu*(1 + K_coop*N),  Pu = max(0, 1 - N - bound)
%       gate = Pu/(1 - bound)
%   Sites held by bound heads cannot switch off (cross-bridge feedback), and at
%   N = 0 the gate is exactly 1, i.e. the calibrated saturating-Ca model.
%
%   PARAMS.tf (thin filament): k_on, k_off, K_coop, and either
%       mode = 'ca'    -> Ca = tf.CaFcn(t)         (Ca-driven)
%       mode = 'clamp' -> gate = tf.gateFcn(t), N frozen (prescribed activation;
%                         measures the cross-bridge bandwidth alone)
%
%   ASSUMPTION (probe only): one binding site per head, so the attached head
%   fraction p1_0 + p2_0 (+ p3_0) is used as the bound-site fraction.
%
%   OUTPUTS = [core outputs, gate, N, Ca].
%
%   See also: dPUdT_CombinedTransitions, caTransient, RunTwitchProbe

n    = numel(PU);
N    = min(1, max(0, PU(n)));
core = PU(1:n-1);
ss   = params.ss;
Ns   = params.NumberOfStates;
bound = params.dS*sum(max(0, core(1:Ns*ss)));
tf   = params.tf;

if strcmp(tf.mode, 'clamp')
    gate = tf.gateFcn(t);  dN = 0;  Ca = NaN;
else
    Ca   = tf.CaFcn(t);
    Pu   = max(0, 1 - N - bound);
    Jon  = tf.k_on*Ca*N*(1 + tf.K_coop*(1 - N));
    Joff = tf.k_off*Pu*(1 + tf.K_coop*N);
    dN   = Joff - Jon;
    gate = Pu/max(1 - bound, eps);
end

p = params;
p.ka = params.ka*gate;
[fc, outputs] = dPUdT_CombinedTransitions(t, core, p);
f = [fc; dN];
outputs = [outputs, gate, N, Ca];
end
