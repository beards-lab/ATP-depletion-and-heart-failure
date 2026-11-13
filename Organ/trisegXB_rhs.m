function [dy, out] = trisegXB_rhs(t, y, C)
%TRISEGXB_RHS  The coupled system: P21 TriSeg heart and circulation with the full
% cross-bridge model as the sarcomere of all three walls. One monolithic DAE
% (decision 0003): rows 1-4 are algebraic (mass-matrix zeros).
%
% PROVENANCE: this file only ASSEMBLES. The equations live in
%   [HOST]  hostP21Geometry, hostP21Ca, hostP21Balance  (copied from the host, marked there)
%   [XB]    Model/dPUdT_CombinedTransitions             (called through xbWall)
%   [GLUE]  xbWall (Ca gate, velocity signal, state packing)
%
%   y = [g (4): xm_LV xm_SEP xm_RV ym;  V (6): V_LV V_RV V_SV V_PV V_SA V_PA;
%        wall LV; wall SEP; wall RV]        (wall layout in xbWall)
%   C from setupTrisegXB
%   out (second output) = [P (6), Q (6), Tm (3), Ls (3), sigma (3), F_active (3),
%        F_passive (3), gate (3), N (3), v (3), velHS (3), LSE (3), bound (3), Ca]

g = y(1:4);  V = y(5:10);
G  = hostP21Geometry(g, C.H);                         % [HOST]
Ca = hostP21Ca(t + C.H.tPhase, C.H);                  % [HOST] (t = 0 is the start of the stored beat)

dy = zeros(C.n, 1);  sigma = zeros(3, 1);
if nargout > 1, aux = zeros(3, 8); end
for k = 1:3
    ik = C.iw{k};
    if nargout > 1
        [dy(ik), sigma(k), aux(k, :)] = xbWall(t, y(ik), G.Ls(k), Ca, C.W(k));   % [XB] + [GLUE]
    else
        [dy(ik), sigma(k)] = xbWall(t, y(ik), G.Ls(k), Ca, C.W(k));
    end
end

if nargout > 1
    [dy(1:4), dy(5:10), ho] = hostP21Balance(g, G, V, sigma, C.H);              % [HOST]
    out = [ho, G.Ls, sigma', aux(:)', Ca];
else
    [dy(1:4), dy(5:10)] = hostP21Balance(g, G, V, sigma, C.H);
end
end
