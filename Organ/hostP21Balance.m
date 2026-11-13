function [res, dV, out] = hostP21Balance(g, G, V, sigma, H)
%HOSTP21BALANCE  TriSeg tension balance, ventricular pressures and the lumped circulation.
%
% PROVENANCE: HOST. Equations copied verbatim from
%   beards-lab/Multiscale-Hx-P21-Model, model/model.m @ 80f86d3 (2025-10-31),
%   lines 296-327 (tension, pressure), 356-362 (flows), 431-434 (the four algebraic
%   equations) and 487-492 (volume derivatives).
% The ONE change against the host: the wall stress `sigma` is an input. The host
% computes it as -Kse*(Lsc - Ls) (l.291-293); here it comes from the cross-bridge
% model's own series element (Organ/xbWall.m).
%
%   g      [xm_LV; xm_SEP; xm_RV; ym] (cm)          G from hostP21Geometry(g, H)
%   V      [V_LV; V_RV; V_SV; V_PV; V_SA; V_PA] (cm^3), the host's states 8-13 in its order
%   sigma  total wall stress (kPa) of LV, SEP, RV
%   res    4x1 residuals of the algebraic equations (host rows 1-4; zero when balanced)
%   dV     6x1 volume derivatives (cm^3/s), same order as V
%   out    [P_LV P_SA P_SV P_RV P_PA P_PV  Q_m Q_a Q_t Q_p Q_SA Q_PA  Tm_LV Tm_SEP Tm_RV] (kPa, cm^3/s, kPa cm)

xm_LV = g(1);  xm_SEP = g(2);  xm_RV = g(3);  ym = g(4);
V_LV = V(1);  V_RV = V(2);  V_SV = V(3);  V_PV = V(4);  V_SA = V(5);  V_PA = V(6);
sigma_LV = sigma(1);  sigma_SEP = sigma(2);  sigma_RV = sigma(3);
Vm_LV = G.Vm(1);  Vm_SEP = G.Vm(2);  Vm_RV = G.Vm(3);
Am_LV = G.Am(1);  Am_SEP = G.Am(2);  Am_RV = G.Am(3);
z_LV  = G.z(1);   z_SEP  = G.z(2);   z_RV  = G.z(3);
Vw_LV = H.Vw(1);  Vw_SEP = H.Vw(2);  Vw_RV = H.Vw(3);
C_SA = H.C.SA;  C_SV = H.C.SV;  C_PA = H.C.PA;  C_PV = H.C.PV;
R_SA = H.R.SA;  R_PA = H.R.PA;
R_m_vlv = H.R.m;  R_a_vlv = H.R.a;  R_t_vlv = H.R.t;  R_p_vlv = H.R.p;

% ---- model.m l.296-327 -------------------------------------------------------
% Representative midwall tension (kPa cm)
Tm_LV  = (Vw_LV  * sigma_LV  / (2 * Am_LV))  * (1 + (z_LV^2)/3  + (z_LV^4)/5);
Tm_SEP = (Vw_SEP * sigma_SEP / (2 * Am_SEP)) * (1 + (z_SEP^2)/3 + (z_SEP^4)/5);
Tm_RV  = (Vw_RV  * sigma_RV  / (2 * Am_RV))  * (1 + (z_RV^2)/3  + (z_RV^4)/5);

% Axial midwall tension component (kPa cm)
Tx_LV  = - Tm_LV  * 2 * xm_LV  * ym / (xm_LV^2  + ym^2);
Tx_SEP = Tm_SEP * 2 * xm_SEP * ym / (xm_SEP^2 + ym^2);
Tx_RV  = Tm_RV  * 2 * xm_RV  * ym / (xm_RV^2  + ym^2);

% Radial midwall tension component (kPa cm)
Ty_LV  = Tm_LV  * (-xm_LV^2  + ym^2) / (xm_LV^2  + ym^2);
Ty_SEP = Tm_SEP * (-xm_SEP^2 + ym^2) / (xm_SEP^2 + ym^2);
Ty_RV  = Tm_RV  * (-xm_RV^2  + ym^2) / (xm_RV^2  + ym^2);

% Ventricular pressure (kPa)
ptrans_LV = 2 * Tx_LV / ym;
ptrans_RV = 2 * Tx_RV / ym;
P_LV      = -ptrans_LV;
P_RV      = ptrans_RV;

% Calculate Pericardial Pressure
P_peri = 0; %exp(s*(Vh/Vh0 - 1))/7.5; % convert mmHg to kPa
P_LV = P_peri + P_LV;
P_RV = P_peri + P_RV;

% Pressure (kPa)
P_SA = V_SA / C_SA;
P_SV = V_SV / C_SV;
P_PA = V_PA / C_PA;
P_PV = V_PV / C_PV;

% ---- model.m l.356-362 -------------------------------------------------------
Q_m = max((P_PV - P_LV) / R_m_vlv, 0);
Q_t = max((P_SV - P_RV) / R_t_vlv, 0);
Q_a = max((P_LV - P_SA) / R_a_vlv, 0);
Q_p = max((P_RV - P_PA) / R_p_vlv, 0);

Q_SA = (P_SA - P_SV) / R_SA;
Q_PA = (P_PA - P_PV) / R_PA;

% ---- model.m l.431-434: the four algebraic equations (mass-matrix rows of zeros) --
dxm_LV  = -V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV;
dxm_SEP = Tx_LV + Tx_SEP + Tx_RV;
dxm_RV  = V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV;
dym     = Ty_LV + Ty_SEP + Ty_RV;

% ---- model.m l.487-492 -------------------------------------------------------
dV_LV = Q_m - Q_a;
dV_SA = Q_a - Q_SA;
dV_SV = Q_SA - Q_t;
dV_RV = Q_t - Q_p;
dV_PA = Q_p - Q_PA;
dV_PV = Q_PA - Q_m;
% ---- end of the copied blocks ------------------------------------------------

res = [dxm_LV; dxm_SEP; dxm_RV; dym];
dV  = [dV_LV; dV_RV; dV_SV; dV_PV; dV_SA; dV_PA];
if nargout > 2
    out = [P_LV, P_SA, P_SV, P_RV, P_PA, P_PV, Q_m, Q_a, Q_t, Q_p, Q_SA, Q_PA, Tm_LV, Tm_SEP, Tm_RV];
end
end
