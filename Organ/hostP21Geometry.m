function G = hostP21Geometry(g, H)
%HOSTP21GEOMETRY  TriSeg wall geometry: midwall shape -> strain -> sarcomere length.
%
% PROVENANCE: HOST. Equations copied verbatim from
%   beards-lab/Multiscale-Hx-P21-Model, model/model.m @ 80f86d3 (2025-10-31), lines 219-246.
% Only the packing of inputs and outputs is new. If the host changes, re-copy those
% lines; Organ/Tests/TestOrganParity.m (part 0) checks this block against the host's
% own stored outputs.
%
%   g = [xm_LV; xm_SEP; xm_RV; ym] (cm)      H from hostP21Load
%   G.Ls  1x3 sarcomere length (um) of LV, SEP, RV; also Vm, Am, Cm, z, eps (1x3 each)

xm_LV = g(1);  xm_SEP = g(2);  xm_RV = g(3);  ym = g(4);
Amref_LV = H.Amref(1);  Amref_SEP = H.Amref(2);  Amref_RV = H.Amref(3);
Vw_LV = H.Vw(1);  Vw_SEP = H.Vw(2);  Vw_RV = H.Vw(3);
Lsref = H.Lsref;

% ---- model.m l.219-246 -------------------------------------------------------
% Volume of spherical cap formed by midwall surface (cm^3)
Vm_LV  = -(pi / 6) * xm_LV  * (xm_LV^2  + 3 * ym^2);
Vm_SEP = (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2);
Vm_RV  = (pi / 6) * xm_RV  * (xm_RV^2  + 3 * ym^2);

% Surface area of midwall surface (cm^2)
Am_LV  = pi * (xm_LV^2  + ym^2);
Am_SEP = pi * (xm_SEP^2 + ym^2);
Am_RV  = pi * (xm_RV^2  + ym^2);

% Curvature of midwall surface (cm^(-1))
Cm_LV  = -2 * xm_LV  / (xm_LV^2  + ym^2);
Cm_SEP =  2 * xm_SEP / (xm_SEP^2 + ym^2);
Cm_RV  =  2 * xm_RV  / (xm_RV^2  + ym^2);

% Ratio of wall thickness to midwall radius of curvature (dimensionless)
z_LV  = 3 * Cm_LV  * Vw_LV  / (2 * Am_LV);
z_SEP = 3 * Cm_SEP * Vw_SEP / (2 * Am_SEP);
z_RV  = 3 * Cm_RV  * Vw_RV  / (2 * Am_RV);

% Myofiber strain (dimensionless)
eps_LV  = 0.5 * log( Am_LV  / Amref_LV  ) - (1/12) * z_LV^2  - 0.019 * z_LV^4;
eps_SEP = 0.5 * log( Am_SEP / Amref_SEP ) - (1/12) * z_SEP^2 - 0.019 * z_SEP^4;
eps_RV  = 0.5 * log( Am_RV  / Amref_RV  ) - (1/12) * z_RV^2  - 0.019 * z_RV^4;

% Sarcomere length (um)
Ls_LV  = Lsref * exp(eps_LV);
Ls_SEP = Lsref * exp(eps_SEP);
Ls_RV  = Lsref * exp(eps_RV);
% ---- end of the copied block -------------------------------------------------

G.Vm  = [Vm_LV, Vm_SEP, Vm_RV];
G.Am  = [Am_LV, Am_SEP, Am_RV];
G.Cm  = [Cm_LV, Cm_SEP, Cm_RV];
G.z   = [z_LV, z_SEP, z_RV];
G.eps = [eps_LV, eps_SEP, eps_RV];
G.Ls  = [Ls_LV, Ls_SEP, Ls_RV];
end
