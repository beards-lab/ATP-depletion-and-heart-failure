function Ca_i = hostP21Ca(t, H)
%HOSTP21CA  The host's analytic cytosolic Ca transient (uM).
%
% PROVENANCE: HOST. Copied from beards-lab/Multiscale-Hx-P21-Model, model/model.m
% @ 80f86d3 (2025-10-31), lines 195-211 (a raised cosine up over k_TS*T and down over
% k_TR*T). The amplitude and the diastolic level are NOT the literals of that file
% (2 and 0.01610): they are read from the stored fit (hostP21Load), whose transient
% runs 0.161 to 2.161 uM, i.e. the stored fits were made with a diastolic 0.161.

T = H.T;
TS_v = H.k_TS * T;                 % Time to maximal systole
TR_v = H.k_TR * T;                 % Time from maximal systole to relaxation
Ca_amplitude = H.Ca.amplitude;
Ca_diastole  = H.Ca.diastole;

% ---- model.m l.204-211 -------------------------------------------------------
tc = mod(t,T);
if tc >= 0 && tc < TS_v
    Ca_i = Ca_diastole + Ca_amplitude*0.5*(1 - cos(pi*(tc)/TS_v));
elseif tc >= TS_v && tc < TR_v + TS_v
    Ca_i = Ca_diastole + Ca_amplitude*0.5*(1 + cos(pi*(tc - TS_v)/TR_v));
else
    Ca_i = Ca_diastole;
end
end
