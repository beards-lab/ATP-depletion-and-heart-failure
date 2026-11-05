function s = saLabel(A, D, pH)
%SALABEL  File label of a standalone TW_ATP condition (ADP in 0.1 uM units, GKATP 0).
s = sprintf('sa_adp_A%g_D%05.0f_pH%.1f', A, round(1e4*D), pH);
end
