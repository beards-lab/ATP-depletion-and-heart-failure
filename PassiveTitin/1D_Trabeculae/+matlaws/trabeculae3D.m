function [S_11, S_22, S_33] = trabeculae3D(par, lambda)
    % S_3D = \sum_{k \in {f,s,n}} \mu_k (W_1 I_{kk} - 1) e_k\otimes e_k
    % W_1 = exp(b*(I_C - 3))
    % For 1D, we simplify with
    %  \mu_f = \mu * (1 - \kappa + \kappa/3)
    %  \mu_s = u_n = \kappa /3 * \mu
    % here \kappa -> 0 is anisotropic, \kappa -> 1 is isotropic
    % par:
    %  par(1) = \mu
    %  par(2) = b
    %  par(3) = kappa
    C_11 = lambda * lambda;
    I_C = C_11 + 2.0 / lambda - 3.0;
    W_1 = exp(par(2) * I_C);
    u_f = par(1) * (1.0 - par(3)* 2.0/3.0);
    u_s = par(1) * (par(3) / 3.0);
    % We export all 3 components for hydrostatic pressure calculations.
    S_11 = u_f * (W_1 * C_11   - 1.0);
    S_22 = u_s * (W_1 / lambda - 1.0);
    S_33 = S_22;
 end