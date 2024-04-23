function [force] = diffeq_sim(f, pars, args, dt, dim, frac_pars)
  % This function simulates the caputo differential equation for the function f on [0,T_f]
  % using the method defined in fracp and the time-step size dt, i.e., it solves
  %
  %    S + \delta * D_t^\alpha S = D_t^\alpha f
  %
  %    f - a function we are to take the fractional derivative of, e.g.
  %        f = f(t, arg), where t is the time
  %
  %    par - parameters of the model Sv, such as modulus
  %
  %    arg - e.g. displacement here
  %
  %    delta - scaling paramter on LHS
  %
  %    alpha - fractional derivative order
  %
  %    N - number of prony terms
  %
  %    T_f - the final time of our simulation ...
  %
  %    dt - the time-step size ...
  %
  %    fracp - the structure that stores how to do fractional derivative
  %            approximation as well as the test structure ...
  %
  % by David Nordsletten, 2018
  % updated by Will Zhang, 2024
  %
    % Some constants, maybe make cross sectional area a parameter later
    reference_length = 0.95;
    x_section_area = 1.0;
    % Calculate strains
    lambda = args / reference_length + 1.0;
    % Initialize Caputo Parameters for fractional derivative
    % LHS D_t^\alpha
    fracpL = frac.caputo_approx_initialize(frac_pars.a, frac_pars.Tf, frac_pars.N);
    fracpL.f_prev = f(pars, lambda(1));
    fracpL.Q = zeros(fracpL.N, dim);
    fracpL.delta = frac_pars.d;
    % RHS D_t^\alpha
    fracpR = frac.caputo_approx_initialize(frac_pars.a, frac_pars.Tf, frac_pars.N);
    fracpR.f_prev = f(pars, lambda(1));
    fracpR.Q = zeros(fracpR.N, dim);
    % Determine dimension of variables
    % Set up output arrays
    N_steps = length(args);
    force   = zeros(N_steps, 1);
    % Set up export counter
    k = 1;
    for i = 1:N_steps % go through the main time loop ...
        % compute the right hand side and update, e.g., solves for v in
        %   \hat{v} = D_t^\alpha S_3D
        [v, fracpR] = frac.caputo_approx_iter(f, pars, lambda(i), fracpR, dt);
        % solve the left hand side, i.e. S_star in the paper, e.g., solves for v in
        %   v + \delta * D_t^\alpha v = \hat{v}
        [v, fracpL] = frac.diffeq_approx1_iter(v, fracpL, dt);
        % Filling out the export data arrays ...
        % Cauchy stress
        % \sigma = F (S_3D - p C^-1) F^T
        cauchy_stress = lambda(i)^2 * (v(1) - v(2) / lambda(i)^3);
        force(i) = cauchy_stress * x_section_area;
    end
end