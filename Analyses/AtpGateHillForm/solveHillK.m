function [K, n, ok] = solveHillK(g2, A, Aref, n)
%SOLVEHILLK  Km that makes a saturable ATP gate reproduce a given low-ATP factor.
%
%   [K, N, OK] = SOLVEHILLK(G2, A, AREF, N) returns the half-saturating [ATP] K such
%   that the ATP_ref-normalised Hill gate with coefficient N equals G2 at [ATP] = A:
%
%       g2 = (A/Aref)^n * (K^n + Aref^n)/(K^n + A^n)
%
%   Closed form, from kappa = K^n:
%
%       kappa = A^n (1 - g2) / (g2 - (A/Aref)^n),      K = kappa^(1/n)
%
%   Feasibility: kappa must be > 0, i.e. G2 must lie strictly between 1 (K^n -> 0,
%   gate switched off because the step is saturated at both concentrations) and
%   (A/Aref)^n (K^n -> Inf, the unsaturated power-law limit). So a suppressive gate
%   (G2 < 1) needs N > log(G2)/log(A/Aref) and an enhancing one (G2 > 1) needs
%   N < that value. OK is false when no positive K exists for the requested N.
%
%   With N omitted the classic Michaelis choice is used: N = +1 for a suppressive
%   gate, N = -1 (ATP-inhibited step) for an enhancing one.
%
%   See also atpGateFactor.

r = A / Aref;
if nargin < 4 || isempty(n)
    if g2 < 1, n = 1; else, n = -1; end
end

kappa = A^n * (1 - g2) / (g2 - r^n);
ok    = isfinite(kappa) && kappa > 0;
if ok
    K = kappa ^ (1/n);
else
    K = NaN;
end
end
