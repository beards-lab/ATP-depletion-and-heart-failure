function f = atpGateFactor(A, Aref, K, n)
%ATPGATEFACTOR  ATP-ref-normalised gate factor, shared by the analysis scripts.
%
%   F = ATPGATEFACTOR(A, AREF, K, N) returns the multiplier applied to a gated rate
%   at [ATP] = A, normalised so that F == 1 exactly at A = AREF:
%
%       f(A) = [A^n/(K^n + A^n)] * [(K^n + Aref^n)/Aref^n]
%            = hill(A) * Corr_ref
%
%   K = Inf (or K <= 0 / NaN) returns the power law (A/Aref)^n, which is the
%   K -> Inf limit of the expression above. N is the Hill coefficient; N < 0 makes
%   the step ATP-INHIBITED (f > 1 below Aref).
%
%   Evaluated in log space so large |N| and large K cannot overflow -- the model
%   itself (evaluateModel) uses the direct expression, which agrees to eps over the
%   parameter range of interest and is checked against this function by RunHillGateConversion.
%
%   A may be a vector; K and N are scalars.

if ~isfinite(K) || K <= 0
    f = (A ./ Aref) .^ n;
    return;
end

lK  = n * log(K);
lA  = n * log(A);
lAr = n * log(Aref);

% log f = lA - logsumexp(lK, lA) + logsumexp(lK, lAr) - lAr
f = exp(lA - iLSE(lK, lA) + iLSE(lK, lAr) - lAr);
end

function y = iLSE(x1, x2)
m = max(x1, x2);
y = m + log(exp(x1 - m) + exp(x2 - m));
end
