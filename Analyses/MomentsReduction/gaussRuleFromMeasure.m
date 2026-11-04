function [x, w] = gaussRuleFromMeasure(s, p, N)
% GAUSSRULEFROMMEASURE  N-node Gauss rule of a discrete measure (s, p).
%
%   The N-node rule reproduces the first 2N moments of the measure exactly, so
%   it is precisely what a QMOM closure with 2N moments would reconstruct
%   (assuming the moment inversion is exact). Built with the Stieltjes
%   procedure on the measure itself, which is stable where a Hankel/Wheeler
%   inversion from raw moments is not. N is reduced if the support is smaller.
s = s(:); p = max(0, p(:));
keep = p > 1e-14 * max(max(p), realmin);
s = s(keep); p = p(keep);
m0 = sum(p);
if m0 <= 0, x = 0; w = 0; return; end
N = min(N, numel(s));
a = zeros(N,1); b = zeros(N,1);
pPrev = zeros(size(s)); pk = ones(size(s)); nPrev = 1;
for k = 1:N
    nk = sum(p .* pk.^2);
    a(k) = sum(p .* s .* pk.^2) / nk;
    if k > 1, b(k) = nk / nPrev; end
    pNext = (s - a(k)) .* pk - b(k) * pPrev;
    pPrev = pk; pk = pNext; nPrev = nk;
end
J = diag(a) + diag(sqrt(b(2:N)), 1) + diag(sqrt(b(2:N)), -1);
[V, D] = eig(J);
x = diag(D); w = m0 * V(1,:)'.^2;
end
