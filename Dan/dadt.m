function f = dadt(~,a,N,ds,r_d,r_a)

L = 0.60; % domain size (micros)
s  = 0:ds:(L-ds);
u = 1 - ds*sum(a);
f = zeros(N,1);
f(1) = + (r_a/ds)*u ; % attach rate
% f = f - (r_d'.*(a.^2)); % de-attach rate
f = f - (r_d'.*a.*(s'/0.4)); % de-attach rate
