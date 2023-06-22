function f = dadt(~,a,N,ds,r_d,r_a, e)
%% passive attachment as a function of space

u = 1 - ds*sum(a);
f = zeros(N,1);
f(1) = + (r_a/ds)*u ; % attach rate
s = (0:length(f)-1)*ds;
f = f - (r_d'.*max(a, 0).*s'.^e); % de-attach rate
