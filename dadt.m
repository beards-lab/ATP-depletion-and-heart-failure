function f = dadt(~,a,N,ds,r_d,r_a, e)
%% passive attachment as a function of space

u = 1 - ds*sum(a);
f = zeros(N,1);
f(1) = + (r_a/ds)*u ; % attach rate
f = f - (r_d'.*min(a, 0).^e); % de-attach rate
