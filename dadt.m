function ddt = dadt(~,a,N,s, ds,r_a, r_d, beta, s0)
%% passive attachment as a function of space
if ~isreal(a)
    breakpointhere = 1;
end
u = 1 - ds*sum(a);
ddt = zeros(N, 1);
ddt(1) = + (r_a/ds)*u ; % attach rate
s = (0:length(ddt)-1)'*ds;
% function so that the beta is 1 as a default modifier
% ddt = ddt - (r_d'.*(max(a, 0).^beta).*((s+1).^(beta-1))); % de-attach rate
% ddt = ddt - (r_d'.*max(a, 0).*((s+1).^(beta-1))); % de-attach rate
ddt = ddt - (r_d'.*(max(a, 0))); % de-attach rate

if ~isreal(ddt)
    breakpointhere = 1;
end
