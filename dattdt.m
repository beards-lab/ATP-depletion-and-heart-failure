function datt = dattdt(t,att, s_ap, r_a, r_d, F2vals)
% passive attachment as a function of space

u = 1 - sum(att);
datt = zeros(length(att), 1);
datt(s_ap) = + (r_a)*u ; % attach rate at the current center
% datt = datt - r_d.*att; % de-attach rate
datt = datt - max(att, 0).*(r_d+abs(F2vals)); % de-attach rate depends on force
