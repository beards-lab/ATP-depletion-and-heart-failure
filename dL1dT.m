function f = dL1dT(~,X,V, ksmu)
L = X(1);
L1 = X(2);
% mu = 10;
% ks = 4;
f(1,:) = V;
% f(2,:) = ks*(L - L1)/mu;
f(2,:) = ksmu*(L - L1);