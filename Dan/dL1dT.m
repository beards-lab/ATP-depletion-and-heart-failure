function f = dL1dT(~,X,V)
L = X(1);
L1 = X(2);

mu = 3;
ks = 6;
f(1,:) = V;
f(2,:) = ks*(L - L1)/mu;
