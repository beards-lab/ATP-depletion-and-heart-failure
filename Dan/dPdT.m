function g = dPdT(~,x,N,V,del_U,alphaL)

alpha0   = 0*1; % spontaneous (not length-dependent) unfolding rate
beta  = 0*1e3;
P = x(1:N); % unfolded number probabilities
L = x(end);

% unfolding
P0 = 1 - sum(P); % probability that whold chain is folded
g(1,:)   = (alphaL*L.^4 + alpha0)*P0;
g(2:N,:) = (alphaL*(max(0,L - (1:N-1)'*del_U ).^4) + alpha0 ).*P(1:N-1) ;
g(1:N-1) = g(1:N-1) - (alphaL*(max(0,L - (1:N-1)'*del_U ).^4) + alpha0).*P(1:N-1) ;
% g(2:N,:) = alpha*exp(1*(L - (1:N-1)'*del_U)./del_U ).*P(1:N-1) ;
% g(1:N-1) = g(1:N-1) - alpha*exp(1*(L - (1:N-1)'*del_U)./del_U ).*P(1:N-1) ;

% refolding
% g(1:N)   = g(1:N) - beta.*P;
% g(1:N-1)  = g(1:N-1) + beta.*P(2:N);
g(1:N)   = g(1:N) - beta*(max(0, (1:N)'*del_U - L ).^1).*P;
g(1:N-1) = g(1:N-1) + beta*(max(0, (2:N)'*del_U - L ).^1).*P(2:N);

g(N+1,:) = V;