function g = dXdT(~,x,Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,V)
s  = (0:1:Nx-1)'.*ds;

pu = reshape( x(1:(Ng+1)*Nx), [Nx,Ng+1]);
% pa = reshape( x((Ng+1)*Nx+1:2*(Ng+1)*Nx), [Nx,Ng+1]);
L = x(end);

% Calculate the un-attached chain velocities for every pu(s,n) entry
deltaF = (kS*max(0,L-s-Ls0).^nS - Fc);
Vc = deltaF/mu;

g = zeros((Ng+1)*Nx + 1,1);

ij = (1:Nx)' + Nx*(0:Ng); % matrix of Ng X Nx indices over all elements

% Attach/dettach rectifier connector
% g(ij)             = g(ij)             - kA*x(ij) + kD*x(ij + (Ng+1)*Nx);
% g(ij + (Ng+1)*Nx) = g(ij + (Ng+1)*Nx) + kA*x(ij) - kD*x(ij + (Ng+1)*Nx);
% g(ij)             = g(ij)             - kA*x(ij) + kD*x(ij + (Ng+1)*Nx).*max(0,deltaF);
% g(ij + (Ng+1)*Nx) = g(ij + (Ng+1)*Nx) + kA*x(ij) - kD*x(ij + (Ng+1)*Nx).*max(0,deltaF);

% UPWIND differencing for sliding (+ direction) for pu 
g(ij(1,:)) = g(ij(1,:)) - (1/ds)*(pu(1,:).*Vc(1,:))';
g(ij(2:Nx,:)) = g(ij(2:Nx,:)) - (1/ds)*pu(ij(2:Nx,:)).*max(0,Vc(2:end,:)) + ...
                (1/ds)*pu(ij(2:Nx,:)-1).*max(0,Vc(1:end-1,:));  

% unfolding for pu states
UR = RU.*pu(ij(:,1:Ng)); % rate of probabilty transitions from n to n+1 states
g(ij(:,2:(Ng+1))) = g(ij(:,2:(Ng+1))) + UR;
g(ij(:,1:Ng))     = g(ij(:,1:Ng))     - UR;

% unfolding for pa states
% NxNg = Nx*(Ng+1);
% UR = RU.*pa(ij(:,1:Ng)); % rate of probabilty transitions from n to n+1 states
% g(ij(:,2:(Ng+1))+NxNg) = g(ij(:,2:(Ng+1))+NxNg) + UR;
% g(ij(:,1:Ng)+NxNg)     = g(ij(:,1:Ng)+NxNg)     - UR;


g(end) = V;
