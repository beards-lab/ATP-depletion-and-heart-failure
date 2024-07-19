function [g outputs]  = dXdT(t,x,Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,L_0,nd,kDf, V)
s  = (0:1:Nx-1)'.*ds;

pu = reshape( x(1:(Ng+1)*Nx), [Nx,Ng+1]);
if size(x, 1) > Nx*(Ng+1) + 1 % we have pa
    pa = reshape( x((Ng+1)*Nx+1:2*(Ng+1)*Nx), [Nx,Ng+1]);
    g = zeros(2*(Ng+1)*Nx + 1,1);
else
    pa = [];
    g = zeros((Ng+1)*Nx + 1,1);
end
L = max(0,x(end));

% Calculate the un-attached chain velocities for every pu(s,n) entry
deltaF = kd*max(0,(L-s)/L_0).^nd - Fp;
Vp = deltaF/mu;

ij = (1:Nx)' + Nx*(0:Ng); % matrix of Ng X Nx indices over all elements

% Attach/dettach rectifier connector - only for Ca
if ~isempty(pa)
    detach = kD*x(ij + (Ng+1)*Nx).*(1 + kDf*max(0, deltaF+Fp));
    g(ij)             = g(ij)             - kA*x(ij) + detach;
    g(ij + (Ng+1)*Nx) = g(ij + (Ng+1)*Nx) + kA*x(ij) - detach;

% g(ij)             = g(ij)             - kA*x(ij) + kD*x(ij + (Ng+1)*Nx).*max(0,deltaF);
% g(ij + (Ng+1)*Nx) = g(ij + (Ng+1)*Nx) + kA*x(ij) - kD*x(ij + (Ng+1)*Nx).*max(0,deltaF);

end

% UPWIND differencing for sliding the midpoint (+ direction) for pu
% like:
% k(1) = k(1) - pu(1)*V(1)/ds;
% k(2:Nx) = k(2:Nx) - (pu(2:Nx)*V(2:Nx) - pu(2:Nx)*Vp(1:Nx-1))/ds

indxs_plus = ij(2:Nx,:); % indexes of look ahead - the ij(...) works like reshape in flattenized array g and pu, which should be s x Ng otherwise
indxs_minus = ij(1:Nx-1,:); % indexes of look behind

% positive velocities - extending
g(indxs_plus) = g(indxs_plus) ...
    - (1/ds)*pu(indxs_plus).*max(0,Vp(2:Nx, :)) + (1/ds)*pu(indxs_minus).*max(0,Vp(1:Nx-1,:));

% negative velocities - shrinking
g(indxs_minus) = g(indxs_minus) ...
    - (1/ds)*pu(indxs_plus).*min(0,Vp(2:Nx,:)) + (1/ds)*pu(indxs_minus).*min(0,Vp(1:Nx-1,:));

% boundary conditions
g(ij(end,:)) = g(ij(end,:)) + (1/ds)*(pu(end,:).*min(0,Vp(end,:)))';
g(ij(1,:)) = g(ij(1,:)) - (1/ds)*(pu(1,:).*max(0,Vp(1,:)))';

% unfolding rate for pu states
UR = RU.*pu(ij(:,1:Ng)); % rate of probabilty transitions from n to n+1 states

% refolding rate - speed up 
if RF > 0
    FR = RF.*pu(ij(:,2:Ng+1)); % rate of probabilty transitions from n+1 to n states
else
    FR = 0;
end

if t > 0.055
    a = 1;
end
g(ij(:,2:(Ng+1))) = g(ij(:,2:(Ng+1))) + UR - FR;
g(ij(:,1:Ng))     = g(ij(:,1:Ng))     - UR + FR;

% quantifying speed
if isa(V, 'function_handle')
    g(end) = V(t);
else
    g(end) = V;
end

% unfolding for pa states
if false || ~isempty(pa)
    NxNg = Nx*(Ng+1);
    UR = RU.*pa(ij(:,1:Ng)); % rate of probabilty transitions from n to n+1 states
    g(ij(:,2:(Ng+1))+NxNg) = g(ij(:,2:(Ng+1))+NxNg) + UR - FR;
    g(ij(:,1:Ng)+NxNg)     = g(ij(:,1:Ng)+NxNg)     - UR + FR;

    % shifting by velocity of the whole strand - UPWIND differencing for attached
    Vs = g(end);
    % positive velocities - extending
    g(indxs_plus + NxNg) = g(indxs_plus + NxNg) ...
        - (1/ds)*pa(indxs_plus)*max(0,Vs) + (1/ds)*pa(indxs_minus)*max(0,Vs);

    % negative velocities - shrinking
    g(indxs_minus + NxNg) = g(indxs_minus + NxNg) ...
        - (1/ds)*pa(indxs_plus)*min(0,Vs) + (1/ds)*pa(indxs_minus)*min(0,Vs);
    
    % boundary conditions
    g(ij(end,:) + NxNg) = g(ij(end,:) + NxNg) + (1/ds)*pa(end,:)'*min(0,Vs);
    g(ij(1,:) + NxNg) = g(ij(1,:) + NxNg) - (1/ds)*pa(1,:)'*max(0,Vs);

    
    Ftit = sum(Fp.*pu*ds, 2)' + sum(Fp.*pa*ds, 2)';
else
    Ftit = sum(Fp.*pu*ds, 2)';
end


try
catch
    wat
end
outputs = [sum(Ftit) Ftit];