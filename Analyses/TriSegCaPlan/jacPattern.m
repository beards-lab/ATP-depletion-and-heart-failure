function [S, info] = jacPattern(P, Y, t)
%JACPATTERN  Union of the finite-difference Jacobian sparsity over several states.
%
%   [S, INFO] = JACPATTERN(P, Y, T): Y is n x m (m states, columns), T 1 x m times.
%   A pattern from one state misses couplings that are numerically zero there (empty
%   bins), so the union over states from a whole beat is taken. INFO: density, the
%   dense rows/cols (> 20 nonzeros) and the number of column groups ode15s needs
%   (colgroup-style greedy count, i.e. RHS calls per Jacobian).
n = size(Y, 1);  S = sparse(n, n);
for k = 1:size(Y, 2)
    y = Y(:, k);  f0 = dPUdT_HostWall(t(k), y, P);
    for j = 1:n
        yj = y;  h = 1e-6*max(1, abs(y(j)));  yj(j) = yj(j) + h;
        S(:, j) = S(:, j) | (dPUdT_HostWall(t(k), yj, P) ~= f0); %#ok<SPRIX>
    end
end
S = S | speye(n);
info.density = nnz(S)/n^2;
info.denseRows = find(sum(S, 2) > 20)';
info.denseCols = find(sum(S, 1) > 20);
info.nGroupsFull = n;
info.nGroups = greedyGroups(S);
% groups if the dense rows were supplied analytically (TriSegCoupling section 3)
Sr = S;  Sr(info.denseRows, :) = 0;  Sr = Sr | speye(n);
info.nGroupsNoDenseRows = greedyGroups(Sr);
end

function g = greedyGroups(S)
% columns that share no row can be perturbed together (Curtis-Powell-Reid)
n = size(S, 2);  grp = zeros(1, n);  g = 0;
rowsUsed = {};
for j = 1:n
    rj = find(S(:, j));
    placed = false;
    for k = 1:g
        if ~any(rowsUsed{k}(rj))
            grp(j) = k;  rowsUsed{k}(rj) = true;  placed = true;  break;
        end
    end
    if ~placed
        g = g + 1;  grp(j) = g;  rowsUsed{g} = false(size(S, 1), 1);  rowsUsed{g}(rj) = true;
    end
end
end
