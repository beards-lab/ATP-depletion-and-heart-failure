%% TEST THE MASS MATRIX AND JACOBIAN FOR SPEEDING UP THE ISMULATION

x0 = out.PU(end, :)';      % a typical state

params.Vums = 0;
rhsFun = @(t,x) dPUdT_CombinedTransitions(t,x, params);       % your complex dxdt function

%%
[Jpat, Mdiag] = analyzeSystemRHSFun(rhsFun, x0);

% [Jpat Mdiag] = safeJPattern(rhsFun, x0);

%% Then feed to ODE solver with the mass matrix
opts = odeset('JPattern', Jpat, 'Mass', Mdiag, 'RelTol',1e-6,'AbsTol',1e-8);
tic
for i = 1:100
    [t,x] = ode15s(rhsFun, [i i+1], x(end, :), opts);
end
toc

%% compare to case without the mass matrix
opts = odeset('AbsTol',1e-4, 'RelTol', 1e-2);
tic
for i = 1:100
    [t,x] = ode15s(rhsFun, [i i+1], x(end, :), opts);
end
toc

%% no difference :((



function [Jpat, Mdiag] = analyzeSystemRHSFun(rhsFun, x0)
% rhsFun = @(t,x) dxdt, your costly function
% x0 = typical state vector (column)
%
% Outputs:
%   Jpat = sparse logical matrix of Jacobian nonzeros
%   Mdiag = diagonal vector for safe mass matrix

    N = length(x0);
    f0 = rhsFun(0, x0);       % evaluate once at typical point
    J = sparse(N,N);

    epsFD = 1e-8;  % finite difference step
    fprintf('Estimating Jacobian sparsity pattern...\n');

    % Finite-difference Jacobian column by column
    for i = 1:N
        if mod(i, round(N/10))==0
            fprintf('  Column %d / %d\n', i, N);
        end
        x_eps = x0;
        x_eps(i) = x_eps(i) + epsFD;
        fi = rhsFun(0, x_eps);
        J(:,i) = fi - f0;
    end

    % Normalize by finite difference step
    J = J / epsFD;

    % Logical sparsity pattern
    Jpat = spones(J);

    fprintf('Nonzero Jacobian entries: %d / %d\n', nnz(Jpat), N*N);

    % --------- Mass Matrix ---------
    % Approximate safe diagonal mass scaling:
    % Take the max absolute row sums
    rowSums = sum(abs(J),2);
    tol = 1e-12;
    Mdiag = max(rowSums, tol);   % avoid zero entries
    % sparse diagonal mass matrix
    Mdiag = sparse(1:N, 1:N, Mdiag, N, N);

    fprintf('Mass matrix: diagonal scaling based on row sums of |J|\n');

end

function [Jpat Mdiag] = safeJPattern(rhsFun, x0)

    Mdiag = speye(length(x0));

    N = length(x0);
    f0 = rhsFun(0, x0);

    % Very small eps can cause numerical noise; use moderate FD step:
    epsFD = 1e-6;

    % Preallocate dense logical matrix (converted to sparse at the end)
    Jmask = false(N, N);

    for j = 1:N
        x_eps = x0;
        x_eps(j) = x_eps(j) + epsFD;
        fj = rhsFun(0, x_eps);

        df = fj - f0;

        % If ANY change is detected above noise, mark nonzero
        Jmask(:, j) = abs(df) > 1e-12;
    end

    % Always include diagonal as a fallback (safe!)
    Jmask = Jmask | eye(N) > 0;

    Jpat = sparse(Jmask);
end
