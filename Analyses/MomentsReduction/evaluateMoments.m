function out = evaluateMoments(T, V, params, closure, K)
% EVALUATEMOMENTS  Integrate dMdT_Moments over a velocity table.
%
%   OUT = EVALUATEMOMENTS(T, V, PARAMS, CLOSURE, K) mirrors evaluateModel's
%   segment loop (Vums = V*rsl0*ML per segment) with no window events -- a
%   moments model has no grid to run off. PARAMS is a grid-model params struct
%   from getParams; its PU0 tail is reused so both models start identically.
mp = params;
mp.MomClosure = closure; mp.ss = K; mp.MomScale = 0.005;
n = 7; J = diag(sqrt(1:n-1), 1) + diag(sqrt(1:n-1), -1); [Vv, D] = eig(J);
mp.cachedGH.x = diag(D); mp.cachedGH.w = Vv(1,:)'.^2;
xbKernels(0, mp, 0); mp.kernelsChecked = true;      % validate the flag set once, not per RHS call
tail = params.PU0(2*params.ss+1:end);
assert(numel(tail) == 9, 'expected 9 tail states (A_reg and A_def on, D0 off)');
PU0 = [zeros(2*K,1); tail(:)];
y0 = PU0;

dt = 5e-4;
out.t = []; out.Y = []; nf = 0; tic;
for vs = 1:numel(T)-1
    mp.Vums = V(vs) * mp.rsl0 * mp.ML;
    tspan = unique([T(vs), T(vs)+dt:dt:T(vs+1), T(vs+1)]);
    if numel(tspan) == 2, tspan = linspace(T(vs), T(vs+1), 3); end
    opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'Stats', 'off');
    sol = ode15s(@(t,y) dMdT_Moments(t, y, mp), [T(vs) T(vs+1)], y0, opts);
    nf = nf + sol.stats.nfevals;
    Y = deval(sol, tspan)';
    if vs > 1, tspan = tspan(2:end); Y = Y(2:end,:); end
    out.t = [out.t; tspan(:)]; out.Y = [out.Y; Y];
    y0 = sol.y(:, end);
end
out.wall = toc; out.nfevals = nf;
out.Force = zeros(size(out.t)); out.nodes = zeros(numel(out.t), 2);
for i = 1:numel(out.t)
    seg = max(1, find(out.t(i) >= T(1:end-1), 1, 'last'));
    mp.Vums = V(seg) * mp.rsl0 * mp.ML;
    [~, o] = dMdT_Moments(out.t(i), out.Y(i,:)', mp);
    out.Force(i) = o(1); out.nodes(i,:) = o(8:9);
end
out.SL = out.Y(:, 2*K+3); out.p1_0 = out.Y(:,1); out.p2_0 = out.Y(:,K+1);
out.closure = closure; out.K = K;
end
