function R = runXbWallOpenLoop(W, LsFcn, CaFcn, tspan, opt)
%RUNXBWALLOPENLOOP  One wall (xbWall) driven by a prescribed Ls(t) and Ca(t).
%
% PROVENANCE: [GLUE]. Used for pre-equilibration and for the parity tests: it runs
% the SAME xbWall the coupled system runs, with the same window slides.
%
%   R = runXbWallOpenLoop(W, LsFcn, CaFcn, [t0 t1], opt)
%   opt.dtOut output step (s, default 0.5e-3); opt.maxStep (default inf)
%   R: t, sigma, aux (columns as xbWall), nShift, massLost, nSteps, nFailed, wall (s),
%      W (with the pivot it ended on), wEnd

if nargin < 5, opt = struct(); end
if ~isfield(opt, 'dtOut'),   opt.dtOut = 0.5e-3; end
if ~isfield(opt, 'maxStep'), opt.maxStep = inf; end
t0 = tspan(1);  t1 = tspan(2);
tOut = (t0:opt.dtOut:t1)';  if tOut(end) < t1, tOut(end+1) = t1; end
w = W.w0;  ts = t0;  tic;
R.t = [];  R.sigma = [];  R.aux = [];  R.nShift = 0;  R.massLost = 0;  R.nSteps = 0;  R.nFailed = 0;
while ts < t1 - 1e-12
    o = odeset('RelTol', W.P.odeRelTol, 'AbsTol', W.P.odeAbsTol, 'MaxStep', opt.maxStep, 'Events', @(t, y) ev(y, W));
    sol = ode15s(@(t, y) xbWall(t, y, LsFcn(t), CaFcn(t), W), [ts t1], w, o);
    R.nSteps = R.nSteps + sol.stats.nsteps;  R.nFailed = R.nFailed + sol.stats.nfailed;
    te = sol.x(end);
    tk = tOut(tOut >= ts & tOut <= te);  if ~isempty(R.t), tk = tk(tk > R.t(end)); end
    if ~isempty(tk)
        Yk = deval(sol, tk);
        s = zeros(numel(tk), 1);  a = zeros(numel(tk), 8);
        for i = 1:numel(tk), [~, s(i), a(i, :)] = xbWall(tk(i), Yk(:, i), LsFcn(tk(i)), CaFcn(tk(i)), W); end
        R.t = [R.t; tk];  R.sigma = [R.sigma; s];  R.aux = [R.aux; a];
    end
    w = sol.y(:, end);  ts = te;
    if ts < t1 - 1e-12
        if ~isfield(sol, 'ie') || isempty(sol.ie)
            error('runXbWallOpenLoop:integration', 'integration stopped at t = %.4g s without a window event', ts);
        end
        [w, W, lost] = shiftXbWall(w, sol.ie(end), W);
        R.nShift = R.nShift + 1;  R.massLost = R.massLost + lost;
    end
end
R.wall = toc;  R.W = W;  R.wEnd = w;
end

function [value, isterminal, direction] = ev(y, W)
% the strain-window event of Model/evaluateModel.m (movingWindow), on the wall's Lsc
pivo = (-y(W.nCore + 1) + W.P.LXBpivot)/2;
value = [-pivo + W.P.s(end) - W.P.A1AttachmentWidth; -pivo + W.P.s(1) + W.P.A1AttachmentWidth];
direction = [-1; 1];  isterminal = [1; 1];
end
