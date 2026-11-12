function R = runHostBeats(P, nBeats, opt)
%RUNHOSTBEATS  Integrate one host-driven wall for nBeats, with evaluateModel's sliding window.
%
%   R = RUNHOSTBEATS(P, NBEATS, OPT) integrates dPUdT_HostWall from P.PU0 over
%   [-OPT.tPre, NBEATS*T]. t < 0 is a diastolic pre-equilibration: Ls frozen at Ls(0)
%   (v = 0) and Ca frozen at Ca(0). The strain window slides exactly as in
%   evaluateModel (same event functions on the pivot vs s(1)/s(end) +-
%   A1AttachmentWidth, same shift of floor(ss/3) bins, vacated bins zeroed), and
%   every slide restarts ode15s. Mass lost at each slide is recorded.
%
%   OPT.tPre (s, default 0.5), OPT.dtOut (s, default 0.5e-3), OPT.JPattern (sparse
%   n x n or [] = full finite-difference Jacobian; see jacPattern).
%   R: t, F (total), Factive, Fpassive, SL, gate, N, Ca, v_host, v_mech, velHS,
%   A_reg, bound, nShift, tShift, massLost, wall (s), wallPerBeat, nSteps, nFailed.

if nargin < 3, opt = struct(); end
if ~isfield(opt, 'tPre'),     opt.tPre = 0.5; end
if ~isfield(opt, 'dtOut'),    opt.dtOut = 0.5e-3; end
if ~isfield(opt, 'JPattern'), opt.JPattern = []; end

T  = P.hostT;
ss = P.ss;  Ns = P.NumberOfStates;
iSL = Ns*ss + 3;  iLSE = Ns*ss + 4;  iA = Ns*ss + 8;
nds = floor(ss/3);

% ---- pre-equilibration: frozen host ----
Pp = P;
L0 = P.host.Ls(0);  C0 = P.tf.CaFcn(0);
Pp.host.Ls = @(t) L0;  Pp.host.v = @(t) 0;  Pp.tf.CaFcn = @(t) C0;
y0 = P.PU0(:);
tic;
[~, Y] = ode15s(@dPUdT_HostWall, [-opt.tPre 0], y0, odeset('RelTol', P.odeRelTol, 'AbsTol', P.odeAbsTol), Pp);
y0 = Y(end, :)';

% ---- beats ----
tEnd = nBeats*T;
tOut = 0:opt.dtOut:tEnd;
ts = 0;  R.t = [];  Yall = [];  R.nShift = 0;  R.tShift = [];  R.dirShift = [];  R.massLost = [];
R.nSteps = 0;  R.nFailed = 0;
while ts < tEnd - 1e-12
    o = odeset('RelTol', P.odeRelTol, 'AbsTol', P.odeAbsTol, 'Events', @win, 'Stats', 'off');
    if ~isempty(opt.JPattern)
        o = odeset(o, 'JPattern', opt.JPattern);
    end
    tspan = unique([ts, tOut(tOut > ts), tEnd]);
    sol = ode15s(@dPUdT_HostWall, [ts tEnd], y0, o, P);
    R.nSteps = R.nSteps + sol.stats.nsteps;  R.nFailed = R.nFailed + sol.stats.nfailed;
    te = sol.x(end);
    tk = tspan(tspan <= te);
    if ~isempty(tk)
        Yk = deval(sol, tk)';
        R.t = [R.t; tk(:)];  Yall = [Yall; Yk]; %#ok<AGROW>
    end
    y0 = sol.y(:, end);  ts = te;
    if isfield(sol, 'ie') && ~isempty(sol.ie) && ts < tEnd - 1e-12
        ie = sol.ie(end);
        m0 = P.dS*sum(y0(1:Ns*ss));
        y0 = shift(y0, ie, ss, Ns, nds);
        if ie == 1, P.LXBpivot = P.LXBpivot - nds*P.dS*2; else, P.LXBpivot = P.LXBpivot + nds*P.dS*2; end
        R.massLost(end+1) = m0 - P.dS*sum(y0(1:Ns*ss));
        R.nShift = R.nShift + 1;  R.tShift(end+1) = ts;  R.dirShift(end+1) = ie;
        % the de-duplicated output grid must not repeat ts
        tOut = tOut(tOut > ts);
    end
end
R.wall = toc;
R.wallPerBeat = R.wall/nBeats;
[R.t, iu] = unique(R.t);  Yall = Yall(iu, :);

% ---- outputs (re-evaluated; pivot changes are piecewise, so recompute per segment) ----
% Force etc. are evaluated with the pivot that was active at each time.
pivots = localPivots(R, P, nds);
n = numel(R.t);
R.F = zeros(n,1); R.gate = R.F; R.N = R.F; R.Ca = R.F; R.vh = R.F; R.vm = R.F; R.velHS = R.F;
R.SL = Yall(:, iSL);  R.LSE = Yall(:, iLSE);  R.A_reg = Yall(:, iA);
R.bound = P.dS*sum(Yall(:, 1:Ns*ss), 2);
for k = 1:n
    Pk = P;  Pk.LXBpivot = pivots(k);
    [~, ok] = dPUdT_HostWall(R.t(k), Yall(k, :)', Pk);
    R.F(k) = ok(1);  R.Factive(k,1) = ok(2);  R.Fpassive(k,1) = ok(3);
    R.gate(k) = ok(end-5);  R.N(k) = ok(end-4);  R.Ca(k) = ok(end-3);
    R.vh(k) = ok(end-2);  R.vm(k) = ok(end-1);  R.velHS(k) = ok(end);
end
R.T = T;  R.ss = ss;  R.velMode = P.host.velMode;
R.Y = single(Yall);  R.pivot = pivots;
end

function [value, isterminal, direction] = win(~, y, P)
ss = P.ss;  Ns = P.NumberOfStates;
SL = y(Ns*ss + 3);  LSE = y(Ns*ss + 4);
pivo = (-(SL - LSE) + P.LXBpivot)/2;
value = [-pivo + P.s(end) - P.A1AttachmentWidth; -pivo + P.s(1) + P.A1AttachmentWidth];
direction = [-1; 1];
isterminal = [1; 1];
end

function y = shift(y, ie, ss, Ns, nds)
for k = 0:Ns-1
    b = k*ss;
    if ie == 1      % move right
        y(b+1:b+ss-nds) = y(b+1+nds:b+ss);  y(b+ss-nds+1:b+ss) = 0;
    else            % move left
        y(b+1+nds:b+ss) = y(b+1:b+ss-nds);  y(b+1:b+nds) = 0;
    end
end
end

function piv = localPivots(R, P, nds)
% pivot active at each output time: start pivot + the shifts that happened before it
p0 = P.LXBpivot;
for j = 1:R.nShift   % undo all shifts to get the starting pivot
    if R.dirShift(j) == 1, p0 = p0 + nds*P.dS*2; else, p0 = p0 - nds*P.dS*2; end
end
piv = p0*ones(numel(R.t), 1);
for j = 1:R.nShift
    d = -nds*P.dS*2;  if R.dirShift(j) == 2, d = -d; end
    piv(R.t >= R.tShift(j)) = piv(R.t >= R.tShift(j)) + d;
end
end
