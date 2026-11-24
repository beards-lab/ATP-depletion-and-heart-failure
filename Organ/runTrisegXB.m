function R = runTrisegXB(C, nBeats, opt)
%RUNTRISEGXB  Integrate the coupled system for nBeats and summarise each beat.
%
% PROVENANCE: [GLUE]. One ode15s call per stretch between strain-window slides; a
% slide in any wall is an event that stops the solver, shifts that wall's window
% (shiftXbWall) and restarts. Per-beat summary as the host's own steady-state test
% uses (peak pressures, volumes).
%
%   R = runTrisegXB(C, nBeats)
%   opt.dtOut (s, 0.5e-3), opt.maxWall (s of wall clock, 3600), opt.useJPattern (true),
%   opt.verbose (true)
%   R.t, R.out (columns: see R.col), R.y (states at the output times, single),
%   R.beat (table: one row per beat), R.C (the system as it ended: pivots, last state in C.y0),
%   R.nSteps, nFailed, nFevals, nPD, nShift (1x3), massLost (1x3), wall (s)

if nargin < 3, opt = struct(); end
if ~isfield(opt, 'dtOut'),       opt.dtOut = 0.5e-3; end
if ~isfield(opt, 'maxWall'),     opt.maxWall = 3600; end
if ~isfield(opt, 'useJPattern'), opt.useJPattern = true; end
if ~isfield(opt, 'verbose'),     opt.verbose = true; end

T = C.H.T;  tEnd = nBeats*T;
tOut = (0:opt.dtOut:tEnd)';
y = C.y0;  ts = 0;  tStart = tic;
R.t = [];  R.out = [];  R.y = [];  R.nShift = zeros(1, 3);  R.massLost = zeros(1, 3);
R.nSteps = 0;  R.nFailed = 0;  R.nFevals = 0;  R.nPD = 0;
stopFcn = @(t, y, flag, varargin) toc(tStart) > opt.maxWall;
while ts < tEnd - 1e-12
    o = odeset('Mass', C.M, 'MassSingular', 'yes', 'RelTol', C.opt.relTol, 'AbsTol', C.opt.absTol, ...
               'MaxStep', T/40, 'Events', @(t, yy) ev(yy, C), 'OutputFcn', stopFcn);
    if opt.useJPattern, o = odeset(o, 'JPattern', C.JPattern); end
    sol = ode15s(@(t, yy) trisegXB_rhs(t, yy, C), [ts tEnd], y, o);
    R.nSteps = R.nSteps + sol.stats.nsteps;  R.nFailed = R.nFailed + sol.stats.nfailed;
    R.nFevals = R.nFevals + sol.stats.nfevals;  R.nPD = R.nPD + sol.stats.npds;
    te = sol.x(end);
    tk = tOut(tOut >= ts & tOut <= te);  if ~isempty(R.t), tk = tk(tk > R.t(end)); end
    if ~isempty(tk)
        Yk = deval(sol, tk);  ok = zeros(numel(tk), 46);
        for i = 1:numel(tk), [~, ok(i, :)] = trisegXB_rhs(tk(i), Yk(:, i), C); end
        R.t = [R.t; tk];  R.out = [R.out; ok];  R.y = [R.y; single(Yk')];
    end
    y = sol.y(:, end);  ts = te;
    if ts < tEnd - 1e-12
        if ~isfield(sol, 'ie') || isempty(sol.ie)
            error('runTrisegXB:integration', 'integration stopped at t = %.4g s (beat %.2f) without a window event', ts, ts/T);
        end
        ie = sol.ie(end);  k = ceil(ie/2);  dirn = ie - 2*(k - 1);
        [y(C.iw{k}), C.W(k), lost] = shiftXbWall(y(C.iw{k}), dirn, C.W(k));
        R.nShift(k) = R.nShift(k) + 1;  R.massLost(k) = R.massLost(k) + lost;
    end
    if opt.verbose && ~isempty(R.t) && floor(ts/T + 1e-9) > floor((ts - (te - sol.x(1)))/T + 1e-9)
        fprintf('   beat %2d done, %.0f s\n', floor(ts/T + 1e-9), toc(tStart));
    end
end
R.wall = toc(tStart);  R.wallPerBeat = R.wall/nBeats;
C.y0 = y;  R.C = C;

% ---- columns of R.out (trisegXB_rhs)
w3 = @(i0) i0 + (0:2);
R.col = struct('P_LV', 1, 'P_SA', 2, 'P_SV', 3, 'P_RV', 4, 'P_PA', 5, 'P_PV', 6, 'Q_m', 7, 'Q_a', 8, 'Q_t', 9, 'Q_p', 10, ...
    'Q_SA', 11, 'Q_PA', 12, 'Tm', w3(13), 'Ls', w3(16), 'sigma', w3(19), 'Fact', w3(22), 'Fpas', w3(25), 'gate', w3(28), ...
    'N', w3(31), 'v', w3(34), 'velHS', w3(37), 'LSE', w3(40), 'bound', w3(43), 'Ca', 46);
R.V_LV = double(R.y(:, 5));  R.V_RV = double(R.y(:, 6));

% ---- per-beat summary
B = table();
for b = 1:nBeats
    k = R.t >= (b - 1)*T - 1e-12 & R.t < b*T - 1e-12;
    VL = R.V_LV(k);  VR = R.V_RV(k);  PL = R.out(k, 1);  PR = R.out(k, 4);  [~, iED] = max(VL);  [~, iEDr] = max(VR);
    B = [B; table(b, max(VL), min(VL), PL(iED), max(PL), max(VR), min(VR), PR(iEDr), max(PR), max(R.out(k, 2)), min(R.out(k, 2)), ...
        max(R.out(k, 16)), min(R.out(k, 16)), max(R.out(k, 19)), min(R.out(k, 19)), ...
        'VariableNames', {'beat', 'EDV_LV', 'ESV_LV', 'EDP_LV', 'Pmax_LV', 'EDV_RV', 'ESV_RV', 'EDP_RV', 'Pmax_RV', 'SBP', 'DBP', ...
                          'LsMax_LV', 'LsMin_LV', 'sigmaMax_LV', 'sigmaMin_LV'})]; %#ok<AGROW>
end
B.SV_LV = B.EDV_LV - B.ESV_LV;
R.beat = B;
end

function [value, isterminal, direction] = ev(y, C)
% strain-window events of the three walls (as Model/evaluateModel.m, movingWindow)
value = zeros(6, 1);
for k = 1:3
    W = C.W(k);
    pivo = (-y(C.iw{k}(W.nCore + 1)) + W.P.LXBpivot)/2;
    value(2*k - 1) = -pivo + W.P.s(end) - W.P.A1AttachmentWidth;
    value(2*k)     = -pivo + W.P.s(1) + W.P.A1AttachmentWidth;
end
direction = repmat([-1; 1], 3, 1);  isterminal = ones(6, 1);
end
