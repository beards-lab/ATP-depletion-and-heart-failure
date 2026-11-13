function C = setupTrisegXB(H, opt)
%SETUPTRISEGXB  Assemble the coupled system: index maps, wall parameters, a
% consistent initial state, mass matrix and Jacobian pattern.
%
% PROVENANCE: [GLUE]. Host parameters and the starting volumes and geometry come
% from hostP21Load (a saved host fit); each wall is a cross-bridge wall from
% setupXbWall, pre-equilibrated at its own end-diastolic length and diastolic Ca.
%
%   C = setupTrisegXB(H)            defaults below
%   C = setupTrisegXB(H, opt)
%   opt.paramSet       cross-bridge snapshot                      'robust1_seed'
%   opt.ss             strain bins (sliding window)               60
%   opt.konMult        calibration factor on the host's k_on      3
%   opt.tauV           time constant of the velocity signal (s)   1e-3
%   opt.hostPassive    host passive law in the model's slot       true  (decision 0004)
%   opt.nativeViscous  keep the model's Maxwell element           false (decision 0004)
%   opt.relTol, opt.absTol   solver tolerances                    1e-3, 1e-6
%   opt.tPre           wall pre-equilibration (s)                 0.5

if nargin < 2, opt = struct(); end
def = struct('paramSet', 'robust1_seed', 'ss', 60, 'konMult', 3, 'tauV', 1e-3, 'hostPassive', true, ...
             'nativeViscous', false, 'relTol', 1e-3, 'absTol', 1e-6, 'tPre', 0.5, 'N0', 1);
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opt, fn{i}), opt.(fn{i}) = def.(fn{i}); end, end

C.H = H;  C.opt = opt;
Ca0 = hostP21Ca(H.tPhase, H);
i0 = 10;  y0 = [H.g0; H.V0];
for k = 1:3
    tf  = struct('k_on', H.tf.k_on*opt.konMult, 'k_off', H.tf.k_off(k), 'K_coop', H.tf.K_coop);
    pas = [];
    if opt.hostPassive, pas = struct('k', H.k_passive(k), 'Lsc0', H.Lsc0, 'gamma', H.gamma); end
    W = setupXbWall(H.Ls0(k), tf, pas, opt, Ca0);
    if k == 1, C.W = W; else, C.W(k) = W; end
    C.iw{k} = i0 + (1:W.nW);  i0 = i0 + W.nW;
    y0 = [y0; W.w0]; %#ok<AGROW>
end
C.n = i0;

% ---- consistent geometry: the four algebraic rows, Newton on g with everything else fixed.
% The stored geometry balanced the HOST's wall stresses; the cross-bridge walls give others.
g = y0(1:4);
for it = 1:30
    r = resG(g, y0, C);
    if norm(r) < 1e-10, break; end
    J = zeros(4);
    for j = 1:4, h = 1e-7*max(1, abs(g(j))); gj = g; gj(j) = gj(j) + h; J(:, j) = (resG(gj, y0, C) - r)/h; end
    dg = -J\r;  lam = 1;
    while lam > 1e-4 && norm(resG(g + lam*dg, y0, C)) > norm(r), lam = lam/2; end
    g = g + lam*dg;
end
C.resG0 = norm(resG(g, y0, C));
y0(1:4) = g;
% the velocity filter starts at rest on the consistent geometry
G = hostP21Geometry(g, H);
for k = 1:3, y0(C.iw{k}(C.W(k).nCore + 2)) = G.Ls(k); end
C.y0 = y0;

% ---- mass matrix: rows 1-4 algebraic (host model_sol.m l.60-64)
C.M = speye(C.n);  C.M(1:4, 1:4) = 0;

% ---- Jacobian pattern. Walls talk to each other only through the geometry:
%   host rows (1-10) depend on g, V and, through the wall stress, on each wall's Lsc and Lf;
%   a wall's rows depend on g (its Ls) and on its own state.
S = sparse(C.n, C.n);  S(1:10, 1:10) = 1;
for k = 1:3
    ik = C.iw{k};  nC = C.W(k).nCore;
    S(1:10, ik(nC + (1:2))) = 1;  S(ik, 1:4) = 1;  S(ik, ik) = 1;
end
C.JPattern = S;
end

function r = resG(g, y, C)
y(1:4) = g;
d = trisegXB_rhs(0, y, C);
r = d(1:4);
end
