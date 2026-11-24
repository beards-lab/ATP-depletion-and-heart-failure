%% TestOrganParity  The three checks that guard the coupled system (decision 0003, point 6).
%
%   0  HOST block:  hostP21Geometry / hostP21Balance / hostP21Ca reproduce the host's
%                   own stored outputs (strain, ventricular pressures, Ca, and the four
%                   algebraic equations balance) on every stored sample of the last beat.
%   A  XB wall:     xbWall, with the Ca gate held open and the model's native passive,
%                   reproduces Model/evaluateModel on a hold-shorten-hold-stretch protocol.
%   B  XB wall:     xbWall driven by the host's recorded Ls(t) and Ca(t) reproduces the
%                   open-loop probe of Analyses/TriSegCaPlan (dPUdT_HostWall) on the same beat.
%
%   matlab -batch "run('Organ/Tests/TestOrganParity.m')"        (about 3 min)
% Tolerances: 0.2 kPa for A and B (the acceptance values of #97).

here = fileparts(mfilename('fullpath'));  root = fileparts(fileparts(here));
cd(root);  addpath(genpath(fullfile(root, 'Model')));  addpath(genpath(fullfile(root, 'Auxiliary')));
addpath(fullfile(root, 'params'));  addpath(fullfile(root, 'Organ'));  addpath(fullfile(root, 'Analyses', 'TriSegCaPlan'));
if ~exist('HOSTROOT', 'var'), HOSTROOT = fullfile(fileparts(root), 'Multiscale-Hx-P21-Model'); end
TOL = 0.2;  pass = true(1, 3);
opt = struct('paramSet', 'robust1_seed', 'ss', 60, 'tauV', 1e-3, 'relTol', 1e-3, 'absTol', 1e-6, 'nativeViscous', true, 'tPre', 0.5, 'N0', 1);

%% 0. host block against the host's stored outputs
H = hostP21Load('Nx51', HOSTROOT);  r = H.ref;  n = numel(r.t);
e = zeros(n, 6);
for i = 1:n
    g = r.g(i, :)';  G = hostP21Geometry(g, H);
    [res, ~, o] = hostP21Balance(g, G, r.V(i, :)', r.sigma(i, :)', H);
    e(i, :) = [max(abs(G.eps - r.eps(i, :))), abs(o(1) - r.P(i, 1)), abs(o(4) - r.P(i, 4)), abs(hostP21Ca(r.t(i) + H.tPhase, H) - r.Ca(i)), ...
               max(abs(res([1 3]))), max(abs(res([2 4])))];
end
fprintf('0  host block, %d stored samples: max |eps| error %.2e; |P_LV| %.2e kPa (range %.1f); |P_RV| %.2e kPa; |Ca| %.2e uM;\n', n, max(e(:, 1)), max(e(:, 2)), max(r.P(:, 1)), max(e(:, 3)), max(e(:, 4)));
fprintf('   algebraic equations on the stored state: volume rows %.2e cm^3, tension rows %.2e kPa cm (wall tension up to %.2f)\n', max(e(:, 5)), max(e(:, 6)), max(r.sigma(:))*max(H.Vw)/2);
pass(1) = max(e(:, 1)) < 1e-8 && max(e(:, 2)) < 2e-3*max(r.P(:, 1)) && max(e(:, 3)) < 2e-3*max(r.P(:, 1)) && max(e(:, 4)) < 1e-4;

%% A. xbWall against evaluateModel (gate open, native passive, fast velocity signal)
Tseg = [0 0.5 0.6 0.8 0.85 1.05];  vML = [0 -1 0 2 0];  SL0 = 2.2;
oA = opt;  oA.tauV = 1e-5;  oA.N0 = 0;  oA.tPre = 0;   % a 10 us velocity signal: as good as the imposed velocity
W = setupXbWall(SL0, struct('k_on', 1, 'k_off', 0, 'K_coop', 0), [], oA);
vum = vML*W.P.ML;  Lk = SL0 + [0, cumsum(vum .* diff(Tseg))];
LsF = @(t) interp1(Tseg, Lk, min(max(t, 0), Tseg(end)));
RA = runXbWallOpenLoop(W, LsF, @(t) 1, [0 Tseg(end)], struct('maxStep', 1e-3));
Pc = W.P;  Pc.Velocity = vML;  Pc.g = [];  Pc.mods = {};
[~, oE] = evaluateModel(@dPUdT_CombinedTransitions, Tseg, Pc);
[tu, iu] = unique(oE.t(:));  Fe = interp1(tu, oE.Force(iu), RA.t);
k = RA.t >= 0.45 & isfinite(Fe);
% the velocity signal lags a step in velocity by tauV; leave 1 ms after each step out
for ts = Tseg(2:end-1), k = k & ~(RA.t >= ts & RA.t < ts + 1e-3); end
dA = RA.sigma(k) - Fe(k);  tk = RA.t(k);  [~, im] = max(abs(dA));
fprintf('A  xbWall vs evaluateModel: max |dF| %.4f kPa at t = %.4f s, rms %.4f kPa (force %.1f..%.1f kPa), gate %.6f..%.6f, %d slides\n', ...
    max(abs(dA)), tk(im), rms(dA), min(Fe(k)), max(Fe(k)), min(RA.aux(:, 3)), max(RA.aux(:, 3)), RA.nShift);
for seg = 1:numel(Tseg) - 1
    ks = tk >= Tseg(seg) & tk < Tseg(seg + 1);
    if any(ks), fprintf('      segment %d (%.2f-%.2f s, %+g ML/s): max |dF| %.4f kPa\n', seg, Tseg(seg), Tseg(seg + 1), vML(seg), max(abs(dA(ks)))); end
end
pass(2) = max(abs(dA)) < TOL;

%% B. xbWall against the open-loop probe on the host beat (velocity signal: 1 ms low-pass of dLs/dt)
Hs = hostSignals('Nx51', HOSTROOT);
Pp = setupWall(Hs, 'LV', struct('ss', 60, 'velMode', 'filtHost', 'tauV', opt.tauV, 'konMult', 3));
Rp = runHostBeats(Pp, 2, struct('tPre', 0.5));
tf = struct('k_on', H.tf.k_on*3, 'k_off', H.tf.k_off(1), 'K_coop', H.tf.K_coop);
W  = setupXbWall(Hs.LV.Ls(0), tf, [], opt, Hs.Ca(0));
RB = runXbWallOpenLoop(W, Hs.LV.Ls, Hs.Ca, [0 2*Hs.T]);
Fp = interp1(Rp.t, Rp.F, RB.t);  k = isfinite(Fp);
dB = RB.sigma(k) - Fp(k);
fprintf('B  xbWall vs open-loop probe, 2 beats: max |dF| %.4f kPa, rms %.4f kPa (force %.1f..%.1f kPa); slides %d vs %d\n', ...
    max(abs(dB)), rms(dB), min(Fp(k)), max(Fp(k)), RB.nShift, Rp.nShift);
pass(3) = max(abs(dB)) < TOL;

fprintf('RESULT: 0 %s, A %s, B %s\n', pf(pass(1)), pf(pass(2)), pf(pass(3)));

function s = pf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
