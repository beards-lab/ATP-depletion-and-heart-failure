%% TestCaActivation.m -- regression + parity checks for the Ca-activation layer
%
% Checks the production path (params.CaModel -> calciumActivation.m, hooked
% into dPUdT_CombinedTransitions) against
%   T1  byte-identity: CaModel = 'none' vs the pre-change RHS/getParams
%       (HEAD copies *_ref.m on the path, if present), on non-trivial states
%   T2  twoState parity vs the probe wrapper dPUdT_CaGate, same integrator:
%       steady force-pCa and a Rice 2008 1 Hz twitch
%   T3  clamp: 10 ms availability pulse (probe E1: 56 % of F_ss)
%   T4  evaluateModel path: beat knots merged, out.Ca / out.CaGate / out.CaState
%   T5  guards: UseCa = 1, MEX + CaModel, unknown CaModel
%   T6  RHS cost: none vs twoState vs ref
%   T7  myosim: = twoState at full overlap (Rice twitch); saturating force vs twoState; MEX guard
% Params: params/fastgate_best.m (the probe's set), isometric, SL fixed ends.

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
res  = fullfile(here, 'results');
opts = odeset('AbsTol', 1e-4, 'RelTol', 1e-3, 'MaxStep', 1e-3);
dt   = 5e-4;
TF   = struct('k_on', 101.185, 'k_off', 723.852, 'K_coop', 9.6846);  % Marzban 2020
base = loadParams('fastgate_best');
P    = getParams(base, [], true, false);  P.Vums = 0;  P.v = 0;
pass = true;

%% T1 -- byte identity of CaModel = 'none' ------------------------------------------
if exist('getParams_ref', 'file') && exist('dPUdT_CombinedTransitions_ref', 'file')
    for snap = {'fastgate_best', fullfile(root, 'params', 'ExampleHighATP.m')}
        b  = loadParams(snap{1});
        Pn = getParams(b, [], true, false);      Pn.Vums = 0;
        Pr = getParams_ref(b, [], true, false);  Pr.Vums = 0;
        assert(isequal(Pn.PU0, Pr.PU0), 'PU0 differs');
        % non-trivial states: a short isometric activation from rest
        [~, Y] = ode15s(@dPUdT_CombinedTransitions, [0 0.02 0.05 0.2], Pn.PU0(:), opts, Pn);
        dmax = 0;
        for k = 1:size(Y, 1)
            for vel = [0 -0.5 0.5]
                Pn.Vums = vel;  Pr.Vums = vel;
                [fn, on, rn] = dPUdT_CombinedTransitions(0.1, Y(k,:)', Pn);
                [fr, orf, rr] = dPUdT_CombinedTransitions_ref(0.1, Y(k,:)', Pr);
                dmax = max([dmax, max(abs(fn - fr)), max(abs(on - orf)), max(abs(rn - rr))]);
                assert(numel(fn) == numel(fr) && numel(on) == numel(orf), 'size differs');
            end
        end
        [~, snapName] = fileparts(snap{1});
        fprintf('T1 %-24s none vs HEAD: max |d(f, outputs, rates)| = %g over %d states x 3 velocities\n', ...
            snapName, dmax, size(Y, 1));
        pass = pass && dmax == 0;
    end
else
    fprintf('T1 skipped (no *_ref.m on the path)\n');
end

%% T2 -- twoState vs probe wrapper ---------------------------------------------------
S   = load(fullfile(res, 'probe_E1.mat'), 'rest');  rest = S.rest;   % [core; N = 1]
nC  = numel(rest) - 1;
mkP = @(extra) setParams(base, extra);
tfx = {'CaModel', 'twoState', 'tf_k_on', TF.k_on, 'tf_k_off', TF.k_off, 'tf_K_coop', TF.K_coop};
Pg  = P;  Pg.tf = TF;  Pg.tf.mode = 'ca';

pCa = [6.0 5.6 5.4 5.0 4.51];
F_new = nan(size(pCa));  F_prb = F_new;
for k = 1:numel(pCa)
    Pk = mkP([tfx, {'CaInput', 'const', 'pCa', pCa(k)}]);
    assert(numel(Pk.PU0) == nC + 1 && Pk.ca.idx == nC + 1, 'Ca state index');
    [~, ~, O] = integrate(@dPUdT_CombinedTransitions, Pk, rest, [0 1.5], dt, opts);
    F_new(k) = O(end, 1);
    Pg.tf.CaFcn = @(t) 10^(6 - pCa(k));
    [~, ~, Og] = integrate(@dPUdT_CaGate, Pg, rest, [0 1.5], dt, opts);
    F_prb(k) = Og(end, 1);
end
rel = max(abs(F_new - F_prb)./abs(F_prb));
fprintf('T2 force-pCa  new: %s kPa\n             probe: %s kPa   max rel diff %.2g\n', ...
    mat2str(F_new, 5), mat2str(F_prb, 5), rel);
pass = pass && rel < 1e-6;

% Rice 2008 1 Hz, 1 s diastolic pre-equilibration, 3 beats; probe E3 timeline
Pk = mkP([tfx, {'CaInput', 'analytic', 'CaTransientKind', 'rice2008', 'CaPeriod', 1, 'CaFirstBeat', 1}]);
knots = unique([0, 1 + (0:3), 1 + (0:2) + 0.005]);
[t, ~, O] = integrate(@dPUdT_CombinedTransitions, Pk, rest, knots, dt, opts);
mN = twitchMetrics(t, O(:,1), 3);
dia = caTransient('rice2008', 1 - 1e-4, struct('T', 1));
Pg.tf.CaFcn = @(t) dia;
[~, Yg] = integrate(@dPUdT_CaGate, Pg, rest, [0 1], dt, opts);
Pg.tf.CaFcn = @(t) caTransient('rice2008', t, struct('T', 1));
[tg, ~, Og] = integrate(@dPUdT_CaGate, Pg, Yg(end,:)', unique([0:3, (0:2) + 0.005]), dt, opts);
mG = twitchMetrics(tg, Og(:,1), 2);
fprintf(['T2 Rice 1 Hz  new: Fdia %.4f Fdev %.4f kPa TTP %.2f RT50 %.2f ms | ' ...
         'probe: Fdia %.4f Fdev %.4f kPa TTP %.2f RT50 %.2f ms\n'], ...
    mN.Fd, mN.Fdev, 1e3*mN.TTP, 1e3*mN.RT50, mG.Fd, mG.Fdev, 1e3*mG.TTP, 1e3*mG.RT50);
pass = pass && abs(mN.Fdev - mG.Fdev) < 1e-3*mG.Fdev && abs(mN.TTP - mG.TTP) < 1e-3;
E3 = load(fullfile(res, 'probe_E3.mat'), 'E3');
fprintf('   saved probe E3 (older working copy): Fdia %.2f Fdev %.2f kPa TTP %.1f RT50 %.1f ms\n', ...
    E3.E3(1).Fd, E3.E3(1).Fdev, 1e3*E3.E3(1).TTP, 1e3*E3.E3(1).RT50);

%% T3 -- clamp pulse -----------------------------------------------------------------
Pc = mkP({'CaModel', 'clamp', 'CaPeriod', 10, 'CaFirstBeat', 0, 'CaClampDuration', 0.010});
assert(isempty(Pc.ca.idx) && numel(Pc.PU0) == nC, 'clamp adds no state');
[t, ~, O] = integrate(@dPUdT_CombinedTransitions, Pc, rest(1:nC), [0 0.010 0.6], dt, opts);
mC = twitchMetrics(t, O(:,1), 0);
Pgc = P;  Pgc.tf = TF;  Pgc.tf.mode = 'clamp';  Pgc.tf.gateFcn = @(t) double(t < 0.010);
[tg, ~, Og] = integrate(@dPUdT_CaGate, Pgc, rest, [0 0.010 0.6], dt, opts);
mCg = twitchMetrics(tg, Og(:,1), 0);
E1 = load(fullfile(res, 'probe_E1.mat'), 'E1');
fprintf('T3 clamp 10 ms  new: Fdev %.4f kPa TTP %.2f ms | probe: Fdev %.4f kPa TTP %.2f ms | saved E1: %.0f %% of F_ss\n', ...
    mC.Fdev, 1e3*mC.TTP, mCg.Fdev, 1e3*mCg.TTP, 100*E1.E1.pulse(1).frac);
pass = pass && abs(mC.Fdev - mCg.Fdev) < 1e-6*mCg.Fdev;

%% T4 -- evaluateModel path (knots + stored outputs) ----------------------------------
Pe = setParams(base, [tfx, {'CaInput', 'analytic', 'CaTransientKind', 'rice2008', ...
    'CaPeriod', 0.5, 'CaFirstBeat', 0.2, 'Velocity', 0}]);
Pe.PU0 = rest;
[~, out] = evaluateModel(@dPUdT_CombinedTransitions, [0 1.2], Pe);
needKn = [0.2 0.205 0.7 0.705];
hit = arrayfun(@(k) any(abs(out.t - k) < 1e-12), needKn);
fprintf('T4 evaluateModel: %d pts, knots hit %s, Ca %.3f..%.3f uM, gate max %.3f, N %.3f..%.3f, peak F %.2f kPa\n', ...
    numel(out.t), mat2str(hit), min(out.Ca), max(out.Ca), max(out.CaGate), ...
    min(out.CaState), max(out.CaState), max(out.Force));
pass = pass && all(hit) && isfield(out, 'CaState');

%% T5 -- guards ----------------------------------------------------------------------
g = {};
b = base;  b.UseCa = 1;
try, getParams(b, [], true, false); g{end+1} = 'UseCa=1 NOT rejected'; catch e, assert(strcmp(e.identifier, 'getParams:UseCa')); end
b = base;  b.UseCa = 0;  q = getParams(b, [], true, false);
if isfield(q, 'UseCa'), g{end+1} = 'UseCa=0 not removed'; end
b = base;  b.CaModel = 'twoState';  b.UseCompiledMex = true;
try, getParams(b, [], true, false); g{end+1} = 'MEX NOT rejected'; catch e, assert(strcmp(e.identifier, 'getParams:CaMex')); end
b = base;  b.CaModel = 'land';
try, getParams(b, [], true, false); g{end+1} = 'unknown model NOT rejected'; catch e, assert(strcmp(e.identifier, 'getParams:CaModel')); end
if isempty(g), fprintf('T5 guards OK\n'); else, fprintf('T5 FAIL: %s\n', strjoin(g, '; ')); pass = false; end

%% T6 -- RHS cost ----------------------------------------------------------------------
Pk = mkP([tfx, {'CaInput', 'analytic', 'CaTransientKind', 'rice2008', 'CaPeriod', 1, 'CaFirstBeat', 1}]);
x0 = rest;  x0(1:P.ss) = 1e-3;  xN = x0(1:nC);
B = 15;  m = 400;  tm = zeros(B, 3);
hasRef = exist('dPUdT_CombinedTransitions_ref', 'file') == 2;
if hasRef, Pr = getParams_ref(base, [], true, false); Pr.Vums = 0; end
for k = 1:200, dPUdT_CombinedTransitions(1.01, xN, P); dPUdT_CombinedTransitions(1.01, x0, Pk); end
for b = 1:B
    tic; for k = 1:m, dPUdT_CombinedTransitions(1.01, xN, P);  end; tm(b,1) = toc/m;
    tic; for k = 1:m, dPUdT_CombinedTransitions(1.01, x0, Pk); end; tm(b,2) = toc/m;
    if hasRef
        tic; for k = 1:m, dPUdT_CombinedTransitions_ref(1.01, xN, Pr); end; tm(b,3) = toc/m;
    end
end
md = 1e6*median(tm);
fprintf('T6 RHS median: none %.1f us | twoState+rice %.1f us (+%.1f %%) | HEAD ref %.1f us\n', ...
    md(1), md(2), 100*(md(2)/md(1) - 1), md(3));

%% T7 -- CaModel 'myosim' (MyoSim form of the same kinetics) -----------------------------------
% (a) with UseOverlap = 0 (overlap 1) it must reproduce twoState; (b) at saturating Ca both keep
% the calibrated force; (c) the MEX guard covers it
b0 = base;  b0.UseOverlap = 0;
rice = {'CaInput', 'analytic', 'CaTransientKind', 'rice2008', 'CaPeriod', 1, 'CaFirstBeat', 0.2};
Pt = setParams(b0, [tfx, rice]);  Pm = setParams(b0, [tfx, rice, {'CaModel', 'myosim'}]);
assert(numel(Pm.PU0) == numel(Pt.PU0) && Pm.PU0(end) == 0 && Pt.PU0(end) == 1, 'myosim state: N_on = 0 at rest');
kn = unique([0, 0.2, 0.205, 1.2]);
[~, ~, Ot] = integrate(@dPUdT_CombinedTransitions, Pt, Pt.PU0(:), kn, dt, opts);
[~, ~, Om] = integrate(@dPUdT_CombinedTransitions, Pm, Pm.PU0(:), kn, dt, opts);
dTM = max(abs(Ot(:, 1) - Om(:, 1)))/max(Ot(:, 1));
Fs = zeros(1, 3);  mods = {'none', 'twoState', 'myosim'};
for k = 1:3
    Pk = setParams(base, [tfx, {'CaModel', mods{k}, 'CaInput', 'const', 'pCa', 4.51}]);
    [~, ~, O] = integrate(@dPUdT_CombinedTransitions, Pk, Pk.PU0(:), [0 1.5], dt, opts);  Fs(k) = O(end, 1);
end
b = base;  b.CaModel = 'myosim';  b.UseCompiledMex = true;  mexOK = false;
try, getParams(b, [], true, false); catch e, mexOK = strcmp(e.identifier, 'getParams:CaMex'); end
okT7 = dTM < 1e-3 && abs(Fs(3)/Fs(2) - 1) < 5e-3 && mexOK;
fprintf('T7 myosim: full overlap vs twoState max rel |dF| %.1e | pCa 4.51: none %.2f, twoState %.2f, myosim %.2f kPa | MEX guard %d -> %s\n', ...
    dTM, Fs, mexOK, string(okT7));
pass = pass && okT7;

fprintf('\nTestCaActivation: %s\n', string(pass));

%% ---- local functions -------------------------------------------------------------
function p = setParams(base, nv)
for i = 1:2:numel(nv), base.(nv{i}) = nv{i+1}; end
p = getParams(base, [], true, false);
p.Vums = 0;  p.v = 0;
end

function [t, Y, O] = integrate(fcn, P, Y0, knots, dt, opts)
% Integrate FCN across KNOTS, restarting at each (as RunTwitchProbe does).
t = [];  Y = [];
for k = 1:numel(knots) - 1
    ts = unique([knots(k):dt:knots(k+1), knots(k+1)]);
    [tk, Yk] = ode15s(fcn, ts, Y0, opts, P);
    if k > 1, tk = tk(2:end);  Yk = Yk(2:end, :); end
    t = [t; tk];  Y = [Y; Yk];  Y0 = Yk(end, :)'; %#ok<AGROW>
end
[~, o] = fcn(t(1), Y(1,:)', P);
O = zeros(numel(t), numel(o));
for i = 1:numel(t)
    [~, O(i,:)] = fcn(t(i), Y(i,:)', P);
end
end

function m = twitchMetrics(t, F, t0)
i0 = find(t >= t0 - 1e-12, 1);
Fd = F(i0);
[Fp, ip] = max(F(i0:end));  ip = ip + i0 - 1;
half = Fd + 0.5*(Fp - Fd);
j = find(F(ip:end) <= half, 1) + ip - 1;
if isempty(j), rt50 = NaN; else, rt50 = interp1(F([j-1 j]), t([j-1 j]), half) - t(ip); end
m = struct('Fd', Fd, 'Fdev', Fp - Fd, 'TTP', t(ip) - t0, 'RT50', rt50);
end
