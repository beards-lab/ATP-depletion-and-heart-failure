function R = humanProbe(p0, opts)
% HUMANPROBE  Species-agnostic readout of the cross-bridge model's headline
% observables, with windows long enough for slow (human beta-MHC) kinetics.
%
%   R = humanProbe(p0[, opts]) returns, at 8 mM MgATP (or opts.MgATP),
%   saturating Ca, SL0 = 2.2 um:
%     F0          isometric steady force [kPa] (active + passive, as Force)
%     FXB0        isometric active (cross-bridge) force [kPa]
%     attached    p1+p2 at isometric steady state;  p1, p2 sub-pools
%     XTOR        gross ATP turnover per head [1/s] (RTD);  XTOR_net = RTD-RDT
%     R2T, R2D    post-stroke exit fluxes to PT / to PD [1/s per head]
%     PT, PD, SRX isometric pool fractions
%     k2eff       effective post-stroke detachment rate [1/s]
%     dwell_ms    attached/XTOR [ms]
%     iso_settled relative force drift over the last 10 % of the isometric run
%     ktr         Brenner-style redevelopment rate (instant detachment, fit on
%                 the 8-92 % band) [1/s];  ktr_t50 = time to 50 % [ms]
%     FV_v, FV_f  ACTIVE (XB) force after 0.1 ML of shortening at each velocity [ML/s, kPa]
%                 (same read-out rule as Model/experiments/runFVExperiment)
%     FV_fnorm    FV_f / FV_f(v=0)
%     Vmax_hill   zero-force intercept of a Hill fit (F+a)(v+b)=(F0+a)b [ML/s]
%     a_over_F0   Hill curvature;  Vmax_lin = linear zero crossing if reached
%     Pmax, v_Pmax  peak normalised power (F/F0 * v) and its velocity [ML/s]
%
%   opts.vels   shortening velocities, ML/s, positive (default [0.1 0.25 0.5 0.75 1 1.5 2 3])
%   opts.Tiso   isometric settling time [s] (default 3)
%   opts.Tktr   redevelopment window [s] (default 2)
%   opts.MgATP  (default 8)
%   opts.fcn    ODE function handle (default @dPUdT_CombinedTransitions)
%
% See also mechProbe, ktrProbe (Auxiliary/diagnostics), which use 0.6 s / 0.15 s
% windows sized for the rodent parametrisation.
if nargin < 2, opts = struct(); end
vels  = getf(opts, 'vels', [0.1 0.25 0.5 0.75 1 1.5 2 3]);
Tiso  = getf(opts, 'Tiso', 3);
Tktr  = getf(opts, 'Tktr', 2);
MgATP = getf(opts, 'MgATP', 8);
fcn   = getf(opts, 'fcn', @dPUdT_CombinedTransitions);

% ---------------- isometric steady state ----------------
p = p0; p.Velocity = 0; p.MgATP = MgATP; p.Ca = 1000; p.ValuesInTime = 1;
p.UseSlack = 0; p.UseKtrProtocol = 0; p.SL0 = 2.2;
p.RunSlack = 0; p.RunForceVelocity = 0; p.RunKtr = 0; p.RunStairs = 0;
if isfield(p, 'PU0'), p = rmfield(p, 'PU0'); end
p = getParams(p, [], true, false);
[~, o] = evaluateModel(fcn, [0 Tiso], p);
t = o.t(:); F = o.Force(:);
R.F0 = F(end); R.FXB0 = o.FXB(end);
R.attached = o.p1_0(end) + o.p2_0(end); R.p1 = o.p1_0(end); R.p2 = o.p2_0(end);
R.XTOR = o.RTD(end); R.XTOR_net = o.RTD(end) - o.RDT(end);
R.R2T = o.R2T(end); R.R2D = o.R2D(end);   % post-stroke exits: to PT (needs ATP) / to PD
R.PT = o.PuATP(end); R.PD = o.PuR(end); R.SRX = o.SR(end) + o.SRD(end);
R.k2eff = (o.R2T(end) + o.R2D(end)) / max(o.p2_0(end), eps);
R.dwell_ms = 1000 * R.attached / max(R.XTOR, eps);
i9 = find(t >= 0.9*Tiso, 1);
R.iso_settled = abs(F(end) - F(i9)) / max(abs(F(end)), eps);   % drift over last 10 %
% activation rise rate (for reference; not Brenner)
w = t > 0.002 & F < 0.92*F(end) & F > 0.08*F(end);
if nnz(w) >= 3, pf = polyfit(t(w), log(1 - F(w)/F(end)), 1); R.k_rise = -pf(1); else, R.k_rise = NaN; end

% ---------------- Brenner-style ktr ----------------
ss = p.ss;
SR = o.SR(end); SRD = o.SRD(end); PD = o.PuR(end); att = R.attached;
PU0 = zeros(1, 2*ss + 7);
PU0(2*ss+1) = SR; PU0(2*ss+2) = 0; PU0(2*ss+3) = o.SL(end); PU0(2*ss+4) = o.LSE(end);
PU0(2*ss+5) = PD + att; PU0(2*ss+6) = SRD; PU0(2*ss+7) = 0;
pk = p; pk.PU0 = PU0; pk.UseSpaceExtension = 0;
[~, o2] = evaluateModel(fcn, [0 Tktr], pk);
t2 = o2.t(:); F2 = o2.Force(:); Fss = F2(end);
w = t2 > 0.001 & F2 < 0.92*Fss & F2 > 0.08*Fss;
if nnz(w) >= 3, pf = polyfit(t2(w), log(1 - F2(w)/Fss), 1); R.ktr = -pf(1); else, R.ktr = NaN; end
i50 = find(F2 >= 0.5*Fss, 1); R.ktr_t50 = 1000*t2(max(i50,1));
R.ktr_Fend_over_F0 = Fss / R.F0;

% ---------------- force-velocity (0.1 ML of shortening at each v) ----------------
pv0 = p; pv0.PU0 = o.PU(end, :); pv0.LXBpivot = o.params.LXBpivot;
N = numel(vels); Fv = nan(1, N);
parple = gcp('nocreate');
if isempty(parple)
    for j = 1:N
        pv = pv0; pv.Velocity = -vels(j);
        try, [~, ov] = evaluateModel(fcn, [0 0.1/vels(j)], pv); Fv(j) = ov.FXB(end); catch, end
    end
else
    fut = cell(1, N);
    for j = 1:N
        pv = pv0; pv.Velocity = -vels(j);
        fut{j} = parfeval(parple, @evaluateModel, 2, fcn, [0 0.1/vels(j)], pv);
    end
    for j = 1:N
        try, [~, ov] = fetchOutputs(fut{j}); Fv(j) = ov.FXB(end); catch, end
    end
end
R.FV_v = [0 vels]; R.FV_f = [R.FXB0 Fv]; R.FV_fnorm = R.FV_f / R.FXB0;   % ACTIVE (XB) force
% Hill fit on positive-force points, F0 pinned at the isometric point
g = isfinite(R.FV_f) & R.FV_f > 0.02*R.FXB0;
v = R.FV_v(g); f = R.FV_fnorm(g);
hill = @(q, v) (1 + q(1)) * q(2) ./ (v + q(2)) - q(1);          % q = [a/F0, b]
if nnz(g) >= 3
    q = fminsearch(@(q) sum((hill(abs(q), v) - f).^2), [0.3 0.5], optimset('Display','off','TolX',1e-6));
    q = abs(q); R.a_over_F0 = q(1); R.b = q(2); R.Vmax_hill = q(2) / q(1);
    vv = linspace(0, R.Vmax_hill, 2000); pw = hill(q, vv) .* vv;
    [R.Pmax, k] = max(pw); R.v_Pmax = vv(k);
else
    R.a_over_F0 = NaN; R.b = NaN; R.Vmax_hill = NaN; R.Pmax = NaN; R.v_Pmax = NaN;
end
R.Vmax_lin = NaN;
ff = R.FV_f; vv = R.FV_v;
for j = 1:N
    if isfinite(ff(j)) && isfinite(ff(j+1)) && ff(j) > 0 && ff(j+1) <= 0
        R.Vmax_lin = interp1([ff(j) ff(j+1)], [vv(j) vv(j+1)], 0, 'linear'); break;
    end
end
end

function v = getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
