function [Fdev, S] = simHumanTwitch(base, th, tb, opts)
%SIMHUMANTWITCH  Developed force of a paced isometric twitch on the grid TB (s after the stimulus).
%
%   [FDEV, S] = SIMHUMANTWITCH(BASE, TH, TB) runs params struct BASE (loaded, scaled; getParams
%   is called here) with a prescribed Ca transient and the twoState thin filament, PRE s of
%   diastole then NB beats of period T, and returns the last beat's developed force (kPa, force
%   minus force at the stimulus) on TB. TH fields (all optional):
%       CaKind, CaOpt        caTransient kind (default 'risedecay') and its option struct
%       tf_k_on, tf_k_off, tf_K_coop     thin filament
%       CaModel              'twoState' (default) | 'myosim'
%       pre, nb, T           protocol (default 2 s, 3 beats, 1 s)
%       shift                model beat onset relative to the data's time zero (s); may be < 0
%   S: diagnostics on TB (F, Ca, gate, bound, SRX, N, PD, PT, Fxb, Fpas, hyd = net ATP hydrolysis
%   per head per s, R2T = post-stroke detachment flux), Fdia, drift between the last two beats, P.
%   Generalises Analyses/CalciumTwitch/simEHTTwitch.m (fixed rice2008 transient).

if nargin < 4
    opts = odeset('AbsTol', 1e-4, 'RelTol', 1e-3, 'MaxStep', 2e-3);
end
pre = getf(th, 'pre', 2);  nb = getf(th, 'nb', 3);  T = getf(th, 'T', 1);
shift = getf(th, 'shift', 0);
if shift ~= 0, tb = mod(tb - shift, T); end
P = base;
P.CaModel = getf(th, 'CaModel', 'twoState');  P.CaInput = 'analytic';
P.CaTransientKind = getf(th, 'CaKind', 'risedecay');
P.CaPeriod = T;  P.CaFirstBeat = pre;  P.CaTransientOpt = getf(th, 'CaOpt', struct());
for f = ["tf_k_on" "tf_k_off" "tf_K_coop"]
    if isfield(th, f), P.(f) = th.(f); end
end
P = getParams(P, [], true, false);
P.Vums = 0;  P.v = 0;
Ns = P.NumberOfStates;  ss = P.ss;
kinks = P.ca.beatKnots(:)';                       % within-beat kinks (onset)
knots = unique([0, pre + (0:nb)*T, reshape(pre + (0:nb-1)'*T + kinks, 1, [])]);
y0 = P.PU0(:);  tl = pre + (nb-1)*T;
G  = unique([tl - T + tb(:); tl + tb(:)])';
tt = [];  Y = [];
for k = 1:numel(knots) - 1
    a = knots(k);  b = knots(k+1);
    ts = unique([a, G(G > a & G < b), b]);
    if numel(ts) == 2, ts = [a, (a+b)/2, b]; end
    [tk, yk] = ode15s(@dPUdT_CombinedTransitions, ts, y0, opts, P);
    if tk(end) < b - 1e-9, error('simHumanTwitch:ode', 'integration stopped at %.4f', tk(end)); end
    y0 = yk(end, :)';
    tt = [tt; tk(2:end)]; Y = [Y; yk(2:end, :)]; %#ok<AGROW>
end
Fl = force(tt, Y, P, tl + tb);  F0 = force(tt, Y, P, tl);
Fdev = Fl - F0;
if nargout > 1
    Fp = force(tt, Y, P, tl - T + tb);
    S.F = Fl;  S.Fdia = F0;  S.drift = max(Fl - F0) - max(Fp - force(tt, Y, P, tl - T));
    z = zeros(size(tb));
    [S.Ca, S.gate, S.bound, S.SRX, S.N, S.PD, S.PT, S.Fxb, S.Fpas, S.hydDRX, S.hydSRX, S.R2T] = deal(z);
    for i = 1:numel(tb)
        j = find(abs(tt - (tl + tb(i))) < 1e-9, 1);
        [~, o, r] = dPUdT_CombinedTransitions(tt(j), Y(j, :)', P);
        S.Ca(i) = o(end-1);  S.gate(i) = o(end);  S.bound(i) = o(5) + o(6);
        S.SRX(i) = Y(j, Ns*ss+1) + Y(j, Ns*ss+6);
        S.N(i)   = Y(j, P.ca.idx(1));
        if strcmp(P.CaModel, 'myosim'), S.N(i) = 1 - S.N(i); end
        S.PD(i)  = Y(j, Ns*ss+5);  S.PT(i) = o(9);
        S.Fxb(i) = o(2);  S.Fpas(i) = o(3);
        S.hydDRX(i) = r(1) - r(2);  S.hydSRX(i) = r(13) - r(14);  S.R2T(i) = r(7);
    end
    S.hyd = S.hydDRX + S.hydSRX;
    S.P = P;  S.tb = tb;
end
end

function F = force(tt, Y, P, tq)
F = zeros(numel(tq), 1);
for i = 1:numel(tq)
    j = find(abs(tt - tq(i)) < 1e-9, 1);
    if isempty(j), error('simHumanTwitch:grid', 'time %.6f not on the output grid', tq(i)); end
    [~, o] = dPUdT_CombinedTransitions(tt(j), Y(j, :)', P);
    F(i) = o(1);
end
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
