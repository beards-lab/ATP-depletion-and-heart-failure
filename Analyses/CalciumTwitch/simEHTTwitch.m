function [Fdev, S] = simEHTTwitch(base, th, tb, opts)
%SIMEHTTWITCH  Developed force of a 1 Hz paced isometric twitch on the grid TB.
%
%   [FDEV, S] = SIMEHTTWITCH(BASE, TH, TB) runs params struct BASE (from
%   loadParams; SL0 already set) with CaModel 'twoState' + 'rice2008' transient,
%   PRE s of diastole and NB beats, and returns the last beat's developed force
%   (kPa, force minus force at the stimulus) at times TB (s after the stimulus).
%   TH is a struct with any of the fields
%       tau1, tau2, t0, Cadia, Camax     Ca transient (s, s, s, uM, uM)
%       tf_k_on, tf_k_off, tf_K_coop     thin filament
%       pre, nb                          protocol (default 1 s, 2 beats)
%       shift                            offset of the model beat from the data's time zero (s)
%       CaModel                          'twoState' (default) or 'myosim'
%   S returns the diagnostics of the last beat on TB: F, Ca, gate, bound, SRX, N, PD, PT, Fxb, Fpas,
%   net ATP hydrolysis hyd = hydDRX (PT->PD minus PD->PT) + hydSRX (SR->SRD minus SRD->SR), heads/s,
%   plus Fdia and the drift of developed force between the last two beats.

if nargin < 4
    opts = odeset('AbsTol', 1e-4, 'RelTol', 1e-3, 'MaxStep', 2e-3);
end
pre = getf(th, 'pre', 1);  nb = getf(th, 'nb', 2);  T = 1;
% shift: the model beat (Ca onset) starts SHIFT s after the recorded stimulus; negative = before
% it (the stimulus record need not be synchronised with force). Sampled periodically.
shift = getf(th, 'shift', 0);
if shift ~= 0, tb = mod(tb - shift, T); end
ca  = struct('tau1', getf(th, 'tau1', 0.020), 'tau2', getf(th, 'tau2', 0.110), ...
             't0', getf(th, 't0', 0.005), 'Cadia', getf(th, 'Cadia', 0.09), 'Camax', getf(th, 'Camax', 1.45));
if abs(ca.tau1 - ca.tau2) < 1e-3                    % the bi-exponential's peak normalisation is singular at tau1 = tau2
    ca.tau2 = ca.tau1 + 1e-3;
end
P = base;
P.CaModel = getf(th, 'CaModel', 'twoState');  P.CaInput = 'analytic';  P.CaTransientKind = 'rice2008';
P.CaPeriod = T;  P.CaFirstBeat = pre;  P.CaTransientOpt = ca;
for f = ["tf_k_on" "tf_k_off" "tf_K_coop"]
    if isfield(th, f), P.(f) = th.(f); end
end
P = getParams(P, [], true, false);
P.Vums = 0;  P.v = 0;
Ns = P.NumberOfStates;  ss = P.ss;

% knots: beat onsets and the transient's onset kink; the last beat is sampled on TB
knots = unique([0, pre + (0:nb)*T, pre + (0:nb-1)*T + ca.t0]);
y0 = P.PU0(:);  tl = pre + (nb-1)*T;
G  = unique([tl - T + tb(:); tl + tb(:)])';          % sample the last two beats on TB
tt = [];  Y = [];
for k = 1:numel(knots) - 1
    a = knots(k);  b = knots(k+1);
    ts = unique([a, G(G > a & G < b), b]);
    if numel(ts) == 2, ts = [a, (a+b)/2, b]; end
    [tk, yk] = ode15s(@dPUdT_CombinedTransitions, ts, y0, opts, P);
    y0 = yk(end, :)';
    tt = [tt; tk(2:end)]; Y = [Y; yk(2:end, :)]; %#ok<AGROW>
end
% force on the last two beats' grids
Fl = force(tt, Y, P, tl + tb);  Fp = force(tt, Y, P, tl - T + tb);
F0 = force(tt, Y, P, tl);
Fdev = Fl - F0;
if nargout > 1
    S.F = Fl;  S.Fdia = F0;  S.drift = max(Fl - F0) - max(Fp - force(tt, Y, P, tl - T));
    [S.Ca, S.gate, S.bound, S.SRX, S.N, S.PD, S.PT, S.Fxb, S.Fpas, S.hydDRX, S.hydSRX] = deal(zeros(size(tb)));
    for i = 1:numel(tb)
        j = find(abs(tt - (tl + tb(i))) < 1e-9, 1);
        [~, o, r] = dPUdT_CombinedTransitions(tt(j), Y(j, :)', P);
        S.Ca(i) = o(end-1);  S.gate(i) = o(end);  S.bound(i) = o(5) + o(6);
        S.SRX(i) = Y(j, Ns*ss+1) + Y(j, Ns*ss+6);
        S.N(i)   = Y(j, P.ca.idx(1));                % non-permissive fraction of sites
        if strcmp(P.CaModel, 'myosim'), S.N(i) = 1 - S.N(i); end   % its state is N_on
        S.PD(i)  = Y(j, Ns*ss+5);  S.PT(i) = o(9);   % primed M.ADP.Pi, M.ATP (2-state outputs)
        S.Fxb(i) = o(2);  S.Fpas(i) = o(3);          % cross-bridge and passive force (kPa)
        % net ATP hydrolysis (heads/s): forward minus reverse, so reversal returns no ATP
        S.hydDRX(i) = r(1) - r(2);                   % PT -> PD  minus PD -> PT
        S.hydSRX(i) = r(13) - r(14);                 % SR -> SRD minus SRD -> SR
    end
    S.hyd = S.hydDRX + S.hydSRX;
    S.P = P;
end
end

function F = force(tt, Y, P, tq)
F = zeros(numel(tq), 1);
for i = 1:numel(tq)
    j = find(abs(tt - tq(i)) < 1e-9, 1);
    if isempty(j), error('simEHTTwitch:grid', 'time %.6f not on the output grid', tq(i)); end
    [~, o] = dPUdT_CombinedTransitions(tt(j), Y(j, :)', P);
    F(i) = o(1);
end
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
