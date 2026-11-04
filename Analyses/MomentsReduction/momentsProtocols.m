function P = momentsProtocols()
% MOMENTSPROTOCOLS  Velocity tables for the PDE-vs-moments comparison.
%
%   P(i).name, P(i).T (segment boundaries, s), P(i).V (ML/s per segment).
%   Every protocol starts with a 1.0 s isometric activation from rest at
%   SL0 = 2.2 um, so the activation transient itself is part of the test.
%   Velocities are in ML/s (1 ML = 2.2 um here: SL0 = 2.2, ML = 2... see
%   evaluateModel: Vums = v*ML).

tAct = 1.0;
P = struct('name', {}, 'T', {}, 'V', {});

% isometric force development (the max-Ca analogue of a twitch rise)
P(end+1) = mk('iso', [0 tAct+0.5], 0);

% ktr: fast release 0.15 ML in 5 ms, 20 ms hold, fast restretch, redevelopment
P(end+1) = mk('ktr', tAct + cumsum([0 0.005 0.020 0.005 0.5]), [-30 0 30 0], tAct);

% slack-restretch: fast release to slack, unloaded hold, SLOW restretch
P(end+1) = mk('slack', tAct + cumsum([0 0.012 0.060 0.120 0.4]), [-10 0 1 0], tAct);

% physiological beat: half-sine ejection 0.15 ML in 0.2 s (peak 1.18 ML/s),
% half-sine refill over 0.4 s, twice
nE = 20; nF = 20; Te = 0.2; Tf = 0.4; dL = 0.15;
te = linspace(0, Te, nE+1); tf = linspace(0, Tf, nF+1);
ve = -dL*pi/(2*Te) * sin(pi*(te(1:end-1)+te(2:end))/2/Te);
vf =  dL*pi/(2*Tf) * sin(pi*(tf(1:end-1)+tf(2:end))/2/Tf);
ve = ve * dL / sum(-ve.*diff(te));  vf = vf * dL / sum(vf.*diff(tf));   % exact length closure
Tb = [te, Te + tf(2:end)];  Vb = [ve, vf];
Tbeat = [Tb, Tb(end) + 0.2, Tb(end) + 0.2 + Tb(2:end)];  Vbeat = [Vb, 0, Vb];
P(end+1) = mk('beat', tAct + [Tbeat, Tbeat(end) + 0.2], [Vbeat, 0], tAct);

% staircase: 4 x (+0.01 ML in 5 ms, 150 ms hold)
Ts = tAct; Vs = [];
for k = 1:4, Ts = [Ts, Ts(end)+0.005, Ts(end)+0.155]; Vs = [Vs, 2, 0]; end %#ok<AGROW>
P(end+1) = mk('stairs', Ts, Vs, tAct);

% force-velocity: isovelocity shortening 0.08 ML at v, then hold
for v = [0.1 0.5 1 2 4]
    P(end+1) = mk(sprintf('fv%g', v), tAct + [0 0.08/v 0.08/v+0.2], [-v 0], tAct); %#ok<AGROW>
end
end

function p = mk(name, T, V, tAct)
if nargin > 3
    T = [0, T]; V = [0, V];   % prepend the isometric activation
end
p.name = name; p.T = T(:)'; p.V = V(:)';
assert(numel(p.V) == numel(p.T) - 1, 'protocol %s: %d segments, %d velocities', name, numel(p.T)-1, numel(p.V));
end
