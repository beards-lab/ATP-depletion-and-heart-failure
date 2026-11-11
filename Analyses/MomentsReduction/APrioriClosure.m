% APrioriClosure.m
% How many moments does the CURRENT model need?
%
% For every sample of every PDE reference run, take the exact distributions
% p1(s), p2(s) and compare every integral the dxdt needs (force, each
% transition flux, SRX-routed flux, A2 hop) against what a closure would
% reconstruct from the distribution's moments:
%   T2  Tewari/Marzban: kernels Taylor-expanded to 2nd order about s = 0,
%       integrated exactly against moments 0..2 (3 moments / state)
%   G3  Gaussian closure from moments 0..2, kernels evaluated exactly by
%       7-point Gauss-Hermite (3 moments / state; arbitrary kernels allowed)
%   QN  N-node Gauss rule matching moments 0..2N-1 (QMOM, 2N moments / state)
% This is the a-priori (frozen-distribution) error: it isolates the closure
% from error accumulation, which RunMomentsModel measures.

here = fileparts(mfilename('fullpath'));
params0 = momentsSetup('fastgate_best');
P = momentsProtocols();
closures = {'T2','G3','Q1','Q2','Q3','Q4'};
kNames = {'F','k1D','k12','k21','k2','SRX','hop'};
[xgh, wgh] = hermiteRule(7);
h = 2e-5;                                     % FD step for the Taylor closure

res = struct();
for ip = 1:numel(P)
    S = load(fullfile(here, 'results', ['ref_' P(ip).name '.mat'])); ref = S.ref;
    n = numel(ref.t); stride = max(1, floor(n/1500)); idx = 1:stride:n;
    E = nan(numel(idx), numel(kNames)); C = nan(numel(idx), numel(kNames), numel(closures));
    shape = nan(numel(idx), 4);
    for jj = 1:numel(idx)
        i = idx(jj);
        seg = find(ref.t(i) >= ref.T(1:end-1), 1, 'last'); vel = ref.V(seg);
        s = params0.s(:) - (-(ref.SL(i) - ref.LSE(i)) + ref.LXBPivot(i))/2;
        p1 = max(0, ref.p1(i,:)') * params0.dS; p2 = max(0, ref.p2(i,:)') * params0.dS;  % masses
        E(jj,:) = fluxes(s, p1, s, p2, params0, vel);
        for ic = 1:numel(closures)
            switch closures{ic}
                case 'T2'
                    C(jj,:,ic) = taylorFluxes(s, p1, p2, params0, vel, h);
                case 'G3'
                    [x1, w1] = gaussClosure(s, p1, xgh, wgh); [x2, w2] = gaussClosure(s, p2, xgh, wgh);
                    C(jj,:,ic) = fluxes(x1, w1, x2, w2, params0, vel);
                otherwise
                    N = str2double(closures{ic}(2));
                    [x1, w1] = gaussRuleFromMeasure(s, p1, N); [x2, w2] = gaussRuleFromMeasure(s, p2, N);
                    C(jj,:,ic) = fluxes(x1, w1, x2, w2, params0, vel);
            end
        end
        m2 = sum(p2); mu = sum(s.*p2)/max(m2,eps); sd = sqrt(max(0, sum((s-mu).^2.*p2)/max(m2,eps)));
        sk = sum((s-mu).^3.*p2)/max(m2,eps)/max(sd,1e-9)^3;
        m1 = sum(p1); mu1 = sum(s.*p1)/max(m1,eps); sd1 = sqrt(max(0, sum((s-mu1).^2.*p1)/max(m1,eps)));
        shape(jj,:) = [sd1, sd, sk, mu];
    end
    scale = max(abs(E), [], 1); scale(scale == 0) = 1;
    scale(6:7) = scale(5);   % SRX-routed and hop fluxes are SHARES of p2 turnover: judge them on the k2 scale
    errMax = squeeze(max(abs(C - E), [], 1)) ./ scale';     % kNames x closures
    res.(matlab.lang.makeValidName(P(ip).name)) = struct('t', ref.t(idx), 'E', E, 'C', C, ...
        'errMax', errMax, 'shape', shape, 'scale', scale);
    fprintf('\n== %s ==  p2 sd max %.2f nm, p1 sd max %.2f nm, |skew p2| max %.1f\n', P(ip).name, ...
        1e3*max(shape(:,2)), 1e3*max(shape(:,1)), max(abs(shape(:,3))));
    fprintf('%-6s', 'kern'); fprintf('%8s', closures{:}); fprintf('   (max |err| / max |flux|, %%)\n');
    for k = 1:numel(kNames)
        fprintf('%-6s', kNames{k}); fprintf('%8.2f', 100*errMax(k,:)); fprintf('\n');
    end
end
save(fullfile(here, 'results', 'apriori.mat'), 'res', 'closures', 'kNames');

%% ---------------------------------------------------------------------------
function F = fluxes(x1, w1, x2, w2, params, vel)
% integrals of every kernel against node masses w at strains x
K1 = xbKernels(x1, params, vel); K2 = xbKernels(x2, params, vel);
F = [sum(w1.*K1.F1) + sum(w2.*K2.F2), sum(w1.*K1.k1D), sum(w1.*K1.k12), ...
     sum(w2.*K2.k21), sum(w2.*K2.k2), sum(w2.*K2.fSRX.*K2.k2), sum(w2.*K2.kHop)];
end

function F = taylorFluxes(s, p1, p2, params, vel, h)
% Tewari-style: every kernel replaced by its 2nd-order Taylor polynomial at 0
xs = [-h; 0; h];
K = xbKernels(xs, params, vel);
fn = fieldnames(K);
for f = 1:numel(fn)
    v = K.(fn{f}); c.(fn{f}) = [v(2), (v(3)-v(1))/(2*h), (v(3)-2*v(2)+v(1))/h^2/2];
end
M1 = [sum(p1), sum(s.*p1), sum(s.^2.*p1)]; M2 = [sum(p2), sum(s.*p2), sum(s.^2.*p2)];
I = @(cc, M) cc*M';
% product kernels (fSRX*k2) expanded as a product of the two Taylor series, truncated at 2nd order
cSR = [c.fSRX(1)*c.k2(1), c.fSRX(1)*c.k2(2)+c.fSRX(2)*c.k2(1), c.fSRX(1)*c.k2(3)+c.fSRX(2)*c.k2(2)+c.fSRX(3)*c.k2(1)];
F = [I(c.F1,M1) + I(c.F2,M2), I(c.k1D,M1), I(c.k12,M1), I(c.k21,M2), I(c.k2,M2), I(cSR,M2), I(c.kHop,M2)];
end

function [x, w] = gaussClosure(s, p, xgh, wgh)
m0 = sum(p);
if m0 <= 0, x = 0; w = 0; return; end
mu = sum(s.*p)/m0; sd = sqrt(max(0, sum((s-mu).^2.*p)/m0));
x = mu + sd*xgh; w = m0*wgh;
end

function [x, w] = hermiteRule(n)
% probabilists' Gauss-Hermite (weight N(0,1)), weights sum to 1
J = diag(sqrt(1:n-1), 1) + diag(sqrt(1:n-1), -1);
[V, D] = eig(J); x = diag(D); w = V(1,:)'.^2;
end
