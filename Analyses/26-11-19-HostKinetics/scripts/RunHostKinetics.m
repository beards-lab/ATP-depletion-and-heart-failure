%% RunHostKinetics  How much rate scaling does the XB model need to follow the P21 host beat? (#98)
%
% Open loop, as in #82: one wall driven by the host's own Ls(t) and Ca(t) (fit Nx51)
% through Analyses/TriSegCaPlan/dPUdT_HostWall. No model file is changed.
%
% Arm families (each arm saved as it finishes; an existing arm is skipped):
%   U<s>  UNIFORM scaling of every rate in the model by s (cross-bridge cycle, SRX,
%         availability gate, viscous passive, thin filament). Done exactly, without
%         listing the rates: the host beat is slowed s-fold (Ls(t/s), Ca(t/s), v/s)
%         and the result is read back in host time. A model whose rates are all s
%         times faster on the real beat is the same system.
%   C<s>  the seven cycle rates only (kah, kamh, ka, kd, k1, k_1, k2), as E7 of #82.
%   N<s>  everything EXCEPT the seven cycle rates by s (U<s> with the cycle rates
%         divided by s again).
%   W<s>_<wall>  U<s> on the septal and RV walls.
%   K<s>_<k>     U<s> with k_on x k instead of x3 (activation amplitude).
% Set ONLY = {'U', ...} to run a subset of families, ARMS to name single arms.
%
%   run('Analyses/26-11-19-HostKinetics/scripts/RunHostKinetics.m')   (~11 min, one core)

here = fileparts(mfilename('fullpath'));  ana = fileparts(here);  root = fileparts(fileparts(ana));
addpath(genpath(fullfile(root, 'Model')));  addpath(genpath(fullfile(root, 'Auxiliary')));
addpath(fullfile(root, 'params'));  addpath(fullfile(root, 'Analyses', 'TriSegCaPlan'));  addpath(here);
if ~exist(fullfile(ana, 'results'), 'dir'), mkdir(fullfile(ana, 'results')); end
if ~exist('HOSTROOT', 'var'), HOSTROOT = fullfile(fileparts(root), 'Multiscale-Hx-P21-Model'); end
if ~exist('FIT', 'var'), FIT = 'Nx51'; end
outF = fullfile(ana, 'results', ['hostKinetics_' FIT '.mat']);
res = struct();  if isfile(outF), L = load(outF);  res = L.res;  end
H = hostSignals(FIT, HOSTROOT);
nB = 4;

list = {};   % name, wall, s (dilation), xrate (cycle rates), konMult
for s = [1 1.5 2 3 4 5 6.5 8 10], list(end+1, :) = {sprintf('U%s', tag(s)), 'LV', s, 1, 3}; end %#ok<SAGROW>
for s = [2 3 5 8 10],             list(end+1, :) = {sprintf('C%s', tag(s)), 'LV', 1, s, 3}; end %#ok<SAGROW>
for s = [5 10],                   list(end+1, :) = {sprintf('N%s', tag(s)), 'LV', s, 1/s, 3}; end %#ok<SAGROW>
for s = [1 5], for w = {'SEP', 'RV'}, list(end+1, :) = {sprintf('W%s_%s', tag(s), w{1}), w{1}, s, 1, 3}; end, end %#ok<SAGROW>
for s = [5 10], for k = [1 10],   list(end+1, :) = {sprintf('K%s_%d', tag(s), k), 'LV', s, 1, k}; end, end %#ok<SAGROW>
if exist('ONLY', 'var'), list = list(ismember(cellfun(@(c) c(1), list(:, 1), 'UniformOutput', false), ONLY), :); end
if exist('ARMS', 'var'), list = list(ismember(list(:, 1), ARMS), :); end

for i = 1:size(list, 1)
    [nm, wall, s, xr, kon] = list{i, :};
    if isfield(res, nm) && ~isfield(res.(nm), 'err'), continue; end
    Hs = dilate(H, s);
    try
        P = setupWall(Hs, wall, struct('ss', 60, 'velMode', 'imposed', 'xrate', xr, 'konMult', kon));
        R = runHostBeats(P, nB, struct('tPre', 0.5, 'dtOut', 0.25e-3*s));
        k = R.t >= (nB - 1)*Hs.T - 1e-9;                 % last beat
        A = struct('s', s, 'xrate', xr, 'konMult', kon, 'wall', wall, ...
            't', (R.t(k) - (nB - 1)*Hs.T)/s, ...         % host time (s)
            'F', R.F(k), 'Fa', R.Factive(k), 'Fp', R.Fpassive(k), 'SL', R.SL(k), 'LSE', R.LSE(k), ...
            'gate', R.gate(k), 'bound', R.bound(k), 'wallPerBeat', R.wallPerBeat, 'nShift', R.nShift, ...
            'massLost', sum(R.massLost), 'peakPerBeat', arrayfun(@(b) max(R.Factive(R.t >= (b-1)*Hs.T & R.t < b*Hs.T)), 1:nB));
        A.m = metrics(A, H, wall);
        res.(nm) = A;
        fprintf('%-9s %-3s s=%4.1f cyc x%5.2f kon x%2d | %5.1f s/beat | active: peak %5.1f min %6.2f kPa, ejection mean %5.2f (host %5.2f, %3.0f %%) | total: slack %5.2f ms, ejection mean %5.2f (host %5.2f)\n', ...
            nm, wall, s, xr, kon, A.wallPerBeat, A.m.FaPeak, A.m.FaMin, A.m.FaEj, A.m.hostXBEj, 100*A.m.FaEj/A.m.hostXBEj, 1e3*A.m.tSlack, A.m.FEj, A.m.hostEj);
    catch ME
        res.(nm) = struct('err', ME.message);  fprintf('%-9s FAILED: %s\n', nm, ME.message);
    end
    save(outF, 'res');
end

function s = tag(x), s = strrep(sprintf('%g', x), '.', 'p'); end

function Hs = dilate(H, s)
% the host beat slowed s-fold: the same lengths and Ca at t/s, velocities divided by s
Hs = H;  Hs.T = s*H.T;
for w = {'LV', 'SEP', 'RV'}
    Ls = H.(w{1}).Ls;  v = H.(w{1}).v;
    Hs.(w{1}).Ls = @(t) Ls(t/s);  Hs.(w{1}).v = @(t) v(t/s)/s;
end
Ca = H.Ca;  Hs.Ca = @(t) Ca(t/s);
end

function m = metrics(A, H, wall)
% all times in host time; ejection window 20-70 ms after the start of the beat
t = A.t;  ej = t >= 0.020 & t <= 0.070;
dt = gradient(t);
m.FaPeak = max(A.Fa);  m.tPeak = t(find(A.Fa == m.FaPeak, 1));
m.FaMin = min(A.Fa(t >= m.tPeak));                    % lowest active stress after the peak
m.FMinEj = min(A.F(ej));
m.tSlack = sum(dt(A.F < 0.05));                       % s with total stress below 0.05 kPa
m.tNeg = sum(dt(A.Fa < 0));                           % s with negative cross-bridge stress
m.FaEj = sum(A.Fa(ej).*dt(ej))/sum(dt(ej));  m.FEj = sum(A.F(ej).*dt(ej))/sum(dt(ej));
hx = interp1(H.raw.t, H.raw.(wall).sigma_XB, t(ej), 'linear', 'extrap');
hs = interp1(H.raw.t, H.raw.(wall).sigma, t(ej), 'linear', 'extrap');
m.hostXBEj = sum(hx.*dt(ej))/sum(dt(ej));  m.hostEj = sum(hs.*dt(ej))/sum(dt(ej));
m.hostXBPeak = max(H.raw.(wall).sigma_XB);
% lowest active stress during the fast-shortening spike, as a fraction of the pre-ejection peak
m.dip = m.FaMin/m.FaPeak;
vCE = gradient(A.SL - A.LSE, t);  m.vCEmax = -min(vCE);   % contractile-element shortening, um/s
end
