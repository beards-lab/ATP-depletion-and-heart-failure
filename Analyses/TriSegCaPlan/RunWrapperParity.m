%% RunWrapperParity  Does the host-wall wrapper reproduce evaluateModel? (#82, Q1)
%
% The coupling plan rests on wrapping the UNCHANGED dPUdT_CombinedTransitions and
% handing it the sarcomere length from outside. This is the parity test for that
% wrapper: the same length protocol is run (a) through evaluateModel, the
% production harness, and (b) through dPUdT_HostWall + runHostBeats with the
% Ca gate held fully open (N = 0, k_off = 0, so gate = 1 exactly). The two force
% traces must agree to integration tolerance.
%
% Protocol (SL 2.2 um, robust1_seed): hold 0.5 s, shorten at 1 ML/s for 0.1 s,
% hold 0.2 s, re-stretch at 2 ML/s for 0.05 s, hold 0.2 s.
%   matlab -batch "run('Analyses/TriSegCaPlan/RunWrapperParity.m')"

here = fileparts(mfilename('fullpath'));  root = fileparts(fileparts(here));
cd(root);  addpath(genpath(fullfile(root, 'Model')));  addpath(genpath(fullfile(root, 'Auxiliary')));
addpath(fullfile(root, 'params'));  addpath(here);

Tseg = [0 0.5 0.6 0.8 0.85 1.05];            % s, segment boundaries
vML  = [0 -1 0 2 0];                         % ML/s per segment (negative = shortening)
SL0  = 2.2;

% --- host-style signals for the same protocol
P0  = loadParams('robust1_seed');  ML = getParams(P0, [], true, false).ML;
vum = vML*ML;                                % um/s
Lk  = SL0 + [0, cumsum(vum .* diff(Tseg))];  % SL at the boundaries
H.T = Tseg(end);  H.HR = 60/H.T;
H.LV.Ls = @(t) interp1(Tseg, Lk, min(max(t, 0), Tseg(end)));
H.LV.v  = @(t) vum(min(max(discretize(min(max(t, 0), Tseg(end) - 1e-12), Tseg), 1), numel(vum)));
H.Ca    = @(t) 1000 + 0*t;                   % irrelevant: the gate is pinned open below

for ss = [60 250]
    P = setupWall(H, 'LV', struct('ss', ss, 'velMode', 'imposed'));
    P.tf.k_off = 0;  P.PU0(end-1) = 0;       % N = 0 for ever  ->  gate = 1
    R = runHostBeats(P, 1, struct('tPre', 1e-3, 'dtOut', 0.5e-3));

    Pc = rmfield(P, {'host', 'tf', 'hostT'});
    Pc.PU0 = P.PU0(1:end-2);  Pc.Velocity = vML;  Pc.g = [];  Pc.mods = {};
    t0 = tic;  [~, o] = evaluateModel(@dPUdT_CombinedTransitions, Tseg, Pc);  tE = toc(t0);
    [tu, iu] = unique(o.t(:));  Fe = interp1(tu, o.Force(iu), R.t);
    k = R.t >= 0.45 & isfinite(Fe);
    d = R.F(k) - Fe(k);
    fprintf('PARITY ss=%3d: max |dF| %.4f kPa, rms %.4f kPa over t >= 0.45 s (force range %.1f..%.1f kPa; %.3f %% of the peak). Wrapper %.1f s (%d slides), evaluateModel %.1f s\n', ...
        ss, max(abs(d)), rms(d), min(Fe(k)), max(Fe(k)), 100*max(abs(d))/max(Fe(k)), R.wall, R.nShift, tE);
    res.(sprintf('ss%d', ss)) = struct('t', R.t, 'Fwrap', R.F, 'Feval', Fe, 'maxAbs', max(abs(d)), 'rms', rms(d));
end
if ~exist(fullfile(here, 'results'), 'dir'), mkdir(fullfile(here, 'results')); end
save(fullfile(here, 'results', 'wrapperParity.mat'), 'res');
