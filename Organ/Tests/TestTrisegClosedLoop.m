%% TestTrisegClosedLoop  Does the coupled system start and run a few beats?
%
% A smoke test of setupTrisegXB + runTrisegXB on the P21 host fit Nx51: consistent
% initial geometry, NBEATS beats without a solver failure, mass conserved in the
% circulation, and the per-beat summary printed.
%   matlab -batch "NBEATS = 2; run('Organ/Tests/TestTrisegClosedLoop.m')"

here = fileparts(mfilename('fullpath'));  root = fileparts(fileparts(here));
cd(root);  addpath(genpath(fullfile(root, 'Model')));  addpath(genpath(fullfile(root, 'Auxiliary')));
addpath(fullfile(root, 'params'));  addpath(fullfile(root, 'Organ'));
if ~exist('HOSTROOT', 'var'), HOSTROOT = fullfile(fileparts(root), 'Multiscale-Hx-P21-Model'); end
if ~exist('NBEATS', 'var'),   NBEATS = 2; end
if ~exist('OPT', 'var'),      OPT = struct(); end

H = hostP21Load('Nx51', HOSTROOT);
tic;  C = setupTrisegXB(H, OPT);  tSetup = toc;
G0 = hostP21Geometry(C.y0(1:4), H);
fprintf('setup %.0f s: %d states (walls %d each); geometry residual %.1e; g0 host %s -> %s; Ls0 host %s -> %s\n', tSetup, C.n, C.W(1).nW, C.resG0, ...
    mat2str(H.g0', 4), mat2str(C.y0(1:4)', 4), mat2str(H.Ls0, 4), mat2str(G0.Ls, 4));
[~, o0] = trisegXB_rhs(0, C.y0, C);
fprintf('   at t = 0: P_LV %.2f, P_RV %.2f kPa (host %.2f, %.2f); wall stress %s kPa (host %s)\n', o0(1), o0(4), H.ref.P(1, 1), H.ref.P(1, 4), mat2str(o0(19:21), 3), mat2str(H.ref.sigma(1, :), 3));

R = runTrisegXB(C, NBEATS, struct());
Vtot0 = sum(C.y0(5:10));  Vtot1 = sum(double(R.y(end, 5:10)));
fprintf('%d beats: %.1f s per beat; steps %d (failed %d), RHS calls %d, Jacobians %d; slides LV/SEP/RV %s; blood volume %.6f -> %.6f cm^3\n', ...
    NBEATS, R.wallPerBeat, R.nSteps, R.nFailed, R.nFevals, R.nPD, mat2str(R.nShift), Vtot0, Vtot1);
disp(R.beat(:, {'beat', 'EDV_LV', 'ESV_LV', 'SV_LV', 'EDP_LV', 'Pmax_LV', 'Pmax_RV', 'SBP', 'DBP', 'LsMax_LV', 'LsMin_LV', 'sigmaMax_LV'}));
if ~exist(fullfile(here, 'results'), 'dir'), mkdir(fullfile(here, 'results')); end
save(fullfile(here, 'results', 'closedLoopSmoke.mat'), 'R', 'H');
