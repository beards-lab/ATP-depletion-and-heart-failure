%% RunHostProbe  Open-loop probes of the full XB model on a TriSeg host beat (#82).
%
% Drives one wall's sarcomere with the P21 host's own Ls(t) and Ca(t) (fit Nx51,
% 325 bpm) through dPUdT_HostWall (the model file is unchanged). Experiments:
%   E1  grid / harness cost: ss = 60 (incumbent sliding window), 120, 200, 4 beats,
%       imposed velocity. Wall clock, window slides, mass lost, force difference.
%   E2  sparsity: Jacobian pattern over a beat, column groups (RHS calls per
%       Jacobian), and the run with JPattern supplied.
%   E3  velocity signal: imposed vs velHS vs zero for the velocity-driven mechanisms.
%   E4  all three walls (LV, SEP, RV) at ss = 60.
%   E5  activation amplitude: k_on x1 (literature) vs x3.
% Results -> results/hostProbe.mat; figures -> figures/.
% Run from the repo root:  run('Analyses/TriSegCaPlan/RunHostProbe.m')

root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root);  addpath(genpath('Model'));  addpath(genpath('Auxiliary'));  addpath('params');
here = fullfile(root, 'Analyses', 'TriSegCaPlan');  addpath(here);
if ~exist(fullfile(here, 'results'), 'dir'), mkdir(fullfile(here, 'results')); end
H = hostSignals('Nx51', 'C:\home\git\Multiscale-Hx-P21-Model');
nB = 4;
res = struct();

%% E1  grid size, sliding window
for ss = [60 120 200]
    P = setupWall(H, 'LV', struct('ss', ss, 'velMode', 'imposed'));
    R = runHostBeats(P, nB, struct('tPre', 0.5));
    res.(sprintf('E1_ss%d', ss)) = R;
    fprintf('E1 ss=%3d n=%3d: %.1f s/beat, slides %d, steps %d (failed %d), mass lost %.2g, F %.1f..%.1f kPa\n', ...
        ss, numel(P.PU0), R.wallPerBeat, R.nShift, R.nSteps, R.nFailed, sum(R.massLost), min(R.F), max(R.F));
end

%% E2  Jacobian sparsity and JPattern
for ss = [60 200]
    P  = setupWall(H, 'LV', struct('ss', ss, 'velMode', 'imposed'));
    R0 = res.(sprintf('E1_ss%d', ss));
    idx = arrayfun(@(f) find(R0.t >= (nB-1)*H.T + f*H.T, 1), [0.02 0.1 0.2 0.35 0.5 0.7 0.9]);
    Pk = P;  Pk.LXBpivot = R0.pivot(idx(1));
    [S, info] = jacPattern(Pk, double(R0.Y(idx, :))', R0.t(idx)');
    y = double(R0.Y(idx(4), :))';
    tic; for k = 1:50, dPUdT_HostWall(R0.t(idx(4)), y, Pk); end; info.rhs_us = toc/50*1e6;
    R = runHostBeats(P, nB, struct('tPre', 0.5, 'JPattern', S));
    res.(sprintf('E2_ss%d', ss)) = struct('info', info, 'R', rmfield(R, 'Y'));
    fprintf('E2 ss=%3d: RHS %.0f us, density %.1f %%, dense rows %s, groups %d of %d (no dense rows: %d); with JPattern %.1f s/beat (steps %d, failed %d)\n', ...
        ss, info.rhs_us, 100*info.density, mat2str(info.denseRows), info.nGroups, numel(P.PU0), info.nGroupsNoDenseRows, R.wallPerBeat, R.nSteps, R.nFailed);
end

%% E3  velocity signal for the velocity-driven mechanisms
for vm = {'velHS', 'zero'}
    P = setupWall(H, 'LV', struct('ss', 60, 'velMode', vm{1}));
    R = runHostBeats(P, nB, struct('tPre', 0.5));
    res.(['E3_' vm{1}]) = R;
    Ri = res.E1_ss60;  Fi = interp1(Ri.t, Ri.F, R.t);  Ai = interp1(Ri.t, Ri.A_reg, R.t);
    lb = R.t >= (nB-1)*H.T;
    fprintf('E3 %-6s: %.1f s/beat | last beat: max|dF| %.2f kPa (peak F %.1f vs %.1f), max|dA_reg| %.3f (A_reg %.3f..%.3f vs %.3f..%.3f)\n', ...
        vm{1}, R.wallPerBeat, max(abs(R.F(lb) - Fi(lb))), max(R.F(lb)), max(Fi(lb)), max(abs(R.A_reg(lb) - Ai(lb))), ...
        min(R.A_reg(lb)), max(R.A_reg(lb)), min(Ai(lb)), max(Ai(lb)));
end

%% E4  the other walls
for w = {'SEP', 'RV'}
    P = setupWall(H, w{1}, struct('ss', 60, 'velMode', 'imposed'));
    R = runHostBeats(P, nB, struct('tPre', 0.5));
    res.(['E4_' w{1}]) = R;
    fprintf('E4 %-3s: %.1f s/beat, slides %d, F %.1f..%.1f kPa\n', w{1}, R.wallPerBeat, R.nShift, min(R.F), max(R.F));
end

%% E5  literature thin filament (k_on x1)
P = setupWall(H, 'LV', struct('ss', 60, 'velMode', 'imposed', 'konMult', 1));
R = runHostBeats(P, nB, struct('tPre', 0.5));
res.E5_kon1 = R;
fprintf('E5 k_on x1: gate max %.3f, F %.1f..%.1f kPa\n', max(R.gate), min(R.F), max(R.F));

res.H = rmfield(H, {'LV', 'SEP', 'RV', 'Ca'});  res.Hraw = H.raw;  res.T = H.T;
fn = fieldnames(res);
for k = 1:numel(fn), if isstruct(res.(fn{k})) && isfield(res.(fn{k}), 'Y'), res.(fn{k}) = rmfield(res.(fn{k}), 'Y'); end, end
save(fullfile(here, 'results', 'hostProbe.mat'), 'res');
disp('saved');
