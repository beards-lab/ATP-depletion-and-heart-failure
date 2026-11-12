%% RunHostProbe  Open-loop probes of the full XB model on a TriSeg host beat (#82).
%
% Drives one wall's sarcomere with the P21 host's own Ls(t) and Ca(t) (fit Nx51,
% 325 bpm) through dPUdT_HostWall (the model file is unchanged). Parameter set:
% robust1_seed, the incumbent. Experiments (each saved as it finishes; a failed
% arm is recorded with its error and the script goes on):
%   E0  host kinematics per wall: Ls range, peak shortening / lengthening speed.
%   E1  grid / harness cost: ss = 60 (incumbent sliding window), 120, 250, 4 beats,
%       imposed velocity. Wall clock, window slides, mass lost, force difference.
%       Plus fastgate_best at ss = 60 (the set of the 2026-11-09 probes).
%   E2  sparsity: Jacobian pattern over a beat, column groups (RHS calls per
%       Jacobian), and the run with JPattern supplied.
%   E3  velocity signal for the velocity-driven mechanisms: imposed (reference),
%       zero, velHS (algebraic loop), low-pass of velHS and of dLs/dt (1, 5 ms).
%   E4  the other walls (SEP, RV) at ss = 60.
%   E5  activation amplitude: k_on x1 (literature), x3, x10.
%   E6  approach to the periodic state: 12 beats.
%   E7  kinetics at the host's rate: cycle rates x2 and x5.
%   E8  kinematic lag: host Ls(t) delayed by 0.1 / 0.5 / 1 / 2 ms (what an explicit
%       co-simulation macro-step does to the wall).
%   E9  the same beat with a C2 spline through the host's stored samples instead of
%       pchip: how much of the step count is the stored signal's kinks.
% Results -> results/hostProbe.mat, log -> results/hostProbe_log.txt.
%   matlab -batch "run('Analyses/TriSegCaPlan/RunHostProbe.m')"     (~25 min)

here = fileparts(mfilename('fullpath'));  root = fileparts(fileparts(here));
cd(root);  addpath(genpath(fullfile(root, 'Model')));  addpath(genpath(fullfile(root, 'Auxiliary')));
addpath(fullfile(root, 'params'));  addpath(here);
if ~exist(fullfile(here, 'results'), 'dir'), mkdir(fullfile(here, 'results')); end
if ~exist('HOSTROOT', 'var'), HOSTROOT = fullfile(fileparts(root), 'Multiscale-Hx-P21-Model'); end
if ~exist('ONLY', 'var'), ONLY = {'E0', 'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9'}; end
outF = fullfile(here, 'results', 'hostProbe.mat');
res = struct();  if isfile(outF), L = load(outF);  res = L.res;  end
diary(fullfile(here, 'results', 'hostProbe_log.txt'));  diary on;
H = hostSignals('Nx51', HOSTROOT);
nB = 4;  T = H.T;
lastBeat = @(R) R.t >= (round(R.t(end)/T) - 1)*T - 1e-9;
base = struct('ss', 60, 'velMode', 'imposed');
fprintf('\n=== RunHostProbe %s, host Nx51, HR %.0f bpm, T %.1f ms ===\n', datestr(now), H.HR, 1e3*T); %#ok<TNOW1,DATST>

%% E0  host kinematics
if any(strcmp(ONLY, 'E0'))
    tt = linspace(0, T, 2001);
    for w = {'LV', 'SEP', 'RV'}
        Ls = H.(w{1}).Ls(tt);  v = H.(w{1}).v(tt);
        res.E0.(w{1}) = struct('LsMin', min(Ls), 'LsMax', max(Ls), 'LsED', Ls(1), 'vShort', -min(v), 'vLength', max(v));
        fprintf('E0 %-3s: Ls %.3f..%.3f um (t = 0: %.3f), peak shortening %.1f um/s = %.2f ML/s, peak lengthening %.1f um/s = %.2f ML/s; host passive %.2f..%.2f kPa\n', ...
            w{1}, min(Ls), max(Ls), Ls(1), -min(v), -min(v)/max(Ls), max(v), max(v)/max(Ls), min(H.raw.(w{1}).sigma_pas), max(H.raw.(w{1}).sigma_pas));
    end
    fprintf('E0 Ca: %.3f..%.3f uM\n', min(H.raw.Ca), max(H.raw.Ca));
    save(outF, 'res');
end

%% E1  grid size, sliding window
if any(strcmp(ONLY, 'E1'))
    for ss = [60 120 250]
        nm = sprintf('E1_ss%d', ss);
        res.(nm) = arm(H, 'LV', setfield(base, 'ss', ss), nB, struct()); %#ok<SFLD>
        R = res.(nm);
        if ~isfield(R, 'err')
            k = lastBeat(R);
            n = size(R.Y, 2);
            fprintf('E1 ss=%3d n=%3d: %.1f s/beat (pre-equilibration %.0f s), slides %d, steps %d (failed %d), RHS calls %d, Jacobians %d (= %.0f %% of the RHS calls), mass lost %.2g, last beat F %.2f..%.2f kPa (all beats max %.1f)\n', ...
                ss, n, R.wallPerBeat, R.wallPre, R.nShift, R.nSteps, R.nFailed, R.nFevals, R.nPD, 100*R.nPD*n/R.nFevals, sum(R.massLost), min(R.F(k)), max(R.F(k)), max(R.F));
        end
        save(outF, 'res');
    end
    for ss = [60 200]            % the set of the 2026-11-09 probes, which read 104.9 kPa at ss = 200
        nm = sprintf('E1_fastgate%d', ss);
        res.(nm) = arm(H, 'LV', struct('ss', ss, 'velMode', 'imposed', 'paramSet', 'fastgate_best'), nB, struct());  R = res.(nm);
        if ~isfield(R, 'err'), k = lastBeat(R); pk = arrayfun(@(b) max(R.F(R.t >= (b-1)*T & R.t < b*T)), 1:nB); fprintf('E1 fastgate_best ss=%d: %.1f s/beat, slides %d, last beat F %.2f..%.2f kPa, peak per beat %s\n', ss, R.wallPerBeat, R.nShift, min(R.F(k)), max(R.F(k)), mat2str(pk, 4)); end
        save(outF, 'res');
    end
end

%% E2  Jacobian sparsity and JPattern
if any(strcmp(ONLY, 'E2'))
    for ss = [60 250]
        R0 = res.(sprintf('E1_ss%d', ss));
        if isfield(R0, 'err'), continue; end
        try
            P  = setupWall(H, 'LV', setfield(base, 'ss', ss)); %#ok<SFLD>
            idx = arrayfun(@(f) find(R0.t >= (nB-1)*T + f*T, 1), [0.02 0.1 0.2 0.35 0.5 0.7 0.9]);
            Pk = P;  Pk.LXBpivot = R0.pivot(idx(1));
            [S, info] = jacPattern(Pk, double(R0.Y(idx, :))', R0.t(idx)');
            y = double(R0.Y(idx(4), :))';
            for k = 1:20, dPUdT_HostWall(R0.t(idx(4)), y, Pk); end      % warm-up
            tic; for k = 1:500, dPUdT_HostWall(R0.t(idx(4)), y, Pk); end; info.rhs_us = toc/500*1e6;
            R = runHostBeats(P, nB, struct('tPre', 0.5, 'JPattern', S));
            res.(sprintf('E2_ss%d', ss)) = struct('info', info, 'wallPerBeat', R.wallPerBeat, 'nSteps', R.nSteps, 'nFailed', R.nFailed, 'Fmax', max(R.F));
            fprintf('E2 ss=%3d: RHS %.0f us, density %.1f %%, %d dense rows, groups %d of %d (without the dense rows: %d); with JPattern %.1f s/beat (steps %d, failed %d)\n', ...
                ss, info.rhs_us, 100*info.density, numel(info.denseRows), info.nGroups, size(S, 1), info.nGroupsNoDenseRows, R.wallPerBeat, R.nSteps, R.nFailed);
        catch ME
            res.(sprintf('E2_ss%d', ss)) = struct('err', ME.message);  fprintf('E2 ss=%3d FAILED: %s\n', ss, ME.message);
        end
        save(outF, 'res');
    end
end

%% E3  velocity signal for the velocity-driven mechanisms
if any(strcmp(ONLY, 'E3'))
    arms = {'zero', 'zero', 5e-3; 'velHS', 'velHS', 5e-3; 'filtHS_1ms', 'filtHS', 1e-3; 'filtHS_5ms', 'filtHS', 5e-3; ...
            'filtHost_1ms', 'filtHost', 1e-3; 'filtHost_5ms', 'filtHost', 5e-3};
    Ri = res.E1_ss60;  ki = lastBeat(Ri);
    for a = 1:size(arms, 1)
        nm = ['E3_' arms{a, 1}];
        res.(nm) = arm(H, 'LV', struct('ss', 60, 'velMode', arms{a, 2}, 'tauV', arms{a, 3}), nB, struct());
        R = res.(nm);
        if ~isfield(R, 'err')
            k = lastBeat(R);  Fi = interp1(Ri.t, Ri.F, R.t(k));  Ai = interp1(Ri.t, Ri.A_reg, R.t(k));
            res.(nm).dFmax = max(abs(R.F(k) - Fi));  res.(nm).dFrms = rms(R.F(k) - Fi);  res.(nm).dAmax = max(abs(R.A_reg(k) - Ai));
            fprintf('E3 %-13s: %.1f s/beat | last beat vs imposed: max|dF| %.2f kPa, rms %.2f (peak F %.1f vs %.1f), max|dA_reg| %.3f (A_reg %.3f..%.3f vs %.3f..%.3f)\n', ...
                arms{a, 1}, R.wallPerBeat, res.(nm).dFmax, res.(nm).dFrms, max(R.F(k)), max(Ri.F(ki)), res.(nm).dAmax, min(R.A_reg(k)), max(R.A_reg(k)), min(Ri.A_reg(ki)), max(Ri.A_reg(ki)));
        end
        save(outF, 'res');
    end
end

%% E4  the other walls
if any(strcmp(ONLY, 'E4'))
    for w = {'SEP', 'RV'}
        nm = ['E4_' w{1}];
        res.(nm) = arm(H, w{1}, base, nB, struct());  R = res.(nm);
        if ~isfield(R, 'err'), k = lastBeat(R); fprintf('E4 %-3s: %.1f s/beat, slides %d, last beat F %.2f..%.2f kPa, passive up to %.1f kPa\n', w{1}, R.wallPerBeat, R.nShift, min(R.F(k)), max(R.F(k)), max(R.Fpassive(k))); end
        save(outF, 'res');
    end
end

%% E5  activation amplitude
if any(strcmp(ONLY, 'E5'))
    for km = [1 10]
        nm = sprintf('E5_kon%d', km);
        res.(nm) = arm(H, 'LV', setfield(base, 'konMult', km), nB, struct());  R = res.(nm); %#ok<SFLD>
        if ~isfield(R, 'err'), k = lastBeat(R); fprintf('E5 k_on x%-2d: gate %.3f..%.3f, last beat F %.2f..%.2f kPa, active %.2f..%.2f kPa\n', km, min(R.gate(k)), max(R.gate(k)), min(R.F(k)), max(R.F(k)), min(R.Factive(k)), max(R.Factive(k))); end
        save(outF, 'res');
    end
    R = res.E1_ss60;  k = lastBeat(R);
    fprintf('E5 k_on x3 : gate %.3f..%.3f, last beat F %.2f..%.2f kPa, active %.2f..%.2f kPa\n', min(R.gate(k)), max(R.gate(k)), min(R.F(k)), max(R.F(k)), min(R.Factive(k)), max(R.Factive(k)));
end

%% E6  approach to the periodic state
if any(strcmp(ONLY, 'E6'))
    res.E6 = arm(H, 'LV', base, 12, struct());  R = res.E6;
    if ~isfield(R, 'err')
        pk = arrayfun(@(b) max(R.F(R.t >= (b-1)*T & R.t < b*T)), 1:12);
        mn = arrayfun(@(b) min(R.F(R.t >= (b-1)*T & R.t < b*T)), 1:12);
        res.E6.peak = pk;  res.E6.trough = mn;
        fprintf('E6 12 beats: %.1f s/beat; peak F per beat %s; change vs beat 12 (%%): %s\n', R.wallPerBeat, mat2str(pk, 4), mat2str(100*(pk/pk(end) - 1), 2));
    end
    save(outF, 'res');
end

%% E7  kinetics at the host's rate
if any(strcmp(ONLY, 'E7'))
    for xr = [2 5]
        nm = sprintf('E7_x%d', xr);
        res.(nm) = arm(H, 'LV', setfield(base, 'xrate', xr), nB, struct());  R = res.(nm); %#ok<SFLD>
        if ~isfield(R, 'err'), k = lastBeat(R); fprintf('E7 rates x%d: %.1f s/beat, last beat active F %.2f..%.2f kPa, total %.2f..%.2f, gate max %.3f\n', xr, R.wallPerBeat, min(R.Factive(k)), max(R.Factive(k)), min(R.F(k)), max(R.F(k)), max(R.gate(k))); end
        save(outF, 'res');
    end
end

%% E8  kinematic lag of the host signal
if any(strcmp(ONLY, 'E8'))
    Ri = res.E1_ss60;
    for dl = [0.1 0.5 1 2]*1e-3
        nm = sprintf('E8_lag%dus', round(dl*1e6));
        res.(nm) = arm(H, 'LV', setfield(base, 'delay', dl), nB, struct());  R = res.(nm); %#ok<SFLD>
        if ~isfield(R, 'err')
            k = lastBeat(R);  Fi = interp1(Ri.t, Ri.F, R.t(k));
            res.(nm).dFmax = max(abs(R.F(k) - Fi));  res.(nm).dFrms = rms(R.F(k) - Fi);
            fprintf('E8 lag %.1f ms: last beat max|dF| %.2f kPa, rms %.2f kPa (peak F %.1f)\n', 1e3*dl, res.(nm).dFmax, res.(nm).dFrms, max(R.F(k)));
        end
        save(outF, 'res');
    end
end

%% E9  smooth host signal: is the step count the host's stored samples?
if any(strcmp(ONLY, 'E9'))
    Hs = hostSignals('Nx51', HOSTROOT, 'spline');
    res.E9_spline = arm(Hs, 'LV', base, nB, struct());  R = res.E9_spline;  Ri = res.E1_ss60;
    if ~isfield(R, 'err')
        k = lastBeat(R);  Fi = interp1(Ri.t, Ri.F, R.t(k));
        fprintf('E9 C2 spline host signal: %.1f s/beat, steps %d (failed %d), RHS calls %d, Jacobians %d; pchip: steps %d (failed %d), RHS calls %d, Jacobians %d; last beat max|dF| %.2f kPa\n', ...
            R.wallPerBeat, R.nSteps, R.nFailed, R.nFevals, R.nPD, Ri.nSteps, Ri.nFailed, Ri.nFevals, Ri.nPD, max(abs(R.F(k) - Fi)));
    end
    save(outF, 'res');
end

res.Hraw = H.raw;  res.T = T;  res.HR = H.HR;
save(outF, 'res');
fprintf('saved %s\n', outF);  diary off;

function R = arm(H, wall, so, nB, ro)
% one probe arm; a failure is returned, not thrown
if ~isfield(ro, 'tPre'), ro.tPre = 0.5; end
try
    P = setupWall(H, wall, so);
    R = runHostBeats(P, nB, ro);
catch ME
    R = struct('err', ME.message);
    fprintf('   arm FAILED (%s, %s): %s\n', wall, jsonencode(so), ME.message);
end
end
