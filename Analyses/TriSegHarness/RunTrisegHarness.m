%% RunTrisegHarness  The closed loop to a periodic state (#97).
%
% P21 host fit Nx51 with the full cross-bridge model (robust1_seed) in all three
% walls: host passive law, host Ca(t), k_on x KONMULT, velocity signal a 1 ms low-pass
% of dLs/dt, 60 strain bins. Runs blocks of 10 beats until the LV peak pressure and
% stroke volume change by less than 1 % from one beat to the next (the host's own
% criterion is of this kind), up to MAXBEATS.
%   matlab -batch "run('Analyses/TriSegHarness/RunTrisegHarness.m')"
% Results -> results/harness_<TAG>.mat; log -> results/harness_<TAG>.log

here = fileparts(mfilename('fullpath'));  root = fileparts(fileparts(here));
cd(root);  addpath(genpath(fullfile(root, 'Model')));  addpath(genpath(fullfile(root, 'Auxiliary')));
addpath(fullfile(root, 'params'));  addpath(fullfile(root, 'Organ'));
if ~exist('HOSTROOT', 'var'), HOSTROOT = fullfile(fileparts(root), 'Multiscale-Hx-P21-Model'); end
if ~exist('KONMULT', 'var'),  KONMULT = 3; end
if ~exist('MAXBEATS', 'var'), MAXBEATS = 60; end
if ~exist('OPT', 'var'),      OPT = struct(); end
if ~exist('TAG', 'var'),      TAG = sprintf('kon%g', KONMULT); end
OPT.konMult = KONMULT;
if ~exist(fullfile(here, 'results'), 'dir'), mkdir(fullfile(here, 'results')); end
diary(fullfile(here, 'results', ['harness_' TAG '.log']));  diary on;

H = hostP21Load('Nx51', HOSTROOT);
C = setupTrisegXB(H, OPT);
fprintf('\n=== RunTrisegHarness %s: %d states, k_on x%g, tauV %.1f ms, host passive %d, native viscous %d ===\n', TAG, C.n, C.opt.konMult, 1e3*C.opt.tauV, C.opt.hostPassive, C.opt.nativeViscous);
B = table();  R = [];  done = false;  nb = 0;  tAll = tic;  S = struct('nSteps', 0, 'nFailed', 0, 'nFevals', 0, 'nPD', 0, 'nShift', zeros(1, 3), 'massLost', zeros(1, 3));
while ~done && nb < MAXBEATS
    R = runTrisegXB(C, 10, struct('verbose', false));
    C = R.C;  b = R.beat;  b.beat = b.beat + nb;  B = [B; b];  nb = nb + 10; %#ok<AGROW>
    for f = {'nSteps', 'nFailed', 'nFevals', 'nPD', 'nShift', 'massLost'}, S.(f{1}) = S.(f{1}) + R.(f{1}); end
    dP = abs(diff(B.Pmax_LV))./B.Pmax_LV(2:end);  dS = abs(diff(B.SV_LV))./B.SV_LV(2:end);
    fprintf('beats %2d-%2d (%.0f s): Pmax_LV %.3f kPa, SV_LV %.2f uL, EDV %.2f uL; last beat-to-beat change: Pmax %.2f %%, SV %.2f %%\n', ...
        nb - 9, nb, toc(tAll), B.Pmax_LV(end), 1e3*B.SV_LV(end), 1e3*B.EDV_LV(end), 100*dP(end), 100*dS(end));
    done = all(dP(end-2:end) < 0.01) && all(dS(end-2:end) < 0.01);
    save(fullfile(here, 'results', ['harness_' TAG '.mat']), 'R', 'B', 'S', 'H', 'nb');
end
wall = toc(tAll);
fprintf('\nperiodic: %d after %d beats; %.1f s per beat; per beat: %.0f steps (%.0f failed), %.0f RHS calls, %.1f Jacobians; slides per beat LV/SEP/RV %s; mass lost at slides %s\n', ...
    done, nb, wall/nb, S.nSteps/nb, S.nFailed/nb, S.nFevals/nb, S.nPD/nb, mat2str(S.nShift/nb, 2), mat2str(S.massLost, 2));
disp(B(max(1, end-4):end, :));
d = H.data;  l = B(end, :);  r = H.ref;
fprintf('last beat against the host fit and the data (kPa, uL):\n');
fprintf('   %-10s %10s %10s %10s\n', '', 'coupled', 'host fit', 'data');
fprintf('   %-10s %10.1f %10.1f %10.1f\n', 'LV EDV', 1e3*l.EDV_LV, 1e3*max(r.V(:, 1)), 1e3*d.EDV_LV);
fprintf('   %-10s %10.1f %10.1f %10.1f\n', 'LV ESV', 1e3*l.ESV_LV, 1e3*min(r.V(:, 1)), 1e3*d.ESV_LV);
fprintf('   %-10s %10.2f %10.2f %10.2f\n', 'LV EDP', l.EDP_LV, r.P(1, 1), d.EDP_LV);
fprintf('   %-10s %10.2f %10.2f %10.2f\n', 'LV Pmax', l.Pmax_LV, max(r.P(:, 1)), d.ESP_LV);
fprintf('   %-10s %10.1f %10.1f %10.1f\n', 'RV EDV', 1e3*l.EDV_RV, 1e3*max(r.V(:, 2)), 1e3*d.EDV_RV);
fprintf('   %-10s %10.1f %10.1f %10.1f\n', 'RV ESV', 1e3*l.ESV_RV, 1e3*min(r.V(:, 2)), 1e3*d.ESV_RV);
fprintf('   %-10s %10.2f %10.2f %10.2f\n', 'RV Pmax', l.Pmax_RV, max(r.P(:, 4)), d.ESP_RV);
T = H.T;  k = R.t >= 9*T - 1e-9;  tk = R.t(k) - 9*T;  dtk = [diff(tk); 0];  c = R.col;
for w = 1:3
    sg = R.out(k, c.sigma(w));
    fprintf('   wall %d: Ls %.3f-%.3f um; stress %.2f-%.1f kPa; slack (< 0.05 kPa) %.1f ms; gate up to %.2f; peak shortening %.1f um/s\n', w, min(R.out(k, c.Ls(w))), max(R.out(k, c.Ls(w))), ...
        min(sg), max(sg), 1e3*sum(dtk(sg < 0.05)), max(R.out(k, c.gate(w))), -min(R.out(k, c.v(w))));
end
Q = R.out(k, c.Q_a);  ej = tk(Q > 0.01*max(Q));
fprintf('   LV ejection: %.1f ms (%.1f-%.1f ms of the beat)\n', 1e3*(max(ej) - min(ej)), 1e3*min(ej), 1e3*max(ej));
diary off;
