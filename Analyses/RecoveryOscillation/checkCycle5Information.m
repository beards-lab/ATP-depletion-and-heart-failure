% checkCycle5Information.m -- the two measurements the `lean5` campaign rests
% on. Read-only; run it before trusting either claim.
%
% CLAIM 1: slack cycle 5 carries real between-preparation information, so
% restricting rsK80 to cycles 1-4 (sec. 3 of conclusions.md) threw away signal
% along with noise.
%
%   The sec. 3 argument was: injecting each run's own measured noise puts
%   rsK80's measurement noise at 32.7 % against a between-preparation spread of
%   16.9 %. A feature cannot be better determined between preparations than its
%   own noise allows, so the 16.9 % must be an n = 3 fluke -- and cycle 5, at
%   29-245 % per-cycle noise against 6-15 % for cycles 1-4, is the cause.
%
%   That is sound about the noise. What it does not test is whether cycle 5
%   carries an ATP signal, and it does: its 2 mM / 8 mM ratio lands INSIDE the
%   range spanned by cycles 1-4. Transducer noise does not reproduce a x0.39
%   ATP effect across three preparations to within a few percent of four other
%   cycles. Both observations can be true only if the noise estimate is
%   incomplete -- and the ATP ratio is the harder evidence, because it is a
%   ratio of two independent measurements rather than a model of the residual.
%
%   Cycle 5 is also not a repeat of cycles 1-4. It is the ONLY re-stretch in
%   the protocol that returns to ~2.0 um SL instead of 2.2, at the smallest
%   amplitude and the slowest ramp. Dropping it removes the protocol's only
%   probe of how the re-stretch recovery depends on the length it recovers at
%   -- and, with the lattice attenuation biting below SL ~ 2.09 um, the only
%   manoeuvre from which d_optimal and sigma_lattice can be identified.
%
% CLAIM 2: `xrate` is a working search coordinate, so MechanismAblation sec.
% 7.4's "add it to the pool" is safe to act on. This is not obvious: Model/
% updateRates.m does rmfield(params,'xrate'), so it is reasonable to fear the
% modifier is discarded before it reaches anything.
%
%   matlab -batch "run('Analyses/RecoveryOscillation/checkCycle5Information.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

%% ---- claim 1 -----------------------------------------------------------
cases = {'1:4', 1:4; '5 only', 5; '1:5', 1:5};
s8 = nan(1,3); zz = nan(1,3);
for i = 1:size(cases,1)
    T = buildPooledTarget(struct('segments', struct('rsK80', cases{i,2})));
    k = find(strcmp(T.name,'rsK80'), 1);
    s8(i) = T.s8(k);
    if isfield(T,'z8'); zz(i) = T.z8(k); elseif isfield(T,'z'); zz(i) = T.z(k); end
end
T  = buildPooledTarget(struct('segments', struct('rsK80', 1:5)));
k  = find(strcmp(T.name,'rsK80'), 1);
kA = find(strcmp(T.name,'rsA'),   1);
kS = find(strcmp(T.name,'steady'),1);

fprintf('\n=== CLAIM 1: is cycle 5 information or noise? ===\n\n');
fprintf('%-10s %10s %10s %10s\n', 'segments', 's8 (%)', 'w=1/s^2', 'z (ATP)');
for i = 1:size(cases,1)
    fprintf('%-10s %10.1f %10.1f %10.2f\n', cases{i,1}, 100*s8(i), 1/s8(i)^2, zz(i));
end
fprintf('\nrsK80 per cycle, 8 mM : '); fprintf('%8.2f', T.m8(k,:)); fprintf('\n');
fprintf('rsK80 per cycle, 2 mM : '); fprintf('%8.2f', T.m2(k,:)); fprintf('\n');
r = T.m2(k,:)./T.m8(k,:);
fprintf('ATP ratio 2/8         : '); fprintf('%8.3f', r); fprintf('\n');
fprintf('  cycle 5 ratio %.3f vs cycles 1-4 range %.3f-%.3f -> %s\n', ...
        r(5), min(r(1:4)), max(r(1:4)), ...
        iVerdict(r(5) >= min(r(1:4)) && r(5) <= max(r(1:4))));

fprintf('\ncycle 5 is a DIFFERENT manoeuvre, not a noisy repeat:\n');
fprintf('  re-stretch amplitude rsA : '); fprintf('%8.3f', T.m8(kA,:)); fprintf('\n');
fprintf('  plateau `steady`         : '); fprintf('%8.3f', T.m8(kS,:)); fprintf('\n');
fprintf('  (cycle 5 returns to ~2.0 um SL; cycles 1-4 to 2.2)\n');

%% ---- claim 2 -----------------------------------------------------------
fprintf('\n=== CLAIM 2: does an xrate modifier reach the rates? ===\n\n');
base = getParams(); base.mods = {'xrate'};
a = getParams(base, 1.0, true, true);
b = getParams(base, 1.5, true, true);
scaled = {'k2','ka','kah','kamh','kd','k1','k_1'};
ok = true;
for i = 1:numel(scaled)
    f = scaled{i}; ratio = b.(f)/a.(f);
    ok = ok && abs(ratio - 1.5) < 1e-9;
    fprintf('  %-6s x%.4f\n', f, ratio);
end
fprintf('  %-6s x%.4f   (SRX: must NOT scale)\n', 'ksr0', b.ksr0/a.ksr0);
fprintf('\n  core cycle scales exactly with the modifier : %s\n', iVerdict(ok));
fprintf('  SRX left alone                              : %s\n', ...
        iVerdict(abs(b.ksr0/a.ksr0 - 1) < 1e-12));
fprintf(['\n  So xrate is a usable common-mode coordinate. Its VALUE is still\n' ...
         '  meaningless in isolation -- it trades back into the individual\n' ...
         '  rates -- so report it as a search coordinate, never as a result.\n\n']);

function s = iVerdict(tf)
if tf; s = 'YES'; else; s = 'NO'; end
end
