% Baseline.m — regression + reference point for the SRX-recoil test.
%
% 1. Proves the new routing is INERT at its defaults (the wiring must not move
%    the current fit by itself).
% 2. Records the reference numbers every variant is compared against.
%
% Base: params/rskR2_w025_opt.m, protocol 03/27 8 mM (same base as
% Analyses/RestretchRecoveryFit §6b, so its numbers are directly comparable).

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if isempty(gcp('nocreate')); parpool('local', 5); end
if ~exist(fullfile(here,'results'), 'dir'); mkdir(fullfile(here,'results')); end

SNAP = 'params/rskR2_w025_opt.m';
pBase = getParams(loadParams(SNAP), [], true, false);
pBase.MaxRunTime = 600;
pBase.RunForceVelocity = 1; pBase.RunSlack = 1; pBase.RunSlackPassive = 1;
pBase.RunKtr = 0; pBase.RunStairs = 0; pBase.RunForceVelocityTime = 0;

% new params must be present and off
fprintf('defaults: SRXFromR2=%g  SRXFromR2HighStrain=%g  sSRXrip=%g  dSRXrip=%g  SRXRecoilToD=%d\n', ...
    pBase.SRXFromR2, pBase.SRXFromR2HighStrain, pBase.sSRXrip, pBase.dSRXrip, pBase.SRXRecoilToD);
assert(pBase.SRXFromR2 == 0 && pBase.SRXFromR2HighStrain == 0, 'routing not off by default');

ds    = load('data/protocol_03_27_2026_8mM_slack.mat');
vt    = ds.velocitytable;
iR    = find(vt(:,2) > 1);
amp   = (vt(iR+1,4) - vt(iR,4))';       % re-stretch amplitude per cycle (ML)
rsK_d = ds.features_data.rsK;

B = srxRecoilCase(pBase, 'baseline', struct(), amp, rsK_d);
if ~B.ok; error('baseline failed: %s', B.err); end

fprintf('\n---- baseline (routing off) ----\n');
fprintf('featTotal  L1 %.3f | L2 %.3f   (%.0f s)\n', B.L1, B.L2, B.time);
fprintf('rsK model %s\n', mat2str(round(B.rsK_m,1)));
fprintf('rsK data  %s   ratio x%.2f\n', mat2str(round(rsK_d,1)), B.rsK_x);
fprintf('ktr x%.2f | steady %.1f kPa | rsA x%.2f\n', B.ktr_x, B.steady, B.rsA_x);
fprintf('rsK vs amplitude: slope %+.0f (data %+.0f)\n', B.slope, ...
    polyfit(amp(:), rsK_d(:), 1)*[1;0]);

% Second run with the routing explicitly zeroed the long way round (fraction 0
% but the high-strain ramp configured) — must be bit-comparable.
Z = srxRecoilCase(pBase, 'inert-check', ...
    struct('SRXFromR2', 0, 'SRXFromR2HighStrain', 0, 'sSRXrip', 0.010, 'SRXRecoilToD', true), ...
    amp, rsK_d);
fprintf('\ninert check: L2 %.6f vs %.6f  (dL2 = %.2g)\n', Z.L2, B.L2, Z.L2 - B.L2);
assert(abs(Z.L2 - B.L2) < 1e-9, 'routing is NOT inert when the fractions are 0');
fprintf('PASS: wiring is inert at zero fraction.\n');

save(fullfile(here,'results','baseline.mat'), 'B', 'amp', 'rsK_d', 'SNAP');
