% RunDestControl.m — is the SRX exit FORCE GATE the whole reason the recoil has
% the wrong sign?
%
% Control experiment. Keep the entry route identical to the B-variants (only the
% high-strain rupture branch is diverted, so the operating point is untouched)
% and change ONLY the destination's exit law:
%
%   dest    exit rate at the post-restretch force (~0.86 F_iso)   force law
%   ----    -------------------------------------------------    ---------
%   PT      kah = 101.6 /s                                        flat  (baseline)
%   P_SRD   kmsrd*exp(F/sigma_srd1) = 342 /s                       ACCELERATED
%   P_SR    ksr2srd = 110 /s then the above                        flat, then acc.
%   D0      k_D0T*MgATP/(K_D0T+MgATP)                              flat, tunable
%
% If the sign flips as soon as the exit is force-flat and slower than kah, then
% the SRX pool fails here for one reason only — it empties FASTER the harder you
% pull, and the post-restretch window is the high-force window.
%
% k_D0T is scanned across and below kah. Run Baseline.m first.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if isempty(gcp('nocreate')); parpool('local', 5); end

L = load(fullfile(here,'results','baseline.mat'));
amp = L.amp; rsK_d = L.rsK_d; B = L.B;
pBase = getParams(loadParams(L.SNAP), [], true, false);
pBase.MaxRunTime = 600;
pBase.RunForceVelocity = 1; pBase.RunSlack = 1; pBase.RunSlackPassive = 1;
pBase.RunKtr = 0; pBase.RunStairs = 0; pBase.RunForceVelocityTime = 0;

d0 = @(s0, kD0T) struct('SRXFromR2HighStrain', 1, 'sSRXrip', s0, 'dSRXrip', 0.005, ...
                        'SRXRecoilToD0', true, 'UseD0State', true, ...
                        'k_D0T', kD0T, 'K_D0T', 0.15);

V = {};
V{end+1} = {'D0 s0=.015 k=500 (fast)', d0(0.015, 500)};
V{end+1} = {'D0 s0=.015 k=100',        d0(0.015, 100)};
V{end+1} = {'D0 s0=.015 k=40',         d0(0.015,  40)};
V{end+1} = {'D0 s0=.015 k=15',         d0(0.015,  15)};
V{end+1} = {'D0 s0=.010 k=40',         d0(0.010,  40)};
V{end+1} = {'D0 s0=.020 k=40',         d0(0.020,  40)};
V{end+1} = {'D0 s0=.010 k=15',         d0(0.010,  15)};

C = cell(1, numel(V));
for i = 1:numel(V)
    fprintf('[%d/%d] %-24s ', i, numel(V), V{i}{1});
    C{i} = srxRecoilCase(pBase, V{i}{1}, V{i}{2}, amp, rsK_d);
    if C{i}.ok
        fprintf('L2 %7.3f  rsK x%.2f  ktr x%.2f  F %.1f  (%.0f s)\n', ...
            C{i}.L2, C{i}.rsK_x, C{i}.ktr_x, C{i}.steady, C{i}.time);
    else
        fprintf('FAILED: %s\n', C{i}.err);
    end
end
save(fullfile(here,'results','destcontrol.mat'), 'C', 'V', 'amp', 'rsK_d');

fprintf('\n%-24s %8s %8s %8s %8s %9s\n', 'variant', 'L2', 'rsK x', 'ktr x', 'steady', 'slope');
fprintf('%-24s %8.3f %8.2f %8.2f %8.1f %+9.0f\n', 'baseline (-> PT)', B.L2, B.rsK_x, ...
        B.ktr_x, B.steady, B.slope);
for i = 1:numel(C)
    if ~C{i}.ok; continue; end
    fprintf('%-24s %8.3f %8.2f %8.2f %8.1f %+9.0f\n', C{i}.name, C{i}.L2, C{i}.rsK_x, ...
            C{i}.ktr_x, C{i}.steady, C{i}.slope);
end
fprintf('%-24s %8s %8.2f %8.2f %8s %+9.0f\n', 'DATA', '-', 1.00, 1.00, '~80', ...
        polyfit(amp(:), rsK_d(:), 1)*[1;0]);

nm = cellfun(@(x) regexprep(strtok(x,'|'), '\[.*$',''), B.fn, 'uni', 0);
for i = 1:numel(C)
    if ~C{i}.ok; continue; end
    d = C{i}.ct2 - B.ct2;
    [~, o] = sort(abs(d), 'descend');
    fprintf('\n%s   (L2 %.3f -> %.3f)   rsK per cycle %s\n', C{i}.name, B.L2, C{i}.L2, ...
            mat2str(round(C{i}.rsK_m,1)));
    for j = o(1:min(6,numel(o)))
        if abs(d(j)) < 1e-3; break; end
        fprintf('    %-20s %7.3f -> %8.3f  %+8.3f\n', nm{j}, B.ct2(j), C{i}.ct2(j), d(j));
    end
end
fprintf('\ndata rsK per cycle %s\n', mat2str(round(rsK_d,1)));
