% RunSRXGate.m — why the SRX detour is a SHORTCUT, and the one parameter that
% could turn it back into a trap.
%
% RunSRXRecoil.m found every routing variant makes the post-restretch recovery
% FASTER. The arithmetic of the transit times says why (all rates from
% params/rskR2_w025_opt.m, F_iso = 76 kPa):
%
%   normal route   p2 -> PT -> PD          1/kah                        = 9.8 ms
%   via SRXD       p2 -> P_SRD -> PD       1/(kmsrd*exp(F/sigma_srd1))
%                                          F=76 kPa: 1/342              = 2.9 ms
%                                          F=0     : 1/19.8             = 50.5 ms
%   via SRXT       p2 -> P_SR -> P_SRD -> PD   1/ksr2srd + the above
%                                          F=76 kPa: 9.1 + 2.9          = 12.0 ms
%                                          F=0     : 9.1 + 50.5         = 59.6 ms
%
% The SRX exit is force-ACCELERATED (exp(+F/sigma_srd1)), so SRX is a fast lane
% at high force and a slow lane at low force. The post-restretch window starts at
% 0.86 F_iso; the post-slack and ktr windows start at zero force. Routing
% detachment through SRX therefore speeds up exactly the window that is already
% 2.35x too fast, and slows the two that are already right.
%
% sigma_srd1 -> large removes that force dependence (exit = kmsrd at every force),
% which is the only way the detour becomes a trap where it is needed. The
% question this script answers is what that costs at the operating point —
% because the same gate carries the NORMAL SRX flow, which is load-bearing for
% isometric force.
%
% C0* are the gate-alone controls (no routing), so the two effects separate.

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

rip = @(s0, varargin) struct('SRXFromR2HighStrain', 1, 'sSRXrip', s0, ...
                             'dSRXrip', 0.005, 'SRXRecoilToD', true, varargin{:});

V = {};
V{end+1} = {'C0a gate only sg=100',   struct('sigma_srd1', 100)};
V{end+1} = {'C0b gate only sg=1e6',   struct('sigma_srd1', 1e6)};
V{end+1} = {'C1 rip.020->SRXD sg=100', rip(0.020, 'sigma_srd1', 100)};
V{end+1} = {'C2 rip.020->SRXD sg=1e6', rip(0.020, 'sigma_srd1', 1e6)};
V{end+1} = {'C3 rip.015->SRXD sg=1e6', rip(0.015, 'sigma_srd1', 1e6)};
V{end+1} = {'C4 rip.010->SRXD sg=1e6', rip(0.010, 'sigma_srd1', 1e6)};
% and with the SRXD return rate itself slowed, once it is force-flat
V{end+1} = {'C5 rip.015 sg=1e6 km=5',  rip(0.015, 'sigma_srd1', 1e6, 'kmsrd', 5)};
V{end+1} = {'C6 all->SRXD  sg=1e6',    struct('SRXFromR2', 1, 'SRXRecoilToD', true, 'sigma_srd1', 1e6)};

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
save(fullfile(here,'results','gate.mat'), 'C', 'V', 'amp', 'rsK_d');

fprintf('\n%-24s %8s %8s %8s %8s %8s\n', 'variant', 'L2', 'rsK x', 'ktr x', 'steady', 'SRX_ss');
sx = @(r) getfielddef(r.features_model, 'SRX_ss');
fprintf('%-24s %8.3f %8.2f %8.2f %8.1f %8.3f\n', 'baseline', B.L2, B.rsK_x, B.ktr_x, ...
        B.steady, sx(B));
for i = 1:numel(C)
    if ~C{i}.ok; continue; end
    fprintf('%-24s %8.3f %8.2f %8.2f %8.1f %8.3f\n', C{i}.name, C{i}.L2, C{i}.rsK_x, ...
            C{i}.ktr_x, C{i}.steady, sx(C{i}));
end

nm = cellfun(@(x) regexprep(strtok(x,'|'), '\[.*$',''), B.fn, 'uni', 0);
for i = 1:numel(C)
    if ~C{i}.ok; continue; end
    d = C{i}.ct2 - B.ct2;
    [~, o] = sort(abs(d), 'descend');
    fprintf('\n%s   (L2 %.3f -> %.3f)\n', C{i}.name, B.L2, C{i}.L2);
    for j = o(1:min(6,numel(o)))
        if abs(d(j)) < 1e-3; break; end
        fprintf('    %-20s %7.3f -> %8.3f  %+8.3f\n', nm{j}, B.ct2(j), C{i}.ct2(j), d(j));
    end
end

function v = getfielddef(s, f)
    if isfield(s, f); v = mean(s.(f), 'omitnan'); else; v = NaN; end
end
