% RunSRXRecoil.m — does routing the strong-bridge detachment into SRX slow the
% post-restretch recovery, and what does it cost?
%
% Two hypotheses, both destination switches (no rate is capped, so the restretch
% and ktr-overstretch force peaks that pin R2(s) are untouched by construction —
% see Analyses/ApoPoolDetachment for why a rate cap is falsified):
%
%   A  SRXFromR2 = f            the WHOLE detachment flux lands in SRX instead of
%                               PT. f = 1 makes SRX an in-series step of every
%                               cycle. Affects every protocol.
%   B  SRXFromR2HighStrain = f  only the high-strain (mechanical-rupture) branch
%                               recoils into SRX. The chemical path at the
%                               operating point is untouched, so the detour only
%                               exists while the fibre is stretched.
%
% Destination is P_SR (SRX.ATP) or, with SRXRecoilToD, P_SRD (SRX.ADP) — the
% nucleotide-consistent landing state for a bridge torn off before ADP release.
%
% Run Baseline.m first. ~20 s per variant.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if isempty(gcp('nocreate')); parpool('local', 5); end

L = load(fullfile(here,'results','baseline.mat'));   % B, amp, rsK_d, SNAP
amp = L.amp; rsK_d = L.rsK_d;
pBase = getParams(loadParams(L.SNAP), [], true, false);
pBase.MaxRunTime = 600;
pBase.RunForceVelocity = 1; pBase.RunSlack = 1; pBase.RunSlackPassive = 1;
pBase.RunKtr = 0; pBase.RunStairs = 0; pBase.RunForceVelocityTime = 0;

V = {};
V{end+1} = {'A1 all->SRXT  f=1.00',  struct('SRXFromR2', 1.00)};
V{end+1} = {'A2 all->SRXT  f=0.50',  struct('SRXFromR2', 0.50)};
V{end+1} = {'A3 all->SRXT  f=0.25',  struct('SRXFromR2', 0.25)};
V{end+1} = {'A4 all->SRXD  f=1.00',  struct('SRXFromR2', 1.00, 'SRXRecoilToD', true)};
V{end+1} = {'A5 all->SRXD  f=0.50',  struct('SRXFromR2', 0.50, 'SRXRecoilToD', true)};

hs = @(f, s0, dst) struct('SRXFromR2HighStrain', f, 'sSRXrip', s0, ...
                          'dSRXrip', 0.005, 'SRXRecoilToD', dst);
V{end+1} = {'B1 rip->SRXT s0=.015', hs(1.00, 0.015, false)};
V{end+1} = {'B2 rip->SRXD s0=.015', hs(1.00, 0.015, true)};
V{end+1} = {'B3 rip->SRXD s0=.010', hs(1.00, 0.010, true)};
V{end+1} = {'B4 rip->SRXD s0=.020', hs(1.00, 0.020, true)};
V{end+1} = {'B5 rip->SRXD s0=.005', hs(1.00, 0.005, true)};
V{end+1} = {'B6 rip->SRXD s0=.010 f=.5', hs(0.50, 0.010, true)};
V{end+1} = {'B7 rip->SRXT s0=.005', hs(1.00, 0.005, false)};

C = cell(1, numel(V));
for i = 1:numel(V)
    fprintf('[%2d/%2d] %-26s ', i, numel(V), V{i}{1});
    C{i} = srxRecoilCase(pBase, V{i}{1}, V{i}{2}, amp, rsK_d);
    if C{i}.ok
        fprintf('L2 %7.3f  rsK x%.2f  ktr x%.2f  F %.1f  (%.0f s)\n', ...
            C{i}.L2, C{i}.rsK_x, C{i}.ktr_x, C{i}.steady, C{i}.time);
    else
        fprintf('FAILED: %s\n', C{i}.err);
    end
end
save(fullfile(here,'results','variants.mat'), 'C', 'V', 'amp', 'rsK_d');

%% ---- summary table -------------------------------------------------------
B = L.B;
fprintf('\n%-26s %8s %8s %8s %8s %8s %8s\n', 'variant', 'L2', 'L1', 'rsK x', 'ktr x', 'steady', 'slope');
fprintf('%-26s %8.3f %8.2f %8.2f %8.2f %8.1f %+8.0f\n', 'baseline', B.L2, B.L1, ...
        B.rsK_x, B.ktr_x, B.steady, B.slope);
for i = 1:numel(C)
    if ~C{i}.ok; continue; end
    fprintf('%-26s %8.3f %8.2f %8.2f %8.2f %8.1f %+8.0f\n', C{i}.name, C{i}.L2, C{i}.L1, ...
            C{i}.rsK_x, C{i}.ktr_x, C{i}.steady, C{i}.slope);
end
fprintf('%-26s %8s %8s %8.2f %8.2f %8.1f %+8.0f\n', 'DATA', '-', '-', 1.00, 1.00, ...
        mean(B.features_model.steady)*0 + 80.0, polyfit(amp(:), rsK_d(:), 1)*[1;0]);

%% ---- which features get better, which worse ------------------------------
% L2 per feature vs baseline, biggest movers first.
nm = cellfun(@(x) regexprep(strtok(x,'|'), '\[.*$',''), B.fn, 'uni', 0);
for i = 1:numel(C)
    if ~C{i}.ok; continue; end
    d = C{i}.ct2 - B.ct2;
    [~, o] = sort(abs(d), 'descend');
    fprintf('\n%s   (L2 %.3f -> %.3f)\n', C{i}.name, B.L2, C{i}.L2);
    for j = o(1:min(7,numel(o)))
        if abs(d(j)) < 1e-3; break; end
        fprintf('    %-22s %7.3f -> %7.3f  %+7.3f  %s\n', nm{j}, B.ct2(j), C{i}.ct2(j), d(j), ...
                repmat('+', 1, min(20, max(0, round(d(j)*4)))));
    end
end
