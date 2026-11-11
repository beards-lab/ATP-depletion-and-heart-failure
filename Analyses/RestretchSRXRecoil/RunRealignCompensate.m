% RunRealignCompensate.m — is the compliant-realignment bill payable?
%
% M1 buys `rsK` at full authority (x2.35 -> x0.76 across the gain scan) with `ktr`,
% `steady` and FV untouched. Its entire bill lands in the restretch TRANSIENT
% SHAPE: vall_y, vall2_dy, ovrsht_dy, peak2, coolDownLS. That is a specific,
% diagnosable failure — the deficit suppresses re-attachment during and just after
% the ramp, so the valley goes too deep and the recovery starts too low.
%
% Those features have their own levers, none of which touch rsK's mechanism:
%   kA2re   PT -> p2 re-attachment during lengthening; adds NEW low-strain heads
%           through the valley. Zero in the baseline, and the one lever previously
%           shown to fill vall_y toward data while holding peak1_y.
%   eta_M   Maxwell dashpot; RestretchRecoveryFit 6b's cheap rsK win (x2.70->x1.58
%           for +9.8) and it is load-bearing for ovrsht_dy / vall2_dy.
% If either pays a useful share, the mechanism is a refit candidate; if neither
% does, the damage is intrinsic and M1 needs a different engagement profile.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if isempty(gcp('nocreate')); parpool('local', 5); end

L = load(fullfile(here,'results','baseline.mat'));
amp = L.amp; rsK_d = L.rsK_d; B = L.B;
pBase = getParams(loadParams(L.SNAP), [], true, false);
pBase.MaxRunTime = 150;
pBase.RunForceVelocity = 1; pBase.RunSlack = 1; pBase.RunSlackPassive = 1;
pBase.RunKtr = 0; pBase.RunStairs = 0; pBase.RunForceVelocityTime = 0;

m1 = @(k, varargin) struct('UseLoadRealign', true, 'k_cr', k, 'tau_cr', 0.020, ...
                           'F_cr', 25, varargin{:});

V = {};
V{end+1} = {'M1 k=3 (reference)',    m1(3)};
V{end+1} = {'M1 k=3 + kA2re=10',     m1(3, 'kA2re', 10)};
V{end+1} = {'M1 k=3 + kA2re=30',     m1(3, 'kA2re', 30)};
V{end+1} = {'M1 k=3 + etaM x0.3',    m1(3, 'eta_M', 0.3*pBase.eta_M)};
V{end+1} = {'M1 k=3 + kA2re30+etaM', m1(3, 'kA2re', 30, 'eta_M', 0.3*pBase.eta_M)};
V{end+1} = {'M1 k=6 + kA2re=30',     m1(6, 'kA2re', 30)};
V{end+1} = {'M1 k=6 + kA2re30+etaM', m1(6, 'kA2re', 30, 'eta_M', 0.3*pBase.eta_M)};
V{end+1} = {'kA2re=30 alone (ctrl)', struct('kA2re', 30)};

matf = fullfile(here,'results','compensate.mat');
if exist(matf,'file'); S = load(matf); C = S.C; else; C = cell(1, numel(V)); end
if numel(C) < numel(V); C(numel(C)+1:numel(V)) = {[]}; end
if ~exist('SEL','var') || isempty(SEL); SEL = 1:numel(V); end

for i = SEL(:)'
    fprintf('[%d/%d] %-22s ', i, numel(V), V{i}{1});
    C{i} = srxRecoilCase(pBase, V{i}{1}, V{i}{2}, amp, rsK_d);
    if C{i}.ok
        fprintf('L2 %8.3f  rsK x%.2f  ktr x%.2f  F %.1f  (%.0f s)\n', ...
            C{i}.L2, C{i}.rsK_x, C{i}.ktr_x, C{i}.steady, C{i}.time);
    else
        fprintf('FAILED: %s\n', C{i}.err);
    end
    save(matf, 'C', 'V', 'amp', 'rsK_d');
end
keepC = ~cellfun(@isempty, C);

nm = cellfun(@(x) regexprep(strtok(x,'|'), '\[.*$',''), B.fn, 'uni', 0);
shape = ismember(nm, {'vall_y','vall2_dy','ovrsht_dy','peak2','coolDownLS'});
fprintf('\n%-22s %9s %8s %8s %8s %10s\n', 'variant', 'L2', 'rsK x', 'ktr x', 'steady', 'shapeL2');
fprintf('%-22s %9.3f %8.2f %8.2f %8.1f %10.3f\n', 'baseline', B.L2, B.rsK_x, B.ktr_x, ...
        B.steady, sum(B.ct2(shape)));
for i = find(keepC)
    if isempty(C{i}) || ~C{i}.ok; continue; end
    fprintf('%-22s %9.3f %8.2f %8.2f %8.1f %10.3f\n', C{i}.name, C{i}.L2, C{i}.rsK_x, ...
            C{i}.ktr_x, C{i}.steady, sum(C{i}.ct2(shape)));
end
fprintf('%-22s %9s %8.2f %8.2f %8s %10s\n', 'DATA', '-', 1.00, 1.00, '~80', '-');

for i = find(keepC)
    if isempty(C{i}) || ~C{i}.ok; continue; end
    d = C{i}.ct2 - B.ct2; [~, o] = sort(abs(d), 'descend');
    fprintf('\n%s   (L2 %.3f)   rsK/cycle %s\n', C{i}.name, C{i}.L2, mat2str(round(C{i}.rsK_m,1)));
    for j = o(1:min(6,numel(o)))
        if abs(d(j)) < 1e-3; break; end
        fprintf('    %-20s %7.3f -> %8.3f  %+8.3f\n', nm{j}, B.ct2(j), C{i}.ct2(j), d(j));
    end
end
fprintf('\ndata rsK/cycle %s\n', mat2str(round(rsK_d,1)));
