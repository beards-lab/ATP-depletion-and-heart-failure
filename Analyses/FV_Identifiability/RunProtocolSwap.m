% RunProtocolSwap.m -- how much of the measured ATP ratio is protocol, in the model?
%
% The 4 and 2 mM sheets use a different length protocol from the 8 mM sheets (a
% 0.01 Lo larger pre-ramp step and a shorter ramp at 0.5-4 ML/s). Three extra
% simulations of the incumbent separate the two effects:
%   8 mM model on the low-ATP protocol     -> protocol effect alone
%   4 and 2 mM model on the 8 mM protocol  -> ATP effect alone
% Compared with the three matched simulations (jac_beard/e001_r*.mat).
%
%   matlab -batch "run('Analyses/FV_Identifiability/RunProtocolSwap.m')"
% Output -> results/protocol_swap.mat and the table printed below.

here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
cd(root); addpath(genpath(root));
wt = strsplit(genpath(fullfile(root, '.claude')), pathsep); wt = wt(~cellfun(@isempty, wt));
if ~isempty(wt), rmpath(wt{:}); end
[~, pBeard, ~, recs] = fvidBase();
base = struct('name', '', 'f', 1, 'kind', 'base');
% {model [ATP], protocol index into recs (1 = 8 mM, 2 = 4 mM, 3 = 2 mM)}
C = [8 2; 4 1; 2 1];
sw = fullfile(here, 'results', 'protocol_swap_sims.mat');
if isfile(sw)
    S = load(sw); Fs = S.Fs;
else
    Fs = cell(1, 3);
    pl = gcp('nocreate'); if isempty(pl), parpool('local', 3); end
    parfor c = 1:3
        r = recs(C(c, 2)); r.atp = C(c, 1); %#ok<PFBNS>
        Fs{c} = simBeard(pBeard, base, r);
    end
    save(sw, 'Fs', 'C');
end
Fm = cell(1, 3);
for a = 1:3, r = load(fullfile(here, 'results', 'jac_beard', sprintf('e001_r%d.mat', a))); Fm{a} = r.r.F; end

rd = @(rec, F) readIsoV(rec.t, rec.L, double(F));
pt = @(R) (R.Fx - R.F0)';
[vS, iv] = sort(round(rd(recs(1), Fm{1}).vMeas * 2) / 2);
P88 = pt(rd(recs(1), Fm{1}));  P88 = P88(iv);        % 8 mM model, 8 mM protocol
P8L = pt(rd(recs(2), Fs{1}));  P8L = P8L(iv);        % 8 mM model, low-ATP protocol
P4L = pt(rd(recs(2), Fm{2}));  P4L = P4L(iv);        % 4 mM model, its own protocol
P2L = pt(rd(recs(3), Fm{3}));  P2L = P2L(iv);        % 2 mM model, its own protocol
P48 = pt(rd(recs(1), Fs{2}));  P48 = P48(iv);        % 4 mM model, 8 mM protocol
P28 = pt(rd(recs(1), Fs{3}));  P28 = P28(iv);        % 2 mM model, 8 mM protocol
nz = @(x) x .* (x > 0.01 * P88(1)) ./ (x > 0.01 * P88(1));   % NaN where the model is slack
T = table(vS(:), nz(P8L)' ./ nz(P88)', nz(P4L)' ./ nz(P88)', nz(P48)' ./ nz(P88)', nz(P4L)' ./ nz(P8L)', ...
    nz(P2L)' ./ nz(P88)', nz(P28)' ./ nz(P88)', nz(P2L)' ./ nz(P8L)', 'VariableNames', ...
    {'v', 'protocolOnly', 'r4_asMeasured', 'r4_on8protocol', 'r4_onLowProtocol', 'r2_asMeasured', 'r2_on8protocol', 'r2_onLowProtocol'});
fprintf('\nModel force ratios at the 1.00 Lo read (NaN = model slack):\n');
fprintf('  protocolOnly     = 8 mM model on the low-ATP protocol / 8 mM model on the 8 mM protocol\n');
fprintf('  r4_asMeasured    = 4 mM on its protocol / 8 mM on its protocol (what the data ratio is)\n');
fprintf('  r4_on8protocol   = 4 mM / 8 mM, both on the 8 mM protocol;  r4_onLowProtocol = both on the low-ATP protocol\n');
disp(T);
save(fullfile(here, 'results', 'protocol_swap.mat'), 'T');
