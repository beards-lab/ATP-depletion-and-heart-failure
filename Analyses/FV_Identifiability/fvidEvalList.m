function E = fvidEvalList(pool, tier)
% fvidEvalList  The list of model evaluations behind the finite-difference Jacobian.
%
%   E(1)            the incumbent (no change)
%   E(2 .. P+1)     each pool parameter x 1.10           (the Jacobian itself)
%   then, for a subset of parameters (the optimiser's compulsory ones plus every
%   eighth of the pool): x 1.05 and x 1/1.10             (step and side check)
%   then three synthetic "truths": every parameter moved at once by a random
%   factor exp(0.10 * N(0,1))                            (recovery test)
%
% tier 'beard' adds pAtp_k2, the exponent of the low-ATP k2 gate, which the 8 mM
% slack battery cannot see. The list is deterministic, so checkpoints stay valid.

STEP = 1.10; NSYNTH = 3; SYNTH_SD = 0.10;
if strcmp(tier, 'beard'), pool = [pool, {'pAtp_k2'}]; end
comp = {'k2', 'ka', 'ksrd', 'sigma_srd2', 'PieceWiseStrainDep2X_logOffset', 'ekSE', 'dr2', 'L_thin', 'k_cr', 'tau_cr'};
sub = unique([find(ismember(pool, comp)), 1:8:numel(pool), find(strcmp(pool, 'pAtp_k2'))]);
E = struct('name', '', 'f', 1, 'kind', 'base');
for k = 1:numel(pool), E(end+1) = struct('name', pool{k}, 'f', STEP, 'kind', 'jac'); end %#ok<AGROW>
for k = sub
    E(end+1) = struct('name', pool{k}, 'f', 1.05, 'kind', 'half'); %#ok<AGROW>
    E(end+1) = struct('name', pool{k}, 'f', 1 / STEP, 'kind', 'minus'); %#ok<AGROW>
end
st = RandStream('mt19937ar', 'Seed', 88);
for d = 1:NSYNTH
    f = exp(SYNTH_SD * randn(st, 1, numel(pool)));
    E(end+1) = struct('name', {pool}, 'f', {f}, 'kind', 'synth'); %#ok<AGROW>
end
end
