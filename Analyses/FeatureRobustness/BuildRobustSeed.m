% BuildRobustSeed.m -- refit9's objective with the robust reader and the ramp-end delay.
%
% Analyses/FeatureRobustness found that most of the objective's numerical noise
% was the READER (the model read at 1 kHz on a grid anchored on its first
% stored time, the data at 4 kHz) and the rest was rsT20, whose anchor is the
% argmin of a flat valley -- a feature whose between-heart spread equals its own
% measurement noise. This seed changes exactly that and nothing else:
%   numerics  SiLU 0.01, clamp, RelTol 1e-3, dense 0.5 ms INSIDE the cycles only
%             (odeOutputMaxSegment 5): 1.16x the old run time, not 2.45x
%   reader    FeatureReadGrid = 'data' (model read at the data's 0.25 ms, on an
%             absolute grid; t0_crossing interpolated, seam artefact gone)
%             SkipOscillationFit = 1 (ktr2_* is scored by nothing)
%   objective rsT20[1:4] -> rsTfitE[1:4], same construction as every ATP term of
%             refit9: fl = nEl*s^2/3, w = LAMBDA/s^2, s = its between-prep spread;
%             guards re-derived at GUARD_X times THIS point's raw cost under the
%             new reader (refit9's rule), since the reader moves their values
% Starting point: refit9_opt, unchanged parameters.
%
%   matlab -batch "run('Analyses/FeatureRobustness/BuildRobustSeed.m')"
here = fileparts(mfilename('fullpath')); root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
res = fullfile(here, 'results');
if ~exist('SEEDTAG', 'var'), SEEDTAG = 'robust1_seed'; end
if ~exist('FROM', 'var'), FROM = 'refit9_opt'; end
if ~exist('NNEIGH', 'var'), NNEIGH = 4; end
LAMBDA = 0.0106294; GUARD_X = 1.2;                 % as BuildRefit9Seed
NUM = struct('SmoothRectifierEps', 0.01, 'ClampStrainDepExtrapolation', 1, 'odeRelTol', 1e-3, 'odeAbsTol', 1e-6, ...
             'odeOutputDt', 5e-4, 'odeOutputMaxSegment', 5, 'FeatureReadGrid', 'data', 'SkipOscillationFit', 1);
GUARD = {'A','peak1_y','vall_y','peak2','vall2_dy','restretchSlopeStart','A_over_steady','vall_t','peak1_dSL'};
ATPC  = {'ktrfit','ktrfit_shape','rsK80','rsTfitE','ktrP_k80','steady','t0_crossing'};
base  = @(e) regexprep(regexprep(e, '\|.*$', ''), '\[.*$', '');
if isempty(gcp('nocreate')), parpool('Threads', 5); end

%% ---- the rsTfitE term, from the same three preparations ------------------
R = poolLoadRuns();
V = nan(3, 4); d = 0;
for i = 1:numel(R)
    if ~strcmp(R(i).cond, '8mM'), continue; end
    d = d + 1; V(d, :) = R(i).feat.rsTfitE(1:4);
end
s = mean(std(V, 0, 1) ./ abs(mean(V, 1)));        % BuildRefit9Seed's cv()
termR = sprintf('rsTfitE[1:4]|fl:%.6g|%.6g', 4*s^2/3, LAMBDA/s^2);
fprintf('rsTfitE: per-prep 8 mM (cycles 1-4)\n'); disp(round(V, 5));
fprintf('  s = %.4f  ->  %s\n', s, termR);

%% ---- objective ------------------------------------------------------------
p = getParams(loadParams(FROM), [], true, false);
for q = fieldnames(NUM)', p.(q{1}) = NUM.(q{1}); end
p.RunForceVelocity = true; p.RunSlack = true; p.RunStairs = false; p.RunKtr = true;
p.RunForceVelocityTime = false; p.RunSlackPassive = true; p.EvalFeatures = true;
p.BreakOnODEUnstable = false; p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0;
p.RunSlackSegments = 'AllPar'; p.MaxRunTime = 600; p.OptimizeOn = 'Feats';
p.ghostSave = ''; p.ghostLoad = ''; p.SaveBest = 0; p.mods = {}; p.g = ones(1,0);
fn = p.fn;
iR = find(strcmp(cellfun(base, fn, 'uni', 0), 'rsT20'));
assert(numel(iR) == 1, 'expected one rsT20 term in %s', FROM);
fprintf('replacing  %s\n     with  %s\n', fn{iR}, termR);
fn{iR} = termR;
iG = find(cellfun(@(e) ~contains(e, ',') && ismember(base(e), GUARD), fn));
fnRaw = fn;
for i = iG, tk = split(fn{i}, '|'); tk = tk(~startsWith(tk, 'fl:')); fnRaw{i} = strjoin(tk, '|'); end

%% ---- guards at GUARD_X of this point's raw cost under the new reader -----
p.fn = fnRaw; p = getParams(p, p.g, true, false);
[~, ~, ctRaw, NM] = evaluateBakersExp(1, p, true, false);
if any(ctRaw >= 100), error('BuildRobustSeed:coverage', 'unscorable term: %s', strjoin(NM(ctRaw >= 100), ', ')); end
w = @(e) str2double(extractAfter(e, find(e == '|', 1, 'last')));
for i = iG
    tk = split(fn{i}, '|'); tk = tk(~startsWith(tk, 'fl:'));
    tk = [tk(1:end-1); {sprintf('fl:%.6g', GUARD_X * ctRaw(i) / w(fn{i}))}; tk(end)];
    fn{i} = strjoin(tk, '|');
end
p.fn = fn; p = getParams(p, p.g, true, false);

%% ---- the seed and NNEIGH neighbours 1e-6 away ------------------------------
t0 = tic; [c1, ~, ct1, NM] = evaluateBakersExp(1, p, true, false); w1 = toc(t0);
rng(7); cN = nan(NNEIGH, 1);
for n = 1:NNEIGH
    pp = p;
    for i = 1:numel(pp.leanPool)
        nm = pp.leanPool{i};
        if ~contains(nm, '__') && isfield(pp, nm) && isscalar(pp.(nm)) && pp.(nm) ~= 0
            pp.(nm) = pp.(nm) * (1 + 1e-6*randn);
        end
    end
    pp = getParams(pp, pp.g, true, false);
    cN(n) = evaluateBakersExp(1, pp, true, false);
end
lab = cellfun(base, NM, 'uni', 0);
C = [c1; cN];
fprintf('\n=== %s = %.5f  (%.0f s)  | %d neighbours 1e-6 away: %s  -> SD %.2f %% ===\n', SEEDTAG, c1, w1, NNEIGH, ...
        mat2str(round(cN', 5)), 100*std(C)/mean(C));
[~, o] = sort(ct1, 'descend');
for i = o(:)'
    if ct1(i) <= 1e-9, continue; end
    tag = ''; if ismember(lab{i}, ATPC), tag = '   ATP'; elseif ismember(lab{i}, GUARD), tag = '   guard'; end
    fprintf('    %-22s %8.4f  (%4.1f %%)%s\n', lab{i}, ct1(i), 100*ct1(i)/c1, tag);
end
fprintf('  ATP-carrying share: %.1f %%\n', 100*sum(ct1(ismember(lab, ATPC)))/c1);

%% ---- write + round-trip ---------------------------------------------------
out = fullfile(root, 'params', [SEEDTAG '.m']);
writeParamsToMFile(out, p, [], sprintf(['%s with the robust reader (FeatureReadGrid data, SkipOscillationFit), ' ...
    'dense output inside the cycles only, rsT20 -> rsTfitE[1:4]; guards %.1fx. cost %.6f'], FROM, GUARD_X, c1));
q = getParams(); s2 = loadParams(out); f2 = fieldnames(s2);
for i = 1:numel(f2), q.(f2{i}) = s2.(f2{i}); end
for fq = fieldnames(NUM)'
    assert(isfield(q, fq{1}) && isequal(q.(fq{1}), NUM.(fq{1})), 'setting %s did not round-trip', fq{1});
end
assert(isequal(q.fn, fn), 'fn did not round-trip');
q = getParams(loadParams(SEEDTAG), [], true, false);         % the way a campaign loads it
c3 = evaluateBakersExp(1, q, true, false);
fprintf('written %s | reloaded cost %.6f (in memory %.6f, diff %.1e)\n', out, c3, c1, c3 - c1);
seed = struct('snap', SEEDTAG, 'from', FROM, 'cost', c1, 'neighbours', cN, 'featCost', ct1(:)', ...
              'featNames', {NM(:)'}, 'fn', {fn}, 'termR', termR, 'when', datestr(now, 'yyyy-mm-dd HH:MM'));
save(fullfile(res, [SEEDTAG '_build.mat']), 'seed', 'ctRaw');
fprintf('DONE\n');
