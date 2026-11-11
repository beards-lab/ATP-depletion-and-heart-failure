%% BuildYieldSeeds.m -- the controlled pair of re-fits for the overstretched-heads question.
%
% Both seeds start at refit9_opt with ONE objective change: the re-stretch
% timing trio (peak1_dSL, vall_t, restretchSlopeStart) is scored on the plain
% chi-2 floor nEl*s^2/3 of its own pooled between-prep spread -- the rule every
% ATP-carrying term already uses -- instead of as a GUARD floored at 1.2x the
% seed's residual (BuildRefit9Seed), which froze the late yield in place.
%
%   yieldA_seed  no gate (control): can the existing parameters move the yield?
%   yieldB_seed  UseYieldGate on the stretched side: p1 rebinds at zero strain
%                (YieldHop1 = 1, weak binding in rapid equilibrium), p2 ruptures
%                to PT. Rates FIXED (2000 / 3000 s^-1, i.e. a 0.3-0.5 ms dwell
%                past the limit); s_yield joins the pool, so the fit chooses the
%                elastic limit and cannot switch the premise off by shrinking k.
%
%   matlab -batch "run('Analyses/OverstretchedHeads/BuildYieldSeeds.m')"
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));

S_POOL = struct('peak1_dSL', 0.191, 'vall_t', 0.194, 'restretchSlopeStart', 0.197);   % DataPooling report, s used
GATE = struct('UseYieldGate', 1, 's_yield', 0.012, 'w_yield', 0.0015, 'k_yield1', 2000, ...
              'YieldHop1', 1, 'k_yield2', 3000);

p = getParams(loadParams('refit9_opt'), [], true, false);
fn = p.fn; base = @(e) regexprep(regexprep(e, '\|.*$', ''), '\[.*$', '');
for q = fieldnames(S_POOL)'
    i = find(cellfun(@(e) ~contains(e, ',') && strcmp(base(e), q{1}), fn));
    assert(numel(i) == 1, 'expected one fn entry for %s', q{1});
    tk = split(fn{i}, '|'); j = find(startsWith(tk, 'fl:'));
    old = tk{j}; nEl = 5; tk{j} = sprintf('fl:%.6g', nEl*S_POOL.(q{1})^2/3);
    fn{i} = strjoin(tk, '|');
    fprintf('%-20s %s -> %s\n', q{1}, old, tk{j});
end
p.fn = fn;

pA = p;
writeParamsToMFile(fullfile(root, 'params', 'yieldA_seed.m'), pA, [], ...
    'refit9_opt; re-stretch timing trio floored at nEl*s^2/3 (pooled s) instead of guarded. Control: no yield gate. Analyses/OverstretchedHeads');
pB = p;
for q = fieldnames(GATE)', pB.(q{1}) = GATE.(q{1}); end
pB.leanPool = [pB.leanPool, {'s_yield'}];
writeParamsToMFile(fullfile(root, 'params', 'yieldB_seed.m'), pB, [], ...
    'refit9_opt + stretched-side yield gate (p1 rebinds, p2 -> PT, rates fixed, s_yield free); trio floored at nEl*s^2/3. Analyses/OverstretchedHeads');

for t = {'yieldA_seed', 'yieldB_seed'}
    s = loadParams(t{1});
    fprintf('%s: UseYieldGate %d | fn trio %s\n', t{1}, s.UseYieldGate, strjoin(s.fn(contains(s.fn, {'peak1_dSL', 'vall_t|', 'restretchSlopeStart'})), '  '));
end
