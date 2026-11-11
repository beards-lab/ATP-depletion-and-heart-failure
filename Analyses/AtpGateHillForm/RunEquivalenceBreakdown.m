% RunEquivalenceBreakdown.m -- post-processes results/equivalence.mat (no simulation).
%
% RunEquivalenceDiagnostic established that (i) the pipeline is bit-reproducible
% (two runs of the same parameters agree exactly) and (ii) the converted gates
% reproduce every rate factor to 3e-16 -- yet the 2 mM cost still moves by ~0.016.
% So the model is amplifying a machine-precision input difference. This script says
% by how much, and in which features, so the size of that noise floor is on record.
%
% Run: matlab -batch "run('Analyses/AtpGateHillForm/RunEquivalenceBreakdown.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root));
S = load(fullfile(here,'results','equivalence.mat'));

pairs = {'fA','fB','power#1 vs power#2 (control)'; ...
         'fA','fC','power#1 vs hill, all gates converted'; ...
         'fA','fD','power#1 vs hill, |p|<1e-3 dropped'};

for ip = 1:size(pairs,1)
    fa = S.(pairs{ip,1}); fb = S.(pairs{ip,2});
    ca = evalFeatureCost(S.fd, fa, S.FN, 2);
    cb = evalFeatureCost(S.fd, fb, S.FN, 2);
    d  = cb - ca;
    fprintf('\n=== %s ===\n', pairs{ip,3});
    fprintf('total cost %.6f -> %.6f   (delta %+.4f)\n', sum(ca), sum(cb), sum(d));
    [~, od] = sort(abs(d), 'descend');
    shown = 0;
    for k = od(:)'
        if abs(d(k)) < 1e-6 || shown >= 6, break; end
        fprintf('   %-26s %9.4f -> %9.4f  (%+.4f)\n', S.FN{k}, ca(k), cb(k), d(k));
        shown = shown + 1;
    end
    if shown == 0, fprintf('   (no feature cost moved by more than 1e-6)\n'); end
end

% absolute, not relative: is ovrsht_dy simply a near-zero quantity?
fprintf('\n=== the amplifying features, in absolute units ===\n');
look = {'ovrsht_dy','vall2_dy','rsK','ktr','steady','peak1_y','vall_y','peak2','A'};
fprintf('%-12s %-34s %-34s %12s\n','feature','power#1','hill (all gates)','max |diff|');
for i = 1:numel(look)
    f = look{i};
    if ~isfield(S.fA,f) || ~isfield(S.fC,f), continue; end
    x = double(S.fA.(f)); y = double(S.fC.(f));
    g = @(v) strtrim(num2str(v(:)','%9.4g'));
    fprintf('%-12s %-34s %-34s %12.3g\n', f, g(x), g(y), max(abs(x(:)-y(:))));
end

fprintf(['\nReading: a 3e-16 change in the rate constants is the ONLY input difference, and the\n' ...
         'pipeline is deterministic, so anything above reflects discrete decisions inside the\n' ...
         'solve/extract chain (event crossings, threshold detections) flipping side. That is the\n' ...
         'numerical noise floor of the 2 mM feature cost.\n']);
fprintf('\nEQUIVALENCE BREAKDOWN DONE\n');
