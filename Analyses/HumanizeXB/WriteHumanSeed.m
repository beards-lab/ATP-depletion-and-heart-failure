% WriteHumanSeed.m -- write params/human_v0_seed.m: the incumbent (robust1_seed)
% with the literature-ratio human (beta-MHC) rate scaling of conclusions.md §4.1.
% No simulation. Re-run after changing the factors below.
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
S = struct('k2',0.1, 'k_2',0.1, 'kah',0.1, 'kamh',0.1, 'k1',0.3, 'k_1',0.3, 'ka',0.1, 'kd',0.1);   % = set v0b_kd0.1 of results/explore2.txt
p = applyScaling(loadParams('robust1_seed'), S);
p.SimTitle = 'human v0 seed';
p.mods = {}; p.g = [];
f = fieldnames(S); txt = '';
for i = 1:numel(f), txt = [txt sprintf('%s x%.3g ', f{i}, S.(f{i}))]; end %#ok<AGROW>
writeParamsToMFile(fullfile(root, 'params', 'human_v0_seed.m'), p, '', ...
    ['human (beta-MHC) v0 seed = robust1_seed with ' strtrim(txt) '; literature ratios, NOT a fit; see Analyses/HumanizeXB/conclusions.md']);
q = getParams(loadParams('human_v0_seed'), [], true, false);
fprintf('human_v0_seed: ka %.1f kd %.1f k1 %.1f k_1 %.2f k2 %.2f k_2 %.3f kah %.1f kamh %.2f (ksr0 %.0f unchanged)\n', ...
    q.ka, q.kd, q.k1, q.k_1, q.k2, q.k_2, q.kah, q.kamh, q.ksr0);
