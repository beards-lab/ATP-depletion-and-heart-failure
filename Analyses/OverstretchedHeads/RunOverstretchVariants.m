%% RunOverstretchVariants.m -- is holding heads beyond ~15 nm central to the fit?
%
% Full objective (evaluateBakersExp, AllPar) on PARAMSET with its own numerics,
% NO re-fit: each variant removes overstretched heads a different way and the
% change in cost and in the ATP-carrying features says how much the fit leans
% on them. The apo-pool (D*) variants are also run at 2 mM MgATP; in the main
% tree nothing else in these runs depends on [ATP], so the 2/8 mM feature
% ratio isolates what the D* route alone would contribute to the ATP effect.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..')));
res = fullfile(here, 'results');

rip = @(K) struct('UseD0State', 1, 'k2rip', 200, 'alphaRip', 100, 'dr3', -0.015, 'k_D0T', 500, 'K_D0T', K);   % 200/s at 15 nm, e-fold 10 nm (the ApoPoolDetachment slope)
if ~exist('VAR', 'var')
VAR = { 'as fitted',                                struct(),                                  8
        'fast detachment beyond +-30 nm',           struct('UseStrictDetachmentAt', 0.030),    8
        'fast detachment beyond +-20 nm',           struct('UseStrictDetachmentAt', 0.020),    8
        'fast detachment beyond +-15 nm',           struct('UseStrictDetachmentAt', 0.015),    8
        'D* rip >15 nm, apo pool, K_D0T 0.15 mM',   rip(0.15),                                 8
        'D* rip >15 nm, apo pool, K_D0T 2 mM',      rip(2),                                    8
        'as fitted @ 2 mM',                         struct(),                                  2
        'D* K_D0T 0.15 mM @ 2 mM',                  rip(0.15),                                 2
        'D* K_D0T 2 mM @ 2 mM',                     rip(2),                                    2 };
end
FEAT = {'ktrfit', 'rsK80', 'rsT20', 'peak1_y', 'vall_y', 'ktrP_k80', 'steady'};

P0 = getParams(loadParams(PARAMSET), [], true, false);
P0.RunSlackSegments = 'AllPar'; P0.MaxRunTime = 900; P0.PlotEachSeparately = 0; P0.PlotFeatureFitting = 0;
P0.ghostSave = ''; P0.ghostLoad = ''; P0.SaveBest = 0; P0.mods = {}; P0.g = ones(1,0);
if isempty(gcp('nocreate')), parpool(5); end
R = struct();
for v = 1:size(VAR, 1)
    P = P0; ov = VAR{v, 2}; f = fieldnames(ov);
    for k = 1:numel(f), P.(f{k}) = ov.(f{k}); end
    P.MgATP = VAR{v, 3};
    P = getParams(P, P.g, true, false);
    t0 = tic;
    try
        params0 = P; evalc('LoadData; RunBakersExp;');
    catch err
        fprintf('%-42s %g mM  FAILED: %s\n', VAR{v,1}, VAR{v,3}, err.message); continue
    end
    [fc, ~, ~] = evalFeatureCost(features_data, features_model, P.fn);
    R(v).name = VAR{v,1}; R(v).MgATP = VAR{v,3}; R(v).cost = sum(fc); R(v).fc = fc(:)'; R(v).wall = toc(t0);
    fe = struct(); for q = FEAT, if isfield(features_model, q{1}), fe.(q{1}) = features_model.(q{1})(:)'; else, fe.(q{1}) = NaN; end, end
    R(v).feat = fe;
    fprintf('%-42s %g mM  cost %8.4f  (%4.0f s) | rsK80 %s | rsT20 ms %s | ktrfit %s | ktrP %.1f\n', VAR{v,1}, VAR{v,3}, sum(fc), R(v).wall, ...
        mat2str(round(fe.rsK80(1:min(4,end)), 1)), mat2str(round(1e3*fe.rsT20(1:min(4,end)), 2)), mat2str(round(fe.ktrfit, 1)), fe.ktrP_k80);
    save(fullfile(res, ['variants_' PARAMSET '.mat']), 'R', 'VAR');
end
lab = cellfun(@(s) regexprep(regexprep(s, '\|.*$', ''), '\[.*$', ''), P0.fn, 'uni', 0);
save(fullfile(res, ['variants_' PARAMSET '.mat']), 'R', 'VAR', 'lab');
