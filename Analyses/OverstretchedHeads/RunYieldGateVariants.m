%% RunYieldGateVariants.m -- earlier detachment on the STRETCHED side only, no re-fit.
%
% ProbeRestretchDwell: at refit9_opt the heads beyond +15 nm are almost all
% p1 (weak, pre-stroke), living ~4 ms there (p2 lives 0.26 ms), and they
% exist ONLY during lengthening. UseYieldGate (default off) adds a bounded
% logistic detachment for s > s_yield, so by construction it cannot touch the
% isometric plateau or shortening. This asks, on the full objective with the
% incumbent's own numerics and nothing re-fitted:
%   1. which features move (isolation: FV / steady / post-slack ktr should not)
%   2. which direction the restretch timing trio moves -- scored RAW, because the
%      guards' deadbands would hide it (peak1_dSL, vall_t, restretchSlopeStart)
%   3. p1-only vs p2-only vs both, and the p2 rupture routed via D* (UseD0State)
%   4. with D*, the model 2 mM / 8 mM ratios -- nothing else in refit9 depends
%      on [MgATP], so any ratio here IS the D* route's contribution
%
%   matlab -batch "run('Analyses/OverstretchedHeads/RunYieldGateVariants.m')"
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
res = fullfile(here, 'results');

g = @(sy, k1, k2) struct('UseYieldGate', 1, 's_yield', sy, 'w_yield', 0.0015, 'k_yield1', k1, 'k_yield2', k2);
d = @(o, K) setfield(setfield(setfield(o, 'UseD0State', 1), 'k_D0T', 500), 'K_D0T', K); %#ok<SFLD>
if ~exist('VAR', 'var')
VAR = { 'as fitted',                         struct(),                        8
        'p1 gate 15 nm, 1000/s',             g(0.015, 1000, 0),               8
        'p1 gate 15 nm, 5000/s',             g(0.015, 5000, 0),               8
        'p1 gate 12 nm, 1000/s',             g(0.012, 1000, 0),               8
        'p1 gate 12 nm, 5000/s',             g(0.012, 5000, 0),               8
        'p1 gate 10 nm, 2000/s',             g(0.010, 2000, 0),               8
        'p2 gate 12 nm, 3000/s -> PT',       g(0.012, 0, 3000),               8
        'p1+p2 gate 12 nm -> PT',            g(0.012, 2000, 3000),            8
        'p1+p2 gate 12 nm, D* K 0.15',       d(g(0.012, 2000, 3000), 0.15),   8
        'p1+p2 gate 12 nm, D* K 0.15 @2mM',  d(g(0.012, 2000, 3000), 0.15),   2
        'p1+p2 gate 12 nm, D* K 5',          d(g(0.012, 2000, 3000), 5),      8
        'p1+p2 gate 12 nm, D* K 5 @2mM',     d(g(0.012, 2000, 3000), 5),      2
        % added after the first sweep (first run by hand, same settings):
        'p1 gate 12 nm 1000/s, A_def off',   setfield(g(0.012, 1000, 0), 'UseLoadRealign', 0), 8 %#ok<SFLD>
        'p1 HOP 12 nm, 2000/s',              setfield(g(0.012, 2000, 0), 'YieldHop1', 1),      8 }; %#ok<SFLD>
end
FEAT = {'peak1_y', 'peak1_dSL', 'vall_y', 'vall_t', 'peak2', 'steady', 'restretchSlopeStart', ...
        'vall2_dy', 'rsK80', 'rsT20', 'ktrfit', 'ktrP_k80', 'A', 'PS_restretchPeak'};
TRIO = {'peak1_dSL', 'vall_t', 'restretchSlopeStart'};

P0 = getParams(loadParams(PARAMSET), [], true, false);
P0.RunSlackSegments = 'AllPar'; P0.MaxRunTime = 900; P0.PlotEachSeparately = 0; P0.PlotFeatureFitting = 0;
P0.ghostSave = ''; P0.ghostLoad = ''; P0.SaveBest = 0; P0.mods = {}; P0.g = ones(1,0);
if isempty(gcp('nocreate')), refreshPool(5); end
lab = cellfun(@(s) regexprep(regexprep(s, '\|.*$', ''), '\[.*$', ''), P0.fn, 'uni', 0);
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
        fprintf('%-34s %g mM  FAILED: %s\n', VAR{v,1}, VAR{v,3}, err.message); continue
    end
    [fc, ~, ~] = evalFeatureCost(features_data, features_model, P.fn);
    R(v).name = VAR{v,1}; R(v).MgATP = VAR{v,3}; R(v).ovr = ov; R(v).cost = sum(fc); R(v).fc = fc(:)'; R(v).wall = toc(t0);
    fe = struct(); fd = struct();
    for q = FEAT
        fe.(q{1}) = NaN; fd.(q{1}) = NaN;
        if isfield(features_model, q{1}), fe.(q{1}) = features_model.(q{1})(:)'; end
        if isfield(features_data, q{1}),  fd.(q{1}) = features_data.(q{1})(:)'; end
    end
    R(v).feat = fe; R(v).data = fd;
    raw = zeros(1, numel(TRIO));             % raw chi-2 of the timing trio, no deadband, signed mean residual
    sgn = zeros(1, numel(TRIO));
    for j = 1:numel(TRIO)
        r = (fe.(TRIO{j}) - fd.(TRIO{j})) / mean(fd.(TRIO{j}));
        raw(j) = sum(r.^2); sgn(j) = mean(r);
    end
    R(v).trioRaw = raw; R(v).trioBias = sgn;
    top = sortrows(table(lab(:), fc(:), 'VariableNames', {'f', 'c'}), 'c', 'descend');
    fprintf(['%-34s %g mM  cost %8.4f (%3.0f s) | pk1 %s | dSL %s | vall %s | vt ms %s | pk2 %s | slope %s | rsK80 %s | ktrfit %.1f | ktrP %.1f\n' ...
             '%38s trio raw chi2 [dSL vt slope] %s  bias %s | top: %s %.3f, %s %.3f, %s %.3f\n'], ...
        VAR{v,1}, VAR{v,3}, sum(fc), R(v).wall, mat2str(round(fe.peak1_y,1)), mat2str(round(1e3*fe.peak1_dSL,1)), ...
        mat2str(round(fe.vall_y,1)), mat2str(round(1e3*fe.vall_t,1)), mat2str(round(fe.peak2,1)), mat2str(round(fe.restretchSlopeStart)), ...
        mat2str(round(fe.rsK80,1)), mean(fe.ktrfit), fe.ktrP_k80, '', mat2str(round(raw,3)), mat2str(round(sgn,3)), ...
        top.f{1}, top.c(1), top.f{2}, top.c(2), top.f{3}, top.c(3));
    save(fullfile(res, ['yieldgate_' PARAMSET '.mat']), 'R', 'VAR', 'lab');
end
