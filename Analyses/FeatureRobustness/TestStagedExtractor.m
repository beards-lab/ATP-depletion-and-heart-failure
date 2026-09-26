%% TestStagedExtractor.m -- a staged extractSlackAttributes must be bit-identical with its switches off.
% STAGE (a folder OUTSIDE the repo, so no genpath can shadow Model/) holds the
% candidate. Checks:
%   1. data path, three 8 mM preparations: every field of Model/'s output equal
%      (isequaln) in the staged output; only new fields may be added
%   2. model path, one real slack run of PARAMSET (dense in cycles): the same,
%      with out.params carrying no new switch
%   3. model path with FeatureReadGrid = 'data' + SkipOscillationFit: values next
%      to robustSlackFeatures' offline 4 kHz read of the same trace
if ~exist('STAGE', 'var'), STAGE = 'C:\Users\fjezek\AppData\Local\Temp\featrob_staging'; end
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
here = fileparts(mfilename('fullpath')); root = fullfile(here, '..', '..');
addpath(genpath(root));
useStage = @(on) stageOn(STAGE, on);
NEWF = {'rsT20e', 'rsTfitE', 'tOn20'};

%% 1. data path
DAYS = {'03_27', '04_03', '04_10'}; okAll = true;
for d = 1:3
    S = load(fullfile(root, 'data', sprintf('protocol_%s_2026_8mM_slack.mat', DAYS{d})));
    a = {S.datatable(:,1), S.datatable(:,3), S.datatable(:,2), S.velocitytable, [], [], false, true};
    useStage(false); f0 = extractSlackAttributes(a{:});
    useStage(true);  f1 = extractSlackAttributes(a{:});
    [ok, msg] = sameFields(f0, f1, NEWF); okAll = okAll && ok;
    fprintf('data %s: %s\n', DAYS{d}, msg);
end

%% 2. model path, switches off
useStage(false);
P = getParams(loadParams(PARAMSET), [], true, false);
P.SmoothRectifierEps = 0.01; P.ClampStrainDepExtrapolation = 1; P.odeRelTol = 1e-3; P.odeAbsTol = 1e-6;
P.odeOutputDt = 5e-4; P.odeOutputMaxSegment = 5;
P.RunSlackSegments = 'AllPar'; P.PlotFeatureFitting = 0; P.PlotEachSeparately = 0; P.mods = {}; P.g = ones(1,0);
P = getParams(P, P.g, true, false);
delete(gcp('nocreate'));
[~, os, fmRun] = runSlackExperiment(P);
D = load(fullfile(root, 'data', P.velocitytableonfile)); vt = D.velocitytable; vt(1,1) = -20;
k = makeMonotonous(os.t); tF = os.t(k); fF = os.Force(k); sF = os.SL(k);
useStage(false); g0 = extractSlackAttributes(tF, fF, sF, vt, [], os, 0, [], os.datatable);
useStage(true);  g1 = extractSlackAttributes(tF, fF, sF, vt, [], os, 0, [], os.datatable);
[ok, msg] = sameFields(g0, g1, NEWF); okAll = okAll && ok;
fprintf('model, switches off: %s\n', msg);
[ok, msg] = sameFields(fmRun, g0, {}); fprintf('model, in-run vs re-extracted (Model/): %s\n', msg);

%% 3. model path, switches on
o2 = os; o2.params.FeatureReadGrid = 'data'; o2.params.SkipOscillationFit = 1;
t0 = tic; g2 = extractSlackAttributes(tF, fF, sF, vt, [], o2, 0, [], os.datatable); w2 = toc(t0);
useStage(true); t0 = tic; extractSlackAttributes(tF, fF, sF, vt, [], os, 0, [], os.datatable); w1 = toc(t0);
[tu, iu] = unique(os.t(:)); gg = (ceil(tu(1)/2.5e-4):floor(tu(end)/2.5e-4))' * 2.5e-4;
LSE = os.SL(:) - os.LXB(:);
rb = robustSlackFeatures(gg, interp1(tu, os.Force(iu)', gg), interp1(tu, os.SL(iu)', gg), vt, false, tu, LSE(iu));
fprintf('\nswitches on: extraction %.1f s (off: %.1f s); ktr2_k all NaN: %d\n', w2, w1, all(isnan(g2.ktr2_k)));
cmp = {'rsK80','rsK80'; 'rsT20','rsT20'; 'rsTfitE','rsTfitE'; 'rsT20e','rsT20e'; 'ktrfit','ktrfit'; ...
       'tOn20','tOn20'; 't0_crossing','t0x'; 'restretchSlopeStart','restretchSlopeStart'; 'peak1_dSL','peak1_dSL'};
for j = 1:size(cmp, 1)
    a1 = g2.(cmp{j,1})(:)'; b1 = rb.(cmp{j,2})(:)';
    fprintf('  %-20s staged-on %-44s offline %-44s max rel diff %.1e\n', cmp{j,1}, mat2str(round(a1, 5)), ...
            mat2str(round(b1, 5)), max(abs(a1 - b1) ./ abs(b1)));
end

%% 4. the parallel path's merge: init hold + cycles (what a campaign on a pool extracts)
Pi = rmfield(P, 'PU0'); Pi.Velocity = 0;
[~, oi] = evaluateModel(str2func(P.modelFcn), [-10 vt(3,1)], Pi);
osc = rmfield(os, 'datatable'); fi = intersect(fieldnames(oi), fieldnames(osc), 'stable');
oi2 = rmfield(oi, setdiff(fieldnames(oi), fi)); osc = rmfield(osc, setdiff(fieldnames(osc), fi));
osc = orderfields(osc, oi2); op = mergeOutStructs([oi2, osc]); op.datatable = os.datatable; op.params = os.params;
k = makeMonotonous(op.t); tP = op.t(k); fP = op.Force(k); sP = op.SL(k);
useStage(false); h0 = extractSlackAttributes(tP, fP, sP, vt, [], op, 0, [], op.datatable);
useStage(true);  h1 = extractSlackAttributes(tP, fP, sP, vt, [], op, 0, [], op.datatable);
[ok, msg] = sameFields(h0, h1, NEWF); okAll = okAll && ok;
fprintf('\nparallel-path merge, switches off: %s\n', msg);
op2 = op; op2.params.FeatureReadGrid = 'data'; op2.params.SkipOscillationFit = 1;
h2 = extractSlackAttributes(tP, fP, sP, vt, [], op2, 0, [], op.datatable);
fprintf('  t0_crossing  serial path %s | parallel path, off %s | parallel path, on %s\n', ...
        mat2str(round(g0.t0_crossing, 4)), mat2str(round(h0.t0_crossing, 4)), mat2str(round(h2.t0_crossing, 4)));
dd = {};
for q = {'rsK80','rsT20','rsTfitE','ktrfit','tOn20','t0_crossing','restretchSlopeStart','peak1_dSL','vall_t','steady','A','ktr'}
    if ~isequaln(h2.(q{1}), g2.(q{1})), dd{end+1} = sprintf('%s (max rel %.1e)', q{1}, max(abs(h2.(q{1}) - g2.(q{1})) ./ abs(g2.(q{1})))); end %#ok<AGROW>
end
fprintf('  switches on, parallel vs serial path, features that differ: %s\n', strjoin([{'none'}, dd], ', '));

%% 5. ktr runner
Pk = P; Pk.UsePooledTarget = true;
useStage(false); clear runKtrExperiment; [~, ~, k0] = runKtrExperiment(Pk);
useStage(true);  clear runKtrExperiment; assert(startsWith(which('runKtrExperiment'), STAGE));
[~, ~, k1] = runKtrExperiment(Pk);
[ok, msg] = sameFields(k0, k1, {}); okAll = okAll && ok;
fprintf('\nktr runner, switch off: %s\n', msg);
Pk.FeatureReadGrid = 'data'; [~, ~, k2] = runKtrExperiment(Pk);
fprintf('  ktrP_k80: off %.4f | on %.4f | ktrP_k63: off %.4f | on %.4f\n', k0.ktrP_k80, k2.ktrP_k80, k0.ktrP_k63, k2.ktrP_k63);

useStage(false);
fprintf('\nBIT-IDENTICAL WITH SWITCHES OFF: %d\n', okAll);

function stageOn(stage, on)
    if on, addpath(stage, '-begin'); elseif contains(path, stage), rmpath(stage); end
    clear extractSlackAttributes runKtrExperiment
    w = which('extractSlackAttributes');
    assert(on == startsWith(w, stage), 'wrong extractSlackAttributes on the path: %s', w);
end

function [ok, msg] = sameFields(a, b, allowedNew)
    fa = fieldnames(a); fb = fieldnames(b);
    miss = setdiff(fa, fb); extra = setdiff(fb, fa);
    diffF = {};
    for i = 1:numel(fa)
        if isfield(b, fa{i}) && ~isequaln(a.(fa{i}), b.(fa{i})), diffF{end+1} = fa{i}; end %#ok<AGROW>
    end
    ok = isempty(miss) && isempty(diffF) && all(ismember(extra, allowedNew));
    msg = sprintf('%d fields equal, %d differ %s, %d missing, new: %s', numel(fa) - numel(diffF) - numel(miss), ...
                  numel(diffF), strjoin(diffF, ','), numel(miss), strjoin(extra', ','));
end
