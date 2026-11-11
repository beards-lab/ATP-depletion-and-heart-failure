%% CompareOutputGrids.m -- which objective terms need dense output, and on which segments?
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..')));
res = fullfile(here, 'results');
fx = @(varargin) struct('SmoothRectifierEps', 0.01, 'ClampStrainDepExtrapolation', 1, 'odeRelTol', 1e-3, 'odeAbsTol', 1e-6, varargin{:});
if ~exist('VAR', 'var')
VAR = { 'dense everywhere', fx('odeOutputDt', 5e-4)
        'dense segs <= 5 s', fx('odeOutputDt', 5e-4, 'odeOutputMaxSegment', 5)
        'no dense',          fx('odeOutputDt', 0) };
end
P0 = getParams(loadParams(PARAMSET), [], true, false);
P0.RunSlackSegments = 'AllPar'; P0.MaxRunTime = 900; P0.PlotEachSeparately = 0; P0.PlotFeatureFitting = 0;
P0.ghostSave = ''; P0.ghostLoad = ''; P0.SaveBest = 0; P0.mods = {}; P0.g = ones(1,0);
if isempty(gcp('nocreate')), parpool(5); end
FC = []; FM = {};
for v = 1:size(VAR, 1)
    P = P0; ov = VAR{v,2}; f = fieldnames(ov); for k = 1:numel(f), P.(f{k}) = ov.(f{k}); end
    P = getParams(P, P.g, true, false);
    params0 = P; evalc('LoadData; RunBakersExp;');
    [fc, ~, ~] = evalFeatureCost(features_data, features_model, P.fn);
    FC(v, :) = fc(:)'; FM{v} = features_model; %#ok<SAGROW>
end
lab = cellfun(@(s) regexprep(regexprep(s, '\|.*$', ''), '\[.*$', ''), P0.fn, 'uni', 0);
fprintf('%-22s', 'term'); fprintf('%20s', VAR{:,1}); fprintf('\n');
for j = 1:numel(lab)
    if all(abs(FC(:, j)) < 1e-9), continue; end
    fprintf('%-22s', lab{j}); fprintf('%20.5f', FC(:, j)); fprintf('\n');
end
fprintf('%-22s', 'TOTAL'); fprintf('%20.5f', sum(FC, 2)); fprintf('\n');
F = {'ktrfit', 'ktrfit_shape', 'rsK80', 'rsT20', 'ktrP_k80', 'peak1_dSL', 'vall_t', 't0_crossing', 'steady'};
for q = F
    fprintf('%-13s', q{1});
    for v = 1:numel(FM), if isfield(FM{v}, q{1}), fprintf('  %-34s', mat2str(round(FM{v}.(q{1})(:)', 4))); end, end
    fprintf('\n');
end
save(fullfile(res, ['outputgrids_' PARAMSET '.mat']), 'FC', 'FM', 'lab', 'VAR');
