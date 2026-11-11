%% RunA1SlipVariants.m -- A1 slip-rebinding to the neighbouring actin site, no re-fit.
%
% yieldB_opt's torn weak heads rebind at ZERO strain (YieldHop1 = 1,
% YieldHopShift = 0), which throws away their extension: the valley and the
% plateau held during the rest of the ramp come out low (Analyses/OverstretchedHeads
% sec. 5, 8). The A2 state has an attachment-shift option (a strained head slips
% one actin monomer, d_actin, and rebinds). Here the same logic for A1: the yield
% flux lands YieldHopShift back along the filament, still strained, and the
% YieldHop1 split sends the rest to PD (unbinding). Full objective at 8 mM on
% the common objective of the campaigns, plus the sub-sample timing features.
%   matlab -batch "run('Analyses/OverstretchedHeads/RunA1SlipVariants.m')"
if ~exist('BASE', 'var'), BASE = 'yieldB_opt'; end
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
res = fullfile(here, 'results');
sl = @(sh, hop) struct('YieldHopShift', sh, 'YieldHop1', hop);
if ~exist('VARS', 'var')
VARS = { 'as fitted (land at 0)',             struct()
         'slip 5.5 nm, all rebind',           sl(0.0055, 1)
         'slip 5.5 nm, 3/4 rebind 1/4 PD',    sl(0.0055, 0.75)
         'slip 5.5 nm, 1/2 rebind 1/2 PD',    sl(0.0055, 0.5)
         'slip 11 nm, all rebind',            sl(0.011, 1)
         'slip 5.5 nm, s_yield 12 nm',        setfield(sl(0.0055, 1), 's_yield', 0.012) }; %#ok<SFLD>
end
fnC = loadParams('yieldA_seed'); fnC = fnC.fn;
labOH = cellfun(@(s) regexprep(regexprep(s, '\|.*$', ''), '\[.*$', ''), fnC, 'uni', 0);
FEATOH = {'peak1_y', 'vall_y', 'peak2', 'restretchSlopeStart', 'rsK80', 'rsT20', 'peak1_dSLq', 'vall_tq'};
if isempty(gcp('nocreate')), refreshPool(5); end
ROH = struct();
for vOH = 1:size(VARS, 1)
    P = getParams(loadParams(BASE), [], true, false);
    ovOH = VARS{vOH, 2}; fOH = fieldnames(ovOH); for j = 1:numel(fOH), P.(fOH{j}) = ovOH.(fOH{j}); end
    P.fn = fnC; P.RunSlackSegments = 'AllPar'; P.MaxRunTime = 900; P.PlotEachSeparately = 0; P.PlotFeatureFitting = 0;
    P.ghostSave = ''; P.ghostLoad = ''; P.SaveBest = 0; P.mods = {}; P.g = ones(1,0);
    P = getParams(P, P.g, true, false);
    tOH = tic; params0 = P; evalc('LoadData; RunBakersExp;');
    fcOH = evalFeatureCost(features_data, features_model, P.fn);
    ROH(vOH).name = VARS{vOH, 1}; ROH(vOH).cost = sum(fcOH); ROH(vOH).fc = fcOH(:)'; ROH(vOH).fm = features_model;
    [~, oOH] = sort(fcOH(:), 'descend');
    fprintf('\n== %-32s cost %.4f (%.0f s) | top: %s %.3f, %s %.3f, %s %.3f\n', VARS{vOH,1}, sum(fcOH), toc(tOH), ...
        labOH{oOH(1)}, fcOH(oOH(1)), labOH{oOH(2)}, fcOH(oOH(2)), labOH{oOH(3)}, fcOH(oOH(3)));
    for q = FEATOH
        if isfield(features_model, q{1}), fprintf('   %-20s %s\n', q{1}, mat2str(features_model.(q{1})(:)', 4)); end
    end
    save(fullfile(res, ['a1slip_' BASE '.mat']), 'ROH', 'VARS', 'labOH');
end
fprintf('\n   data vall_y %s | peak2 %s | peak1_y %s\n', mat2str(features_data.vall_y(:)', 4), mat2str(features_data.peak2(:)', 4), mat2str(features_data.peak1_y(:)', 4));
