%% EvaluateYieldCampaigns.m -- score the controlled pair on one instrument.
%
% For refit9_opt (incumbent), yieldA_opt (control re-fit, no gate) and
% yieldB_opt (re-fit with the stretched-side yield gate), all on the SAME
% objective (yieldA_seed's fn: trio floored at nEl*s^2/3):
%   1. total, per-feature breakdown, restretch features vs data, trio raw chi-2
%   2. the overstretched-head probe (ProbeRestretchDwell) -- the physiology
%   3. the held-out staircase (RunStairsCheck)
% and on yieldB_opt only, the low-ATP question: D* (UseD0State) on the p2
% rupture flux at 8 and 2 mM, with physiological and phenomenological K_D0T.
%   matlab -batch "run('Analyses/OverstretchedHeads/EvaluateYieldCampaigns.m')"
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
res = fullfile(here, 'results');
if ~exist('SNAPS', 'var'), SNAPS = {'refit9_opt', 'yieldA_opt', 'yieldB_opt', 'yieldC_opt'}; end
fnCommon = loadParams('yieldA_seed'); fnCommon = fnCommon.fn;
FEAT = {'peak1_y', 'peak1_dSL', 'vall_y', 'vall_t', 'peak2', 'vall2_dy', 'restretchSlopeStart', 'steady', 'A', 'rsK80', 'rsT20', 'ktrfit', 'ktrP_k80'};
if isempty(gcp('nocreate')), refreshPool(5); end
featLabOH = cellfun(@(s) regexprep(regexprep(s, '\|.*$', ''), '\[.*$', ''), fnCommon, 'uni', 0);

if ~exist('SKIP_OBJ', 'var'), SKIP_OBJ = false; end   % section 1 already saved in results/campaign_eval.mat
EVOH = struct();
for kOH = 1:numel(SNAPS) * ~SKIP_OBJ
    P = getParams(loadParams(SNAPS{kOH}), [], true, false);
    P.fn = fnCommon; P.RunSlackSegments = 'AllPar'; P.MaxRunTime = 900; P.PlotEachSeparately = 0; P.PlotFeatureFitting = 0;
    P.ghostSave = ''; P.ghostLoad = ''; P.SaveBest = 0; P.mods = {}; P.g = ones(1,0);
    P = getParams(P, P.g, true, false);
    params0 = P; evalc('LoadData; RunBakersExp;');
    fc = evalFeatureCost(features_data, features_model, P.fn);
    EVOH(kOH).snap = SNAPS{kOH}; EVOH(kOH).cost = sum(fc); EVOH(kOH).fc = fc(:)'; EVOH(kOH).fm = features_model; EVOH(kOH).fd = features_data;
    fprintf('\n===== %s : common-objective cost %.4f =====\n', SNAPS{kOH}, sum(fc));
    [~, oOH] = sort(fc(:), 'descend');
    for i = oOH(1:8)', if fc(i) > 1e-4, fprintf('   %-22s %.4f\n', featLabOH{i}, fc(i)); end, end
    for q = FEAT
        if isfield(features_model, q{1})
            fprintf('   %-20s model %-40s data %s\n', q{1}, mat2str(features_model.(q{1})(:)', 4), mat2str(features_data.(q{1})(:)', 4));
        end
    end
    if P.UseYieldGate, fprintf('   gate: s_yield %.2f nm, w %.2f nm, k1 %.0f (hop %.2f), k2 %.0f\n', 1e3*P.s_yield, 1e3*P.w_yield, P.k_yield1, P.YieldHop1, P.k_yield2); end
end
if ~SKIP_OBJ, save(fullfile(res, 'campaign_eval.mat'), 'EVOH', 'featLabOH', 'SNAPS'); end

%% physiology probe + staircase
for kOH = 1:numel(SNAPS)
    PARAMSET = SNAPS{kOH}; TAG = SNAPS{kOH}; OVR = struct(); RUNKTR = true;
    run(fullfile(here, 'ProbeRestretchDwell.m'));
end
run(fullfile(here, 'RunStairsCheck.m'));

%% low ATP on the gated fit: does D* carry a restretch-specific [ATP] effect?
if any(strcmp(SNAPS, 'yieldB_opt'))
    PB = getParams(loadParams('yieldB_opt'), [], true, false);
    PB.fn = fnCommon; PB.RunSlackSegments = 'AllPar'; PB.MaxRunTime = 900; PB.PlotEachSeparately = 0; PB.PlotFeatureFitting = 0;
    PB.ghostSave = ''; PB.ghostLoad = ''; PB.SaveBest = 0; PB.mods = {}; PB.g = ones(1,0);
    VOH = {'no D*', 8, struct(); 'no D*', 2, struct(); ...
         'D* K 0.15', 8, struct('UseD0State', 1, 'K_D0T', 0.15); 'D* K 0.15', 2, struct('UseD0State', 1, 'K_D0T', 0.15); ...
         'D* K 5', 8, struct('UseD0State', 1, 'K_D0T', 5); 'D* K 5', 2, struct('UseD0State', 1, 'K_D0T', 5)};
    AOH = struct();
    for vOH = 1:size(VOH, 1)
        P = PB; ovOH = VOH{vOH, 3}; fOH = fieldnames(ovOH); for j = 1:numel(fOH), P.(fOH{j}) = ovOH.(fOH{j}); end
        P.MgATP = VOH{vOH, 2}; P = getParams(P, P.g, true, false);
        params0 = P; evalc('LoadData; RunBakersExp;');
        AOH(vOH).name = VOH{vOH,1}; AOH(vOH).MgATP = VOH{vOH,2}; AOH(vOH).fm = features_model;
    end
    fprintf('\n===== yieldB_opt: model 2 mM / 8 mM ratios (data: rsK80 x0.38, ktr x0.53-0.57, peak1_dSL x1.17, vall_t x0.79) =====\n');
    for vOH = 1:2:size(VOH, 1)
        r = @(q) mean(AOH(vOH+1).fm.(q)(1:min(4,end)) ./ AOH(vOH).fm.(q)(1:min(4,end)));
        fprintf('   %-10s rsK80 x%.3f | ktrfit x%.3f | peak1_dSL x%.3f | vall_t x%.3f | steady x%.3f\n', ...
            VOH{vOH,1}, r('rsK80'), r('ktrfit'), r('peak1_dSL'), r('vall_t'), r('steady'));
    end
    save(fullfile(res, 'campaign_atp.mat'), 'AOH', 'VOH');
end
