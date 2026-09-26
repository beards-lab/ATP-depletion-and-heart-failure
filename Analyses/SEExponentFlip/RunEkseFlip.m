%% RunEkseFlip.m -- does flipping the series-elastic exponent above 1 land in a better basin?
%
% The fitted series element is Force = (kSE*LSE)^ekSE with ekSE = 0.845 < 1:
% infinitely stiff at zero load and SOFTENING with load, the opposite of a
% tendon-like toe. Each variant sets ekSE and rescales kSE so that ONE property
% at the isometric operating force F_REF is unchanged:
%   'stiff'  the SE stiffness dF/dLSE = ekSE*F/LSE (what the transients see)
%   'ext'    the SE extension LSE at F_REF (what sets the CE length)
% Full objective, no re-fit: this only ranks seeds for a proper re-fit.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
if ~exist('F_REF', 'var'), F_REF = 71.45; end          % kPa, model steady at SL 2.2
if ~exist('EK', 'var'), EK = [1.0 1.2 1.5]; end
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..')));
res = fullfile(here, 'results');

P0 = getParams(loadParams(PARAMSET), [], true, false);
P0.RunSlackSegments = 'AllPar'; P0.MaxRunTime = 900; P0.PlotEachSeparately = 0; P0.PlotFeatureFitting = 0;
P0.ghostSave = ''; P0.ghostLoad = ''; P0.SaveBest = 0; P0.mods = {}; P0.g = ones(1,0);
e0 = P0.ekSE; k0 = P0.kSE;
L0  = F_REF^(1/e0) / k0;               % SE extension at F_REF (um)
Kop = e0 * F_REF / L0;                 % SE stiffness at F_REF (kPa/um)
fprintf('as fitted: ekSE %.4f, kSE %.1f -> LSE(F_REF) %.2f nm, stiffness %.0f kPa/um\n', e0, k0, 1e3*L0, Kop);

V = {'as fitted', e0, k0};
for e = EK
    Ls = e * F_REF / Kop;  V(end+1, :) = {sprintf('ekSE %.1f, same stiffness', e), e, F_REF^(1/e)/Ls}; %#ok<SAGROW>
    V(end+1, :) = {sprintf('ekSE %.1f, same extension', e), e, F_REF^(1/e)/L0};                      %#ok<SAGROW>
end
if isempty(gcp('nocreate')), parpool(5); end
R = struct();
for v = 1:size(V, 1)
    P = P0; P.ekSE = V{v,2}; P.kSE = V{v,3}; P = getParams(P, P.g, true, false);
    t0 = tic; params0 = P; evalc('LoadData; RunBakersExp;');
    [fc, ~, ~] = evalFeatureCost(features_data, features_model, P.fn);
    R(v).name = V{v,1}; R(v).ekSE = V{v,2}; R(v).kSE = V{v,3}; R(v).cost = sum(fc); R(v).fc = fc(:)'; R(v).wall = toc(t0);
    R(v).rsK80 = features_model.rsK80(:)'; R(v).rsT20 = features_model.rsT20(:)'; R(v).ktrfit = features_model.ktrfit(:)';
    R(v).ktrP = features_model.ktrP_k80; R(v).steady = features_model.steady(:)';
    fprintf('%-28s ekSE %.2f kSE %8.1f  cost %.4f (%3.0f s) | rsK80 %s | rsT20 ms %s | ktrfit %s | ktrP %.1f | steady %.1f\n', V{v,1}, V{v,2}, V{v,3}, ...
        sum(fc), R(v).wall, mat2str(round(R(v).rsK80(1:4),1)), mat2str(round(1e3*R(v).rsT20(1:4),2)), mat2str(round(R(v).ktrfit,1)), R(v).ktrP, R(v).steady(1));
    save(fullfile(res, ['ekseflip_' PARAMSET '.mat']), 'R', 'V');
end
lab = cellfun(@(s) regexprep(regexprep(s, '\|.*$', ''), '\[.*$', ''), P0.fn, 'uni', 0);
save(fullfile(res, ['ekseflip_' PARAMSET '.mat']), 'R', 'V', 'lab');
