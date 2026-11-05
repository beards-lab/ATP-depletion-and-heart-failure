function [K, S] = adpCases()
%ADPCASES  Cases of the ADP study. pH 7.2 (unchanged) / 6.9 (mildly acidic) / 6.6 (heavily acidic);
%   ATP 8/6/4/2 mM; ADP from adpScenario (CK equilibrium) on the 'mild' and 'severe' PCr paths.
%   D8(pH) = control ADP at 8 mM and that pH. Every case sets the Ca side (TW_ATP: ATP, ADP, pH) and
%   the XB side (MgATP, MgADP; both ATP gates + ADP competition in the R2 gate) separately:
%     ref    (8, D8 | 8, D8)                         per pH
%     atp    (A, D8 | A, D8)  ATP alone, ADP held     per pH, A < 8
%     full   (A, D  | A, D )  ATP + ADP               per pH, A < 8, path
%     adpCa  (A, D  | A, D8)  ADP on the Ca side only per pH, A < 8, path
%   so  ATP effect = atp - ref,  Ca-side ADP = adpCa - atp,  XB-side ADP = full - adpCa.
%   S lists the standalone T-World conditions (A, ADP, pH) the cases need.
pHs = [7.2 6.9 6.6];  As = [6 4 2];  paths = {'mild', 'severe'};
K = struct('label', {}, 'kind', {}, 'path', {}, 'pH', {}, 'A', {}, 'DCa', {}, 'DXB', {});
add = @(K, kind, path, p, a, dca, dxb) [K, struct('label', sprintf('ad_%s_%s_A%g_pH%.1f', kind, path, a, p), ...
    'kind', kind, 'path', path, 'pH', p, 'A', a, 'DCa', dca, 'DXB', dxb)];
for p = pHs
    D8 = adpScenario(8, p, 'control');
    K = add(K, 'ref', 'none', p, 8, D8, D8);
    for a = As
        K = add(K, 'atp', 'none', p, a, D8, D8);
        for q = paths
            D = adpScenario(a, p, q{1});
            K = add(K, 'full', q{1}, p, a, D, D);
            K = add(K, 'adpCa', q{1}, p, a, D, D8);
        end
    end
end
S = unique([[K.A]', [K.DCa]', [K.pH]'], 'rows', 'stable');
end
