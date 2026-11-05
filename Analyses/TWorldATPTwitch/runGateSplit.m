function runGateSplit(idx)
%RUNGATESPLIT  Which XB ATP gate carries the XB-side effect: Ca side at 8 mM, XB at MgATP 2 mM with
%   only UseAtpDetachR2 (ATP-binding detachment, K_T2) or only UseAtpKmsrd (SRX-ADP return, K_srx)
%   on, at pH 7.2 (low activation) and 6.6 (high activation). Same protocol as runCoupled.
here = fileparts(mfilename('fullpath'));
twPaths('twatp');
[ATP, ADP] = atpLadder();
G = struct('gate', {'R2', 'Kmsrd', 'R2', 'Kmsrd'}, 'pH', {7.2, 7.2, 6.6, 6.6});
if nargin < 1, idx = 1:numel(G); end
tb = (0:0.001:0.999)';
for i = idx
    g = G(i);  f = fullfile(here, 'results', sprintf('gs_%s_pH%.1f.mat', g.gate, g.pH));
    if exist(f, 'file'), continue; end
    sa = load(fullfile(here, 'results', sprintf('sa_twatp_A8_pH%.1f.mat', g.pH)));  sa = sa.r;
    base = loadParams('human_twitch_t1_isoT');
    base.UseAtpDetachR2 = strcmp(g.gate, 'R2');  base.UseAtpKmsrd = strcmp(g.gate, 'Kmsrd');  base.MgATP = 2;
    th = struct('CaInput', 'tworld', 'CaTworld', struct('coupling', 'land', 'model', @model_TWorld_shock_ATP, ...
        'X0', sa.Xend, 'param', struct('extraParams', [0 310 g.pH 8 ADP(ATP == 8) 0])));
    [Fdev, S] = simHumanTwitch(base, th, tb);  S = rmfield(S, 'P');
    r = struct('g', g, 'tb', tb, 'Fdev', Fdev, 'S', S);  save(f, 'r');
    fprintf('gate %-5s pH %.1f: Fdev %.3f kPa, SRX min %.3f, bound max %.4f, ATP/beat %.3f\n', g.gate, g.pH, max(Fdev), ...
        min(S.SRX), max(S.bound), trapz(tb, S.hyd));
end
end
