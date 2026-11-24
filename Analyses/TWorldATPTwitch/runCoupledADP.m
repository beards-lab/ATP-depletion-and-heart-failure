function runCoupledADP(idx)
%RUNCOUPLEDADP  Coupled TW_ATP -> XB twitch for the ADP study (adpCases), as runCoupled but with ADP
%   on both sides. XB side: both MgATP gates on plus UseAdpCompR2 (MgADP competes with ATP in the R2
%   gate, K_D 0.194 mM), normalised at (8 mM, control ADP at pH 7.2) so the fitted twitch is the
%   reference. MgADP = ADP * Mg/(Kd_MgADP + Mg) with TW_ATP's Mg 0.5 mM and Kd 0.68 mM; MgATP = total
%   ATP as in runCoupled. T-World starts from runStandalone('adp') steady states.
here = fileparts(mfilename('fullpath'));
twPaths('twatp');
K = adpCases();
if nargin < 1 || isempty(idx), idx = 1:numel(K); end
fMg = 0.5/(0.68 + 0.5);
tb = (0:0.001:0.999)';
for i = idx
    k = K(i);  f = fullfile(here, 'results', [k.label '.mat']);
    if exist(f, 'file'), fprintf('%s cached\n', k.label); continue; end
    fs = fullfile(here, 'results', [saLabel(k.A, k.DCa, k.pH) '.mat']);
    t0 = tic;  while ~exist(fs, 'file') && toc(t0) < 3600, pause(10); end   % standalone may still be running
    pause(2);  sa = load(fs);  sa = sa.r;
    base = loadParams('human_twitch_t1_isoT');
    base.UseAtpDetachR2 = 1;  base.UseAtpKmsrd = 1;  base.UseAdpCompR2 = 1;
    base.MgATP = k.A;  base.MgADP = fMg*k.DXB;  base.MgADP_ref = fMg*adpScenario(8, 7.2, 'control');
    th = struct('CaInput', 'tworld', 'CaTworld', struct('coupling', 'land', 'model', @model_TWorld_shock_ATP, ...
        'X0', sa.Xend, 'param', struct('extraParams', [0 310 k.pH k.A k.DCa 0])));
    tic;
    try
        [Fdev, S] = simHumanTwitch(base, th, tb);  err = '';
    catch ME
        Fdev = [];  S = [];  err = ME.message;
    end
    r = struct('k', k, 'wall', toc, 'tb', tb, 'Fdev', Fdev, 'err', err);
    if ~isempty(S)
        S = rmfield(S, 'P');  r.S = S;
        r.CaSA = interp1(sa.beat(2).t/1e3, sa.beat(2).Ca, tb, 'linear', 'extrap');
    end
    save(f, 'r');
    if isempty(err)
        fprintf('%-30s %5.0f s | Fdev %.2f kPa, Fdia %.2f, drift %.2e | max|Ca-CaSA| %.2e uM\n', k.label, r.wall, ...
            max(Fdev), S.Fdia, S.drift, max(abs(S.Ca(:) - r.CaSA(:))));
    else
        fprintf('%-30s FAILED: %s\n', k.label, err);
    end
end
end
