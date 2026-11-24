function runCoupled(idx)
%RUNCOUPLED  Paced isometric twitch of the XB model driven by the TW_ATP myocyte (CaInput 'tworld',
%   coupling 'land': the 93 T-World states inside the XB ODE, open loop) for the cases IDX of
%   coupledCases. T-World starts from the 150-beat 1 Hz steady state of runStandalone (sa_*.mat)
%   at the same ATP/ADP/pH, so the coupled beats need no further T-World pre-pacing; the XB model
%   gets simHumanTwitch's default 2 s diastole + 3 beats. Cached as results/<label>.mat.
here = fileparts(mfilename('fullpath'));
twPaths('twatp');
K = coupledCases();
if nargin < 1 || isempty(idx), idx = 1:numel(K); end
[ATP, ADP] = atpLadder();
tb = (0:0.001:0.999)';
for i = idx
    k = K(i);  f = fullfile(here, 'results', [k.label '.mat']);
    if exist(f, 'file'), fprintf('%s cached\n', k.label); continue; end
    sa = load(fullfile(here, 'results', sprintf('sa_twatp_A%g_pH%.1f.mat', k.aCa, k.pH)));  sa = sa.r;
    base = loadParams(['human_twitch_' k.core]);
    base.UseAtpDetachR2 = 1;  base.UseAtpKmsrd = 1;       % both ATP gates; exactly 1 at MgATP_ref = 8
    base.MgATP = k.aXB;
    th = struct('CaInput', 'tworld', 'CaTworld', struct('coupling', 'land', 'model', @model_TWorld_shock_ATP, ...
        'X0', sa.Xend, 'param', struct('extraParams', [0 310 k.pH k.aCa ADP(ATP == k.aCa) 0])));
    tic;
    try
        [Fdev, S] = simHumanTwitch(base, th, tb);  err = '';
    catch ME
        Fdev = [];  S = [];  err = ME.message;
    end
    r = struct('k', k, 'wall', toc, 'tb', tb, 'Fdev', Fdev, 'err', err);
    if ~isempty(S)
        S = rmfield(S, 'P');  r.S = S;
        r.CaSA = interp1(sa.beat(2).t/1e3, sa.beat(2).Ca, tb, 'linear', 'extrap');   % standalone beat (uM)
    end
    save(f, 'r');
    if isempty(err)
        fprintf('%-34s %5.0f s | Fdev %.2f kPa, Fdia %.2f, drift %.2e | max|Ca-CaSA| %.2e uM\n', k.label, r.wall, ...
            max(Fdev), S.Fdia, S.drift, max(abs(S.Ca(:) - r.CaSA(:))));
    else
        fprintf('%-34s FAILED: %s\n', k.label, err);
    end
end
end
