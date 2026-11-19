function W = setupXbWall(Ls0, tf, pas, opt, Ca0)
%SETUPXBWALL  Parameters and initial state of one cross-bridge wall (see xbWall).
%
% PROVENANCE: [GLUE]. The cross-bridge parameters are a snapshot from params/
% (opt.paramSet), passed through the model's own getParams. Two things are set here
% on top of the snapshot, both by decision 0004:
%   - pas (optional): the HOST's passive law, k*(L - Lsc0)^gamma (host model.m
%     l.262-264), placed in the model's OWN passive slot (k_pas, Lsc0, gamma). The
%     model evaluates it at the contractile length Lsc where the host uses Ls; the
%     two differ by the series-element length, 0.01 um at peak force.
%   - opt.nativeViscous = false switches off the model's Maxwell element, the
%     velocity-dependent part of its passive, which was fitted at high Ca.
%
%   W = setupXbWall(Ls0, tf, pas, opt)        state at the snapshot's initial condition
%   W = setupXbWall(Ls0, tf, pas, opt, Ca0)   pre-equilibrated opt.tPre s at Ls0 and Ca0
%   tf   struct k_on, k_off, K_coop (thin filament; k_on already carries opt.konMult)
%   pas  [] (native passive) or struct k, Lsc0, gamma
%   opt  paramSet, ss, tauV, relTol, absTol, nativeViscous, tPre, N0

P = loadParams(opt.paramSet);
P.SL0 = Ls0;  P.Slim_l = Ls0 - 1.2;  P.Slim_r = Ls0 + 1.2;   % wide limits; the cap sets the grid
P.MaxStrainArraySize = opt.ss;
P.Velocity = 0;  P.odeRelTol = opt.relTol;  P.odeAbsTol = opt.absTol;
P.MaxRunTime = 3600;  P.BreakOnODEUnstable = false;
if ~isempty(pas)
    P.k_pas = pas.k;  P.Lsc0 = pas.Lsc0;  P.gamma = pas.gamma;
end
if ~opt.nativeViscous, P.UseMaxwellDashpot = false; end
P = getParams(P, [], true, false);

W.P = P;
W.nCore = numel(P.PU0);  W.nBins = P.NumberOfStates*P.ss;
W.iSL = W.nBins + 3;  W.iLSE = W.nBins + 4;
W.nW  = W.nCore + 3;  W.nds = floor(P.ss/3);
W.tf  = tf;  W.tauV = opt.tauV;
core0 = P.PU0(:);
W.w0  = [core0; Ls0 - core0(W.iLSE); Ls0; opt.N0];

if nargin > 4 && opt.tPre > 0
    % frozen length and Ca; the strain window may still slide while force develops
    R = runXbWallOpenLoop(W, @(t) Ls0, @(t) Ca0, [-opt.tPre 0], struct('dtOut', opt.tPre));
    W = R.W;  W.w0 = R.wEnd;
end
end
