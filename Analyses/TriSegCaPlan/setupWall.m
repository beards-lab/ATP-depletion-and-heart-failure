function P = setupWall(H, wall, opt)
%SETUPWALL  Params for one host-driven wall: incumbent + Ca gate + host kinematics + grid.
%
%   P = SETUPWALL(H, WALL, OPT)
%   OPT.paramSet  parameter snapshot (default 'robust1_seed', the incumbent; the first
%                 probes of 2026-11-09 used 'fastgate_best', the TriSegCoupling set)
%   OPT.xrate     multiplier on the cross-bridge cycle rates (kah, kamh, ka, kd, k1,
%                 k_1, k2; updateRates), default 1
%   OPT.tauV      time constant (s) of the filtered velocity signal, default 5e-3
%   OPT.delay     kinematic lag (s) applied to the host's Ls(t) and v(t), default 0
%   OPT.ss        strain bins: 60 = incumbent sliding window; larger = wider grid
%   OPT.velMode   'imposed' | 'velHS' | 'zero'   (see dPUdT_HostWall)
%   OPT.konMult   k_on multiplier on the Marzban 2020 thin filament (default 3, which
%                 gives a ~37 % twitch on this core, Analyses/CalciumTwitch E5)
%   OPT.relTol    ode15s RelTol (default 1e-3, the CalciumTwitch probe value)
%   The grid is symmetric about the starting (end-diastolic) SL; with ss = 60 the
%   window is the incumbent's (+-0.06 um half-sarcomere) and evaluateModel slides it.

if ~isfield(opt, 'paramSet'), opt.paramSet = 'robust1_seed'; end
if ~isfield(opt, 'xrate'),    opt.xrate = 1; end
if ~isfield(opt, 'tauV'),     opt.tauV = 5e-3; end
if ~isfield(opt, 'delay'),    opt.delay = 0; end
if ~isfield(opt, 'ss'),       opt.ss = 60; end
if ~isfield(opt, 'velMode'),  opt.velMode = 'imposed'; end
if ~isfield(opt, 'konMult'),  opt.konMult = 3; end
if ~isfield(opt, 'relTol'),   opt.relTol = 1e-3; end

P = loadParams(opt.paramSet);
if opt.xrate ~= 1          % scaled here: getParams applies params.xrate only with its modifiers
    for rn = {'kah', 'kamh', 'ka', 'kd', 'k1', 'k_1', 'k2'}
        if isfield(P, rn{1}) && isnumeric(P.(rn{1})), P.(rn{1}) = P.(rn{1})*opt.xrate; end
    end
end
P.SL0 = H.(wall).Ls(0);
P.Slim_l = P.SL0 - 1.2;  P.Slim_r = P.SL0 + 1.2;   % wide enough; the cap sets ss
P.MaxStrainArraySize = opt.ss;
P.Velocity = 0;
P.odeRelTol = opt.relTol;  P.odeAbsTol = 1e-6;
P.MaxRunTime = 3600;
P.BreakOnODEUnstable = false;
P = getParams(P, [], true, false);

P.host = H.(wall);  P.host.velMode = opt.velMode;  P.host.tauV = opt.tauV;
if opt.delay > 0
    Ls0 = H.(wall).Ls;  v0 = H.(wall).v;  dl = opt.delay;
    P.host.Ls = @(t) Ls0(t - dl);  P.host.v = @(t) v0(t - dl);
end
P.tf.k_on   = 101.185*opt.konMult;   % uM^-1 s^-1  (Marzban 2020 rat)
P.tf.k_off  = 723.852;               % s^-1
P.tf.K_coop = 9.6846;
P.tf.CaFcn  = H.Ca;
P.PU0 = [P.PU0(:); 1; 0];            % N = 1: fully non-permissive before the first beat; vf = 0
P.hostT = H.T;
end
