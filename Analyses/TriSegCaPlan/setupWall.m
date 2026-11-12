function P = setupWall(H, wall, opt)
%SETUPWALL  Params for one host-driven wall: incumbent + Ca gate + host kinematics + grid.
%
%   P = SETUPWALL(H, WALL, OPT)
%   OPT.paramSet  incumbent snapshot (default 'fastgate_best', the TriSegCoupling set)
%   OPT.ss        strain bins: 60 = incumbent sliding window; larger = wider grid
%   OPT.velMode   'imposed' | 'velHS' | 'zero'   (see dPUdT_HostWall)
%   OPT.konMult   k_on multiplier on the Marzban 2020 thin filament (default 3, which
%                 gives a ~37 % twitch on this core, Analyses/CalciumTwitch E5)
%   OPT.relTol    ode15s RelTol (default 1e-3, the CalciumTwitch probe value)
%   The grid is symmetric about the starting (end-diastolic) SL; with ss = 60 the
%   window is the incumbent's (+-0.06 um half-sarcomere) and evaluateModel slides it.

if ~isfield(opt, 'paramSet'), opt.paramSet = 'fastgate_best'; end
if ~isfield(opt, 'ss'),       opt.ss = 60; end
if ~isfield(opt, 'velMode'),  opt.velMode = 'imposed'; end
if ~isfield(opt, 'konMult'),  opt.konMult = 3; end
if ~isfield(opt, 'relTol'),   opt.relTol = 1e-3; end

P = loadParams(opt.paramSet);
P.SL0 = H.(wall).Ls(0);
P.Slim_l = P.SL0 - 1.2;  P.Slim_r = P.SL0 + 1.2;   % wide enough; the cap sets ss
P.MaxStrainArraySize = opt.ss;
P.Velocity = 0;
P.odeRelTol = opt.relTol;  P.odeAbsTol = 1e-6;
P.MaxRunTime = 3600;
P.BreakOnODEUnstable = false;
P = getParams(P, [], true, false);

P.host = H.(wall);  P.host.velMode = opt.velMode;
P.tf.k_on   = 101.185*opt.konMult;   % uM^-1 s^-1  (Marzban 2020 rat)
P.tf.k_off  = 723.852;               % s^-1
P.tf.K_coop = 9.6846;
P.tf.CaFcn  = H.Ca;
P.PU0 = [P.PU0(:); 1];               % N = 1: fully non-permissive before the first beat
P.hostT = H.T;
end
