function [Ca, opt] = caTransient(kind, t, opt)
%CATRANSIENT  Candidate prescribed cytosolic [Ca2+](t) driving signals, in uM.
%
%   CA = CATRANSIENT(KIND, T) evaluates transient KIND at times T (s). Every
%   transient is periodic with period OPT.T and a beat starts at each multiple
%   of OPT.T. CA = CATRANSIENT(KIND, T, OPT) overrides the defaults below.
%
%   [CA, OPT] = CATRANSIENT(KIND, T, OPT) also returns OPT "prepared" for the
%   ODE hot path (getParams does this once): every default and derived
%   coefficient is resolved into OPT.c (the 'marzban2020' pchip included),
%   OPT.id replaces the string switch, and OPT.kinks lists the times within a
%   beat (s after onset) where the signal has a derivative kink, so the
%   integrator can be restarted there (evaluateModel merges them as knots).
%   Passing a prepared OPT back skips all of that: ~1 us instead of ~5 us.
%
%   KIND
%     'rice2008'     bi-exponential; rat, 24 C (Rice et al. 2008, CellML defaults):
%                    rise tau1 = 20 ms, decay tau2 = 110 ms, 0.09 -> 1.45 uM, onset 5 ms.
%     'marzban2020'  log-normal in cycle phase phi = mod(t,T)/T; rat, 2-10 Hz table,
%                    pchip in rate (Marzban 2020 PhysiomeS000007; identical table in
%                    Multiscale-Hx-P21-Model/model/compare_ca_curve.m).
%     'tewari_tanh'  tanh - exponential; human, Piacentino 1 Hz shape time-scaled by
%                    the Kang restitution Td(T) (cross-bridge-model XB.mo, block Ca2).
%     'fibersim2c'   FiberSim two-compartment SR/myofilament pulse, closed form
%                    (defaults of demo_files/electrical_stimulation/twitch).
%     'pCa'          constant 10^(6-OPT.pCa) uM.
%     'risedecay'    exponential rise x sigmoidal (Hill) decay, with beat-to-beat superposition:
%                    h(x) = (1 - exp(-x/tau_r))^m / (1 + (x/t50)^hdec), x = time since onset t0,
%                    summed over the current and the NSUP previous beats, then mapped so
%                    that the periodic signal's minimum is Cadia and its maximum Camax.
%                    Rise (tau_r, m) and decay (t50, hdec) are independent, so a transient
%                    can be prescribed from literature metrics (TTP, peak->50 %, peak->90 %).
%                    Defaults are the hiPSC-EHT set of Analyses/HumanTwitch/
%                    literature_ca_transient.md (de Lange 2021 anchor): Cadia 0.15,
%                    Camax 0.60 uM, tau_r 0.025 s, m 1, t50 0.265 s, hdec 3, t0 0.005 s,
%                    nsup 3 -> TTP ~88 ms, peak->50 % ~188 ms, peak->90 % ~480 ms.
%
%   OPT fields: T (period, s; default 1), pCa, and per-kind overrides.

if nargin < 3, opt = struct(); end
if ~isfield(opt, 'prepared')
    opt = prepare(kind, opt);
end
c  = opt.c;
tt = mod(t, c.T);                                 % time since beat onset

switch opt.id
    case 1  % rice2008
        x  = max(0, tt - c.t0);
        Ca = c.Cad + (c.Cam - c.Cad)/c.beta*(exp(-x/c.tau2) - exp(-x/c.tau1));

    case 2  % marzban2020
        phi = max(tt/c.T, 1e-9);
        Ca  = c.c4(1)./phi.*exp(-c.c4(2)*(log(phi) + c.c4(3)).^2) + c.c4(4);

    case 3  % tewari_tanh
        Ca  = max(0, c.a*(0.5*(1 - tanh(7.66/c.r*(tt - 0.289*c.r))) - exp(-162.2/c.r*tt))) ...
              + c.Ca0;

    case 4  % fibersim2c
        % Ca_myofil (M): first-order approach to C*ko/(ks+ko), ko = kl+ka during the
        % stimulus and kl after it (inlined: no anonymous function in the hot path)
        kon = c.kl + c.ka;  qon = c.C*kon/(c.ks+kon);  qoff = c.C*c.kl/(c.ks+c.kl);
        Ca  = 1e6*((tt <= c.D).*(qon + (c.m0 - qon).*exp(-(c.ks+kon)*tt)) ...
                 + (tt >  c.D).*(qoff + (c.mD - qoff).*exp(-(c.ks+c.kl)*(tt - c.D))));

    case 5  % pCa
        Ca = c.Ca*ones(size(t));

    case 6  % risedecay (periodic superposition of the last nsup beats)
        g = zeros(size(tt));
        for k = 0:c.nsup
            x = tt - c.t0 + k*c.T;
            g = g + (x > 0).*(1 - exp(-max(x, 0)/c.tau_r)).^c.m./(1 + (max(x, 0)/c.t50).^c.hdec);
        end
        Ca = c.Cad + (c.Cam - c.Cad)*(g - c.gmin)/(c.gmax - c.gmin);
end
end

function opt = prepare(kind, opt)
% Resolve defaults and derived constants once; see the header for KIND.
c.T = getOpt(opt, 'T', 1);
switch kind
    case 'rice2008'
        id = 1;
        c.tau1 = getOpt(opt, 'tau1', 0.020);  c.tau2 = getOpt(opt, 'tau2', 0.110);
        c.Cad  = getOpt(opt, 'Cadia', 0.09);  c.Cam  = getOpt(opt, 'Camax', 1.45);
        c.t0   = getOpt(opt, 't0', 0.005);
        r      = c.tau1/c.tau2;
        c.beta = r^(c.tau1/(c.tau2-c.tau1)) - r^(c.tau2/(c.tau2-c.tau1));  % normalises the peak to Camax
        kinks  = c.t0;                                                  % onset of the rise

    case 'marzban2020'
        id = 2;
        tab = [2      3      4      5      6      7      8      9      10;      % Hz
               0.0838 0.1306 0.1802 0.2557 0.3099 0.3613 0.4080 0.4539 0.4602; % a
               0.7513 0.8756 1.0274 1.4988 1.6107 1.6741 1.7902 2.1398 1.9832; % b
               2.2370 2.0486 1.9480 1.8520 1.6932 1.6773 1.5988 1.4952 1.4524; % c
               0.1865 0.1815 0.1709 0.1693 0.1610 0.1661 0.1425 0.1458 0.1222]; % Ca0
        c.c4  = pchip(tab(1,:), tab(2:5,:), 1/c.T);
        kinks = [];

    case 'tewari_tanh'
        id = 3;
        Td1   = 0.36/(0.39^2.2 + 1);                        % Kang 80 % relaxation time at 1 Hz
        c.r   = (0.36*c.T^2.2/(0.39^2.2 + c.T^2.2))/Td1;    % duration relative to 1 Hz
        c.a   = 0.56 + 0.1655*(1/c.T - 1);
        c.Ca0 = getOpt(opt, 'Ca0', 0.040);
        kinks = [];

    case 'fibersim2c'
        id = 4;
        c.C  = getOpt(opt, 'Ca_content', 1e-3);  c.ks = getOpt(opt, 'k_serca', 20);
        c.kl = getOpt(opt, 'k_leak', 6e-4);      c.ka = getOpt(opt, 'k_act', 8.2e-2);
        c.D  = getOpt(opt, 'stim', 0.010);
        myo  = @(x0, ko, s) c.C*ko/(c.ks+ko) + (x0 - c.C*ko/(c.ks+ko)).*exp(-(c.ks+ko)*s);
        c.m0 = c.C*c.kl/(c.ks+c.kl);                                         % resting
        c.mD = myo(c.m0, c.kl+c.ka, c.D);                                    % end of stimulus
        kinks = c.D;                                                         % end of the stimulus

    case 'pCa'
        id = 5;
        c.Ca  = 10^(6 - opt.pCa);
        kinks = [];

    case 'risedecay'
        id = 6;
        c.tau_r = getOpt(opt, 'tau_r', 0.025);  c.m     = getOpt(opt, 'm', 1);
        c.t50   = getOpt(opt, 't50', 0.265);    c.hdec  = getOpt(opt, 'hdec', 3);
        c.t0    = getOpt(opt, 't0', 0.005);
        c.Cad   = getOpt(opt, 'Cadia', 0.15);   c.Cam   = getOpt(opt, 'Camax', 0.60);
        c.nsup  = getOpt(opt, 'nsup', 3);
        % periodic min/max of the raw superposed shape, once, on a fine grid
        tg = linspace(0, c.T, 4001);  tg = tg(1:end-1);  g = zeros(size(tg));
        for k = 0:c.nsup
            x = tg - c.t0 + k*c.T;
            g = g + (x > 0).*(1 - exp(-max(x, 0)/c.tau_r)).^c.m./(1 + (max(x, 0)/c.t50).^c.hdec);
        end
        c.gmin = min(g);  c.gmax = max(g);
        kinks = c.t0;                                                        % onset

    otherwise
        error('caTransient:kind', 'Unknown transient ''%s''.', kind);
end
opt.c = c;  opt.id = id;  opt.kinks = kinks;  opt.prepared = true;
end

function v = getOpt(opt, name, default)
if isfield(opt, name), v = opt.(name); else, v = default; end
end
